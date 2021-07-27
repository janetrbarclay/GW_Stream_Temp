import os, sys
import pandas as pd
import xarray as xr

sys.path.append("../river-dl")

from gw_stream_temp.fetch_data import get_NHM_gis_data
from gw_stream_temp.preprocess_modflow_data import compile_catchment_discharge, get_catchment_nodes, compile_model_outputs, make_model_shapefile
from gw_stream_temp.preprocess_modpath_data import make_modpath_model
from gw_stream_temp.visualize import plot_by_perlocal
from river_dl.preproc_utils import prep_data
from river_dl.evaluate import combined_metrics
from river_dl.postproc_utils import plot_obs
from river_dl.predict import predict_from_io_data
from river_dl.train import train_model
from river_dl import loss_functions as lf
from river_dl.gw_utils import prep_annual_signal_data, calc_pred_ann_temp,calc_gw_metrics

modelName = config["modelName"]
flowModelName = config["flowModelName"]

modelDir = config['modelDir']
outDir = config['out_dir']
rasterPath = config['rasterPath']
workingDir = os.getcwd()

loss_function = lf.multitask_rmse(config['lambdas'])

module other_workflow:
    snakefile: "../river-dl/Snakefile"
    config: config
       
use rule * from other_workflow as other_*

#this allows us to import all the rules from Snakefile but write a custom train_model_local_or_cpu rule
use rule train_model_local_or_cpu from other_workflow as other_train_model_local_or_cpu with:
    output:
        ""
use rule prep_io_data from other_workflow as other_prep_io_data with:
    output:
        ""
use rule prep_ann_temp from other_workflow as other_prep_ann_temp with:
    output:
        ""
use rule plot_prepped_data from other_workflow as other_plot_prepped_data with:
    output:
        ""
        
#modify rule all to include the additional gw output files        
use rule all from other_workflow as other_all with:
    input:
        expand("{outdir}/GW_summary.csv",outdir=outDir),
        expand("{outdir}/GW_stats_{partition}.csv",
                outdir=outDir,
                partition=['trn', 'tst','val']
        ),
        expand("{outdir}/{metric_type}_metrics.csv",
                outdir=outDir,
                metric_type=['overall', 'month', 'reach', 'month_reach'],
        ),
        expand("{outdir}/GW_{model_metric}.png", outdir=outDir, model_metric=['Per_Local','q_all','q_std','q_std_per']),
        '{}/{}.mpend'.format(modelDir,flowModelName)
 
rule get_NHM_data:
    output:
        directory('data_NHM/GFv1.1.gdb')
    run:
        get_NHM_gis_data(item='5e29d1a0e4b0a79317cf7f63',filenameLst=['GFv1.1.gdb.zip'], destination='data_NHM')
        
rule make_model_shapefile:
    output:
        '{}/modelGrid.shp'.format(outDir)
    run:
        make_model_shapefile(modelDir,modelName, flowModelName, output[0], rasterPath)
        
rule get_model_outputs:
    output:
        '{}/Model_Outputs.csv'.format(outDir)
    run:
        compile_model_outputs(modelDir,modelName,output[0])
        
rule create_catchment_dictionaries:
    input:
        'data_NHM/GFv1.1.gdb',
        '{}/modelGrid.shp'.format(outDir),
    output:
        '{}/local_catch_dict.npy'.format(outDir),
        '{}/upstream_catch_dict.npy'.format(outDir),
    run:
        get_catchment_nodes(input[0], input[1],5070,output[0],output[1])
        
rule compile_discharge:
    input:
        '{}/Model_Outputs.csv'.format(outDir),
        '{}/local_catch_dict.npy'.format(outDir),
        '{}/upstream_catch_dict.npy'.format(outDir),
    output:
        '{}/CatchmentDischarge.csv'.format(outDir)
    run:
        compile_catchment_discharge(input[0],input[1],input[2],output[0])
        
rule write_modpath_files:
    output:
        '{}/{}.mpsim'.format(modelDir,flowModelName)
    run:
        make_modpath_model(modelDir, modelName, flowModelName)

rule run_modpath:
    input:
        '{}/{}.mpsim'.format(modelDir,flowModelName)
    output:
        '{}/{}.mpend'.format(modelDir,flowModelName)
    shell:
        'cd {modelDir}; mpath7 {input}; cd {workingDir}' 
    
        
def get_segment_list(gw_file_in, pretrain_file_in, seg_col="seg_id_nat"):
    gwDF = pd.read_csv(gw_file_in)
    ds_pre = xr.open_zarr(pretrain_file_in)
    seg_list = [x for x in gwDF[seg_col] if x in ds_pre[seg_col]]
    return seg_list


rule prep_io_data:
    input:
         config['obs_temp'],
         config['obs_flow'],
         config['sntemp_file'],
         config['dist_matrix'],
         '{}/CatchmentDischarge.csv'.format(outDir),
    output:
        "{outdir}/prepped.npz"
    run:
        prep_data(input[0], input[1], input[2], input[3],
                  x_vars=config['x_vars'],
                  catch_prop_file=None,
                  exclude_file=None,
                  train_start_date=config['train_start_date'],
                  train_end_date=config['train_end_date'],
                  val_start_date=config['val_start_date'],
                  val_end_date=config['val_end_date'],
                  test_start_date=config['test_start_date'],
                  test_end_date=config['test_end_date'],
                  primary_variable=config['primary_variable'],
                  log_q=False, segs=get_segment_list(input[4],input[2]),
                  out_file=output[0])

rule prep_ann_temp:
    input:
         config['obs_temp'],
         config['sntemp_file'],
         "{outdir}/prepped.npz",
         '{}/CatchmentDischarge.csv'.format(outDir),
    output:
        "{outdir}/prepped_withGW.npz",
    run:
        prep_annual_signal_data(input[0], input[1], input[2],
                  train_start_date=config['train_start_date'],
                  train_end_date=config['train_end_date'],
                  val_start_date=config['val_start_date'],
                  val_end_date=config['val_end_date'],
                  test_start_date=config['test_start_date'],
                  test_end_date=config['test_end_date'], 
                  gwVarList = config['gw_vars'],
                   segs=get_segment_list(input[3],input[1]),
                  out_file=output[0])
                  
# use "train" if wanting to use GPU on HPC
#rule train:
#    input:
#        "{outdir}/prepped_withGW.npz"
#    output:
#        directory("{outdir}/trained_weights/"),
#        directory("{outdir}/pretrained_weights/"),
#    params:
#        # getting the base path to put the training outputs in
#        # I omit the last slash (hence '[:-1]' so the split works properly
#        run_dir=lambda wildcards, output: os.path.split(output[0][:-1])[0],
#        pt_epochs=config['pt_epochs'],
#        ft_epochs=config['ft_epochs'],
#        lamb=config['lamb'],
#        lamb2=config['lamb2'],
#        lamb3=config['lamb3'],
#        loss = config['loss_type'],
#    shell:
#        """
#        module load analytics cuda10.1/toolkit/10.1.105 
#        run_training -e /home/jbarclay/.conda/envs/rgcn --no-node-list "python {code_dir}/train_model_cli.py -o {params.run_dir} -i {input[0]} -p {params.pt_epochs} -f {params.ft_epochs} --lamb {params.lamb} --lamb2 {params.lamb2} --lamb3 {params.lamb3} --model rgcn --loss {params.loss} -s 135"
#        """
 
 
# use "train_model" if wanting to use CPU or local GPU
rule train_model_local_or_cpu:
    input:
        "{outdir}/prepped_withGW.npz"
    output:
        directory("{outdir}/trained_weights/"),
        directory("{outdir}/pretrained_weights/"),
    params:
        # getting the base path to put the training outputs in
        # I omit the last slash (hence '[:-1]' so the split works properly
        run_dir=lambda wildcards, output: os.path.split(output[0][:-1])[0],
    run:
        train_model(input[0], config['pt_epochs'], config['ft_epochs'], config['hidden_size'],
                    loss_func=loss_function, out_dir=params.run_dir, model_type='rgcn', num_tasks=2, loss_type=config['loss_type'], lamb2=config['lamb2'],lamb3=config['lamb3'])



 
rule plot_discharge:
    input:
        "{outdir}/CatchmentDischarge.csv",
        "{outdir}/reach_metrics.csv",
        "{outdir}/prepped_withGW.npz",
    output:
        "{outdir}/GW_{model_metric}.png"
    run:
        plot_by_perlocal(input[0],input[1],input[2],output[0], plotCol = wildcards.model_metric, axisTitle=wildcards.model_metric)
        
rule compile_pred_GW_stats:
    input:
        "{outdir}/prepped_withGW.npz",
        "{outdir}/trn_preds.feather",
        "{outdir}/tst_preds.feather",
        "{outdir}/val_preds.feather"
    output:
        "{outdir}/GW_stats_trn.csv",
        "{outdir}/GW_stats_tst.csv",
        "{outdir}/GW_stats_val.csv",
    run: 
        calc_pred_ann_temp(input[0],input[1],input[2], input[3], output[0], output[1], output[2])
        
rule calc_gw_summary_metrics:
    input:
        "{outdir}/GW_stats_trn.csv",
        "{outdir}/GW_stats_tst.csv",
        "{outdir}/GW_stats_val.csv",
    output:
        "{outdir}/GW_summary.csv",
        "{outdir}/GW_scatter.png",
        "{outdir}/GW_boxplot.png",
    run:
        calc_gw_metrics(input[0],input[1],input[2],output[0], output[1], output[2])
