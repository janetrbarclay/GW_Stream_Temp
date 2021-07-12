import os, sys
import pandas as pd
import xarray as xr

sys.path.append("../river-dl")

from gw_stream_temp.fetch_data import get_NHM_gis_data
from gw_stream_temp.preprocess_data import compile_catchment_discharge, get_catchment_nodes, compile_model_outputs, make_model_shapefile
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
outDir = config['outDir']
rasterPath = config['rasterPath']

rule all:
    input:
        "{}/prepped.npz".format(outDir)
        
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
