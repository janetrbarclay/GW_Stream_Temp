import os

from gw_stream_temp.fetch_data import get_NHM_gis_data
from gw_stream_temp.preprocess_data import compile_catchment_discharge, get_catchment_nodes, compile_model_outputs, make_model_shapefile

modelName = "MONTAGUE_drb1.04_mf6_250"
flowModelName = "drb1.04"
modelDir = "/home/jbarclay/GW_Models/DRB/MONTAGUE/model/MONTAGUE_drb1.04_mf6_250_SY05"
outDir = "out_{}".format(modelName)
rasterPath = '{}/MONTAGUE_250_idomain.tif'.format(modelDir)

rule all:
    input:
        "{}/CatchmentDischarge.csv".format(outDir)
        
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
