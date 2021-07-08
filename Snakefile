import os

from gw_stream_temp.fetch_data import get_NHM_gis_data
from gw_stream_temp.preprocess_data import compile_catchment_discharge, get_catchment_nodes, compile_model_outputs, make_model_shapefile

modelName = "CoastalCT_ext"
modelDir = "data_model_CoastalCT_ext/layers_GHB_sea_septic_DRN0_mf6_SFR_PEST_layerCOF"
outDir = "out_{}".format(modelName)
rasterPath = '{}/ibound.tif'.format(modelDir)

rule all:
    input:
        "{}/CatchmentDischarge.csv".format(outDir)
        
rule get_NHM_data:
    output:
        directory('{}/GFv1.1.gdb'.format(outDir))
    run:
        get_NHM_gis_data(item='5e29d1a0e4b0a79317cf7f63',filenameLst=['GFv1.1.gdb.zip'], destination=outDir)
        
rule make_model_shapefile:
    output:
        '{}/modelGrid.shp'.format(outDir)
    run:
        make_model_shapefile(modelDir,modelName, output[0], rasterPath)
        
rule get_model_outputs:
    output:
        '{}/Model_Outputs.csv'.format(outDir)
    run:
        compile_model_outputs(modelDir,modelName,output[0])
        
rule create_catchment_dictionaries:
    input:
        '{}/GFv1.1.gdb'.format(outDir),
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