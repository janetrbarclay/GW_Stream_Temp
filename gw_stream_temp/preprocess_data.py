# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 08:29:46 2021

@author: jbarclay
"""

import os, math
import flopy as fp
from osgeo import gdal, osr
import geopandas as gpd
import numpy as np
from geopandas.tools import sjoin
import pandas as pd
from copy import deepcopy

  

def make_model_shapefile(modelpth="./", thisModelName="MONTAGUE_drb1.04_mf6_250", flow_model_name = "MONTAGUE_drb1.04_mf6_250", out_file ="MONTAGUE_drb1.04_mf6_250_grid.shp", rasterPath = None):
    """
    creates a shapefile of a groundwater model grid
    :param modelpth: [str] file path to model files
    :param thisModelName: [str] base name of the model files
    :param flow_model_name: [str] name of the groundwater flow model within the modflow simulation
    :param out_file: [str] name of the resulting shapefile
    :param rasterPath: [str] path to a geotif file that can be used for geolocating the model shapefile
    """
    
    #################################
    ## get the modflow version
    nam_file = os.path.join(modelpth,'{}.nam'.format(thisModelName))
    
    o_file = open(nam_file)
    line = o_file.readlines()
    o_file.close()
    if any([x.startswith("BEGIN") for x in line]):
        mf_version = "mf6"
        mfpth = "mf6"
    elif any(["NWT" in x for x in line]):
        mf_version = "mfnwt"
        mfpth = "mfnwt"

    
    #load the model
    if mf_version=="mfnwt":
        ml = fp.modflow.Modflow.load(nam_file, version='mfnwt', exe_name=mfpth,  verbose=False, model_ws=modelpth, load_only=None)
    elif mf_version=="mf6":
        sim = fp.mf6.MFSimulation.load(thisModelName, version='mf6', exe_name=mfpth, sim_ws=modelpth)
        ml = sim.get_model(flow_model_name) 
        
    #get the spatial reference if needed
    if ml.modelgrid.xoffset==0.0 and rasterPath is not None:
        set_model_spatial_reference(ml, rasterPath)
    
    #export the basic grid as a shapefile
    fp.export.shapefile_utils.model_attributes_to_shapefile(out_file, ml,package_names=['dis'])
    

def set_model_spatial_reference(ml, rasterPath):
    """
    sets the spatial reference for the groundwater model if it isn't included in the model files
    :param ml: [modflow model object] flopy modflow model object
    :param rasterPath: [str] filepath to a geotif for geolocating the model grid
    """
    tiffDS = gdal.Open(rasterPath)
    gt = tiffDS.GetGeoTransform()
    
    # =============================================================================
    #     gt : 6-element geotransform list [C, A, B, F, E, D]. Gives the coordinates of one pixel
    #         (the upper left pixel). If there is no rotation, B=D=0. If cells are square, A=-E.   
    #         Letter designations come from the original documentation.
    #         
    #         C = x coordinate in map units of the upper left corner of the upper left pixel
    #         A = distance from C along x axis to upper right pixel corner of the upper left pixel
    #         B = distance from C along x axis to lower left pixel corner of the upper left pixel,
    #         F = y coordinate in map units of the upper left corner of the upper left pixel
    #         E = distance from C along y axis to lower left pixel corner of the upper left pixel
    #         D = distance from C along y axis to upper right pixel corner of the upper left pixel
    # =============================================================================
    
    #x coordinate of the lower left corner
    XLL = gt[0]+ml.modelgrid.nrow*gt[2]
    #y coordinate of the lower left corner
    YLL = gt[3]+ml.modelgrid.nrow*gt[5]
    #angle of rotation, in degrees counter-clockwise around the lower left corder
    angrot = math.atan(-1*gt[2]/gt[5])*360/2/math.pi
    
    # xul = gt[0]
    # yul = gt[3]
    # XLL = xul-math.sin(-1*angrot/180*math.pi)*ml.modelgrid.delr[0]*ml.modelgrid.nrow
    # YLL = yul-math.cos(-1*angrot/180*math.pi)*ml.modelgrid.delr[0]*ml.modelgrid.nrow
    
    #get the projection
    prj = tiffDS.GetProjection()
    srs = osr.SpatialReference()
    srs.ImportFromWkt(prj)
    prj4 = (srs.ExportToProj4()).split('+')
    prj4 = dict([item.split('=') for item in prj4 if len(item.split('=')) == 2])
    
    
    ml.modelgrid.set_coord_info(xoff=XLL,yoff=YLL,angrot=angrot,proj4=prj4,merge_coord_info=False)
    
    
    
def compile_model_outputs(modelpth="./", thisModelName="MONTAGUE_drb1.04_mf6_250", out_file = "resultsAgg.csv"):
    """
    creates a csv of groundwater discharge for each node in the groundwater model
    :param modelpth: [str] file path to model files
    :param thisModelName: [str] base name of the model files
    :param out_file: [str] name of the resulting csv file
    """
    #################################
    ## get the budget file path
    nam_file = os.path.join(modelpth,'{}.nam'.format(thisModelName))
    
    o_file = open(nam_file)
    line = o_file.readlines()
    o_file.close()
    
    if any([x.startswith("BEGIN") for x in line]):
        mf_version = "mf6"
        #open the OC file to get the budget path
        oc_file = [x for x in line if x.strip().startswith("OC6")][0].split()[1]
        o_file = open(os.path.join(modelpth,oc_file))
        line_oc = o_file.readlines()
        o_file.close()
        budget_file = [x for x in line_oc if "BUDGET" in x and "FILEOUT" in x][0].split()[-1]
        
    elif any(["NWT" in x for x in line]):
        mf_version = "mfnwt"
        budget_file = [x for x in line if x.startswith("DATA(BINARY)") and ".cb" in x][0].split()[2]
    
    
    #open the budget file
    cell_budget = fp.utils.CellBudgetFile(os.path.join(modelpth,budget_file))
    
    #get the list of gw discharge related records
    gw_dis_records = [x for x in cell_budget.get_unique_record_names() if any([y in x.decode() for y in ['DRN','RIV','GHB','SFR','HEAD DEP BOUNDS'] ])]
    
    #get the last timestep for each period
    kstpkperLst = [[x for x in cell_budget.get_kstpkper() if x[1] == y][-1] for y in range(cell_budget.nper)]
    
    resultsDF = pd.DataFrame(columns=['node','record','per','q'])
    for thisRec in gw_dis_records:
        for x in kstpkperLst:
            tempDF = pd.DataFrame(cell_budget.get_data(text=thisRec,kstpkper=x)[0])
            tempDF['per']=x[1]
            tempDF['record']=thisRec.decode().strip()
            resultsDF=resultsDF.append(tempDF, ignore_index=True)
            
            
    #aggregate the results by type and then overall
    resultsAgg = resultsDF[['node','q','record']].groupby(by=['record','node'],as_index=False).mean().groupby(by='node',as_index=False).sum()
    
    resultsAgg.to_csv(out_file)

def get_catchment_nodes(NHM_gdb=r'C:\Users\jbarclay\OneDrive - DOI\StreamTemp\Analysis\Data\GFv1.1.gdb',model_shapefile="modelshapefile.shp", model_epsg=None, local_out_file = 'localCatchDict.npy', upstream_out_file = 'upstreamCatchDict.npy'):
    """
    creates a numpy dictionaries (saved out to files) of 1) the model nodes within each catchment and 2) the catchments in or upstream of each catchment
    :param NHM_gdb: [str] file path to geodatabase of the national hydrologic model framework files
    :param model_shapefile: [str] file path of the groundwater model grid shapefile
    :param model_eps: [str] epsg code for the model grid shapefile
    :param local_out_file: [str] path to the output file for saving the dictionary of gw model nodes within each catchment
    :param upstream_out_file: [str] path to the output file for saving the dictionary of catchments in or upstream of each catchment
    """    
    
    #open files
    catchmentGDF = gpd.read_file(NHM_gdb, layer='nhru_v1_1_simp')
    reachGDF = gpd.read_file(NHM_gdb,layer='nsegment_v1_1')
    modelGDF = gpd.read_file(model_shapefile)
    
    #create a uniform ibound column
    modelGDF['ibound']= modelGDF[[x for x in modelGDF.columns if "ibound" in x or "idomain" in x]].max(axis=1)

    
    #make a dictionary matching the v1 and v1.1 id's
    matchDict = {row.nsegment_v1_1:row.seg_id_nhm for row in reachGDF.itertuples()}
    matchDict[0]=0
    
    #add the seg_id_nat
    catchmentGDF['seg_id_nhm'] = [matchDict[x] for x in catchmentGDF.hru_segment_v1_1]
    reachGDF['toseg_id_nhm'] = [matchDict[x] for x in reachGDF.tosegment_v1_1]
      
    
    if modelGDF.crs is None and model_epsg is not None:
        modelGDF.set_crs(epsg=model_epsg, inplace=True)
    
    #reproject the model grid to the catchments
    modelGDF.to_crs(crs=catchmentGDF.crs, inplace=True)
    
    #join the catchments and the model grid
    joinDF = sjoin(catchmentGDF,modelGDF.loc[modelGDF.idomain_1==1],how="inner")
    
    #calculate the inbound fraction of area to filter our catchments that are largely outside the model area
    fracArea = joinDF[['seg_id_nhm','ibound']].groupby('seg_id_nhm',as_index=False).mean()
    joinDF = joinDF.loc[joinDF.seg_id_nhm.isin(fracArea.seg_id_nhm.loc[fracArea.ibound>0.33])]
    
    #dictionary matching seg_id_nhm : node
    catchDict = {x:joinDF.loc[(joinDF.seg_id_nhm==x),"node"].values for x in np.unique(joinDF.seg_id_nhm)}
    np.save(local_out_file,catchDict)    
           
    get_upstream_catchments(reachGDF,catchDict, upstream_out_file)
    

    

def get_upstream_catchments(reachGDF,catchDict, out_file = 'upstreamCatchDict.npy'):
    """
    creates a numpy dictionarys (saved out to files) of the catchments upstream of each catchment
    :param reachGDF: [str] geodataframe of reaches within the national hydrologic model framework
    :param catchDict: [str] dictionary of model nodes for each catchment
    :param out_file: [str] output file for saving the dictionary of catchments in or upstream of each catchment
    """   
    
    reachDF_scratch = deepcopy(reachGDF.loc[reachGDF.seg_id_nhm.isin(catchDict.keys()),['nsegment_v1_1','tosegment_v1_1','seg_id_nhm','toseg_id_nhm','Version']])

    
    #dictionary matching seg_id_nhm : all upstream seg_id_nhm (including the current seg_id_nhm)
    upStreamDict = {x:[x] for x in np.unique([reachDF_scratch.seg_id_nhm,reachDF_scratch.toseg_id_nhm])}

    
    i=0
    while reachDF_scratch.shape[0]>0:
        i = i+1
        thisGroup = reachDF_scratch.loc[~(reachDF_scratch.seg_id_nhm.isin(np.unique(reachDF_scratch.toseg_id_nhm))),['seg_id_nhm','toseg_id_nhm']]
        for row in thisGroup.itertuples():
            if row.toseg_id_nhm!=0:

                    upStreamDict[row.toseg_id_nhm].extend(upStreamDict[row.seg_id_nhm])

        reachDF_scratch = reachDF_scratch.loc[~(reachDF_scratch.seg_id_nhm.isin(thisGroup.seg_id_nhm))]

    np.save(out_file,upStreamDict)


def compile_catchment_discharge(node_discharge_file="resultsAgg.csv", catchDictFile = 'localCatchDict.npy', upStreamDictFile = 'upstreamCatchDict.npy', out_file = "CatchmentDischarge.csv"):
    """
    compiles the groundwater discharge (total and as a percent of upstream - including the local catchment - discharge) for each catchment
    :param node_discharge_file: [str] csv file of groundwater discharge for each model node
    :param catchDictFile: [str] file path to dictionary of model nodes for each catchment
    :param upStreamDictFile: [str] file path to the dictionary of catchments in or upstream of each catchment
    :param out_file: [str] output csv file of the compiled discharge
    """       
    node_discharge = pd.read_csv(node_discharge_file)

    catchDict = np.load(catchDictFile, allow_pickle=True)[()]
    upStreamDict = np.load(upStreamDictFile, allow_pickle=True)[()]
    
    localDis = [(x,np.sum(node_discharge.loc[node_discharge.node.isin(catchDict[x]),"q"])) for x in catchDict.keys()if x !=0]
    dischargeDF = pd.DataFrame(localDis,columns=['seg_id_nat','q_local'])
    dischargeDF['q_all']=[np.sum(dischargeDF.q_local.loc[dischargeDF.seg_id_nat.isin(upStreamDict[x])]) for x in dischargeDF.seg_id_nat]
    dischargeDF['Per_Local']=dischargeDF['q_local']/dischargeDF['q_all']*100
    
    dischargeDF.to_csv(out_file)
