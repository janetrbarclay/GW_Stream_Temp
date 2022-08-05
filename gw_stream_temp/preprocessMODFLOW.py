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
import sys                
                    

def load_modflow_model(modelpth="./", thisModelName="MONTAGUE_drb1.04_mf6_250", flow_model_name = "MONTAGUE_drb1.04_mf6_250"):
    """
    loads a modflow groundwater model
    :param modelpth: [str] file path to model files
    :param thisModelName: [str] base name of the model files
    :param flow_model_name: [str] name of the groundwater flow model within the modflow simulation
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
        if flow_model_name in sim.model_names:
            ml = sim.get_model(flow_model_name)
        else:
            flow_model_name = [x for x in sim.model_names if sim.model_dict[x].model_type=='gwf']
            assert len(flow_model_name)>1, "Multiple flow models are present and none match the given flow model name"
            ml = sim.get_model(flow_model_name[0])
            
        
    return ml

def make_model_shapefile(modelpth="./", thisModelName="MONTAGUE_drb1.04_mf6_250", flow_model_name = "MONTAGUE_drb1.04_mf6_250", out_file ="MONTAGUE_drb1.04_mf6_250_grid.shp", rasterPath = None):
    """
    creates a shapefile of a groundwater model grid
    :param modelpth: [str] file path to model files
    :param thisModelName: [str] base name of the model files
    :param flow_model_name: [str] name of the groundwater flow model within the modflow simulation
    :param out_file: [str] name of the resulting shapefile
    :param rasterPath: [str] path to a geotif file that can be used for geolocating the model shapefile
    """
    
    ml = load_modflow_model(modelpth,thisModelName,flow_model_name)
        
    #get the spatial reference if needed
    rasterPrj = None
    if ml.modelgrid.xoffset==0.0 and rasterPath is not None:
        rasterPrj = set_model_spatial_reference(ml, rasterPath)
    
    #export the basic grid as a shapefile
    fp.export.shapefile_utils.model_attributes_to_shapefile(out_file, ml,package_names=['dis'])
    
    #apply the spatial reference to the shapefile, if needed
    if not os.path.exists(out_file.replace("shp","prj")):
        modelGDF = gpd.read_file(out_file)
        modelGDF.set_crs(rasterPrj,inplace=True)
        modelGDF.to_file(out_file)
    

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
    
    return prj
    
    
    
def compile_model_outputs(modelpth="./", outputpth=None, thisModelName="MONTAGUE_drb1.04_mf6_250", out_file = "resultsAgg.feather"):
    """
    creates a feather of groundwater discharge for each node in the groundwater model
    :param modelpth: [str] file path to model files
    :param thisModelName: [str] base name of the model files
    :param out_file: [str] name of the resulting feather file
    """
    #################################
    ## get the budget file path
    if not outputpth:
        outputpth = modelpth
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
    cell_budget = fp.utils.CellBudgetFile(os.path.join(outputpth,budget_file))
    
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
    
    #also get the month-to-month variation in q
    resultsAggSTD = resultsDF[['node','q','per']].groupby(by=['per','node'],as_index=False).sum()[['node','q']].groupby(by="node",as_index=False).std().rename(columns={'q':'q_std'})
    
    resultsAgg=resultsAgg.merge(resultsAggSTD)
    
    #and normalize the std by the q
    resultsAgg['q_std_per']=resultsAgg['q_std']/resultsAgg['q']
    
    resultsAgg.to_feather(out_file)

def compile_catchment_params(modelpth="./", thisModelName="MONTAGUE_drb1.04_mf6_250", flow_model_name = "MONTAGUE_drb1.04_mf6_250",paramDict={'hk':['mean',[0]]}, reachIdx = "seg_id_nat", catchDictFile = 'localCatchDict.npy', out_file = "CatchmentDischarge.feather"):
    """
    compiles model parameters for each catchment
    :param paramDict: dictionary of model parameters to extract
    :param catchDictFile: [str] file path to dictionary of model nodes for each catchment
    :param reachIdx: [str] reach identifier
    :param out_file: [str] output feather file of the compiled discharge
    """  
    ml = load_modflow_model(modelpth,thisModelName,flow_model_name)
    mf_version = ml.version
    
    catchDict = np.load(catchDictFile, allow_pickle=True)[()]
    resultsDF = pd.DataFrame({reachIdx:[int(x) for x in catchDict.keys()if x !=0]})
    
    if mf_version=="mf6":
        paramsDF = pd.DataFrame({'node':[x+1 for x in range(ml.dis.nrow.array*ml.dis.ncol.array)],'ibound':ml.dis.idomain.array[0].ravel()})
    
        for thisVar in paramDict.keys():
            if thisVar=="hk":
                varArray = ml.npf.k.array
                for thisLayer in paramDict[thisVar][1]:
                    paramsDF[thisVar+"_"+str(thisLayer)]=np.ravel(varArray[int(thisLayer)])
                    if paramDict[thisVar][0]=="mean":
                        resultsDF[thisVar+"_"+str(thisLayer)]=[np.mean(paramsDF.loc[paramsDF.node.isin(catchDict[x]),thisVar+"_"+str(thisLayer)]) for x in catchDict.keys()if x !=0]
    
    resultsDF.to_feather(out_file)