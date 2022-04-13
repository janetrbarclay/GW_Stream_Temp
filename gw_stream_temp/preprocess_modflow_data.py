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
import gc
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

# =============================================================================
# def get_NHM_catchment_nodes(NHM_gdb=r'C:\Users\jbarclay\OneDrive - DOI\StreamTemp\Analysis\Data\GFv1.1.gdb',model_shapefile="modelshapefile.shp", local_out_file = 'localCatchDict.npy', upstream_out_file = 'upstreamCatchDict.npy', model_epsg=None):
#     """
#     creates a numpy dictionaries (saved out to files) of 1) the model nodes within each catchment and 2) the catchments in or upstream of each catchment
#     :param NHM_gdb: [str] file path to geodatabase of the national hydrologic model framework files
#     :param model_shapefile: [str] file path of the groundwater model grid shapefile
#     :param model_eps: [str] epsg code for the model grid shapefile
#     :param local_out_file: [str] path to the output file for saving the dictionary of gw model nodes within each catchment
#     :param upstream_out_file: [str] path to the output file for saving the dictionary of catchments in or upstream of each catchment
#     """    
#     
#     #open files
#     catchmentGDF = gpd.read_file(NHM_gdb, layer='nhru_v1_1_simp')
#     reachGDF = gpd.read_file(NHM_gdb,layer='nsegment_v1_1')
#     modelGDF = gpd.read_file(model_shapefile)
#     
#     #create a uniform ibound column
#     modelGDF['ibound']= modelGDF[[x for x in modelGDF.columns if "ibound" in x or "idomain" in x]].max(axis=1)
# 
#     
#     #make a dictionary matching the v1 and v1.1 id's
#     matchDict = {row.nsegment_v1_1:row.seg_id_nhm for row in reachGDF.itertuples()}
#     matchDict[0]=0
#     
#     #add the seg_id_nat
#     catchmentGDF['seg_id_nhm'] = [matchDict[x] for x in catchmentGDF.hru_segment_v1_1]
#     reachGDF['toseg_id_nhm'] = [matchDict[x] for x in reachGDF.tosegment_v1_1]
#       
#     
#     if modelGDF.crs is None and model_epsg is not None:
#         modelGDF.set_crs(epsg=model_epsg, inplace=True)
#     
#     #reproject the model grid to the catchments
#     modelGDF.to_crs(crs=catchmentGDF.crs, inplace=True)
#     
#     #join the catchments and the model grid
#     joinDF = sjoin(catchmentGDF,modelGDF.loc[modelGDF.idomain_1==1],how="inner")
#     
#     #calculate the inbound fraction of area to filter our catchments that are largely outside the model area
#     fracArea = joinDF[['seg_id_nhm','ibound']].groupby('seg_id_nhm',as_index=False).mean()
#     joinDF = joinDF.loc[joinDF.seg_id_nhm.isin(fracArea.seg_id_nhm.loc[fracArea.ibound>0.33])]
#     
#     #dictionary matching seg_id_nhm : node
#     catchDict = {x:joinDF.loc[(joinDF.seg_id_nhm==x),"node"].values for x in np.unique(joinDF.seg_id_nhm)}
#     np.save(local_out_file,catchDict)    
#            
#     get_NHM_upstream_catchments(reachGDF,catchDict, upstream_out_file)
#     
# 
# =============================================================================
    

# =============================================================================
# def get_NHM_upstream_catchments(reachGDF,catchDict, out_file = 'upstreamCatchDict.npy'):
#     """
#     creates a numpy dictionarys (saved out to files) of the catchments upstream of each catchment
#     :param reachGDF: [str] geodataframe of reaches within the national hydrologic model framework
#     :param catchDict: [str] dictionary of model nodes for each catchment
#     :param out_file: [str] output file for saving the dictionary of catchments in or upstream of each catchment
#     """   
#     
#     reachDF_scratch = deepcopy(reachGDF.loc[reachGDF.seg_id_nhm.isin(catchDict.keys()),['nsegment_v1_1','tosegment_v1_1','seg_id_nhm','toseg_id_nhm','Version']])
# 
#     
#     #dictionary matching seg_id_nhm : all upstream seg_id_nhm (including the current seg_id_nhm)
#     upStreamDict = {x:[x] for x in np.unique([reachDF_scratch.seg_id_nhm,reachDF_scratch.toseg_id_nhm])}
# 
#     
#     i=0
#     while reachDF_scratch.shape[0]>0:
#         i = i+1
#         thisGroup = reachDF_scratch.loc[~(reachDF_scratch.seg_id_nhm.isin(np.unique(reachDF_scratch.toseg_id_nhm))),['seg_id_nhm','toseg_id_nhm']]
#         for row in thisGroup.itertuples():
#             if row.toseg_id_nhm!=0:
# 
#                     upStreamDict[row.toseg_id_nhm].extend(upStreamDict[row.seg_id_nhm])
# 
#         reachDF_scratch = reachDF_scratch.loc[~(reachDF_scratch.seg_id_nhm.isin(thisGroup.seg_id_nhm))]
# 
#     np.save(out_file,upStreamDict)
# =============================================================================

def check_fix_geometries(gdf):
    if not np.all(gdf['geometry'].is_valid):
        gdf.geometry[~gdf['geometry'].is_valid] = gdf.geometry[~gdf['geometry'].is_valid].buffer(0)
        
    return gdf

def actualsize(input_obj):
    memory_size = 0
    ids = set()
    objects = [input_obj]
    while objects:
        new = []
        for obj in objects:
            if id(obj) not in ids:
                ids.add(id(obj))
                memory_size += sys.getsizeof(obj)
                new.append(obj)
        objects = gc.get_referents(*new)
    return memory_size


def get_catchment_nodes(gdb = None, reach_files=None, catchment_files=None, attribute_files = None, networkCode = "NHM", reachIdx = "seg_id_nhm",model_shapefile="modelshapefile.shp", local_out_file = 'localCatchDict.npy', upstream_out_file = 'upstreamCatchDict.npy', model_crs=None, network_crs=None):
    """
    creates a numpy dictionaries (saved out to files) of 1) the model nodes within each catchment and 2) the catchments in or upstream of each catchment
    :param model_shapefile: [str] file path of the groundwater model grid shapefile
    :param model_eps: [str] epsg code for the model grid shapefile
    :param local_out_file: [str] path to the output file for saving the dictionary of gw model nodes within each catchment
    :param upstream_out_file: [str] path to the output file for saving the dictionary of catchments in or upstream of each catchment
    """    
    
    #open files
    if gdb:
        catchmentGDF = gpd.read_file(gdb, layer=catchment_files)
        reachGDF = gpd.read_file(gdb,layer=reach_files)
        if network_crs:
            catchmentGDF.set_crs(crs=network_crs,inplace=True, allow_override=True)
            reachGDF.set_crs(crs=network_crs,inplace=True, allow_override=True)
    elif reach_files and catchment_files:
        reachGDF = gpd.GeoDataFrame(pd.concat([gpd.read_file(os.path.join(i,"Hydrography","NHDFlowline.shp")) for i in reach_files],ignore_index=True), crs=gpd.read_file(reach_files[0]).crs)
        catchmentGDF = gpd.GeoDataFrame(pd.concat([gpd.read_file(os.path.join(i,"Catchment.shp")) for i in catchment_files], ignore_index=True), crs=gpd.read_file(catchment_files[0]).crs)
        if network_crs:
            catchmentGDF.set_crs(crs=network_crs,inplace=True, allow_override=True)
            reachGDF.set_crs(crs=network_crs,inplace=True, allow_override=True)
    
    modelGDF = gpd.read_file(model_shapefile)
    
    if modelGDF.crs.srs!= model_crs and model_crs is not None:
        modelGDF.set_crs(crs=model_crs, inplace=True, allow_override=True)
    
    #reproject the catchments to the model grid
    if catchmentGDF.crs!=modelGDF.crs:
        catchmentGDF.to_crs(crs=modelGDF.crs, inplace=True)
    if reachGDF.crs!=modelGDF.crs:
        reachGDF.to_crs(crs=modelGDF.crs, inplace=True)

    
    #check for invalid geometries
    reachGDF = check_fix_geometries(reachGDF)
    catchmentGDF = check_fix_geometries(catchmentGDF)
    modelGDF = check_fix_geometries(modelGDF) 
    
    
    #network specific processing
    if networkCode=="NHDPlus":
        #add a COMID column to the catchments
        catchmentGDF['COMID']=catchmentGDF['FEATUREID']
        #read in the attribute data
        attributeDF = pd.concat([gpd.read_file(os.path.join(i,"PlusFlowlineVAA.dbf")) for i in attribute_files],ignore_index=True)
        #clean up the COMID heading
        attributeDF.rename(columns={attributeDF.columns[[x.lower()=="comid" for x in attributeDF.columns]].values[0]:"COMID"}, inplace=True)
        attributeDF = attributeDF.drop(columns="geometry")
        reachGDF = pd.merge(left=reachGDF, right=attributeDF,on="COMID")
    elif networkCode=="NHM":
        #make a dictionary matching the v1 and v1.1 id's
        matchDict = {row.nsegment_v1_1:row.seg_id_nhm for row in reachGDF.itertuples()}
        matchDict[0]=0
        
        #add the seg_id_nat
        catchmentGDF['seg_id_nhm'] = [matchDict[x] if x in matchDict.keys() else np.nan for x in catchmentGDF.hru_segment_v1_1]
        reachGDF['toseg_id_nhm'] = [matchDict[x] if x in matchDict.keys() else np.nan for x in reachGDF.tosegment_v1_1]
          
    with open("log.txt","w+") as f:
        f.write("gets to here\n")   
        
    catchmentAreaDF = pd.DataFrame(catchmentGDF[reachIdx]).assign(Catch_Area = catchmentGDF.area)
    #clip the catchments and the reaches to the model grid
    catchmentGDF = gpd.clip(catchmentGDF,modelGDF, keep_geom_type=True)
    reachGDF = gpd.clip(reachGDF,modelGDF, keep_geom_type=True)

    
    #create a uniform ibound column
    modelGDF['ibound']= modelGDF[[x for x in modelGDF.columns if "ibound" in x or "idomain" in x]].max(axis=1)
    
    
    #join the catchments and the model grid centroids (the catchments are joined to the whole model grid to enable filtering by the percent included)

    joinDF = sjoin(catchmentGDF,gpd.GeoDataFrame(data=modelGDF[['node','ibound']], geometry = modelGDF.centroid),how="inner")
    
    #join the reaches and the model grid
    reachGDF = sjoin(reachGDF,modelGDF,how="inner")
    
    #calculate the inbound fraction of area to filter our catchments that are largely outside the model area
    #fracArea = joinDF[[reachIdx,'ibound']].groupby(reachIdx,as_index=False).mean()
    fracArea = joinDF[[reachIdx,'node']].merge(catchmentAreaDF).merge(pd.DataFrame(modelGDF[['node','ibound']]).assign(Model_Area = modelGDF.area))
    fracArea = fracArea.loc[fracArea.ibound==1].groupby([reachIdx,'Catch_Area'],as_index=False).sum()
    fracArea['PerActive']=fracArea['Model_Area']/fracArea['Catch_Area']
    fracArea = fracArea.merge(joinDF[[reachIdx,'ibound']].groupby(reachIdx,as_index=False).mean().rename(columns={'ibound':'mean_ibound'}))
    #this requires > 33% of the catchment to be in the active model area and > 33% of the overlapped model cells to be active 
    joinDF = joinDF.loc[joinDF[reachIdx].isin(fracArea[reachIdx].loc[(fracArea.mean_ibound>0.33)&(fracArea.PerActive>0.33)])]
    
    #filter the reaches to those in active model cells or with catchments within the active area
    activeReaches = [x for x in reachGDF.loc[(reachGDF.ibound==1)|(reachGDF[reachIdx].isin(joinDF[reachIdx])),reachIdx].drop_duplicates().values]
    
    #get the reaches downstream of active reaches. this will likely add reaches outside the primary model area (ie reaches in adjacent basins with boundary mismatches in the heeadwaters), but will be important for getting reaches in large rivers / bays that aren't included in the modflow model but are in the stream temp work
    i = 0
    while i < 10:
        i = i + 1
        if networkCode=="NHM":
            newReaches = reachGDF.loc[(~reachGDF[reachIdx].isin(activeReaches)) & (reachGDF[reachIdx].isin(reachGDF.toseg_id_nhm[reachGDF[reachIdx].isin(activeReaches)])),reachIdx].drop_duplicates()
        if networkCode=="NHDPlus":
            newReaches = reachGDF.loc[(~reachGDF[reachIdx].isin(activeReaches)) & (reachGDF['Hydroseq'].isin(reachGDF.DnHydroseq[reachGDF[reachIdx].isin(activeReaches)])),reachIdx].drop_duplicates()
        if len(newReaches)==0:
            break
        activeReaches.extend(newReaches.values)
    reachGDF = reachGDF.loc[reachGDF[reachIdx].isin(activeReaches)]
    
    #add the segments that are missing hrus, this only adds the model cells that overlap those segments, not the full watersheds. It removes those cells from their other segments to prevent double counting
    reachesToAdd = reachGDF.loc[~reachGDF[reachIdx].isin(joinDF[reachIdx]),[reachIdx,"node"]]
    joinDF = pd.concat([joinDF.loc[~(joinDF.node.isin(reachesToAdd.node))],reachesToAdd])
                        
    #dictionary matching COMID : node
    catchDict = {x:joinDF.loc[(joinDF[reachIdx]==x),"node"].values for x in np.unique(joinDF[reachIdx])}
    np.save(local_out_file,catchDict)    
    
    print("and to here")

    with open("log.txt","a") as f:
        f.write("and to here\n")

    del joinDF, catchmentGDF, modelGDF
    gc.collect(generation=2)


    if networkCode=="NHM":       
        get_NHM_upstream_catchments(reachGDF.loc[reachGDF[reachIdx].isin(catchDict.keys()),['seg_id_nhm','toseg_id_nhm']].drop_duplicates(), upstream_out_file)
    elif networkCode=="NHDPlus":
        get_NHDPlus_upstream_catchments(reachGDF.loc[reachGDF[reachIdx].isin(catchDict.keys()),['COMID','Hydroseq','DnHydroseq']].drop_duplicates(), upstream_out_file)


def get_NHM_upstream_catchments(reachDF, out_file = 'upstreamCatchDict.npy'):
    """
    creates a numpy dictionarys (saved out to files) of the catchments upstream of each catchment
    :param reachGDF: [str] geodataframe of reaches within the national hydrologic model framework
    :param catchDict: [str] dictionary of model nodes for each catchment
    :param out_file: [str] output file for saving the dictionary of catchments in or upstream of each catchment
    """   


    
    #dictionary matching seg_id_nhm : all upstream seg_id_nhm (including the current seg_id_nhm)
    upStreamDict = {x:[x] for x in np.unique([reachDF.seg_id_nhm,reachDF.toseg_id_nhm])}

    print("and also to here")
    i=0
    while reachDF.shape[0]>0:
        i = i+1
        if i%100==0:
             print(i)
        thisGroup = reachDF.loc[~(reachDF.seg_id_nhm.isin(np.unique(reachDF.toseg_id_nhm))),['seg_id_nhm','toseg_id_nhm']]
        for row in thisGroup.itertuples():
            if row.toseg_id_nhm!=0:
                try:
                    upStreamDict[row.toseg_id_nhm].extend(upStreamDict[row.seg_id_nhm])
                except:
                    pass

        reachDF = reachDF.loc[~(reachDF.seg_id_nhm.isin(thisGroup.seg_id_nhm))]

    np.save(out_file,upStreamDict)

def get_NHDPlus_upstream_catchments(reachDF,out_file = 'upstreamCatchDict.npy'):
    """
    creates a numpy dictionarys (saved out to files) of the catchments upstream of each catchment
    :param reachGDF: [str] geodataframe of reaches within the national hydrologic model framework
    :param catchDict: [str] dictionary of model nodes for each catchment
    :param out_file: [str] output file for saving the dictionary of catchments in or upstream of each catchment
    """   
    
    
    gc.collect(generation=2)

    with open("log.txt","a") as f:
        f.write("and also to here\n")
    
    #dictionary matching Hydroseq : all upstream COMID (including the current COMID)
    upStreamDict = {x[0]:[x[1]] for x in reachDF[['Hydroseq','COMID']].values}


    i=0
    while reachDF.shape[0]>0:
        i = i+1

        #this works from the most upstream reaches (those that aren't downstream of other reaches)
        thisGroup = reachDF.loc[~(reachDF.Hydroseq.isin(np.unique(reachDF.DnHydroseq)))]
        for row in thisGroup.itertuples():
            if row.DnHydroseq!=0:
                try:
                    upStreamDict[row.DnHydroseq].extend(upStreamDict[row.Hydroseq])
                except:
                    pass

        reachDF = reachDF.loc[~(reachDF.COMID.isin(thisGroup.COMID))]
    #switch the upstream dictionary from Hydroseq keys to COMID keys
    keyList = [x for x in upStreamDict.keys()]
    for oldKey in keyList:
        newKey = upStreamDict[oldKey][0]
        upStreamDict[newKey]=upStreamDict.pop(oldKey)
    
    np.save(out_file,upStreamDict)


def compile_catchment_discharge(node_discharge_file="resultsAgg.feather", reachIdx = "seg_id_nat", catchDictFile = 'localCatchDict.npy', upStreamDictFile = 'upstreamCatchDict.npy', out_file = "CatchmentDischarge.feather"):
    """
    compiles the groundwater discharge (total and as a percent of upstream - including the local catchment - discharge) for each catchment
    :param node_discharge_file: [str] feather file of groundwater discharge for each model node
    :param catchDictFile: [str] file path to dictionary of model nodes for each catchment
    :param upStreamDictFile: [str] file path to the dictionary of catchments in or upstream of each catchment
    :param out_file: [str] output feather file of the compiled discharge
    """       
    node_discharge = pd.read_feather(node_discharge_file)

    catchDict = np.load(catchDictFile, allow_pickle=True)[()]
    upStreamDict = np.load(upStreamDictFile, allow_pickle=True)[()]
    
    localDis = [(x,np.sum(node_discharge.loc[node_discharge.node.isin(catchDict[x]),"q"]),np.mean(node_discharge.loc[node_discharge.node.isin(catchDict[x]),"q_std"]),np.mean(node_discharge.loc[node_discharge.node.isin(catchDict[x]),"q_std_per"])) for x in catchDict.keys()if x !=0]
    dischargeDF = pd.DataFrame(localDis,columns=[reachIdx,'q_local','q_std','q_std_per'])
    dischargeDF['q_all']=[np.sum(dischargeDF.q_local.loc[dischargeDF[reachIdx].isin(upStreamDict[x])]) for x in dischargeDF[reachIdx]]
    dischargeDF['Per_Local']=dischargeDF['q_local']/dischargeDF['q_all']*100
    dischargeDF['nDown'] = [np.sum([x in upStreamDict[y] for y in upStreamDict.keys()]) for x in dischargeDF[reachIdx]] #this counts the number of reaches within the upStreamDict for which the current reach is upstream. it is used to filter out streams near the boundaries of models (where the downstream network is more fully represented in a different model)
    
    dischargeDF.to_feather(out_file)
    

def aggregate_catchment_discharge (dischargeFiles, out_file, spatial_idx_name):
    """
    combines compiled discharge from multiple models into 1 dataframe
    :param out_file: [str] list of discharge files to aggregate
    :param out_file: [str] output feather file of the compiled discharge
    """
    
    dischargeDF = pd.read_feather(dischargeFiles[0])
    for i in range(1,len(dischargeFiles)):
        dischargeDF = pd.concat([dischargeDF,pd.read_feather(dischargeFiles[i])], ignore_index=True)
        
    #keep the rows with the greatest # of downstream segments for each segment (filters out segments that are included in model edges but their downstream segments are not)
    idx = dischargeDF.groupby([spatial_idx_name])['nDown'].transform('max') == dischargeDF['nDown']
    dischargeDF = dischargeDF[idx]
    
    #aggregates the discharge by segment
    dischargeDF_agg = dischargeDF.groupby(spatial_idx_name,as_index=False).mean()
    dischargeDF_agg.to_feather(out_file)
    
    #get the number of discharge values, mean and std values for each segment (some segments are represented in multiple models)
    dischargeDF_agg_stats = dischargeDF[[spatial_idx_name,"q_local"]].groupby(spatial_idx_name,as_index=False).agg(['mean','std','count']).droplevel(0,axis=1).reset_index()
    dischargeDF_agg_stats.to_feather(out_file.replace(".feather","_stats.feather"))
