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



def check_fix_geometries(gdf):
    '''
    checks and repairs (by adding a buffer) geopandas geometries
    '''
    if not np.all(gdf['geometry'].is_valid):
        gdf.geometry[~gdf['geometry'].is_valid] = gdf.geometry[~gdf['geometry'].is_valid].buffer(0)
        
    return gdf


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

        
    catchmentAreaDF = pd.DataFrame(catchmentGDF[reachIdx]).assign(Catch_Area = catchmentGDF.area)
    #clip the catchments and the reaches to the model grid
    catchmentGDF = gpd.clip(catchmentGDF,modelGDF, keep_geom_type=True)
    reachGDF = gpd.clip(reachGDF,modelGDF, keep_geom_type=True)

    
    #create a uniform ibound column
    modelGDF['ibound']= modelGDF[[x for x in modelGDF.columns if "ibound" in x or "idomain" in x]].max(axis=1)
    
    
    #join the catchments and the model grid centroids (the catchments are joined to the whole model grid to enable filtering by the percent included)
    
    ## modelGDF.node is 1 based
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
    
    
def compile_catchment_discharge_csv(node_discharge_file = "drn_obs.csv", reachIdx = "seg_id_nat", upStreamDictFile = 'upstreamCatchDict.npy', out_file = "CatchmentDischarge.feather"):
    """
    compiles the groundwater discharge (total and as a percent of upstream - including the local catchment - discharge) for each catchment
    :param node_discharge_file: [str] csv file of groundwater discharge for each model node for each model timestep (reachs as row names, timesteps in column 1)
    :param upStreamDictFile: [str] file path to the dictionary of catchments in or upstream of each catchment
    :param out_file: [str] output feather file of the compiled discharge
    """  
    
    rawDischargeDF = pd.read_csv(node_discharge_file)
    
    upStreamDict = np.load(upStreamDictFile, allow_pickle=True)[()]
    
    dischargeDF = pd.DataFrame(data = {reachIdx: rawDischargeDF.columns[1:],'q_local':rawDischargeDF.iloc[:,1:].mean(),'q_std':rawDischargeDF.iloc[:,1:].std()})
    dischargeDF['q_std_per']=dischargeDF['q_std']/dischargeDF['q_local']
    dischargeDF['seg_id_nat']=[int(x) for x in dischargeDF['seg_id_nat']]
    dischargeDF['q_all']=[np.sum(dischargeDF.q_local.loc[dischargeDF[reachIdx].isin(upStreamDict[int(x)])]) for x in dischargeDF[reachIdx]]
    dischargeDF['Per_Local']=dischargeDF['q_local']/dischargeDF['q_all']*100
    dischargeDF['nDown'] = [np.sum([x in upStreamDict[y] for y in upStreamDict.keys()]) for x in dischargeDF[reachIdx]] #this counts the number of reaches within the upStreamDict for which the current reach is upstream. it is used to filter out streams near the boundaries of models (where the downstream network is more fully represented in a different model)

    dischargeDF.reset_index(inplace=True,drop=True) 
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
    
    
def aggregate_catchment_params (paramFiles, out_file, spatial_idx_name):
    """
    combines compiled parameter values from multiple models into 1 dataframe
    :param out_file: [str] list of parameter files to aggregate
    :param out_file: [str] output feather file of the compiled parameters
    """
    
    paramDF = pd.read_feather(paramFiles[0])
    for i in range(1,len(paramFiles)):
        paramDF = pd.concat([paramDF,pd.read_feather(paramFiles[i])], ignore_index=True)
        
    
    #aggregates the param values by segment
    paramDF_agg = paramDF.groupby(spatial_idx_name,as_index=False).mean()
    paramDF_agg.to_feather(out_file)
    
    #get the number of discharge values, mean and std values for each segment (some segments are represented in multiple models)
    paramDF_agg_stats = paramDF[[spatial_idx_name,"hk_0"]].groupby(spatial_idx_name,as_index=False).agg(['mean','std','count']).droplevel(0,axis=1).reset_index()
    paramDF_agg_stats.to_feather(out_file.replace(".feather","_stats.feather"))
    
    

