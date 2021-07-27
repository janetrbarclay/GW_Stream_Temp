import os
import flopy as fp
import numpy as np
import pandas as pd
from flopy.utils import Util2d


from gw_stream_temp.preprocess_modflow_data import load_modflow_model



def make_modpath_model(modelpth="./", thisModelName="MONTAGUE_drb1.04_mf6_250", flow_model_name = "MONTAGUE_drb1.04_mf6_250", totalParts = 3e6):
    """
    creates and runs a modpath model based on an existing modflow model
    :param modelpth: [str] file path to model files
    :param thisModelName: [str] base name of the model files
    :param flow_model_name: [str] name of the groundwater flow model within the modflow simulation
    :param totalParts: [int] approximate number of particles for each flow period
    """
    
    mf = load_modflow_model(modelpth,thisModelName,flow_model_name)
    
        #############################################
    ## Get the temporal specs of the model
    nper = mf.nper
    ssPer = list(np.where(mf.modeltime.steady_state)[0])
    #this will cause problems if there are multiple steady state periods with no transient periods
    if nper>1:
        perLen = int(np.mean(mf.modeltime.perlen[~(mf.modeltime.steady_state)]) )
        #increase the number of particles
        
    else:
        perLen = int(mf.modeltime.perlen)
        
    
    if nper!=len(ssPer):
        perPerYear = int(365/perLen)
        partStartPer = 0
    else:
        partStartPer = 0
        perPerYear = 1
        
    totalParts = int(totalParts * perPerYear) 
    
    mp = fp.modpath.Modpath7(modelname=flow_model_name, 
                         simfile_ext='mpsim', 
                         namefile_ext='mpnam', 
                         version='modpath7',
                         flowmodel=mf, 
                         headfilename=mf.oc.head_filerecord.array[0][0],
                         budgetfilename=mf.oc.budget_filerecord.array[0][0],
                         model_ws=modelpth,
                         verbose=True)
    
    rech = mf.get_package('RCH').recharge.array[:perPerYear]
    option_flags = [2,1,1,1,2,2,1,3,1,2,1,1] #the second digit controls the tracking direction (1=forward, 2=backward)
    
    mp7_pgList = []
    partID = 0
    ibound = mf.dis.idomain
    rechTot = np.sum([np.multiply(rech[thisPer],ibound.array[0]) for thisPer in range(len(rech))])
    
    for thisPer in range(len(rech)):
        rechFract = rech[thisPer]*ibound.array[0]/rechTot
        rech_int = (rechFract*totalParts)
        rech_int[(rech_int>0)&(rech_int<1)]=1
        rech_int = rech_int.astype(int)
        rech_sqr = np.round(np.sqrt(rech_int))**2
        rech_int = Util2d(mp,(rech.shape[2],rech.shape[3]),np.int,rech_int,name="rech_int")
        sqrList = np.unique(rech_sqr)
        sqrList = sqrList[sqrList!=0]

        for thisSqr in sqrList:
            thisSqrRt = int(np.sqrt(thisSqr))
            nameField = "Part_"+str(int(thisSqr))+"_Per"+str(thisPer)
            temp = np.zeros_like(rech_sqr)
            temp [rech_sqr==thisSqr]=1
            v = [tuple(coord) for coord in np.argwhere(np.stack([temp,np.zeros(temp.shape)])==1).tolist()]
            partIDs = [x for x in range(partID, len(v)+partID)]
            partID = max(partIDs)+1 #this is the starting partID for the next set
            sd = fp.modpath.FaceDataType(drape=1, rowdivisions6=thisSqrRt, columndivisions6=thisSqrRt, verticaldivisions1=0, verticaldivisions2=0, verticaldivisions3=0, verticaldivisions4=0, horizontaldivisions1=0, horizontaldivisions2=0, horizontaldivisions3=0, horizontaldivisions4=0, rowdivisions5=0, columndivisions5=0)
            v2 = [x for x in list(np.where(temp.ravel()==1))]
            p = fp.modpath.NodeParticleData(subdivisiondata=sd, nodes=v2)
            pg1 = fp.modpath.ParticleGroupNodeTemplate(particlegroupname=nameField, particledata=p, filename=os.path.join(modelpth,"arrays","%s.txt"%nameField), releasedata=[thisPer*perLen+perLen/2])
            mp7_pgList.append(pg1)
    timePtList_base = [0,1,2,5,10,20,50]
    timePtList=[]
    for startPer in range(len(rech)):
        timePtList.extend([x +startPer * perLen +perLen/2 for x in timePtList_base])
        timePtList.extend([365+x*30 for x in range(36)])
        timePtList.extend([np.max(timePtList)+x*90 for x in range(40)])
        timePtList.extend([np.max(timePtList)+x*365 for x in range(25)])
    mpsim = fp.modpath.Modpath7Sim(mp,
                               trackingdirection="forward",
                               weaksinkoption="pass_through",
                               weaksourceoption="pass_through",
                               budgetoutputoption='summary',
                               simulationtype = 'timeseries',
                               stoptimeoption='extend',
                               zonedataoption='on',
                               timepointdata = [len(timePtList),timePtList],
                               stopzone=-1,
                               referencetime=(partStartPer,0,0),
                               particlegroups=mp7_pgList,
                               extension='mpsim')
    mpbas = fp.modpath.Modpath7Bas(mp,
                               porosity=0.3,
                               defaultiface={'DRN': 6, 'RCH': 6, 'SFR': 6, 'GHB': 6, 'RIV':6})
    #make a directory for the particle arrays if one doesn't already exist
    if not "arrays" in os.listdir(modelpth):
        os.mkdir(os.path.join(modelpth,"arrays"))
    mp.write_input()
