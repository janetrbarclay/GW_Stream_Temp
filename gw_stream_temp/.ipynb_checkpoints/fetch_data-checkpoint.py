# -*- coding: utf-8 -*-
"""
Created on Fri Jun 25 11:38:06 2021

@author: jbarclay
"""

import sciencebasepy
import zipfile
import os



def get_NHM_gis_data(item='5e29d1a0e4b0a79317cf7f63',filenameLst=['GFv1.1.gdb.zip'], destination='data/NHM'):
    """
    fetches files from a science base item and unzips them (if appropriate)
    :param item: [str] science base item number
    :param filenameLst: [str] list of files to download
    :param destination: [str] path to save the downloaded files
    """
    sb = sciencebasepy.SbSession()
    item_json = sb.get_item(item)
    fileList = sb.get_item_file_info(item_json)
    for filename in filenameLst:
        fileDict = [x for x in fileList if x['name']==filename][0]
        sb.download_file(fileDict['url'], os.path.join(destination,filename))
        
        if filename.endswith("zip"):
            with zipfile.ZipFile(os.path.join(destination,filename)) as z:
                z.extractall(destination)
        


    







    
    
    
    
    