# -*- coding: utf-8 -*-
"""
Created on Wed Sep 15 12:05:59 2021

@author: jbarclay
"""

import sciencebasepy
import os
import zipfile

def get_sciencebase_data(item, fileName, outPath,outFile):
    sb = sciencebasepy.SbSession()
    
    item_json = sb.get_item(item)
    fileList = sb.get_item_file_info(item_json)
        
    fileDict = [x for x in fileList if x['name']==fileName][0]
    sb.download_file(fileDict['url'], os.path.join(outPath,fileName))
    
    if fileName.endswith("zip"):
        with zipfile.ZipFile(os.path.join(outPath,fileName)) as z:
            z.extractall(os.path.dirname(outFile))