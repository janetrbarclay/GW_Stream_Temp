# -*- coding: utf-8 -*-
"""
Created on Fri Jun 25 11:38:06 2021

@author: jbarclay
"""

import sciencebasepy
import zipfile
import os
import requests
import py7zr
import time
from urllib.request import urlopen, Request
from urllib.error import HTTPError, URLError



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
                
def make_vpu_dict():
    '''
    creates a dictionary of vpu's, etc for downloading the NHDPlus network. 
    this is based on the script found here: https://github.com/Msawtelle/PyNHDPlus/blob/master/NHDPlus_Extractor/NHDPlus_Extractor.py
    '''
    vpu_dict = {'01': {'vpu':['01'],'DA':'NE'},
                     '02': {'vpu':['02'],'DA':'MA'},
                     '03': {'vpu':['03N','03S','03W'],'DA':'SA'},
                     '04': {'vpu':['04'],'DA':'GL'},
                     '05': {'vpu':['05'],'DA':'MS'},
                     '06': {'vpu':['06'],'DA':'MS'},
                     '07': {'vpu':['07'],'DA':'MS'},
                     '08': {'vpu':['08'],'DA':'MS'},
                     '09': {'vpu':['09'],'DA':'SR'},
                     '10': {'vpu':['10U','10L'],'DA':'MS'},
                     '11': {'vpu':['11'],'DA':'MS'},
                     '12': {'vpu':['12'],'DA':'TX'},
                     '13': {'vpu':['13'],'DA':'RG'},
                     '14': {'vpu':['14'],'DA':'CO'},
                     '15': {'vpu':['15'],'DA':'CO'},
                     '16': {'vpu':['16'],'DA':'GB'},
                     '17': {'vpu':['17'],'DA':'PN'},
                     '18': {'vpu':['18'],'DA':'CA'},
                     '20': {'vpu':['20'],'DA':'HI'},
                     '21': {'vpu':['21'],'DA':'CI'},
                     '22': {'vpu':['22A','22G','22M'],'DA':'PI'}
                     }
                     
    return vpu_dict

    
def get_NHDPlus_gis_data(vpu,DA,file = 'NHDSnapshot', destination="data_NHDPlus"):
    '''
    this downloads NHDPlus data
    #this is based on the script found here: https://github.com/Msawtelle/PyNHDPlus/blob/master/NHDPlus_Extractor/NHDPlus_Extractor.py
    :param vpu: [str] vector processing unit for the NHDPlus
    :param DA: [str] drainage area abbreviation for the NHDPlus
    :param file: [str] NHDPlus file to download
    :param destination: [str] location to save the downloaded file
    '''
    if not os.path.isdir(destination):
        os.mkdir(destination)
    
    #make the url
    base_url = 'https://edap-ow-data-commons.s3.amazonaws.com/NHDPlusV21/Data/NHDPlus'
    alternate_url_list = ['03N', '03S', '03W','05', '06', '07', '08',
                                '10U', '14', '15', '10L', '11','22AS', '22GU', '22MP']
    if vpu in alternate_url_list:   
        url = '{0}{1}/NHDPlus{2}/NHDPlusV21_{1}_{2}_{3}'.format(base_url,DA,vpu,file)
    else:
        url = '{0}{1}/NHDPlusV21_{1}_{2}_{3}'.format(base_url,DA,vpu,file)
        
    #get the version number
    i = 1
    while i < 25:

        final_url = '{}_{:02d}.7z'.format(url,i)

        try:
            time.sleep(0.3)
            req = Request(final_url, data=None,headers={'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/58.0.3029.81 Safari/537.36'})
            response = urlopen(req)
            print('FOUND '+ final_url)
            break
        except HTTPError or URLError:
            pass
        i += 1
    output_path = os.path.join(destination,"NHDPlusV21_{}_{}_{}_{:02d}.7z".format(vpu,DA,file,i))
    req = requests.get(final_url, stream=True)

    with open(output_path,'wb') as output_file:
        for chunk in req.iter_content(chunk_size=1025):
            if chunk:   
                output_file.write(chunk)
                
    #write the path for the zip file to a txt file for later processing
    with open(os.path.join(destination,"NHDPlusV21_{}_{}_{}.pth".format(vpu,DA,file)),"w+") as f:
        f.write(output_path)

def unzip_NHDPlus_gis_data(filePath):
    '''
    unzips downloaded NHDPlus data
    '''
    with open(filePath) as f:
        archivePath = f.readline()
        
    with py7zr.SevenZipFile(archivePath, mode='r') as z:
        z.extractall(os.path.dirname(archivePath))

    
    
    
    
    
