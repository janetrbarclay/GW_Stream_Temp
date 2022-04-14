import requests, zipfile, os

def get_gw_archive_data(archiveCode, basefile, destination):
    #make the url
    url = "https://water.usgs.gov/GIS/dsdl/gwmodels/{}/{}".format(archiveCode,basefile)
    req = requests.get(url, stream=True)

    with open(destination,'wb') as output_file:
        for chunk in req.iter_content(chunk_size=1025):
            if chunk:
                output_file.write(chunk)

def unzip_gw_archive_data(zippedFile,output):
    #get the destination directory
    archiveDir = os.path.dirname(zippedFile)
    subDirList = ['ancillary','bin','georef','model','output','source']
    subDir = [x for x in subDirList if os.path.basename(zippedFile).startswith(x)]
    subDir = subDir[0] if len(subDir)==1 else "ancillary"
        
    if subDir =="ancillary" and os.path.basename(zippedFile)!="ancillary.zip":
        subDir = os.path.join(subDir,os.path.basename(zippedFile).split(".")[0])
    

    outPath = os.path.join(archiveDir,subDir)
    
    #some directories shouldn't have subdirectories
    outPath = archiveDir if subDir in ['bin','source','georef'] else outPath

    with zipfile.ZipFile(zippedFile) as z:
        z.extractall(outPath)

    with open(output, "w+") as f:
        f.write(outPath)
    
def get_model_list(archiveDir, outFile):
    model_list = [x for x in os.listdir(os.path.join(archiveDir,"model")) if os.path.isdir(os.path.join(archiveDir,"model",x))]
    with open(outFile,"w") as f:
        [f.write("{}\n".format(x)) for x in model_list]
    
    return outFile
	
