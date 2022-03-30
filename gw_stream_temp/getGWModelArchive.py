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
    subDir = "model" if "model" in zippedFile else "output" if "output" in zippedFile else "ancillary"

    outPath = os.path.join(archiveDir,subDir)

    with zipfile.ZipFile(zippedFile) as z:
        z.extractall(outPath)

    with open(output, "w+") as f:
        f.write(outPath)
    

