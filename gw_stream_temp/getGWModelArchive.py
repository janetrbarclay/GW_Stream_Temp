import urllib.request

def get_gw_archive_data(archiveCode, basefile, destination):
    #make the url
    url = "https://water.usgs.gov/GIS/dsdl/gwmodels/{}/{}".format(archiveCode,basefile)
    urllib.request.urlretrieve(url, destination)
    