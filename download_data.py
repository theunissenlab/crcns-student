from download import download
import os.path as op

# NOTE: IF YOU DON'T HAVE PYTHON OR THE DOWNLOAD MODULE INSTALLED
# SIMPLY OPEN THE URL BELOW AND UNZIP TO A FOLDER CALLED `data` IN
# THE ROOT OF THE CLASS REPOSITORY.

# Url to data for the class
url = 'https://www.dropbox.com/s/hkrfyj2poizx9br/data.zip?dl=1'
# XXX another URL for course presentation material?

# Download and unzip the data for the course
path_download = op.join(op.dirname(op.abspath(__file__)), 'data')
download(url, path_download, zipfile=True)
