from download import download
import os.path as op

# NOTE: IF YOU DON'T HAVE PYTHON OR THE DOWNLOAD MODULE INSTALLED
# SIMPLY OPEN THE URL BELOW AND UNZIP TO A FOLDER CALLED `data` IN
# THE ROOT OF THE CLASS REPOSITORY.

# Url to data for the class
url = 'https://www.dropbox.com/s/m9riuepg5xfdyoi/data.zip?dl=0'

# Download and unzip the data for the course
path_download = op.dirname(op.abspath(__file__))
download(url, path_download, kind='zip', replace=True)
