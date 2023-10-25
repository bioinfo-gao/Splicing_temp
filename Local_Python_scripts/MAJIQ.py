# create temporary working directory for this notebook
from tempfile import TemporaryDirectory
tempdir = TemporaryDirectory()
%cd $tempdir.name
tempdir = "/Users/zgao1/Documents/Biogen_Project/Python_scripts"
%cd $tempdir.name
# download archive with example data
!curl -LO http://majiq.biociphers.org/data/majiq_het_vignette.zip
# unzip and show what files are in the archive
!unzip -qn majiq_het_vignette.zip

# reticulate::repl_python()
# Python 3.8.12 (/Users/zgao1/Library/r-miniconda/envs/r-reticulate/bin/python)
# Reticulate 1.24 REPL -- A Python interpreter in R.
# Enter 'exit' or 'quit' to exit the REPL and return to R.
# >>> from tempfile import TemporaryDirectory
# >>> tempdir = TemporaryDirectory()
# >>> 
ls()
