# https://support.rstudio.com/hc/en-us/articles/200486138-Changing-R-versions-for-the-RStudio-Desktop-IDE#:~:text=Run%20the%20installer%20from%20CRAN,alias%20directly%20using%20ln%20%2Ds
# change the R version in mac

ls -l /Library/Frameworks/R.framework/Versions/
# total 0
# drwxrwxr-x  6 root  admin  192 Jan 11 12:12 3.3
# lrwxr-xr-x  1 root  admin    3 Jan 11 12:12 Current -> 3.3
# 
# 
# If you want to override the version of R selected by RStudio's default behavior then you can set the RSTUDIO_WHICH_R environment variable to the R executable that you want to run against. For example, to force RStudio to use the R executable located at /usr/local/bin:


conda create --name R41 
conda activate R41
conda install -c conda-forge r-base # R4.1.2

which R
# /opt/anaconda3/envs/R41/bin/R

export RSTUDIO_WHICH_R=/opt/anaconda3/envs/R41/bin/R    #export RSTUDIO_WHICH_R=/usr/local/bin/R