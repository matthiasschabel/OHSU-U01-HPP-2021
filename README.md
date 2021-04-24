# OHSU-U01-HPP-2021
MATLAB source code and pre-processed data for analysis of in vivo placental T2* mapping in pregnant human study participants

Research was supported by the NICHD Human Placenta Project U01 HD087182 (Frias)

This code was developed and tested on both Mac OS 10.15.7 and Ubuntu Linux using MATLAB 2019a with the Statistics and Machine Learning Toolbox

## Installation

Add all .m files to the MATLAB path and save.

## Instructions for use

(1) Download the post-processed study data file "U01_postprocessed_data.mat" (https://www.dropbox.com/s/n5g3q6lo1r1nh88/U01_postprocessed_data.mat?dl=0) 

(2) Load study data into the MATLAB workspace (```load U01_postprocessed_data```). The data file is approximately 1.3GB in size, so this may take a couple of minutes to complete. 

(3) Open "analyzeHumanT2starData.m" in the MATLAB code editor

(4) Set ```rootDirectory``` to an appropriate path on your system

(5) Comment out the line ```exportFiles = false;``` if you want EPS figures and XLSX output to be saved

(6) Run ```analyzeHumanT2starData``` in the MATLAB command window. This will take a few minutes and open a large number of figure windows, reproducing all figures in the manuscript as well as a number of figures containing supporting information. 

