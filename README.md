# OHSU-U01-HPP-2021
MATLAB source code and post-processed MRI data for analysis of in vivo placental T2* mapping in pregnant human study participants as described
in the following prepring, currently in peer-review: [Assessing placental function across gestation: a multi-institutional study of BOLD-MRI for the prediction of adverse pregnancy outcomes](https://www.researchsquare.com/article/rs-406266/v1). A presentation on this project at the 2021 NIH Human Placenta Project Meeting can be found at 1:36:00 [here](https://videocast.nih.gov/watch=42033?jwsource=cl) and the PowerPoint slides for that presentation can be downloaded [here](https://www.dropbox.com/s/0bf2u1jhlm2yalr/202105%20HPP.pptx?dl=0).

Research was supported by the NICHD Human Placenta Project U01 HD087182 (Frias)

This code was developed and tested on both Mac OS 10.15.7 and Ubuntu Linux using MATLAB 2019a with the Statistics and Machine Learning Toolbox

## Installation

Add all .m files in the ```src``` directory to the MATLAB path and save.

## Instructions for use

(1) Download the post-processed study data file [U01_postprocessed_data.mat](https://www.dropbox.com/s/01pastzz5qw7elt/U01_postprocessed_data.mat?dl=0). This data set is copyrighted and may not be used for academic research or commercial purposes without written permission. Please contact [Matthias Schabel](schabelm@ohsu.edu) with questions or requests.

(2) Load study data into the MATLAB workspace (```load U01_postprocessed_data```). The data file is approximately 1.3GB in size, so this may take a couple of minutes to complete. 

(3) Open "analyzeHumanT2starData.m" in the MATLAB code editor

(4) Set ```rootDirectory``` to an appropriate path on your system

(5) Comment out the line ```exportFiles = false;``` if you want EPS figures and XLSX output to be saved

(6) Run ```analyzeHumanT2starData``` in the MATLAB command window. This will take a few minutes and open 57 figure windows, reproducing figures 3-6 from the manuscript as well as a number of figures containing supporting information. Specifically (denoting manuscript figure number x, panel y as MFx.y and output figure z as Fz):
  
  * F7 -> MF1.A and MF1.C
  * F8 -> MF1.B
  * F31 -> MF1.D and MF1.F
  * F32 -> MF1.E
  * F51 -> MF4
  * F54 -> MF5
  * F34-F37, F39-F42, F44-F47 -> MF6

