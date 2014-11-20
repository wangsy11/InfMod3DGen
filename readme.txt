------------------------------------------------------------------------------------
Implementation of algorithm for 'Inferential modeling of 3D chromatin structure'
Commandline Readme
November 19, 2014
------------------------------------------------------------------------------------
   
INSTALLATION
Put all the MATLAB script files in your MATLAB path. 

USAGE
Directly run the m-file ** ChrMod_main.m ** in the directory. 
Example: do the following at the MATLAB command window:  

    >> ChrMod_main(5,100)

Here we have two parameters: the first one decides which chromosome you want to 
model and the second one sets the volume of the ensemble. In the above example, 
we are generating a 100-conformation ensemble for chromosome #5. The concept of 
ensemble is described in our paper.

INPUT
We used **interactions_HindIII_fdr0.0001_intra.txt** file as input. These data 
was downloaded from Duan et al. 
(http://noble.gs.washington.edu/proj/yeast-architecture/)


OUTPUT
The output is automatically saved in mat-file format, named by two input parame-
ters linked together. All variables of the modeling process is stored at the same 
time. The variable 'Ensemble' is the spatial structure information, with its fir-
st dimension corresponds to each single conformation in ensemble, and the next d-
imension presents the 3D coordinates of the end points of segments in the chromo-
some.


NOTES
This software was developed and tested on MATLAB R2012b and Windows7 operating sy-
stems. We have also tested it on linux servers as well.

PROBLEMS
If you encounter any problem, please do not hesitate to contact us.

CITATION
Siyu Wang, Jinbo Xu and Jianyang Zeng*. Inferential modeling of 3D chromatin stru-
cture.

------------------------------------------------------------------------------------
Comments and bug-reports are higly appreciated.

Siyu Wang, Tsinghua University
wangsy11@mails.tsinghua.edu.cn
------------------------------------------------------------------------------------