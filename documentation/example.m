clear

t=clock;

%FindPath path
FindPathPath='/myPath/FindPath/functions';
addpath(FindPathPath);

%ModelPath
modelPath='/myPath/FindPath/data';
addpath(modelPath);

%csv to model path
cvs2cbm='/myPath/FindPath/MakeModel';
addpath(cvs2cbm);

%name of the host CBM
modelName='iMM904_flux.xml';

% name of the SAR database 
dbCSVName='dbtest.csv';

%path to EFMTools
EFMTools='/myPath/efmtool';

% Load models
oriPath=pwd;
cd(modelPath)
 model=readCbModel(modelName);
dbModel=makeCBModelFromCSV(dbCSVName,';','dbModel');
cd(oriPath);

%Define options
options=struct();
options.metIn='xyl-D[e]';
options.metOut={'xu5p-D[c]','glx[c]','akg[c]'};
options.maxLength=10;
options.maxPwy=15;
options.resPath='/myPath/FindPath/results';
options.namePwy='test';
options.weightEfficiency=1;
options.ignoreRxns={'EX_atp(e)','EX_adp(e)','EX_xu5p-D(e)','EX_h(e)','EX_h(e)','EX_akg(e)','EX_pyr(e)','EX_oea(e)','EX_rea(e)','EX_glx(e)','EX_nad(e)','EX_nadp(e)','EX_nadph(e)','EX_nadh(e)'};
%%
%run
changeCobraSolver('glpk');
pause(2)
[pwys,score]=FindPath(model,dbModel,EFMTools,options);

%Return elpsed time from t (i.e. sart of the code)
Time=etime(clock, t)
