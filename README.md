This repo contains the code and associated files for "", currently found on biorxiv (https://www.biorxiv.org/).

## Repo organization and simulation pipeline
The root directory contains two folders: "experimentX (e.g. experiment1)" and "publicSimInput". "publicSimInput" contains input information that will be shared among the simulations for all experiments. Different experiments may differ in regard to the parameter settings for the source program, types of output metrics generated from the source program, types of analyses and so on. The folders are organized to allow for maximized user customization.

There are three main phases to the experiments and the folders are structured based on these phases:
1. Simulations - run simulations and generate raw output data.
2. H5 data processing - process the raw output data into Hierarchical Data Format (HDF5).  
3. Analyses - from HDF5 file, generate calculated metrics (see Supplementary Methods for a detailed list) for the downstream data analyses and visualizations.

In order to recreate the analysis, you simply need to clone the repository, and create a local directory structure that should match the following. Please note, some directories (e.g. "experiment1/analyses/calculatedMetrics") will be missing from the repository and therefore need to be created manually:
```

├── experiment1
│   ├── publicVariables
│   ├── simulations
|   |   ├── ac_amarelCode
|   |   ├── al_amarelLog
|   |   ├── sourceCode
|   |   ├── createSimInputCode
|   |   ├── simInput
|   |   └── simRawOutput
|   |   └── trialmRNAcountData
|   ├── dataProcessing
|   |   ├── createH5Code
|   |   └── simOutputH5data
│   └── analyses
|   |   ├── analysesCode
|   |   ├── calculatedMetrics
|   |   ├── figures
|   |   └── externalData
└── publicSimInput

```

Then, download the simulation output h5 data (https://www.) and place it in /experiment1/dataProcessing/simOutputH5Data/. Upon downloading the data, you should change the filenames to the following new file name: "experiment1_output.h5". 

Next, download the trial mRNA count data (https://www.) and place it in /experiment1/simulations/trialmRNAcountData/. Upon downloading the data, you should change the filenames to the following new file name: "trial_mRNAcount.feather". This data was generated from previous trial simulations, and is used to calculate scaled mRNA synthesis rates for experiment1 (see Methods), as well as to compare the scaled VS non-scaled mRNA synthesis rates in supplementary figures. 

After that, it should just be a matter of running the code in the specified order. **To ensure smooth running of the code, start with a clean R environment for each Rmd**. Include step 1 and 2 if the simulations are to be conducted with different parameters. To analyze current data sets, skip step 1 and 2. 

#### Note: no need to run serperately, but the following code will be sourced before each of other RMD files are run. 

`/experiment1/publicVariables/createPublicVariables.Rmd`


#### Step 1: deploy simulations on Rutgers Amarel cluster and collect raw output

1. `experiment1/simulations/createSimInputCode/experiment1_createSimInput.Rmd`  
2. `experiment1/simulations/ac_amarelCode/experiment1_createJobs.Rmd` 


#### Step 2: process data from raw output to a HDF5 file

1. `experiment1/dataProcessing/createH5Code/experiment1_createH5.Rmd`  


#### Step 3: analyze and visualize the data

1. `experiment1/analyses/analysesCode/experiment1_calculateMetrics.Rmd`  
2. `experiment1/analyses/analysesCode/experiment1_makeFig1.Rmd` 
3. `experiment1/analyses/analysesCode/experiment1_makeFig2.Rmd`
4. `experiment1/analyses/analysesCode/experiment1_makeFig3.Rmd` 
5. `experiment1/analyses/analysesCode/experiment1_makeFig4.Rmd` 
6. `experiment1/analyses/analysesCode/experiment1_makeFig5.Rmd` 
7. `experiment1/analyses/analysesCode/experiment1_makeSuppFigs.Rmd` 




