# Mixed-Effect Model for Brain-Genetic Trajectories

This repository contains an R-based analysis for studying mixed-effect models applied to brain-genetic trajectories. It includes R scripts for data analysis, the required data files, and results generated from the models.

## Directory Structure

- **Rscripts/**  
  This directory contains all R scripts used for data analysis and modeling.

- **data/**  
  Contains the datasets required for analysis. The data files can be loaded by using relative paths with the `here` package.

- **results/**  
  This directory stores the output of the analysis, including model summaries, figures, and tables.

- **.Rproj**  
  The R project file (`Analisi.Rproj`) which can be opened directly in RStudio to set up the project environment.

## Setup Instructions

Follow the steps below to set up the R project and run the scripts for analysis:

### 1. Open the R project

To begin, we advise the user download the project folder, place it in the desired directory and set up the project environment by opening directly in RStudio the file 'Analisi.Rproj'. This will ensure that all relative paths to data and directories are correctly defined. The 'here' package is used for path management in this project.

### 2. Install Dependencies

After environment set up, make sure that all required libraries in this project are installed in your working station. Please refer to the 'libraries.R' file for a list of these required dependencies.

### 3. Run and explore individual scripts

R scripts for this project are organized in three sequential steps:

- **1_importAndClean**: scripts in this directory process the raw data from the 'data/rawdata' directory, to obtain cleaned datasets stored in 'data/cleandata' ready for downstream analysis.

- **2_exploratory**: the aim of this set of scripts is to provide descriptive statistics and exploratory plots of the data, either of baseline characteristics and longitudinal data.

- **3_modeling**: main cross-sectional and longitudinal analysis to assess brain-genetic trajectories.

Output for each analysis are saved in the respective directory from 'results/'.
