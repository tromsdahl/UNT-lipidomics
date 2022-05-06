# UNT-lipidomics
This repo is a collection of R and Python scripts used for lipidomics analysis for the lipidomics methods developed at the University of North Texas.  

The lipidomics methods are part of an upcoming manuscript to be published. If you use any of these scripts or the lipidomics methods, please cite this article:  

> TBD

If you have any questions on the use of these scripts or on the lipidomics methods themselves, feel free to reach out through GitHub, here, or by email: [romsdahl.trevor@gmail.com](romsdahl.trevor@gmail.com).

## Description of Contents
### 1. R scripts
  - This folder contains two R scripts:  
    - neutral_lipids_data_process_seed.R  
    - polar_lipids_data_process_seed.R  
  - Each script is used with data collected from the neutral lipidomics method (C30) or the polar lipidomics method (NH~2~), respectively. 
  - Note: each script says it is for seeds, but they can be used with any tissue desired. For example, to use with leaf tissue, simply use the find and replace function within RStudio to remove 'seed' and replace with 'leaf'.  
  - Instructions on what other changes need to be made in the scripts is described by the file "steps_for_using_R_scripts.docx" found under "03_resources".  

### 2. Python scripts
  - This folder contains 8 Python scripts:  
    - dag_fa_mrm_generator.py
    - lipid_calc_app.py
    - pc_fa_mrm_generator.py
    - pe_fa_mrm_generator.py
    - pg_fa_mrm_generator.py
    - pi_fa_mrm_generator.py
    - ps_fa_mrm_generator.py
    - tag_fa_mrm_generator.py
  - The "mrm generator" scripts are useful for generating MRM tables (Q1 and Q3 masses) to easily load into Analyst or whichever mass spectrometry vendor software is being used. To use these scripts you will need to have Python3 installed, along with the Numpy and Pandas libraries. These are easily installed by downloading [Anaconda](https://www.anaconda.com/). The only changes required to use these scripts are to change what fatty acids are relevant to the biological samples being analyzed on line 5. After that you can run the script and it should output .csv files with the MRM transitions for that lipid class.
  - An additional script located here is the "lipid_calc_app.py". This is useful to quickly look up the mass of a lipid molecular species along with its most likely or potential adducts. This requires to enter in a lipid molecular species in the correct format (e.g. PC-34:2 -- the letters for the lipid class, followed by a '-', number of carbons, then ':', and lastly the number of unsaturations).  
  
### 3. Resources
  - This folder contains three files:  
    - 2021_lipidomics_workshop.pdf
    - lipidomics_workshop_instructions.pdf
    - steps_for_using_R_scripts.docx
  - The first two files are a presentation and set of instructions for a workshop given in the fall of 2021 on how to use the R scripts with the lipidomics methods developed here.
  - The last file are a set of step-by-step instructions taken to prepare the data table loaded into RStudio and used with the R scripts to process the collected lipidomics data.
