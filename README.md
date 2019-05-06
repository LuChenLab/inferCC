# lung_cancer_10x
Lung cancer 10X analysis

Scripts used in this project

1. `0_python`:
    - Python scripts to call CellRanger, velocyto and scripts to check the sequencing quality
 
2. `0_01_seurat`:
    - R scripts to run Seurat with default parameters, used for conditional exploration

3. `1_R_markdown`:
    - Rmarkdown files, the real R code for the whole data process
        - batch effect correction
        - cluster identification
        - cell identification
        - etc
    
4. `2_R_script_each_cell`:
    - R scripts to analysis different cells

5. `3_python_construct_loom`:
    - extract data from loom file, and integrate into one
    
99. `99_files`:
    - additional files contains meta info of all the cells etc.