![image](https://github.com/user-attachments/assets/376a1791-5f66-4b20-a81e-65c2a8e24f80)# DEMINING

DEMINING: a stepwise computational framework to directly <ins>d</ins>etect <ins>e</ins>xpressed DNA and RNA <ins>m</ins>utations <ins>i</ins>n R<ins>N</ins>A deep sequenc<ins>ing</ins> data (DEMINING).

Version: 1.0 2023/12/15

Author: Zhi-Can Fu (fuzhican@picb.ac.cn) 


## 1. Installation
### 1.1 Download DEMINING git repo and decompress DEMINING.zip with password
        git clone git@github.com:fuzhican/DEMINING.git
        cd DEMINING; unzip DEMINING.zip
### 1.2 Create conda environment and activate 
        conda env create -f environment.yml 
        conda activate ~/DEMINING_env


## 2. Usage
        Usage: DEMINING -f Function -1 Path_of_fastq1 -2 Path_of_fastq2 -o Output_path -n Output_name -c Path_of_config_file -t Maximum_threads -g Genome_build_version [ -s Site_list_file -b BAM_file ]

        Arguments:
                [-f Function, "Read_mapping", "Variant_classification" or "All_steps"(default "All_steps")]
                [-1 Path of fastq1]
                [-2 Path of fastq2]
                [-o Output directory (default current directory)]
                [-n Output name]
                [-c Path of config file]
                [-t Maximum_threads]
                [-g Genome build version, "hg38" or "mm10"]
                [-s Site list file (No column name, only required when Function == "Variant_classification")]
                [-b BAM file (Only required when Function == "Variant_classification")]
                [-h show this help message and exit]
                [-v version for DEMINING]
      
## 3. Example
* Distinguish DNA and RNA mutations directly from RNA-seq fastq files

        ./DEMINING -f "All_steps" -1 Test_data/HG00145_chr22_R1.fastq.gz -2 Test_data/HG00145_chr22_R2.fastq.gz -o DEMINING_test -c Test_data/DEMINING_test.conf -n HG00145_chr22 -t 2 -g hg38
        
 
   
## 4. Output

    DEMINING_test/
    ├── HG00145_chr22_DEMINING.tsv      # Predict result
    └── HG00145_chr22_DEMINING_tmp      # Internal files of DEMINING 

        
- HG00145_chr22_DEMINING.tsv

    | **Field**      | **Discription**      | 
    | ---------- | :-----------:  |
    | Site     | Genomic coordinates |
    | DeepDDR_rawScore | Probability of RNA mutation from DEMINING|
    | DeepDDR_DM(DNA_mutation)_or_RM(RNA_mutation)| Predict label of site (DM: DNA mutation; RM: RNA mutation) |

    
    

## 5. Citation
Fu, Z.C.#, Gao, B.Q.#, Nan, F., Ma, X.K., and Yang, L.* (2024). DEMINING: a deep learning model embedded framework to distinguish RNA editing from DNA mutations in RNA sequencing data.


## 6. License
Copyright ©2023 Fudan University. All Rights Reserved.

Licensed GPLv3 for open source use or contact YangLab (liyang_fudan@fudan.edu.cn) for commercial use.

Permission to use, copy, modify, and distribute this software and its documentation for educational, research, and not-for-profit purposes, without fee and without a signed licensing agreement, is hereby granted, provided that the above copyright notice in all copies, modifications, and distributions.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
