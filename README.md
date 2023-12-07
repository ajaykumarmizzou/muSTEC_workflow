# muSTEC: muSTEC : An Integrated Bioinformatics Analysis Suite for the Identification, Characterization, and Outbreak Investigation of Shiga Toxin-Producing E. coli (STEC) strains using Whole Genome Sequencing (WGS) data
muSTEC workflow takes E-coli accessions list and performs SNP variant calling using gatk or bcftools for single/paired end fastq reads. 

## Prerequisites:
-  Create conda environment using given conda_env.yaml file, 
        conda env create -f conda_env.yaml
-  Activate conda environment using,
        conda activate conda_env


## Usage:
Step: 1 Clone the repository on your computer/server using,
        git clone https://github.com/ajaykumarmizzou/muSTEC_workflow

Step: 2 Pass your accession list of SNP variant calling and call the setup.sh,
        ./setup.sh accession_list.txt

Step: 3 Your data will download in /data directory and results will store in /results directory.

## Rights and Permissions
MIT License

Copyright (c) 2023 Ajay Kumar

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


