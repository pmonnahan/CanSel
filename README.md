# CanSel
Automating calculation of numerous population-genetic metrics within and across populations.  Input is a set of VCF files (separated by chromosome) and a population key (e.g. _accessory/igsry_samples.txt_).  Data is compressed/converted into 'zarr' format (https://zarr.readthedocs.io), which allows for efficient access of very large datasets (e.g. 1000 Genomes).  Sliding-window calculation of metrics is accomplished via scikit-allel.  The entire workflow is automated and configured for computation on an HPC system via Snakemake.

## Requirements

Assumes indels have been removed from VCF

#### Anaconda

If on cluster do,
    
    
    module load conda
    
 otherwise if you have Homebrew,
 
    brew install conda
    
    
#### Pip
This will load with Anaconda if on cluster.  Otherwise, use conda or brew.
 
#### Snakemake

    pip3 install --user snakemake pyaml

#### Setup Virtual Environment

    conda create -n <env_name> python=2.7 anaconda
    source activate <env_name>
    
    
Make sure that you provide <env_name> in the config.yml file at the entry "python_virtual_env"  

#### _scikit-allel_    
    
    pip install scikit-allel[full]
    
## Running the workflow
    
#### Clone CanSel repo

    git clone https://github.com/pmonnahan/CanSel.git <project_directory>

The critical files responsible for executing the pipeline are contained in the ./workflow subdirectory contained within the cloned repo. They are:
* Snakefile
* config.yml
* cluster.yaml

The Snakefile is the primary workhouse of snakemake, which specifies the dependencies of various parts of the pipeline and coordinates their submission as jobs to the MSI cluster. No modifications to the Snakefile are necessary. However, in order for the Snakefile to locate all of the necessary input and correctly submit jobs to the cluster, both the config.yaml and cluster.yaml need to be modified. Open these files and change the entries that are indicated with 'MODIFY'.
Once these files have been modified, the entire pipeline can be run from within the cloned folder via:

    snakemake --cluster "qsub -l {cluster.l} -M {cluster.M} -A {cluster.A} -m {cluster.m} -o {cluster.o} -e {cluster.e} -r {cluster.r}" --cluster-config workflow/cluster.yaml -j 32
where -j specifies the number of jobs that can be submitted at once.

The pipeline is currently set up to run on the small queue on the mesabi cluster, which has a per-user submission limit of 500 jobs. This is more than enough for the entire pipeline, so running with -n 500 will submit all necessary jobs as soon as possible. If -j is small (e.g. 32), snakemake will submit the first 32 jobs and then submit subsequent jobs as these first jobs finish.

The attractive feature of _snakemake_ is its ability to keep track of the progress and dependencies of the different stages of the pipeline. Specifically, if an error is encountered or the pipeline otherwise stops before the final step, snakemake can resume the pipeline where it left off, avoiding redundant computation for previously completed tasks.
To run a specific part of the pipeline, do:
    
    snakemake -R <rule_name> --cluster "qsub -l {cluster.l} -M {cluster.M} -A {cluster.A} -m {cluster.m} -o {cluster.o} -e {cluster.e} -r {cluster.r}" --cluster-config workflow/cluster.yaml -j 20 --rerun-incomplete
where rule_name indicates the 'rule' (i.e. job) in the Snakefile that you wish to run.

If running locally, simply remove all text associated with --cluster.