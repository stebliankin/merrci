# Database Preparation

This step is designed to prepare the database on a single machine and than transfer the files to the cluster. 
We don't run this step on a cluster, because the network on distrubuted system is usually slower, but bush files in the folder "preprcessingModules"
can be modified for job submission by providing the Slurm specific headings.

### Requirments
Before preprocessing step, make sure you have installed the following dependancies on the local machine you run the scripts:

* [Python](https://www.python.org/)
* [BioPython](https://biopython.org/)
* [Pandas](https://pandas.pydata.org/)
* [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* [Kraken2](https://ccb.jhu.edu/software/kraken2/)


## Steps

First, you need to set up the config file "configLocal.sh" for your local machine:
* SCRIPTS_DIR_LOCAL - path to RAPToR scripts on a local machine
* DB_LOCATION - location of RAPToR database on local machine
* TREADS - number of threads to use in the preprocessing step

All other paths will be set automatically.

Second, run the preprocessing:

    ./run_preprocessing.sh

Note, that kraken2 version 2.0.8 has a bag in downloading reference genomes:
 
    rsync_from_ncbi.pl: unexpected FTP path (new server?) for na
 
 To fix it, simply add the following after line 38 in the rsync_from_ncbi.pl file of Kraken2 repository:

    next if($ftp_path eq "na");

Once the database has built, you need to transfer it to the cluster:

    scp -r $DB_LOCATION $userName"@"$hostName":"$DB_CLUSTER
    
Where $DB_CLUSTER is the location where you want to put RAPToR database in the cluster environment
