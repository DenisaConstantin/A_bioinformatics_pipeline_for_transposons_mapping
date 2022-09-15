# A_bioinformatics_pipeline_for_transposons_mapping

Before using this pipeline, make sure you have SEQKIT and BLAST installed.

BLAST
  sudo apt-get install ncbi-blast+
  
SEQKIT
  https://bioinf.shenwei.me/seqkit/download/

1_genome.sh
This script processes the genomes (all chromosome from FlyBase) for the following scripts. It must be run without D. melanogaster in the folder. Subsequently, it is added and the following command is run:
makeblastdb -in dmel.fasta -input_type fasta -dbtype nucl -title dmel -out dmel

2_remnants_processing.sh
Following this script, several folders are obtained and named after the remnants of the transposable element that is analyzed. Each folder contains two files; a file with the remnant and a file with the BLAST result (remnant vs D. melanogaster).

Folder: ./folder/
Remnant: H{}, HB{}, Tc1{}, Tc1-2{}, etc.

3_upstream_extraction.sh/3_downstream_extraction.sh
In each folder there is the coord.txt file, the file with the extracted sequence and a folder dedicated to BLAST results.

To use this pipeline:
  1- download the archive (all chromosome from FlyBase);
  2- unzip archive;
  3- move the genomes to the "genome" folder;
  4- make sure the scripts are marked as executable (Properties -> Permissions -> Allow executing file as program);
  5- run the scripts in the following order:
    1_genome.sh;
    2_remnants_processing.sh;
    3_upstream_extraction.sh or 3_downstream_extraction.sh.
