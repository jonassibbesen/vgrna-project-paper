set -e

echo "Setting up"

# args: VARS POPULATION BAM POPULATION SAMPLE CPU

# get dependencies
sudo apt update
sudo apt install -qq -y python3
python3 -c 'print("successfully installed python")'
sudo apt install -qq -y time awscli git-all python3 python3-pip samtools
sudo pip install numpy pysam tables

# get WASP
git clone https://github.com/bmvdgeijn/WASP.git

echo "Downloading files"

# Download variant index
aws s3 cp s3://vg-k8s/users/jeizenga/wasp/h5${VARS}/ . --recursive --no-progress
aws s3 cp s3://vg-k8s/users/jeizenga/wasp/samples/${POPULATION}_samples.txt . --no-progress

# get reads
aws s3 cp s3://vg-k8s/users/jeizenga/wasp/${BAM} mapped_reads.bam --no-progress
samtools sort --threads ${CPU} mapped_reads.bam > mapped_reads.sorted.bam
samtools index mapped_reads.sorted.bam

echo "Computing ASE"

/usr/bin/time -v python WASP/mapping/get_as_counts.py. --snp_tab snp_tab.h5 --snp_index snp_index.h5  --haplotype haplotypes.h5 --samples ${POPULATION}_samples.txt --genotype_sample ${SAMPLE} mapped_reads.sorted.bam > ase_tab.txt


aws cp ase_tab.txt s3://vg-k8s/users/jeizenga/wasp/ase/${BAM}_ase_`date +"%N"`.txt --no-progress

