set -e

echo "Setting up"

# Set number of threads
CPU=1

# get dependencies
#echo $PATH
#ls -lab /usr/bin/
sudo apt update
sudo apt install -qq -y python2
#echo "After installing python2"
#ls -lab /usr/bin/
python2 -c 'print "successfully installed python2" '
sudo ln -s `which python2` /usr/bin/python
python -c 'print "successfully aliased python" '
sudo apt install -qq -y time awscli git-all samtools gzip samtools make wget bedtools default-jre zip curl sed bcftools
curl https://bootstrap.pypa.io/pip/2.7/get-pip.py --output get-pip.py
sudo python2 get-pip.py
sudo pip2 install numpy pandas scipy

echo "Downloading STAR"

# get STAR
wget https://github.com/alexdobin/STAR/archive/refs/tags/2.7.9a.tar.gz
tar -xzf 2.7.9a.tar.gz

echo "Downloading AlleleSeq2"

# get AlleleSeq2
git clone https://github.com/trgaleev/AlleleSeq2.git
sed 's/$(STAR)/$(STAR) --genomeSuffixLengthMax 150/' < AlleleSeq2/makePersonalGenome.mk > AlleleSeq2/subPersonalGenome.mk

echo "Downloading liftOver"

# get liftOver
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver
chmod +x liftOver

echo "Downloading vcf2diploid"

# get vcf2diploid
wget http://alleleseq.gersteinlab.org/vcf2diploid_v0.2.6a.zip
unzip vcf2diploid_v0.2.6a.zip

echo "Downloading inputs"

# Download input data
aws s3 cp s3://vg-k8s/users/jeizenga/alleleseq/inputs/header_primary.bam . --no-progress
aws s3 cp s3://vg-data/1kg_GRCh38/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa . --no-progress
aws s3 cp  s3://vg-k8s/users/jsibbesen/vgrna/data/transcripts/gencode29/gencode.v29.primary_assembly.annotation.gtf.gz . --no-progress
gunzip ./*.gz
sed 's/^chr//' < gencode.v29.primary_assembly.annotation.gtf > gencode.v29.primary_assembly.annotation.rename.gtf
rm gencode.v29.primary_assembly.annotation.gtf

samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa
cat <(grep KI Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai) <(grep GL Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai) | awk '{ printf "%s:1-%s\n", $1, $2 }' > scaffolds_and_decoy_regions.txt
samtools faidx -r scaffolds_and_decoy_regions.txt Homo_sapiens.GRCh38.dna.primary_assembly.fa | sed 's/:1-[0-9]\+//' > Homo_sapiens.GRCh38.dna.scaffolds_and_decoys.fa

# note: not doing Y to avoid a column-less VCF
for c in `seq 1 22; echo X`; do
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/variants/1kg_NA12878/${c}/1kg_NA12878_${c}.vcf.gz .
done
bcftools concat `ls *.vcf.gz` > variants_all.vcf
rm ./*.vcf.gz

echo "Beginning AlleleSeq indexing pipeline"
make -f ./AlleleSeq2/subPersonalGenome.mk \
        PL=./AlleleSeq2/ \
        VCF2DIPLOID_DIR=./vcf2diploid_v0.2.6a/ \
        LIFTOVER=./liftOver \
        BEDTOOLS_intersectBed=`which intersectBed` \
        STAR=./STAR-2.7.9a/bin/Linux_x86_64_static/STAR \
        SAMTOOLS=`which samtools` \
        N_THREADS=${CPU} \
        VCF_SAMPLE_ID=NA12878 \
        REFGENOME=Homo_sapiens.GRCh38.dna.primary_assembly.fa \
        ANNOTATION=gencode.v29.primary_assembly.annotation.rename.gtf \
        FILE_PATH_BAM=header_primary.bam \
        FILE_PATH_VCF=variants_all.vcf \
        OUTPUT_DIR=output \
        output/maternal.gencode.v29.primary_assembly.annotation.rename.gtf output/paternal.gencode.v29.primary_assembly.annotation.rename.gtf


for parental in paternal maternal;
do
mkdir ${parental}_indexes
./STAR-2.7.9a/bin/Linux_x86_64_static/STAR --runThreadN ${CPU} --runMode genomeGenerate --genomeDir ./${parental}_indexes/ --genomeFastaFiles `ls output/*_${parental}.fa` Homo_sapiens.GRCh38.dna.scaffolds_and_decoys.fa --sjdbGTFfile output/${parental}.gencode.v29.primary_assembly.annotation.rename.gtf --sjdbOverhang 150
done
        

# Upload index outputs
aws s3 cp --recursive ./output/ s3://vg-k8s/users/jeizenga/alleleseq/diploid_outputs/ --no-progress
aws s3 cp --recursive paternal_indexes s3://vg-k8s/users/jeizenga/alleleseq/paternal_indexes/ --no-progress
aws s3 cp --recursive maternal_indexes s3://vg-k8s/users/jeizenga/alleleseq/maternal_indexes/ --no-progress

