set -e

# Set indexing mode. 
# 	"def": Default
# 	"decoy": Use genome as decoy
MODE="def"

# Set genome (FASTA) prefix
GENOME_PREFIX="Homo_sapiens.GRCh38.dna.primary_assembly"

# Set transcripts (FASTA) prefix
TRANSCRIPTS_PREFIX="1kg_nonCEU_af001_gencode100"

# Set output name prefix
OUT_PREFIX="salmon_index_1kg_nonCEU_af001_gencode100"

# Set number of threads
CPU=1

if [ "${MODE}" = "def" ]; then

	# Construct salmon index
	/usr/bin/time -v bash -c "salmon index -p ${CPU} --keepDuplicates -i ${OUT_PREFIX} -t ${TRANSCRIPTS_PREFIX}.fa.gz"

elif [ "${MODE}" = "decoy" ]; then

	# Generate decoy list
	/usr/bin/time -v bash -c "grep '"'>'"' ${GENOME_PREFIX}.fa | cut -d ' ' -f 1 | sed -e 's/>//g' > decoys.txt; wc -l decoys.txt"

	# Construct salmon index
	/usr/bin/time -v bash -c "salmon index -p ${CPU} --keepDuplicates -d decoys.txt -i ${OUT_PREFIX} -t ${TRANSCRIPTS_PREFIX}.fa.gz; rm "
fi
