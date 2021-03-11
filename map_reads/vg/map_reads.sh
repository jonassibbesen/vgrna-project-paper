set -e

# Set mapping mode. 
# 	"map": Use vg map
#	"mpmap": Use vg mpmap
#	"mpmap_nosplice": Use vg mpmap without novel splice-junction detection (used when mapping to gene-only graphs)
MODE="mpmap"

# Set read files
READ_1="sim_1kg_NA12878_gencode100_ENCSR000AED_rep1_uni_vg_h1_1.fq.gz"
READ_2="sim_1kg_NA12878_gencode100_ENCSR000AED_rep1_uni_vg_h1_2.fq.gz"

# Set graph prefix
GRAPH_PREFIX="1kg_nonCEU_af001_gencode100"

# Set index prefix
INDEX_PREFIX="1kg_nonCEU_af001_gencode100_index"

# Set output name prefix
OUT_PREFIX="mpmap_1kg_nonCEU_af001_gencode100_sim_vg_ENCSR000AED_rep1_uni_h1"

# Set number of threads
CPU=1

# Use faster vg map
if [ "${MODE}" = "map" ]; then

	# Map reads
	/usr/bin/time -v bash -c "vg map -t ${CPU} --try-up-to 16 --mate-rescues 32 -x ${GRAPH_PREFIX}.xg -g ${INDEX_PREFIX}.gcsa -f ${READ_1} -f ${READ_2} > ${OUT_PREFIX}.gam"

# Use default vg mpmap
elif [ "${MODE}" = "mpmap" ]; then

	# Map reads
	/usr/bin/time -v bash -c "vg mpmap -t ${CPU} -n rna -x ${GRAPH_PREFIX}.xg -g ${INDEX_PREFIX}.gcsa -d {INDEX_PREFIX}.dist -f ${READ_1} -f ${READ_2} > ${OUT_PREFIX}.gamp"

# Use vg mpmap in non-splicing mode
elif [ "${MODE}" = "mpmap_nosplice" ]; then

	# Map reads
	/usr/bin/time -v bash -c "vg mpmap -t ${CPU} -n rna --not-spliced -x ${GRAPH_PREFIX}.xg -g ${INDEX_PREFIX}.gcsa -d {INDEX_PREFIX}.dist -f ${READ_1} -f ${READ_2} > ${OUT_PREFIX}.gamp"
fi	
