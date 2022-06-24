set -e

# Set mapper name
MAPPER="mpmap"

# Set alignment (GAM or GAMP) prefix
ALIGN_PREFIX="mpmap_1kg_nonCEU_af001_gencode100_sim_vg_ENCSR000AED_rep1_uni_h1"

# Set graph (XG) prefix
GRAPH_PREFIX="1kg_nonCEU_af001_gencode100"

# Set simulated alignment (GAM) prefix
SIM_PREFIX="sim_1kg_NA12878_gencode100_ENCSR000AED_rep1_uni_vg_h1"

# Set output name prefix
OUT_PREFIX="mpmap_1kg_nonCEU_af001_gencode100_sim_vg_ENCSR000AED_rep1_uni_h1"

# Set number of threads
CPU=1

if [ -f "${ALIGN_PREFIX}.gamp" ]; then

	# Convert to gamp to gam
	/usr/bin/time -v bash -c "vg view -K -G ${ALIGN_PREFIX}.gamp > ${ALIGN_PREFIX}.gam"
fi

# Calculate alignment statistics
/usr/bin/time -v bash -c "vg stats -a ${ALIGN_PREFIX}.gam"

# Calculate distance statistics
/usr/bin/time -v bash -c "vg gampcompare -t ${CPU} -d -a ${MAPPER} -G ${GRAPH_PREFIX}.xg ${ALIGN_PREFIX}.gam ${SIM_PREFIX}.gam | grep -v ^distance | sort -nr -k4 | sort -ns -k1 | sort -su -k6 | cut -f1-5 | sort -k1n -k2n -k3n -k4n | uniq -c > ${OUT_PREFIX}_dist_gam.txt; gzip ${OUT_PREFIX}_dist_gam.txt"

if [ -f "${ALIGN_PREFIX}.gamp" ]; then

	# Calculate distance statistics
	/usr/bin/time -v bash -c "vg gampcompare -t ${CPU} -d -a ${MAPPER} ${GRAPH_PREFIX}.xg ${ALIGN_PREFIX}.gamp ${SIM_PREFIX}.gam | grep -v ^distance | sort -nr -k4 | sort -ns -k1 | sort -su -k6 | cut -f1-5 | sort -k1n -k2n -k3n -k4n | uniq -c > ${OUT_PREFIX}_dist_gamp.txt; gzip ${OUT_PREFIX}_dist_gamp.txt"
fi
