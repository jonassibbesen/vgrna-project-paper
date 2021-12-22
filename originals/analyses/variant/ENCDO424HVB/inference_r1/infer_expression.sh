set -e

# Set file name prefixes
OUT_PREFIX="rpvg_mpmap_r1_${NAME}_${TRANSCRIPTS}"

IDX=0

# Download all first reads
for READ in $(echo ${READS_1} | sed -e 's/,/ /g'); do 

	# Download reads
	IDX=$((IDX+1))
	aws s3 cp ${READ} reads_${IDX}_1.fq.gz --no-sign-request --no-progress
done

IDX=0

# Download all second reads
for READ in $(echo ${READS_2} | sed -e 's/,/ /g'); do 
	
	# Download reads
	IDX=$((IDX+1))
	aws s3 cp ${READ} reads_${IDX}_2.fq.gz --no-sign-request --no-progress
done

# Concatenate reads
/usr/bin/time -v bash -c 'cat $(for i in $(seq 1 '"${IDX}"'); do echo reads_${i}_1.fq.gz; done) > reads_1.fq.gz; zcat reads_1.fq.gz | wc -l; cat $(for i in $(seq 1 '"${IDX}"'); do echo reads_${i}_2.fq.gz; done) > reads_2.fq.gz; zcat reads_2.fq.gz | wc -l'

# Download graph
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${GRAPH}/${GRAPH}.xg . --no-progress

# Download index
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/mapping/vg/indexes/${GRAPH}/ . --recursive --no-progress

# Map reads
/usr/bin/time -v bash -c "vg mpmap -t ${CPU} -n rna -x ${GRAPH}.xg -g ${GRAPH}_index.gcsa -d ${GRAPH}_index.dist -f reads_1.fq.gz -f reads_2.fq.gz > alignments.gamp"

# Download transcripts
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${TRANSCRIPTS}/${TRANSCRIPTS}.gbwt . --no-progress
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${TRANSCRIPTS}/${TRANSCRIPTS}.gbwt.ri . --no-progress

# Download transcript info
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/${TRANSCRIPTS}/ . --recursive --exclude "*" --include "*.txt.gz" --no-progress

# Concatenate transcript info
/usr/bin/time -v bash -c "zcat */*.txt.gz | head -n 1 > transcript_info.txt; zcat */*.txt.gz | grep -v ^Name >> transcript_info.txt; wc -l transcript_info.txt"

# Infer expression
/usr/bin/time -v bash -c "rpvg -t ${CPU} -r ${SEED} -e rf -n 1000 -i haplotype-transcripts -g ${GRAPH}.xg -p ${TRANSCRIPTS}.gbwt -a alignments.gamp -f transcript_info.txt -o ${OUT_PREFIX}"

# Compress expression values
/usr/bin/time -v bash -c "gzip ${OUT_PREFIX}.txt; gzip ${OUT_PREFIX}_joint.txt"

# Upload expression 
aws s3 sync . s3://vg-k8s/users/jsibbesen/vgrna/analyses/variant/ENCDO424HVB/inference_r1/${NAME}/${TRANSCRIPTS}/ --exclude "*" --include "${OUT_PREFIX}*" --no-progress
