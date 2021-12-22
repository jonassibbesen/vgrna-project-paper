aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/reads/sim/1kg_NA12878_gencode100/ENCSR000AED_rep1/vg_r1/ . --recursive --exclude "*" --include "*txt.gz" --profile vg-dev
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/expression/1kg_NA12878_gencode100/ENCSR000AED_rep1/1kg_NA12878_gencode100_ENCSR000AED_rep1_rsem.isoforms.results . --profile vg-dev
