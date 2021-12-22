aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/reads/sim/1kg_NA12878_gencode100/SRR1153470/vg_r1/ . --recursive --exclude "*" --include "*txt.gz" --profile vg-dev
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/expression/1kg_NA12878_gencode100/SRR1153470/1kg_NA12878_gencode100_SRR1153470_rsem.isoforms.results . --profile vg-dev
