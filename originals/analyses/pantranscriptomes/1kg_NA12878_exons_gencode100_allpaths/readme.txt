aws s3 sync s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/variants/1kg_NA12878_exons/ . --profile vg-dev --exclude "*" --include "*vcf.gz*"
aws s3 sync s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/1kg_NA12878_exons_gencode100_allpaths/ . --profile vg-dev --exclude "*" --include "*txt.gz"
