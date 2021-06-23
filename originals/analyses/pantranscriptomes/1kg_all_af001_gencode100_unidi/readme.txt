aws s3 sync s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/variants/1kg_all_exons/ . --profile vg-dev --exclude "*" --include "*vcf.gz*"
aws s3 sync s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/1kg_all_af001_gencode100_unidi/ . --profile vg-dev --exclude "*" --include "*txt.gz"
