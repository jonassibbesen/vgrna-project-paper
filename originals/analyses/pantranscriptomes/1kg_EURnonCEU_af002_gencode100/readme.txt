aws s3 sync s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/variants/1kg_EURnonCEU_exons/ . --profile vg-dev --exclude "*" --include "*vcf.gz*"
aws s3 sync s3://vg-k8s/users/jsibbesen/vgrna/benchmark/whole_genome/data/graphs/1kg_EURnonCEU_af002_gencode100/ . --profile vg-dev --exclude "*" --include "*txt.gz"
