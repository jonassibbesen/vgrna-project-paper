aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/analyses/variant/ENCDO424HVB/ . --recursive --exclude "*" --include "*_unidi.txt.gz" --include "*_unidi_haps.txt.gz" --profile vg-dev
aws s3 cp s3://vg-k8s/users/jsibbesen/vgrna/analyses/variant/ENCDO424HVB/ . --recursive --exclude "*" --include "*r1*_unidi.txt.gz" --include "*r1*_unidi_joint.txt.gz" --exclude "*history*" --profile vg-dev
