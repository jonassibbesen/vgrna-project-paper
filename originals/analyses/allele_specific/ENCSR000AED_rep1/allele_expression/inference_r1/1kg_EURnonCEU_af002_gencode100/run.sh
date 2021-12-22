for c in $(seq 1 22; echo X; echo Y); do cd ${c}; screen -d -m bash -c "bash ../calc_allele_exp.sh ${c} > calc_allele_exp_log.txt 2>&1"; cd ..; done
