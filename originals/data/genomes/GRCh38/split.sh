GENOME_DIR=/public/groups/vg/jsibbese/vg_data/1kg_GRCh38/genome/

/public/groups/vg/jsibbese/vgrna/tools/samtools-1.9/samtools faidx ${GENOME_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa $(grep -v "KI\|GL" ${GENOME_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai | cut -f1) > Homo_sapiens.GRCh38.dna.primary_assembly_chromosomes.fa

/public/groups/vg/jsibbese/vgrna/tools/samtools-1.9/samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly_chromosomes.fa
/public/groups/vg/jsibbese/vgrna/tools/samtools-1.9/samtools dict Homo_sapiens.GRCh38.dna.primary_assembly_chromosomes.fa > Homo_sapiens.GRCh38.dna.primary_assembly_chromosomes.dict

/public/groups/vg/jsibbese/vgrna/tools/samtools-1.9/samtools faidx ${GENOME_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa $(grep "KI\|GL" ${GENOME_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai | cut -f1) > Homo_sapiens.GRCh38.dna.primary_assembly_scaffolds.fa

/public/groups/vg/jsibbese/vgrna/tools/samtools-1.9/samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly_scaffolds.fa
/public/groups/vg/jsibbese/vgrna/tools/samtools-1.9/samtools dict Homo_sapiens.GRCh38.dna.primary_assembly_scaffolds.fa > Homo_sapiens.GRCh38.dna.primary_assembly_scaffolds.dict
