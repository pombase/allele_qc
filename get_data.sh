mkdir -p data
mkdir -p results

# Get phenotype annotations from PomBase
curl -k https://www.pombase.org/data/annotations/Phenotype_annotations/phenotype_annotations.pombase.phaf.gz --output  data/phenotype_annotations.phaf.gz
gzip -d data/phenotype_annotations.phaf.gz

# Get unique lines with allele types, and remove deletion and wild-type alleles
cut -d $'\t' -f 2,4,9,10,11,12,18 data/phenotype_annotations.phaf|sort|uniq|grep -v $'\t'deletion|grep -v wild_type > data/alleles_pre_format.tsv
python format_alleles.py data/alleles_pre_format.tsv > data/alleles.tsv

# Get peptide sequences from PomBase
curl -k https://www.pombase.org/data/genome_sequence_and_features/artemis_files/chromosome1.contig --output  data/chromosome1.contig
curl -k https://www.pombase.org/data/genome_sequence_and_features/artemis_files/chromosome2.contig --output  data/chromosome2.contig
curl -k https://www.pombase.org/data/genome_sequence_and_features/artemis_files/chromosome3.contig --output  data/chromosome3.contig
curl -k https://www.pombase.org/data/genome_sequence_and_features/artemis_files/mating_type_region.contig --output  data/mating_type_region.contig
curl -k https://www.pombase.org/data/genome_sequence_and_features/artemis_files/pMIT.contig --output  data/pMIT.contig
curl -k https://www.pombase.org/data/genome_sequence_and_features/artemis_files/telomeric.contig --output  data/telomeric.contig
