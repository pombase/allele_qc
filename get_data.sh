mkdir -p data
mkdir -p results

# Get phenotype annotations from PomBase
curl -k https://www.pombase.org/data/annotations/Phenotype_annotations/phenotype_annotations.pombase.phaf.gz --output  data/phenotype_annotations.phaf.gz
gzip -d data/phenotype_annotations.phaf.gz

# Get unique lines with allele types, and remove deletion and wild-type alleles
cut -d $'\t' -f 2,4,9,10,11,12,18 data/phenotype_annotations.phaf|sort|uniq|grep -v $'\t'deletion|grep -v wild_type > data/alleles_pre_format.tsv
python format_alleles.py data/alleles_pre_format.tsv > data/alleles.tsv

# Get genome sequence with annotations
curl -k https://www.pombase.org/releases/latest/pombe-embl/chromosome1.contig --output  data/chromosome1.contig
curl -k https://www.pombase.org/releases/latest/pombe-embl/chromosome2.contig --output  data/chromosome2.contig
curl -k https://www.pombase.org/releases/latest/pombe-embl/chromosome3.contig --output  data/chromosome3.contig
curl -k https://www.pombase.org/releases/latest/pombe-embl/mating_type_region.contig --output  data/mating_type_region.contig
curl -k https://www.pombase.org/releases/latest/pombe-embl/pMIT.contig --output  data/pMIT.contig
curl -k https://www.pombase.org/releases/latest/pombe-embl/telomeric.contig --output  data/telomeric.contig

# Store the genome as a dictionary using pickle
python load_genome.py