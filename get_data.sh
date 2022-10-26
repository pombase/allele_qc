mkdir -p data
mkdir -p results
GREEN='\033[0;32m'
NC='\033[0m' # No Color

echo -e "${GREEN}getting allele data${NC}"
# Get phenotype annotations from PomBase
curl -k https://www.pombase.org/data/annotations/Phenotype_annotations/phenotype_annotations.pombase.phaf.gz --output  data/phenotype_annotations.phaf.gz
gzip -fd data/phenotype_annotations.phaf.gz

# Get unique lines with allele types, and remove deletion and wild-type alleles
cut -d $'\t' -f 2,4,9,10,11,12,18 data/phenotype_annotations.phaf|sort|uniq|grep -v $'\t'deletion|grep -v wild_type > data/alleles_pre_format.tsv
python format_alleles.py data/alleles_pre_format.tsv > data/alleles.tsv

echo -e "${GREEN}Getting contig files${NC}"
# Get genome sequence with annotations
curl -k https://www.pombase.org/releases/latest/pombe-embl/chromosome1.contig --output  data/chromosome1.contig
curl -k https://www.pombase.org/releases/latest/pombe-embl/chromosome2.contig --output  data/chromosome2.contig
curl -k https://www.pombase.org/releases/latest/pombe-embl/chromosome3.contig --output  data/chromosome3.contig
curl -k https://www.pombase.org/releases/latest/pombe-embl/mating_type_region.contig --output  data/mating_type_region.contig
curl -k https://www.pombase.org/releases/latest/pombe-embl/pMIT.contig --output  data/pMIT.contig
curl -k https://www.pombase.org/releases/latest/pombe-embl/telomeric.contig --output  data/telomeric.contig

echo -e "${GREEN}Getting fasta files${NC}"
# Get individual sequences in fasta format
curl -k https://www.pombase.org/releases/latest/fasta/feature_sequences/cds+introns+utrs.fa.gz --output data/cds+introns+utrs.fa.gz
curl -k https://www.pombase.org/releases/latest/fasta/feature_sequences/cds+introns.fa.gz --output data/cds+introns.fa.gz
curl -k https://www.pombase.org/releases/latest/fasta/feature_sequences/cds.fa.gz --output data/cds.fa.gz
curl -k https://www.pombase.org/releases/latest/fasta/feature_sequences/five_prime_utrs.fa.gz --output data/five_prime_utrs.fa.gz
curl -k https://www.pombase.org/releases/latest/fasta/feature_sequences/three_prime_utrs.fa.gz --output data/three_prime_utrs.fa.gz
curl -k https://www.pombase.org/releases/latest/fasta/feature_sequences/introns_within_cds.fa.gz --output data/introns_within_cds.fa.gz
curl -k https://www.pombase.org/releases/latest/fasta/feature_sequences/peptide.fa.gz --output data/peptide.fa.gz
gzip -fd data/*.fa.gz

# Store the genome as a dictionary using pickle
echo -e "${GREEN}Loading the genome${NC}"
python load_sequences.py data/*.fa
python load_genome.py data/*.contig

# Get updates to genome coordinates from PomBase, these can be used to update alleles that used
# previous gene feature coordinates.
echo -e "${GREEN}Getting coordinate changes${NC}"
curl -k https://raw.githubusercontent.com/pombase/genome_changelog/master/all_coordinate_changes_file.tsv --output data/all_coordinate_changes_file.tsv
