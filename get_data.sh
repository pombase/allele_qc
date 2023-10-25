set -e

mkdir -p data
mkdir -p results
GREEN='\033[0;32m'
NC='\033[0m' # No Color


echo -e "${GREEN}getting allele data${NC}"
# Get phenotype annotations from PomBase
curl -k https://curation.pombase.org/dumps/latest_build/exports/all-allele-details.tsv --output  data/alleles_raw.tsv

# Replace the first line by the right column names
echo -ne 'systematic_id\tgene_name\tallele_id\tallele_name\tallele_description\tallele_type\treference\n' > data/alleles.tsv
tail -n +2 data/alleles_raw.tsv|sort|uniq|grep -v $'\t'deletion|grep -v wild_type >> data/alleles.tsv

# Download phaf and canto alleles to filter out the ones without annotations
curl -k https://curation.pombase.org/data/pombe-allele-table.tsv --output data/alleles_canto.tsv
curl -k https://www.pombase.org/data/annotations/Phenotype_annotations/phenotype_annotations.pombase.phaf.gz --output  data/phenotype_annotations.phaf.gz
gzip -fd data/phenotype_annotations.phaf.gz
cut -f 2,4,9,10,11,12,18 data/phenotype_annotations.phaf|sort|uniq|grep -v $'\t'deletion|grep -v wild_type > data/alleles_phaf.tsv
python filter_alleles_pombase.py data/alleles.tsv data/alleles_canto.tsv data/alleles_phaf.tsv
rm data/alleles_raw.tsv data/alleles_canto.tsv data/alleles_phaf.tsv

echo -e "${GREEN}Getting contig files${NC}"
# Get genome sequence with annotations
curl -k https://curation.pombase.org/dumps/latest_build/pombe-embl/chromosome1.contig --output  data/chromosome1.contig
curl -k https://curation.pombase.org/dumps/latest_build/pombe-embl/chromosome2.contig --output  data/chromosome2.contig
curl -k https://curation.pombase.org/dumps/latest_build/pombe-embl/chromosome3.contig --output  data/chromosome3.contig
curl -k https://curation.pombase.org/dumps/latest_build/pombe-embl/mating_type_region.contig --output  data/mating_type_region.contig
curl -k https://curation.pombase.org/dumps/latest_build/pombe-embl/pMIT.contig --output  data/pMIT.contig
curl -k https://curation.pombase.org/dumps/latest_build/pombe-embl/telomeric.contig --output  data/telomeric.contig

# Store the genome as a dictionary using pickle
echo -e "${GREEN}Loading the genome${NC}"
python load_genome.py data/*.contig

# Get updates to genome coordinates from PomBase, these can be used to update alleles that used
# previous gene feature coordinates.
echo -e "${GREEN}Getting coordinate changes${NC}"
curl -k https://raw.githubusercontent.com/pombase/genome_changelog/master/results/only_modified_coordinates.tsv --output data/only_modified_coordinates.tsv

# Get a genome versions where genome sequences changed. This is important to retrieve the
# sequences of features for which the coordinates were defined in old genomes.
echo -e "${GREEN}Getting genome versions where genome sequence changed${NC}"

# Table of when changes in sequence happenned
curl -k https://raw.githubusercontent.com/pombase/genome_changelog/master/results/genome_sequence_changes.tsv --output data/genome_sequence_changes.tsv

mkdir -p data/old_genome_versions/

# Download those versions
# We only want to download the old versions, so we skip the most recent ones
skip_recent_version=""
while read row; do
    year=$(echo "$row"|cut -f2 -d$'\t'|cut -f1 -d'-')
    revision=$(echo "$row"|cut -f1 -d$'\t')
    contig=$(echo "$row"|cut -f3 -d$'\t')
    output_folder="data/old_genome_versions/$contig"
    output_file="${output_folder}/${revision}.contig"
    if [[ $skip_recent_version != *"$contig"* ]]; then
        skip_recent_version="$skip_recent_version $contig"
        continue
    fi
    mkdir -p $output_folder
    if test -s $output_file;then
        continue
    fi
    # We know last change in pre_svn was in 2007
    if (( $year < 2008)); then
        curl -k https://www.pombase.org/data/genome_sequence_and_features/artemis_files/OLD/${revision}/${contig}.contig > $output_file
    else
        svn export --force -r ${revision} https://curation.pombase.org/pombe-embl-repo/trunk/${contig}.contig $output_file
    fi

done < <(tail -n +2 data/genome_sequence_changes.tsv)

echo -e "${GREEN}Downloading protein modification data${NC}"
curl -k https://www.pombase.org/data/annotations/modifications/pombase-chado.modifications.gz --output data/pombase-chado.modifications.gz
gzip -fd data/pombase-chado.modifications.gz

# Download mod ontology to check which amino_acids are allowed for which MOD terms
curl -kL http://purl.obolibrary.org/obo/mod.obo -o data/mod.obo
python make_mod_dict.py

# Download the pombase genome fasta and gtf files (for transvar)
curl -k https://curation.pombase.org/dumps/latest_build/fasta/chromosomes/Schizosaccharomyces_pombe_all_chromosomes.fa.gz -o data/pombe_genome.fa.gz
gzip -fd data/pombe_genome.fa.gz

# TODO
# curl -k path/to/gtf_file/in/nightly/release -o data/pombe_genome.gtf
