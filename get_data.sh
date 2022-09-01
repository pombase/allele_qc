mkdir -p data
mkdir -p results
# Get phenotype annotations from PomBase
curl -k https://www.pombase.org/data/annotations/Phenotype_annotations/phenotype_annotations.pombase.phaf.gz --output  data/phenotype_annotations.phaf.gz
gzip -d data/phenotype_annotations.phaf.gz

# Get unique lines with allele types, and remove deletion and wild-type alleles
cut -d $'\t' -f 2,4,9,10,11,12 data/phenotype_annotations.phaf|uniq|grep -v $'\t'deletion|grep -v wild_type > data/alleles.tsv
