python format_alleles_sgd.py

# python load_genome.py --output data/sgd/genome.pickle --config data/sgd/config.sgd.json data/sgd/genome_embl_files/*.embl

# Remove unknown ids (not in gff), or pseudogene (YLL016W), no main feature (YJL018W)

missing_genes="R0010W YSC0029 R0040C YLL016W YSC0032 YJL018W"

for missing_gene in $missing_genes; do
    grep -v $missing_gene data/sgd/alleles_description_name.tsv > data/sgd/alleles_description_name.tsv.tmp
    mv data/sgd/alleles_description_name.tsv.tmp data/sgd/alleles_description_name.tsv

    grep -v $missing_gene data/sgd/alleles_description_semicolon.tsv > data/sgd/alleles_description_semicolon.tsv.tmp
    mv data/sgd/alleles_description_semicolon.tsv.tmp data/sgd/alleles_description_semicolon.tsv
done


python allele_qc.py --genome data/sgd/genome.pickle\
                    --alleles data/sgd/alleles_description_name.tsv\
                    --output    results/sgd/allele_description_name_qc.tsv

python allele_qc.py --genome data/sgd/genome.pickle\
                    --alleles data/sgd/alleles_description_semicolon.tsv\
                    --output    results/sgd/allele_description_semicolon_qc.tsv


# TODO
# - Sort out mitochondrial translation
# - Error with YSC0029, present in alleles, but missing in the genome gff
# - Include references for alleles
# - Mapping of allele types
# - Deal with examples like `YBR038W	CHS2	chs2-S4E	S4E		S4	phosphomimetic mutant S14E S60E S69E S100E`
#   where it seems like S4E is the description.