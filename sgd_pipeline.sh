# set -e

# cd data/sgd/

# python get_sgd_alleles.py alleles_sgd_raw.tsv
# bash convert_sgd_gff2embl.sh
# curl -kL https://raw.githubusercontent.com/pombase/genome_changelog/master/sgd/all_previous_seqs.tsv -o all_previous_seqs.tsv
# curl http://sgd-archive.yeastgenome.org/sequence/S288C_reference/orf_protein/orf_trans_all.fasta.gz -o current_protein_seqs.fasta.gz
# gzip -fd current_protein_seqs.fasta.gz

# cd ../..

# python format_alleles_sgd.py

# python load_genome.py --output data/sgd/genome.pickle --config data/sgd/config.sgd.json data/sgd/genome_embl_files/*.embl

# # Remove unknown ids (not in gff), or pseudogene (YLL016W), no main feature (YJL018W)

# missing_genes="R0010W YSC0029 R0040C YLL016W YSC0032 YJL018W"

# for missing_gene in $missing_genes; do
#     grep -v $missing_gene data/sgd/alleles_description_name.tsv > data/sgd/alleles_description_name.tsv.tmp
#     mv data/sgd/alleles_description_name.tsv.tmp data/sgd/alleles_description_name.tsv

#     grep -v $missing_gene data/sgd/alleles_description_semicolon.tsv > data/sgd/alleles_description_semicolon.tsv.tmp
#     mv data/sgd/alleles_description_semicolon.tsv.tmp data/sgd/alleles_description_semicolon.tsv
# done


# python allele_qc.py --genome data/sgd/genome.pickle\
#                     --alleles data/sgd/alleles_description_name.tsv\
#                     --output    results/sgd/allele_description_name_qc.tsv

# python allele_qc.py --genome data/sgd/genome.pickle\
#                     --alleles data/sgd/alleles_description_semicolon.tsv\
#                     --output    results/sgd/allele_description_semicolon_qc.tsv


python allele_auto_fix.py --genome data/sgd/genome.pickle\
                          --coordinate_changes_dict data/sgd/coordinate_changes_dict.json\
                          --allele_results_errors results/sgd/allele_description_name_qc_errors.tsv\
                          --output_dir    results/sgd/description_name

# TODO
# - Error with YSC0029, present in alleles, but missing in the genome gff
# - Mapping of allele types
# - Deal with examples like `YBR038W	CHS2	chs2-S4E	S4E		S4	phosphomimetic mutant S14E S60E S69E S100E`
#   where it seems like S4E is the description.
# - Use the allele unique identifiers
