set -e

# python allele_qc.py --genome data/sgd/genome.pickle\
#                     --alleles data/sgd/alleles_description_name.tsv\
#                     --output    results/sgd/allele_description_name_qc.tsv

# python allele_qc.py --genome data/sgd/genome.pickle\
#                     --alleles data/sgd/alleles_description_semicolon.tsv\
#                     --output    results/sgd/allele_description_semicolon_qc.tsv


python allele_auto_fix.py --genome data/sgd/genome.pickle\
                          --coordinate_changes_dict data/sgd/coordinate_changes_dict.json\
                          --allele_results results/sgd/allele_description_name_qc.tsv\
                          --output_dir    results/sgd/description_name

python allele_auto_fix.py --genome data/sgd/genome.pickle\
                        --coordinate_changes_dict data/sgd/coordinate_changes_dict.json\
                        --allele_results results/sgd/allele_description_semicolon_qc.tsv\
                        --output_dir    results/sgd/description_semicolon


# python allele_transvar.py\
#     --genome data/sgd/genome.pickle\
#     --allele_results results/sgd/allele_description_name_qc.tsv\
#     --exclude_transcripts data/frame_shifted_transcripts.tsv\
#     --output results/sgd/description_name/allele_description_name_transvar.tsv\
#     --genome_fasta data/sgd/genome_sequence.fsa\
#     --transvardb data/sgd/features.gtf.transvardb\
#     --sgd_mode True

# python allele_transvar.py\
#     --genome data/sgd/genome.pickle\
#     --allele_results results/sgd/allele_description_semicolon_qc.tsv\
#     --exclude_transcripts data/frame_shifted_transcripts.tsv\
#     --output results/sgd/description_semicolon/allele_description_semicolon_transvar.tsv\
#     --genome_fasta data/sgd/genome_sequence.fsa\
#     --transvardb data/sgd/features.gtf.transvardb\
#     --sgd_mode True

# Temporary removal of type_fix cases
cat results/sgd/description_name/allele_auto_fix.tsv|grep -v type_error > results/sgd/description_name/allele_auto_fix.tsv.tmp
mv results/sgd/description_name/allele_auto_fix.tsv.tmp results/sgd/description_name/allele_auto_fix.tsv

cat results/sgd/description_semicolon/allele_auto_fix.tsv|grep -v type_error > results/sgd/description_semicolon/allele_auto_fix.tsv.tmp
mv results/sgd/description_semicolon/allele_auto_fix.tsv.tmp results/sgd/description_semicolon/allele_auto_fix.tsv


# TODO
# - Error with YSC0029, present in alleles, but missing in the genome gff
# - Mapping of allele types
# - Deal with examples like `YBR038W	CHS2	chs2-S4E	S4E		S4	phosphomimetic mutant S14E S60E S69E S100E`
#   where it seems like S4E is the description.
# - Use the allele unique identifiers
