set -e

cd data/sgd/

# Use intermine to get the latest SGD alleles.
# TODO: include unique identifier
python get_sgd_alleles.py alleles_sgd_raw.tsv

# TODO: download the latest genome
# convert the gff to embl, using emblmygff3 docker image (see translation*.json), which are used for the transformation
bash convert_sgd_gff2embl.sh

# Download all previous protein sequence (visit the repo), as well as the current protein sequences, to make the dictionary.
curl -kL https://raw.githubusercontent.com/pombase/all_previous_sgd_peptide_sequences/master/all_previous_seqs.tsv -o all_previous_seqs.tsv
curl http://sgd-archive.yeastgenome.org/sequence/S288C_reference/orf_protein/orf_trans_all.fasta.gz -o current_protein_seqs.fasta.gz
gzip -fd current_protein_seqs.fasta.gz

cd ../..

# Extract allele descriptions from their name or description field.
# TODO: use a description field provided by SGD, this applies also to all commands below
# that use _description_name or description_semicolon
python format_alleles_sgd.py

# Load the genome to a pickle file
python load_genome.py --output data/sgd/genome.pickle --config data/sgd/config.sgd.json data/sgd/genome_embl_files/*.embl

# Remove unknown ids (not in gff), or pseudogene (YLL016W), no main feature (YJL018W)
# TODO: Check why these are missing
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


python allele_auto_fix.py --genome data/sgd/genome.pickle\
                          --coordinate_changes_dict data/sgd/coordinate_changes_dict.json\
                          --allele_results results/sgd/allele_description_name_qc.tsv\
                          --output_dir    results/sgd/description_name

python allele_auto_fix.py --genome data/sgd/genome.pickle\
                        --coordinate_changes_dict data/sgd/coordinate_changes_dict.json\
                        --allele_results results/sgd/allele_description_semicolon_qc.tsv\
                        --output_dir    results/sgd/description_semicolon

cat results/sgd/*/allele_auto_fix.tsv|grep _fix > results/sgd/allele_qc_fixed.tsv


python allele_transvar.py\
    --genome data/sgd/genome.pickle\
    --allele_results results/sgd/allele_description_name_qc.tsv\
    --exclude_transcripts data/frame_shifted_transcripts.tsv\
    --output results/sgd/description_name/allele_description_name_transvar.tsv\
    --genome_fasta data/sgd/genome_sequence.fsa\
    --transvardb data/sgd/features.gtf.transvardb

python allele_transvar.py\
    --genome data/sgd/genome.pickle\
    --allele_results results/sgd/allele_description_semicolon_qc.tsv\
    --exclude_transcripts data/frame_shifted_transcripts.tsv\
    --output results/sgd/description_semicolon/allele_description_semicolon_transvar.tsv\
    --genome_fasta data/sgd/genome_sequence.fsa\
    --transvardb data/sgd/features.gtf.transvardb

# TODO
# - Error with YSC0029, present in alleles, but missing in the genome gff
# - Mapping of allele types
# - Deal with examples like `YBR038W	CHS2	chs2-S4E	S4E		S4	phosphomimetic mutant S14E S60E S69E S100E`
#   where it seems like S4E is the description.
# - Use the allele unique identifiers
