set -e

cd data/sgd/

# Download and extract gff, delete the rest
curl -kL http://sgd-archive.yeastgenome.org/sequence/S288C_reference/genome_releases/S288C_reference_genome_Current_Release.tgz -o sgd_genome.tgz
tar -xvzf sgd_genome.tgz -C .
rm sgd_genome.tgz
mv S288C*/*.gff.gz sgd_genome.gff.gz
rm -rf S288C*
gzip -fd sgd_genome.gff.gz

# Split into fasta and gff
perl -ne 'if ($found) { print; } elsif (m/##FASTA/) { $found = 1; }' sgd_genome.gff > genome_sequence.fsa
perl -ne 'print; last if /##FASTA/' sgd_genome.gff > features.gff

# Remove problematic ARS lines that are not used
grep -v 'SGD	ARS' features.gff > features.gff.tmp
mv features.gff.tmp features.gff

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
