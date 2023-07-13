# This script was used to check whether some errors had been introduced due to a bug
# See issue: https://github.com/pombase/allele_qc/issues/60

# Use the original data (pre-changes) from the changelog, by piping it into data/pombase-chado.modifications
cat change_log/protein_modification_auto_fix_external_data_21032023.tsv|cut -f1-9 > data/pombase-chado.modifications

# Run the usual pipeline
python protein_modification_qc.py
python protein_modification_auto_fix.py

# Compare the results with the original data
cat results/protein_modification_auto_fix.tsv |cut -f1-9,11-12|sort> new_file.tsv
cat change_log/protein_modification_auto_fix_external_data_21032023.tsv|cut -f1-9,11-12|sort>old_file.tsv

# See that the only difference is the row that had two solution K1001 and other.
diff -u new_file.tsv old_file.tsv

# Thankfully no errors had been introduced by the bug