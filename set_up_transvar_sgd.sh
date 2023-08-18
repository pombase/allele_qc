set -e

# Set up transvar database
transvar config -k reference -v "data/sgd/genome_sequence.fsa" --refversion sgd_genome
transvar index --ensembl data/sgd/features.gtf --reference data/sgd/genome_sequence.fsa

# Hacky way to use the functions inside another script
cp $(which transvar) ./transvar_main_script.py

# See whether it works
transvar panno -i 'Q0050:p.V2A' --ensembl data/sgd/features.gtf.transvardb --reference data/sgd/genome_sequence.fsa

# transvar ganno -i 'I:g.2832796A>T' --ensembl data/sgd/features.gtf.transvardb --reference data/sgd/genome_sequence.fsa
# transvar ganno -i 'I:g.2832795T>A' --ensembl data/sgd/features.gtf.transvardb --reference data/sgd/genome_sequence.fsa
# transvar panno -i 'SPBC1198.04c:p.T566S' --ensembl data/sgd/features.gtf.transvardb --reference data/sgd/genome_sequence.fsa
# transvar panno -i 'SPAPB1A10.09:p.S372_N374delinsAAA' --ensembl data/sgd/features.gtf.transvardb --reference data/sgd/genome_sequence.fsa  --gseq
