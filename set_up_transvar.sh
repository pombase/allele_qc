set -e

# Set up transvar database
transvar config -k reference -v "data/pombe_genome.fa" --refversion pombe_genome
transvar index --ensembl data/pombe_genome.gtf --reference data/pombe_genome.fa

# At this point, we could delete the lib folder, but it's very small

# See whether it works
transvar panno -i 'SPAC3F10.09:p.E2A' --ensembl data/pombe_genome.gtf.transvardb --reference data/pombe_genome.fa
transvar ganno -i 'I:g.2832796A>T' --ensembl data/pombe_genome.gtf.transvardb --reference data/pombe_genome.fa
transvar ganno -i 'I:g.2832795T>A' --ensembl data/pombe_genome.gtf.transvardb --reference data/pombe_genome.fa

# Hacky way to use the functions inside another script
cp $(which transvar) ./transvar_main_script.py

transvar panno -i 'SPBC1198.04c:p.T566S' --ensembl data/pombe_genome.gtf.transvardb --reference data/pombe_genome.fa
transvar panno -i 'SPAPB1A10.09:p.S372_N374delinsAAA' --ensembl data/pombe_genome.gtf.transvardb --reference data/pombe_genome.fa  --gseq
