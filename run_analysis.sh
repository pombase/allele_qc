set -e
bash get_data.sh
python build_alignment_dict_from_genome.py

# Check and fix protein modification
python protein_modification_qc.py
python protein_modification_auto_fix.py
python protein_modification_transvar.py

# Check and fix allele descriptions, types and names
python allele_qc.py
python allele_auto_fix.py
python allele_transvar.py
