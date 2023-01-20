set -e
bash get_data.sh
python build_alignment_dict.py
python protein_modification_qc.py
python protein_modification_auto_fix.py
python allele_qc.py
python allele_auto_fix.py
