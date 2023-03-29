set -e
# No longer needed, since data folder is committed to git
# bash get_data.sh
# python build_alignment_dict.py
uvicorn api:app --host 0.0.0.0 --port 80