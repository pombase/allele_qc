set -e
bash get_data.sh
python build_alignment_dict.py
uvicorn api:app --host 0.0.0.0 --port 80