export TRANSVAR_CFG="$(pwd)/data/transvar.cfg"
export TRANSVAR_DOWNLOAD_DIR="$(pwd)/data/transvar_download"

# These are only needed if you install transvar locally, but it does not harm to have them
export PATH="$(pwd)/lib/htslib/bin:$PATH"
export PATH="$(pwd)/lib/samtools/bin:$PATH"