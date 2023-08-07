mkdir -p lib_src

GREEN='\033[0;32m'
NC='\033[0m' # No Color

root_dir=$(pwd)

echo -e "${GREEN}Downloading libraries${NC}"

curl -kL https://github.com/samtools/samtools/releases/download/1.18/samtools-1.18.tar.bz2 -o lib_src/samtools-1.18.tar.bz2
tar -xjf lib_src/samtools-1.18.tar.bz2 -C lib_src/
rm lib_src/samtools-1.18.tar.bz2

curl -kL https://github.com/samtools/htslib/releases/download/1.18/htslib-1.18.tar.bz2 -o lib_src/htslib-1.18.tar.bz2
tar -xjf lib_src/htslib-1.18.tar.bz2 -C lib_src/
rm lib_src/htslib-1.18.tar.bz2

mkdir -p lib
mkdir -p lib/samtools
mkdir -p lib/htslib

echo -e "${GREEN}Installing samtools${NC}"

cd ${root_dir}/lib_src/samtools-1.18
./configure --prefix="${root_dir}/lib/samtools"
make
make install
rm -rf "${root_dir}/lib/samtools/share"

echo -e "${GREEN}Installing htslib${NC}"

cd ${root_dir}/lib_src/htslib-1.18
./configure --prefix="${root_dir}/lib/htslib"
make
make install

rm -rf "${root_dir}/lib/htslib/lib"
rm -rf "${root_dir}/lib/htslib/include"
rm -rf "${root_dir}/lib/htslib/share"

cd "${root_dir}"
rm -rf lib_src
