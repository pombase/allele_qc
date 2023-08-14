#!/bin/bash

set -e
# Pull the docker image
docker pull quay.io/biocontainers/emblmygff3:2.2--pyhdfd78af_1

# If you want to explore the directory structure:
# docker run -it -v $PWD:/workdir --workdir /workdir --name emblmygff3 quay.io/biocontainers/emblmygff3:2.2--pyhdfd78af_1

docker run --rm --workdir /workdir -v $PWD:/workdir quay.io/biocontainers/emblmygff3:2.2--pyhdfd78af_1 /bin/bash -c \
        "EMBLmyGFF3 features.gff genome_sequence.fsa \
        --topology linear \
        --molecule_type 'genomic DNA' \
        --transl_table 1  \
        --species 'Saccharomyces cerevisiae' \
        --locus_tag LOCUSTAG\
        --project_id PRJXXXXXXX \
        --use_attribute_value_as_locus_tag ID \
        -o result.embl
        "

# Split the embl file into individual chromosomes

cp "result copy.embl" result.embl

input_file="result.embl"
output_files="chrI chrII chrIII chrIV chrIX chrV chrVI chrVII chrVIII chrX chrXI chrXII chrXIII chrXIV chrXV chrXVI chrmt"
mkdri -p genome_embl_files

for output_file in $output_files; do
    perl -ne 'print; if (m#//#) {last}' result.embl >"genome_embl_files/${output_file}.embl"
    perl -i -ne 'if ($found) { print; } elsif (m#//#) { $found = 1; }' result.embl
done

# Clean up
rm result.embl
