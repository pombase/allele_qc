# Pull the docker image
docker pull quay.io/biocontainers/emblmygff3:2.2--pyhdfd78af_1
docker run --rm -v $PWD/scripts:/scripts:ro -v $PWD/data:/data quay.io/biocontainers/bioconvert:1.1.1--pyhdfd78af_0 /bin/bash -c \
        "bioconvert gff32gtf --force data/pombe_genome.gff3 data/pombe_genome_unprocessed.gtf;\