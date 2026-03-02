#!/bin/bash

mkdir -p align quant

process_read() {
    file=$1
    basename="$(basename $file | sed 's/.trim.fastq.gz//')"
    echo "Processing $basename"
    #mkdir quant/${basename}
    #minimap2 -ax sr -t 6 --eqx -N 100 ../../references/gencode_v39.mmi $file \
    #    | samtools view -u - \
    #    | samtools sort -n@ 4 - \
    #    > align/${basename}.bam
    oarfish --model-coverage -j 6 -a align/${basename}.bam -o quant/${basename}/${basename}
}

export -f process_read

find ../../reads/long/trim/*.fastq.gz | xargs -I % -P 6 bash -c 'process_read "$@"' _ %
