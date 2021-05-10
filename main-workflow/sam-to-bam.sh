#! /usr/bin/env bash

echo "Converting to BAM..."

samtools view -b output/alignment.sam > output/alignment.bam

echo "Done."
