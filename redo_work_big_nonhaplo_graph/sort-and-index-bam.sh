#! /usr/bin/env bash

# The first command-line arg to this function is a BAM file. ls should
# complain if this file does not exist. NOTE that this is an extremely
# shoddy way of doing this, in part because it'll print out the filename
# to the console even upon success. A better solution would be e.g.
# https://stackoverflow.com/a/638980
ls $1

# # Based on http://quinlanlab.org/tutorials/samtools/samtools.html
samtools sort $1 -T _tmpfile -o $2

echo "sorted the BAM file $1: the sorted BAM file is named $2"

samtools index $2

echo "indexed the sorted BAM file $2"
