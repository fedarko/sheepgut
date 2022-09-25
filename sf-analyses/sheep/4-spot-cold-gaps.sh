#! /usr/bin/env bash

strainFlye spot cold-gaps \
    --bcf output/call-p50/naive-calls.bcf \
    --min-length 5000 \
    --output-coldspots output/coldspots-p50calls-minlen5000-nocirc-exactpvals.tsv
