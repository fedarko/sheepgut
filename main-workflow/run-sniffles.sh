#! /usr/bin/env bash
# NOTE: this is currently unused in the paper. if you'd like to run this
# script, you will need to install sniffles first.

OUTVCF=output/sniffles-minimap2-results.vcf

# Stop the script when one of the steps breaks.
# Apparently there is Discourse (tm) about whether or not this is good
# practice, but I'm gonna go ahead and use it anyway. See
# https://stackoverflow.com/q/19622198 for details.
set -e

# NOTE: right now i'm getting a weird error about one of the reads not having
# the same number of bases as quality values from NGMLR. I extracted this read
# with zcat / grep and I don't see a problem; haven't been able to replicate
# the error with just a small subset of the reads. For now, I'm just gonna use
# the minimap2 alignment with Sniffles.

## See align-reads-to-edges.sh (which uses minimap2 instead of ngmlr) for
## details about this. basically, it's a way of avoiding having to repeatedly
## specify the read file locations: just yank the locations from the config file
#read READS_FILES < ../config/input-reads
#echo "Input reads: `echo $READS_FILES`"
#echo "`echo $READS_FILES | wc -w` reads files, total."
#
#echo "Running ngmlr..."
## ngmlr doesn't like it when I specify -q (the reads file) multiple times.
## However, ngmlr does accept read files if they're just piped in via stdin;
## so we can use multiple reads files as input by just catting all of them
## at ngmlr. Hack for this c/o https://unix.stackexchange.com/a/20286.
#cat $READS_FILES | ngmlr \
#    -t 4 \
#    -r output/all_edges.fasta \
#    -o output/ngmlr-alignment.sam
#echo "Done running ngmlr."
#exit
#
#echo "Converting ngmlr alignment to BAM..."
#samtools view -b output/ngmlr-alignment.sam > output/ngmlr-alignment.bam
#echo "Done."
#
#./sort-and-index-bam.sh output/ngmlr-alignment.bam output/ngmlr-aln-sorted.bam

echo "Running sniffles on the sorted/indexed BAM file..."
sniffles \
    --ccs_reads \
    -m output/fully-filtered-and-sorted-aln.bam \
    -v $OUTVCF
echo "Done running sniffles."

# Extract SVs relevant to the three selected edges
# https://stackoverflow.com/a/5694596
cat $OUTVCF | grep $'^edge_6104\t' > output/edge_6104_sv.vcf
cat $OUTVCF | grep $'^edge_1671\t' > output/edge_1671_sv.vcf
cat $OUTVCF | grep $'^edge_2358\t' > output/edge_2358_sv.vcf
