#! /usr/bin/env bash
# The purpose of running CheckM is to get a sense of edge sequence quality --
# e.g. are any of them really obvious chimeras, etc.
#
# NOTE that this is a bit inefficient since (at least per its wiki, as of
# reading) CheckM calls Prodigal anyway, so in theory we could call Prodigal
# first and then pass its predicted gene info to CheckM. However, running
# Prodigal is fast and I'm not sure the time savings is worth the increase
# in complexity (being able to use mostly default parameters is nice)...
# May be worth addressing if this is extended to lots of FASTA files at once.
#
# Use of lineage_wf based on
# https://github.com/Ecogenomics/CheckM/wiki/Quick-Start
# -- may be overkill since I really just want QA, but ok with me.
#
# Also, use of -t 10 based on Misha's use of CheckM on the sheep gut dataset;
# you may want to adjust depending on your system's hardware.

checkm lineage_wf -t 10 -x fasta ../seqs/ ../seqs/checkm/
