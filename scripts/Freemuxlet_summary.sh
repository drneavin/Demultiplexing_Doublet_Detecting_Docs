#!/bin/bash



gunzip -c $1/*.clust1.samples.gz | awk 'BEGIN{OFS=FS="\t"}{print $5,$6}' | awk 'BEGIN{FS=OFS="\t"} $1=="DBL" {$2="DBL"}1' | awk 'BEGIN{FS="\t"}{print $2}' | sed -E 's/,[0-9]+//g'|  tail -n+2 | sort | uniq -c | sed -E 's/^ +//g' | sed -E 's/ /\t/g' | awk 'BEGIN{FS=OFS="\t"}{print($2,$1)}' | sed '1 i\Classification\tAssignment N' > $1/freemuxlet_summary.tsv
