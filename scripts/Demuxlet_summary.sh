#!/bin/bash


awk 'BEGIN{OFS=FS="\\t"}{print $5,$6}' $1/*OUT_impute_vars.best  | awk 'BEGIN{FS=OFS="\t"} $1=="DBL" {$2="doublet"}1' | awk 'BEGIN{OFS=FS="\t"}{print $2}' | sed -E "s/,[0-9]+_[0-9]+,[0-9].[0-9]+\t/\t/g" | 
awk 'BEGIN{OFS=FS="\t"}{print $5,$6}' *OUT_impute_vars.best  | awk 'BEGIN{FS=OFS="\t"} $1=="DBL" {$2="doublet"}1' | awk 'BEGIN{OFS=FS="\t"}{print $2}' | sed -E "s/,.+,[0-9]\.[0-9]+//g" | tail -n+2 | sort | uniq -c | sed -E 's/^ +//g' | sed -E 's/ /\t/g' | awk 'BEGIN{FS=OFS="\t"}{print($2,$1)}' | sed '1 i\Classification\tAssignment N' > $1/Demuxlet_summary.tsv
