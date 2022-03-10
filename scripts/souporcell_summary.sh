#!/bin/bash


awk 'BEGIN{FS=OFS="\t"} $2=="unassigned" {$3="unassigned"}1' $1 | awk 'BEGIN{FS=OFS="\t"}{print $3}' | sed -E 's|[0-9]+/[0-9]+|doublet|g' | tail -n+2 | sort | uniq -c | sed -E 's/^ +//g' | sed 's/ /\t/g' | sed '1 i\Assignment N\tClassification' | awk 'BEGIN{FS=OFS="\t"}{print($2,$1)}'