#!/bin/bash


awk 'BEGIN{FS=OFS="\t"}{print $2}' $1 | tail -n+2 | sed -E 's|DBL-[0-9]+|DBL|g' | sort | uniq -c | sed -E 's/^ +//g' | sed 's/ /\t/g' | sed '1 i\Assignment N\tClassification' | awk 'BEGIN{FS=OFS="\t"}{print($2,$1)}'