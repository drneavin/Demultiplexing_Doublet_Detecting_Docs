#!/bin/bash

gunzip -c $1 | tail -n +2 | awk '{print $3}' | sort | uniq -c | sed -E 's/^ +//g' | sed -E 's/ /\t/g' | awk 'BEGIN{FS=OFS="\t"}{print($2,$1)}' | sed '1 i\Classification\tAssignment N' 