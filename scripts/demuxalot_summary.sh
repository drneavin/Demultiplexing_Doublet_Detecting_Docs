#!/bin/bash

gunzip -c $1 | tail -n +2 | awk '{print $2}' | sed -E "s/^.+\+.+$/doublet/g" | sort | uniq -c