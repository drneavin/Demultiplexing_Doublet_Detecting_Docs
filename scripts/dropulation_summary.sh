#!/bin/bash

gunzip -c $1 | tail -n +2 | awk '{print $3}' | sort | uniq -c