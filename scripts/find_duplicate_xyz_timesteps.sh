#!/bin/bash
echo "The file $1 has the following timesteps duplicated"
grep i $1 | awk '{print $3}' | sed 's/,//g' | awk '{lsum[$1]++} END { for (i in lsum) if (lsum[i] > 1) print i, "(" lsum[i] ")"}'
