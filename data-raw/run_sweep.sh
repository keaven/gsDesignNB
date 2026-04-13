#!/bin/bash
cd /Users/Anderkea/Documents/GitHub/gsDesignNB
nohup Rscript data-raw/jensen_correction_broad_sweep.R > data-raw/jensen_sweep.log 2>&1 &
echo "PID: $!"
