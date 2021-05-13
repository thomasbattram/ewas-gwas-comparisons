#!/bin/bash
# ------------------------------------------------
# movement of data to rdsf
# ------------------------------------------------

rdsf_space="/projects/MRC-IEU/research/projects/ieu2/p1/042"

rsync -r data/. $rdsf_space/working/data/
