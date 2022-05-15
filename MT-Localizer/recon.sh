#!/bin/bash

echo $PWD
recon-all -i ../data/fsfast/AnatDir/t1.nii.gz -s S1 -all -openmp 4
