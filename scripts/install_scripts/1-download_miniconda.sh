#!/bin/bash
# Download miniconda insatllation script

outdir=$HOME
link="https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"

cd $outdir

wget $link
