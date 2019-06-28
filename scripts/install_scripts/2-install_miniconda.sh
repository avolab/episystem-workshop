#!/bin/bash
# Install miniconda

installdir="$HOME"
installpath="${installdir}/Miniconda3-latest-Linux-x86_64.sh"
[[ ! -e $installpath ]] && echo "$installpath not found, exiting" && exit 1

cd $installdir

bash $installpath

