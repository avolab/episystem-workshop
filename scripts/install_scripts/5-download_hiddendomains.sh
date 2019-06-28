#!/bin/bash

cd $HOME
link="https://sourceforge.net/projects/hiddendomains/files/hiddenDomains.3.0.tar.gz/download"
outfile="$HOME/hiddenDomains.3.0.tar.gz"

wget $link -O $outfile

