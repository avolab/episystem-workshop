#!/bin/sh
# Jake Yeung
# 9-run_hiddendomains.sh
#  
# 2019-06-26

# export PATH="$PATH:/Users/yeung/data/scchic/for_episys/hiddenDomains"
# hd=
# export PATH="$HOME/hiddenDomains$PATH:$HOME/hiddenDomains"
export PATH="$HOME/hiddenDomains:$PATH"
# export PATH="$PATH:/gnu/store/ks8rldr7skvj9x7qaa7fgds8a7rpz5pl-hiddendomains-3.0/bin"
maindir="/projects/cdbfc3e9-b653-448c-b278-d583cb89eba1/episystem-workshop/data"
outmain="/projects/cdbfc3e9-b653-448c-b278-d583cb89eba1/episystem-workshop/outputs"

chromsizes="$maindir/chromsizes/chromsizes.mm10.filt.txt"
[[ ! -e $chromsizes ]] && echo "$chromsizes not found, exiting" && exit 1

bamdir1="$maindir/sorted_bams_filtered"
# bamdir2="$maindir/sorted_bams_filtered_H3K4me1"

minlength=1000

minpost=0.3
outdir="$outmain/hiddenDomains_output"
mkdir $outdir
mapq=30

for bamdir in $bamdir1; do
	for inf in $(ls -d $bamdir/*.bam); do
		bname=$(basename $inf)
		bname=${bname%%.*}
		bname=${bname}.minpost_${minpost}.mapq_${mapq}
		# echo "hiddenDomains -g $chromsizes -b $minlength -t $inf -o $outdir/$bname -p $minpost -q $mapq"
		hiddenDomains -g $chromsizes -b $minlength -t $inf -o $outdir/$bname -p $minpost -q $mapq
		# hiddenDomains -g $chromsizes -b $minlength -t $inf -o $outdir
		# exit 0
	done
done
