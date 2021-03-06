#!/usr/bin/env bash
# merge a set of CRAM files
# cram_merge.sh file_with_cram_filenames out_file
# * one cram filename per line

if [ "$2-" == "-" ]; then
    echo "ERROR: usage cram_merge.sh file_with_cram_filenames out_file" 
    exit 1
fi

if [ ! -e $1 ]; then
    echo "ERROR: file $1 not found"
    exit 1
fi

# ensure that the files are one per line
# assumes no spaces in the path
cat $1| tr " " "\n" | grep -v "^$" > $2.lst

# check if all files exist
set -e
cat $2.lst |  ( while read f; do
		    echo -n "Checking if $f exists..."
		    if [ ! -e "$f" ]; then
			echo 
			echo "ERROR: file $f not found"
			exit 1
		    fi
		    echo "ok."
		    done
)
set +e
# get reference
cram1=`head -n 1 $2.lst`
echo "Getting reference from $cram1..."
ref=`samtools view -H $cram1 |grep UR|sed "s/.*UR://g;s/ .*//"|head -n 1`
echo "Getting reference from $cram1...done."
echo "ref=$ref"
if [ "$ref-" == "-" ]; then
    samtools view -H $cram1 
    echo "ERROR:Unable to find a reference in the cram file?!"
    exit 1
fi
if [ ! -e "$ref" ]; then
    echo "Reference file $ref not found"
    exit 1
fi
set -e
samtools merge -f -b $2.lst $2.tmp.bam
samtools view -C -T $ref $2.tmp.bam -o $2.tmp.cram && mv $2.tmp.cram $2

rm -f $2.tmp.bam $2.lst
exit 0

