#!/bin/env bash
# =========================================================
# Copyright 2012-2018,  Nuno A. Fonseca (nuno dot fonseca at gmail dot com)
#
# This file is part of iRAP.
#
# This is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with iRAP.  If not, see <http://www.gnu.org/licenses/>.
#
#
# =========================================================
# Script to merge a set of quantification matrices
FILES=$*
echo $*  1>&2

## Read the list of files from stdin
NAMES=$FILES
if [ "$*-" == "-stdin-" ]; then
    set -e
    read -e -t  60 NAMES
    set -- $NAMES
    set +e
    FILES=$NAMES
else
    if [ "$1-" == "-" ]; then
	echo "ERROR: no file names provided" 
	exit 1
    fi
fi

#set -e
f1=$1
shift 1

if [ ! -e $f1 ]; then
    echo "ERROR: file $f1 not found"
    exit 1
fi

if [ "$MAX_THREADS-" == "-" ]; then
    MAX_THREADS=1
fi
if [ "$FILES_PER_THREAD-" == "-" ]; then
    FILES_PER_THREAD=500
fi

function my_wait {
    MAX_THREADS=$1
    set +e
    let jobs_run=0
    let jobs_run=$(jobs -r | wc -l)
    set -e		
    if [ $jobs_run -gt $MAX_THREADS ]; then
	echo -n "waiting for slot..."  1>&2
	builtin wait -n
	ret=$?
	if [ "$ret-" != "0-" ] && [ "$ret-" != "127-" ] ; then
	    echo "Failed to process batch."  1>&2
	    exit 1
	fi
	echo "done."  1>&2
    fi    
}

if [ "$MAX_THREADS" != "1" ]; then
    # split the list of files into chunks of $FILES_PER_THREAD
    echo "Using $MAX_THREADS threads..."  1>&2
#    set -eux
    set -e
    ## check if version of bash is recent enough
    bash_version=$(bash --version|head -n 1|cut -f 4 -d\ |cut -f 1 -d\()
    major=$(echo $bash_version|cut -f 1 -d\.)
    minor=$(echo $bash_version|cut -f 2 -d\.)
    rev=$(echo $bash_version|cut -f 3 -d\.)
    if [ $major -lt 4 ] ||
	   ( [ $major -eq 4 ] && [ $minor -lt 4 ] ) ||
	   ( [ $major -eq 4 ] && [ $minor -eq 4 ] &&  [ $rev -lt 18 ]) ; then
	echo "Invalid version of bash: expected 4.4.18+ and got $bash_version"  1>&2
	exit 1
    fi
    ## bask ok...carry on
    let chunk=1
    fname_prefix=`mktemp .tmp.$$.XXXXXXXXXX`
    set +e
    let files_chunk=0
    set -e
    sel_files=""
    ofiles=""
    for f in $FILES; do
	sel_files="$sel_files $f"
	let files_chunk=$files_chunk+1
	if [ $files_chunk == $FILES_PER_THREAD ]; then
	    ofile=$fname_prefix.$chunk.tsv
	    my_wait $MAX_THREADS
	    irap_merge_tsv.sh -stdin <<EOF > $ofile &
$sel_files
EOF
	    ofiles="$ofiles $ofile"
	    sel_files=""
	    let chunk=$chunk+1
	    set +e
	    let files_chunk=0
	    set -e
	fi
	#echo $files_chunk  1>&2
    done
    # last chunk
    if [ $files_chunk -gt 0 ]; then
	my_wait $MAX_THREADS
	ofile=$fname_prefix.$chunk.tsv
	irap_merge_tsv.sh -stdin <<EOF > $ofile &
$sel_files
EOF
	ofiles="$ofiles $ofile"
	set +e
	let files_chunk=0
	set -e
    fi
    echo "Waiting for all jobs to finish..."   1>&2    
    set +e
    builtin wait
    if [ "$?-" != "0-" ]; then
	rm -f $fname_prefix $ofiles
	echo "Failed to merge."  1>&2
	exit 1
    fi
    set -e
    echo Merging $chunk files   1>&2
    # merge all tmp files (in parallel if necessary)
    ofile=$fname_prefix.tsv
    if [ $chunk == 1 ]; then
	## avoid infinite recursion
	export MAX_THREADS=1
    fi
    FILES_PER_THREAD=50 irap_merge_tsv_p.sh -stdin <<EOF > $ofile 
$ofiles
EOF
    cat $fname_prefix.tsv
    #rm -f $fname_prefix $fname_prefix.*    
else
    ## just call irap_merge_tsv.sh
    set -e
irap_merge_tsv.sh -stdin <<EOF
$FILES
EOF
fi
exit 0
