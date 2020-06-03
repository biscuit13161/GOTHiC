#!/usr/env bash

filename=$1
if [ `echo $filename | grep "\.bam$"` ]; then
	echo bam
	samtools view -h $1 | grep -v "^@" | cut -f -4 > $1.txt
else
	echo sam
	grep -v "^@" $1 | cut -f -4 > $1.txt
fi

