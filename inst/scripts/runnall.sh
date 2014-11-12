#!/bin/sh

for f in $( ls *R); 
do
	Rscript $f;
done

