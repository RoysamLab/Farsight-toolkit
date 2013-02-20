#!/bin/bash
IFS='
'
j=`find $1 -printf "%d\n" | sort -u | tail -n 1` 
echo "Max dir depth:" $j
for (( i=0; i<=j ; i++ ))
do
  for name in `find -mindepth $i -maxdepth $i -iname "* *" -printf "%p\n"`
    do
        newname=`echo "$name" | tr " " "_"`
	echo "$name" "$newname"
	mv "$name" "$newname"
    done
done
