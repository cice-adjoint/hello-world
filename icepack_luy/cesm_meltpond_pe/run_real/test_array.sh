#!/bin/bash 

#num[0]="0"
#echo "${num[0]}"

numlist=(0 1 "2" 3 4)
echo "numlist"
echo "numlist(1)"
echo "numlist[1]"
echo "{numlist[1]}"

echo "$numlist"
echo "$numlist(1)"
echo "$numlist[1]"
echo "$numlist(1)"
echo "$(numlist[1])"
echo '${numlist[1]}'
echo ${numlist[2]}
echo "${numlist[2]}"
