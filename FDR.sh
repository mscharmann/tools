#!/bin/sh

### from:
## http://blog.mcbryan.co.uk/2013/01/false-discovery-rates-and-large-files.html

if [ $# -ne 2 ]
then
echo "usage: PtoFDR file column"
exit 1
fi
 
export LC_ALL=C
sort -S 10G -k"$2","$2"gr $1 | \
  awk -v col=$2 -v numrows=`wc -l $1 | awk '{print $1}'` '
  function min(a,b)
  {
  if (a <= b)
    return a
  else
    return b
  }
 
  BEGIN { cummin = 1.0; OFS="\t"; }
  {
    cummin = min(cummin,$col*(numrows/(numrows - NR + 1)));
    for (i = 1; i <= col ; i++){
      printf $i"\t"
    }
    printf cummin"\t";
 
    for (i = col+1; i <= NF; i++){
      printf $i"\t"
    }
    printf "\n"
  }'

