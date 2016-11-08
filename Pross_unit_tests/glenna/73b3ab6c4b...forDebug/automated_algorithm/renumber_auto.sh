#!/bin/sh

#renumber a pdb file from 1->nres, ignoring chainID, TER statements and anything else taht doesn't have an ATOM card.

pdb=$1

awk 'BEGIN {num=0}
  {
  if( $1=="ATOM" ){
    if( prev_num != substr($0,23,4) ) ++num;
    prev_num = substr($0,23,4);
    printf substr( $0, 1, 22 );
    printf "%4d",num
    print substr( $0,27,10000);
  }
  else print $0
  }' $pdb #> $pdb.TMP

#mv $pdb.TMP $pdb
