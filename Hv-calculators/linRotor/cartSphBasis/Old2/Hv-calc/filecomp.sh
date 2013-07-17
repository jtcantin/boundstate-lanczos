#!/bin/sh

#if diff linAlgCompCNvect.cube linAlgCompMatVec.cube > /dev/null ; then
#  echo Same
#else
#  echo Different
#fi

if diff HvFullVectorOutput2013-06-26-1526 test2 ; then
  echo Same
else
  echo Different
fi