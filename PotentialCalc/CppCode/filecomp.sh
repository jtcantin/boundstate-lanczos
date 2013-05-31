#!/bin/sh

#if diff linAlgCompCNvect.cube linAlgCompMatVec.cube > /dev/null ; then
#  echo Same
#else
#  echo Different
#fi

if diff linAlgCompCNvect.cube linAlgCompMatVec.cube ; then
  echo Same
else
  echo Different
fi