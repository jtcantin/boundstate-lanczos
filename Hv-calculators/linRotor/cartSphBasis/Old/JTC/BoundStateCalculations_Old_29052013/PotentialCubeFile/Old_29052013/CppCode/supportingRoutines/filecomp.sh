#!/bin/sh

#if diff linAlgCompCNvect.cube linAlgCompMatVec.cube > /dev/null ; then
#  echo Same
#else
#  echo Different
#fi

if diff sII_unitCelltest_CM_EA_ZYZ.xyz test ; then
  echo Same
else
  echo Different
fi