#!/bin/bash

stm=$PWD
#echo $stm

dir="$stm"
echo "Set nuSIM path:"
nuSIMPATH="$dir"
echo "    " $nuSIMPATH
export nuSIMPATH

add="/01-Code"
#echo $add
dir="$stm$add"
#echo $dir
echo "Set PYTHON path:"
PYTHONPATH="${PYTHONPATH}:$dir"
echo "    " $PYTHONPATH
export PYTHONPATH
