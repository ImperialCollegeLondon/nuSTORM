#!/bin/bash

# set up the PYTHONPAATH must be sourced before running nuSIM

stm=$PWD
#echo $stm

dir="$stm"
echo "Set nuSIM path:"
nuSIMPATH="$dir"
echo "    " $nuSIMPATH
export nuSIMPATH

add="/01-Code"
add2="/12-examples"
#echo $add
dir="$stm$add"
dir2="$stm$add2"
#echo $dir
echo "Set PYTHON path:"
PYTHONPATH="${PYTHONPATH}:$dir:$dir2"
echo "    " $PYTHONPATH
export PYTHONPATH
