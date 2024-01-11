#!/bin/bash

# set up the PYTHONPAATH must be sourced before running nuSIM

# Last Update: 2023-01-23				Author:  Paul Kyberd
# add 04-Studies to the python path

stm=$PWD
#echo $stm

dir="$stm"
echo "Set nuSIM path:"
nuSIMPATH="$dir"
echo "    " $nuSIMPATH
export nuSIMPATH

add1="/01-Code"
add2="/12-examples"
add3="/04-Studies"
#echo $add
dir1="$stm$add1"
dir2="$stm$add2"
dir3="$stm$add3"

#echo $dir
echo "Set PYTHON path:"
PYTHONPATH="${PYTHONPATH}:$dir1:$dir2:$dir3"
echo "    " $PYTHONPATH
export PYTHONPATH
