ps
#!/bin/bash

stm=$PWD
#echo $stm

dir="$stm"
echo "Set nuAnalysis path:"
nuAnalysisPATH="$dir"
echo "    " $nuAnalysisPATH
export nuAnalysisPATH

add="/01-Code"
#echo $add
dir="$stm$add"
#echo $dir
echo "Set nuAnalysisCPATH:"
nuAnalysisCPATH="$dir"
echo "    " $nuAnalysisCPATH
export nuAnalysisCPATH
