#!/bin/bash

g++ 01-Code/RunControl.cpp 01-Code/nuAnalysis.cpp 03-Skeleton/nuAnalysis-skeleton.cpp -I01-Code `root-config --cflags --libs` -o 12-Bin/nuAnalysis-skeleton.exe
