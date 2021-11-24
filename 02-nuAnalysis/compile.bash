#!/bin/bash

g++ 01-Code/RunControl.cpp 02-Tests/RunControlTst.cpp -I01-Code `root-config --cflags --libs` -o 12-Bin/RunControlTst.cpp
