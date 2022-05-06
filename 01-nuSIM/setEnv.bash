#!/bin/bash

#
#	Initialises things to run a nuSIM study. 
#	First time it is run it sets up the whole system and
#	creates an empty file nuSIM.init.
#	On subsequent runs it checks for the existance of the
#	file and then only sets up the required logical names
#
#	Version 1.0							14 April 2022
#	Paul Kyberd
#
#	first run startup.bash from the nuSIM directory. This script
#	assumes nuSIMPATH is set up correctly


#   checking if this script has been run in this directory
initFile='nuSIM.init'
if test -f "$initFile"; then
	echo "      $initFile exists ... just setting up logical names"
	# still need to set the StudyDir
	StudyDir="$PWD"
	export StudyDir
	echo What is the name of the study 
	read StudyName
	if test -d "$StudyName"; then
		export StudyName
	else
		echo "StudyName directory does not exist - try again"
		return
	fi
else
	echo "Need to iniitalise"
	#  first create the initialisation flag file
	touch nuSIM.init
	# Set the StudyDir
	StudyDir="$PWD"
	export StudyDir

	# now copy across the runNumber file from nuSIM
	ln -sf ${nuSIMPATH}/runNumber runNumber

	# now get the name of the study

	echo What is the name of the study 
	read StudyName
	mkdir $StudyName
	export StudyName

	# copy across a control file for Pi decay in the production straight, Pion Flash
	cp $nuSIMPATH/04-Studies/PSPiFlash.dict .
	sed "s/<studyname>/$StudyName/" PSPiFlash.dict > PSPiFlash.dict.tmp
	mv PSPiFlash.dict.tmp PSPiFlash.dict
	# use the # delimiter to avoid clash with the directory /
	sed "s#<studyDir>#$StudyDir#" PSPiFlash.dict > PSPiFlash.dict.tmp
	mv PSPiFlash.dict.tmp PSPiFlash.dict
	# copy across a plots control file for Pi and mu in the production straight
	cp $nuSIMPATH/04-Studies/plotsPSPiFlash.dict $StudyName/
	
	# copy across a control file for Pi decay in the production straight and muon decay in the ring: 
	#		muon signal
	cp $nuSIMPATH/04-Studies/MuRingDcy.dict .
	sed "s/<studyname>/$StudyName/" MuRingDcy.dict > MuRingDcy.dict.tmp
	mv MuRingDcy.dict.tmp MuRingDcy.dict
	# use the # delimiter to avoid clash with the directory /
	sed "s#<studyDir>#$StudyDir#" MuRingDcy.dict > MuRingDcy.dict.tmp
	mv MuRingDcy.dict.tmp MuRingDcy.dict
	# copy across a plots control file for Pi and mu in the production straight
	cp $nuSIMPATH/04-Studies/plotsMuRingDcy.dict $StudyName/

	
	# copy across a control file for muon decay in the production straight. pion flash neutrnios ignored: 
	#		muon background
	cp $nuSIMPATH/04-Studies/PSMuDcy.dict .
	sed "s/<studyname>/$StudyName/" PSMuDcy.dict > PSMuDcy.dict.tmp
	mv PSMuDcy.dict.tmp PSMuDcy.dict
	# use the # delimiter to avoid clash with the directory /
	sed "s#<studyDir>#$StudyDir#" PSMuDcy.dict > PSMuDcy.dict.tmp
	mv PSMuDcy.dict.tmp PSMuDcy.dict
	# copy across a plots control file for Pi and mu in the production straight
	cp $nuSIMPATH/04-Studies/plotsPSMuDcy.dict $StudyName/


fi



echo "Checking the environment variables for nuSIM\n"

echo "nuSIMPATH is $nuSIMPATH\n"

echo "PYTHONPATH is $PYTHONPATH\n"

echo "StudyDir is $StudyDir\n"

echo "StudyName is $StudyName\n"
