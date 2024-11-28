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
	if [[ ${1} == '--studyname' ]]; then
		StudyName=${2}
		echo "Study name parsed... Study Name: ${StudyName}"
	else
		echo What is the name of the study?
		read StudyName
	fi
	if test -d "$StudyName"; then
		export StudyName

	else
		echo "StudyName directory does not exist."
		echo Would you like to create directory $StudyName? Yes or No? Warning: Capitalization important!
		read CreateStudy
		if [[ "$CreateStudy" == "Yes" ]]; then
			echo "StudyName directory will be created."
			mkdir $StudyName
			export StudyName
		elif [[ "$CreateStudy" == "No" ]]; then
			echo "StudyName directory will not be created. Try again."
			return
		else
			echo "Answer not recognized. Try again."
			return
		fi
	fi
	if test -f "$StudyName/nuSIMstudy.init"; then
		echo "Reentering previous study."
		echo "Setup done."
	else
		echo "First time this study is used. Need to initialise and copy control files over."
		touch $StudyName/nuSIMstudy.init

		# copy across a control file for Pi decay in the production straight, Pion Flash
		cp $nuSIMPATH/04-Studies/PSPiFlash.dict .
		sed "s/<studyname>/$StudyName/" PSPiFlash.dict > PSPiFlash.dict.tmp
		mv PSPiFlash.dict.tmp PSPiFlash.dict
		# use the # delimiter to avoid clash with the directory /
		sed "s#<studyDir>#$StudyDir#" PSPiFlash.dict > PSPiFlash.dict.tmp
		mv PSPiFlash.dict.tmp $StudyName/PSPiFlash.dict
		# copy across a plots control file for Pi and mu in the production straight
		cp $nuSIMPATH/04-Studies/plotsPSPiFlash.dict $StudyName/

		# copy across a control file for Pi decay in the production straight and muon decay in the ring:
		#		muon signal
		cp $nuSIMPATH/04-Studies/MuRingDcy.dict .
		sed "s/<studyname>/$StudyName/" MuRingDcy.dict > MuRingDcy.dict.tmp
		mv MuRingDcy.dict.tmp MuRingDcy.dict
		# use the # delimiter to avoid clash with the directory /
		sed "s#<studyDir>#$StudyDir#" MuRingDcy.dict > MuRingDcy.dict.tmp
		mv MuRingDcy.dict.tmp $StudyName/MuRingDcy.dict
		# copy across a plots control file for Pi and mu in the production straight
		cp $nuSIMPATH/04-Studies/plotsMuRingDcy.dict $StudyName/

		# copy across a control file for muon decay in the production straight. pion flash neutrnios ignored:
		#		muon background
		cp $nuSIMPATH/04-Studies/PSMuDcy.dict .
		sed "s/<studyname>/$StudyName/" PSMuDcy.dict > PSMuDcy.dict.tmp
		mv PSMuDcy.dict.tmp PSMuDcy.dict
		# use the # delimiter to avoid clash with the directory /
		sed "s#<studyDir>#$StudyDir#" PSMuDcy.dict > PSMuDcy.dict.tmp
		mv PSMuDcy.dict.tmp $StudyName/PSMuDcy.dict
		# copy across a plots control file for Pi and mu in the production straight
		cp $nuSIMPATH/04-Studies/plotsPSMuDcy.dict $StudyName/
	fi
else
	echo "Need to initalise"
	#  first create the initialisation flag file
	touch nuSIM.init
	# Set the StudyDir
	StudyDir="$PWD"
	export StudyDir
	runNumberPath="${nuSIMPATH}/runNumber"

	if [ -e "$runNumberPath" ] && [[ $(<"$runNumberPath") =~ ^-?[0-9]+$ ]]; then
		echo "runNumber exists at $runNumberPath and is an integer. Doing nothing."
	elif [ -L "$runNumberPath" ] || ! [[ $(<"$runNumberPath") =~ ^-?[0-9]+$ ]]; then
		echo "runNumber is a symbolic link or is not an integer. Deleting it and setting to 0."
		rm -f "$runNumberPath"
		echo 0 > "$runNumberPath"
	else
		echo "runNumber does not exist. Creating it and setting to 0."
		echo 0 > "$runNumberPath"
	fi
	# now copy across the runNumber file from nuSIM
	ln -sf ${nuSIMPATH}/runNumber runNumber

	# now get the name of the study
	if [[ ${1} == '--studyname' ]]; then
		StudyName=${2}
		echo "Study name parsed... Study Name: ${StudyName}"
	else
		echo What is the name of the study?
		read StudyName
	fi
	mkdir $StudyName
	export StudyName

	#  first create the initialisation flag file in study
	touch $StudyName/nuSIMstudy.init

	# copy across a control file for Pi decay in the production straight, Pion Flash
	cp $nuSIMPATH/04-Studies/PSPiFlash.dict $StudyName/
	sed "s/<studyname>/$StudyName/" $StudyName/PSPiFlash.dict > $StudyName/PSPiFlash.dict.tmp
	mv $StudyName/PSPiFlash.dict.tmp $StudyName/PSPiFlash.dict
	# use the # delimiter to avoid clash with the directory /
	sed "s#<studyDir>#$StudyDir#" $StudyName/PSPiFlash.dict > $StudyName/PSPiFlash.dict.tmp
	mv $StudyName/PSPiFlash.dict.tmp $StudyName/PSPiFlash.dict
	# copy across a plots control file for Pi and mu in the production straight
	cp $nuSIMPATH/04-Studies/plotsPSPiFlash.dict $StudyName/

	# copy across a control file for Pi decay in the production straight and muon decay in the ring:
	#		muon signal
	cp $nuSIMPATH/04-Studies/MuRingDcy.dict $StudyName/
	sed "s/<studyname>/$StudyName/" $StudyName/MuRingDcy.dict > $StudyName/MuRingDcy.dict.tmp
	mv $StudyName/MuRingDcy.dict.tmp $StudyName/MuRingDcy.dict
	# use the # delimiter to avoid clash with the directory /
	sed "s#<studyDir>#$StudyDir#" $StudyName/MuRingDcy.dict > $StudyName/MuRingDcy.dict.tmp
	mv $StudyName/MuRingDcy.dict.tmp $StudyName/MuRingDcy.dict
	# copy across a plots control file for Pi and mu in the production straight
	cp $nuSIMPATH/04-Studies/plotsMuRingDcy.dict $StudyName/


	# copy across a control file for muon decay in the production straight. pion flash neutrnios ignored:
	#		muon background
	cp $nuSIMPATH/04-Studies/PSMuDcy.dict $StudyName/
	sed "s/<studyname>/$StudyName/" $StudyName/PSMuDcy.dict > $StudyName/PSMuDcy.dict.tmp
	mv $StudyName/PSMuDcy.dict.tmp $StudyName/PSMuDcy.dict
	# use the # delimiter to avoid clash with the directory /
	sed "s#<studyDir>#$StudyDir#" $StudyName/PSMuDcy.dict > $StudyName/PSMuDcy.dict.tmp
	mv $StudyName/PSMuDcy.dict.tmp $StudyName/PSMuDcy.dict
	# copy across a plots control file for Pi and mu in the production straight
	cp $nuSIMPATH/04-Studies/plotsPSMuDcy.dict $StudyName/


fi



echo "Checking the environment variables for nuSIM"

echo "nuSIMPATH is $nuSIMPATH"

echo "PYTHONPATH is $PYTHONPATH"

echo "StudyDir is $StudyDir"

echo "StudyName is $StudyName"
