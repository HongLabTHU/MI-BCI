#!/bin/bash

# set
set -e

# parse args
POSITIONAL=()
# help\
if [ $# -lt 2 ]
then
    echo "-p|--project-dir project root directory \n
          --all run all process \n
          --surf run surface analysis \n
          --volume run native volume space analysis.\n"
    exit 1;
fi

run_surf=false;
run_volume=false;

while [[ $# -gt 0 ]]
do 
	key="$1"

    case $key in
        -p|--project-dir)
		project="$2"
		shift
		shift
		;;
	    --all)
	    run_surf=true;
	    run_volume=true;
	    shift
	    ;;
	    --surf)
	    run_surf=true;
	    shift
	    ;;
	    --volume)
	    run_volume=true;
	    shift
	    ;;
		*)
		POSITIONAL+=("$1")
		shift
		;;
	esac
done
# set positional arguments
set -- "${POSITIONAL[@]}"

# cd root dir
cd $project
echo $PWD

sess_f="./sessid"

# preprocess
preproc-sess -d $PWD -sf $sess_f -surface self lhrh -mni305 -per-run -sliceorder odd -fwhm 5 -fsd bold

if [ "${run_surf}" = true ]; then
    # lh
    mkanalysis-sess -surface self lh -fwhm 5 -analysis mt.lh -fsd bold -stc odd -TR 3 -event-related -paradigm mt.par -nconditions 1 -fslhrf 1 -refeventdur 21 -polyfit 2 -nskip 4 -mcextreg -per-run -force
    mkcontrast-sess -analysis mt.lh -contrast mt -a 1 -c 0
    # rh
    mkanalysis-sess -surface self rh -fwhm 5 -analysis mt.rh -fsd bold -stc odd -TR 3 -event-related -paradigm mt.par -nconditions 1 -fslhrf 1 -refeventdur 21 -polyfit 2 -nskip 4 -mcextreg -per-run -force
    mkcontrast-sess -analysis mt.rh -contrast mt -a 1 -c 0

    # run
    # lh
    selxavg3-sess -d $PWD -sf $sess_f -a mt.lh -max-threads
    # rh
    selxavg3-sess -d $PWD -sf $sess_f -a mt.rh -max-threads
fi

if [ "${run_volume}" = true ]; then
    # native space
    # native
    mkanalysis-sess -native -fwhm 5 -analysis mt.native -fsd bold -stc odd -TR 3 -event-related -paradigm mt.par -nconditions 1 -fslhrf 1 -refeventdur 21 -polyfit 2 -nskip 4 -mcextreg -per-session -force
    mkcontrast-sess -analysis mt.native -contrast mt -a 1 -c 0

    # native
    selxavg3-sess -d $PWD -sf $sess_f -a mt.native -max-threads
fi

exit 0;
