#!/bin/bash

# This script is used to intiate multiple STRUCTURE runs, all of them using the same specified K value,
# specified number of replicates, and mainparams.txt/extraparams.txt parameter files.
# It does this by reading in all of the STRUCTURE input files in the specified directory,
# and for each, generating an output directory, to which STRUCTURE output files are written.

# %%% VARIABLE DECLARATIONS %%%
# Specified variables: use a while statement and the getopts package to specify which letters correspond to which variables
while getopts "i:r:k:" arg; do
# case/esac syntax, for specified variables
  case $arg in
	# Argument passed to the BASH script, specifying the directory containing the STRUCTURE input files (.str; only these files!)
	i) structInFolder=$OPTARG;;
	# Number of replicates. This is the number of times (for EACH RUN) that STRUCTURE is run.
	# Note that this is not the number of STRUCTURE runs to perform. That's determined by the number of files in the inputFiles directory specified above
	r) repNumber=$OPTARG;;
	# K value to used in STRUCTURE model
	k) kValue=$OPTARG;;
  esac
done
# Fixed varialbes
# Building the list of all the input files in the directory, by using wildcard expansion
inputFiles="$structInFolder/*.str"
# Declare COUNTER variable (used below to create output directories)
COUNTER=1
# Specify filepath to STRUCTURE executable
STR_PATH=/usr/local/bin/structure

# %%% PRINT VARIABLE SPECIFICATIONS %%%
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
echo "   STRUCTURE input files to be processed: $inputFiles"
echo "   Number of STRUCTURE replicate runs, for each input file: $repNumber"
echo "   K value to use in STRUCTURE run: $kValue"
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"

# %%% PRINT START TIME %%%
# To get the start processing time, get the current system time, taken at the beginning of the run, and print it
now=$(date)
printf "\n"
echo "%%% Current date/time: $now"
# End processing time can be taken from the output file timestamps (i.e. nohup.out)

# %%% START STRUCTURE RUNS %%%
# OUTER LOOP: for each file within the input file directory (where inFile is the input file being processed, in that directory), do the following
for inFile in $inputFiles
do
	printf "\n"
	echo "%%% Processing input file $inFile"
	# Make an output directory, to store STRUCTURE outputs
	out_dir=Run${COUNTER}_Outputs
	mkdir $out_dir
	# INNER LOOP: iterate based on the specified number of replicates. Loop uses C++ syntax (3 parameter loop control expression)
	for((r=1; r<=repNumber; r++)); do
		echo "Running STRUCTURE. Output file: $out_dir/parentAndHybrid_${kValue}_${r}"
		# Run STRUCTURE using current input file, with mainparams/extraparams in the specified directory.
		# -o argument: specifies current output file (iterates with inner loop)
		# > outfile: specifies the file to take printed text and output it (not used for downstream processing)
		# &: tells BASH to run STRUCTURE in the background
		$STR_PATH -i $inFile -m mainparams.txt -e extraparams.txt -K $kValue -o $out_dir/parentAndHybrid_${kValue}_${r} > outfile_${COUNTER}_${r} &
                sleep 5s
        done
	# Increment outer loop COUNTER variable, for output names
	((COUNTER++))
done
printf "\n"
echo "%%% All STRUCTURE runs have been initiated!"
