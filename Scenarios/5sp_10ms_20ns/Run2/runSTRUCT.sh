#!/bin/bash/

# To get the start processing time, get the current system time, taken at the beginning of the run, and print it
now=$(date)
echo "Current date/time: $now"
# End processing time can be taken from the output file timestamps (i.e. nohup.out)

# STRUCTURE command
# K = 8, five replicates
for r in {1..5} ;
do
   /usr/local/bin/structure -i parentandhybrid2.str -m mainparams.txt -e extraparams.txt -K 5 -o parentandhybrid5Deme_$r > Outputs/parentandhybrid5Deme_$r &
   sleep 3s
done

echo "All runs started!"
