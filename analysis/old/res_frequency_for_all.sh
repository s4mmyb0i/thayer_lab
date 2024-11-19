#!/bin/bash

# A script to run which_captured.py for all four constructs

list_of_constructs=("P" "PL" "AP" "APL")

for string in "${list_of_constructs[@]}"
do
    echo "Processing the $string construct"
    python3 res_frequency.py "$string"
    echo -e "Moving on from the $string construct\n"
done

echo "Done!"