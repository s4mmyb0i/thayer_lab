#!/bin/bash

# A script to generate random walks for all four constructs

source header.sh

# Extract the number of walks from the project directory name
for project_dir in ../../project_1000000; do
    if [[ -d "$project_dir" ]]; then
        # Get the project name and extract the number of walks
        project_name=$(basename "$project_dir")
        random_walk_number_of_walks="${project_name#project_}"  # Use 'project_' instead of 'projects_'

        # Ensure the extracted value is a valid integer
        if ! [[ "$random_walk_number_of_walks" =~ ^[0-9]+$ ]]; then
            echo "Error: Invalid random walk number of walks extracted from project name: $project_name"
            continue
        fi

        echo "Number of steps in each walk: $random_walk_steps"
        echo "Number of walks: $random_walk_number_of_walks"
        echo -e "Lower threshold value: $threshold_value_lower\n"

        list_of_constructs=("P" "PL" "AP" "APL")

        for string in "${list_of_constructs[@]}"
        do
            echo "Processing the $string construct in project $project_name"
            python3 random_walk.py "$string" "$random_walk_steps" "$random_walk_number_of_walks" "$threshold_value_lower" "project_$random_walk_number_of_walks"
            echo -e "Moving on from the $string construct\n"
        done
    fi
done

echo "Done!"