#!/bin/bash

# Script to copy a Python script to multiple directories, run it, and delete it afterward.

# Arguments: path to script location and the script name
location=$1
script=$2
echo "Script that will be run: $script"

# Define the projects and constructs
projects=(project_100 project_1000 project_10000 project_100000 project_1000000)
constructs=(P PL AP APL)

# Function to run the script in the target directories
run() {
    project=$1
    echo "Processing project: $project"
    
    for construct in "${constructs[@]}"; do
        # Define the target directory for each construct
        target_dir="../$project/analysis/$construct"
        
        # Check if the target directory exists
        if [ -d "$target_dir" ]; then
            # Copy the script to the target directory
            # cp "./"
            cp "./$location/$script" "$target_dir/"
            
            # Move to the target directory
            cd "$target_dir" || { echo "Failed to enter directory $target_dir"; exit 1; }

            # Run the script
            python3 "$script" "$construct"

            # Delete the copied script
            rm "$script"

            # Go back to the parent directory
            cd - || exit 1
        else
            echo "Directory $target_dir does not exist, skipping."
        fi
    done
}

# Loop through each project directory and run the function
for project in "${projects[@]}"; do
    run "$project"
done