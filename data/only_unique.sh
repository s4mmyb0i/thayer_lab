#!/bin/bash

# A script to save unique paths (from the random_walks files) in new files

list_of_constructs=("P" "PL" "AP" "APL")

# Define the function to process each file
process_file() {
    local project="$1"
    local construct="$2"
    local file="$3"

    output_dir="$project/data/$construct/unique_walk"
    base_name=$(basename "$file")

    # Create the output directory if it doesn't exist
    mkdir -p "$output_dir"

    # Sort the file and remove duplicate lines
    (head -n 1 "$file" && tail -n +2 "$file" | sort | uniq) > "$output_dir/$base_name"

    echo "Processed random walk file $file and saved unique paths to $output_dir/$base_name"
}

export -f process_file

# Loop through each project directory
for project_dir in ../../project_*; do
    if [[ -d "$project_dir" ]]; then
        project=$(basename "$project_dir")

        # Find all CSV files for each construct and pass to parallel for processing
        for construct in "${list_of_constructs[@]}"; do
            input_dir="$project_dir/data/$construct/random_walk"
            
            # Use find to list all .csv files and process each in parallel
            find "$input_dir" -name "*.csv" -type f | parallel process_file "$project_dir" "$construct" {}
        done

        echo "Processed all files for all constructs in project $project"
    fi
done

echo "Done!"