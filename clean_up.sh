#!/bin/bash

# A script to delete all the files in specified subdirectories within project directories. Select which subfolder to delete using the SUBFOLDER_NAME variable.

# Store the subfolder name as a variable
SUBFOLDER_NAME="random_walk"
# SUBFOLDER_NAME="."

SUBDIRECTORY_TYPE="data"
# SUBDIRECTORY_TYPE="analysis"

# Navigate to the base directory where your project directories are located
BASE_DIR="../"

# Loop through each project directory
for PROJECT_DIR in "$BASE_DIR"project_*; do
    # Check if it's a directory
    if [[ -d $PROJECT_DIR ]]; then
        # Loop through each construct directory
        for CONSTRUCT_DIR in "$PROJECT_DIR"/"$SUBDIRECTORY_TYPE"/*; do
            # Check if it's a directory
            if [[ -d $CONSTRUCT_DIR ]]; then
                # Define the subfolder directory path
                TARGET_SUBFOLDER="$CONSTRUCT_DIR/$SUBFOLDER_NAME"
                
                # Check if the target subfolder exists
                if [[ -d $TARGET_SUBFOLDER ]]; then
                    # Remove all files in the target subfolder
                    echo "Deleting files in $TARGET_SUBFOLDER..."
                    rm -f "$TARGET_SUBFOLDER"/*
                fi
            fi
        done
    fi
done

echo "Cleanup completed."