#!/bin/bash

# A script to run everything at once

# echo "Running create_network_for_all.sh..."
# ./create_network_for_all.sh

echo "Running random_walk_for_all.sh..."
./random_walk_for_all.sh

echo "Running remove_duplicate_paths_for_all.sh..."
./only_unique.sh

# cd analysis

# echo "Running how_end_for_all.sh..."
# ./how_end_for_all.sh

# ./mean_median_for_all.sh

# echo "Running res_frequency_for_all.sh..."
# ./res_frequency_for_all.sh

# echo "Running which_captured_for_all.sh..."
# ./which_captured_for_all.sh