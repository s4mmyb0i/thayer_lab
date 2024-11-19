#!/bin/bash

# A script to run create_network.py for all four constructs

source header.sh

# exec &> create_networks.out

echo "Lower threshold value: $threshold_value_lower"
echo "Random Walk Steps: $random_walk_steps"
echo "Random Walk Number Of Walks: $random_walk_number_of_walks"

echo "Processing the P construct"
python3 create_network.py P 827 "$threshold_value_lower"
echo "Finished processing the P construct"\n

echo "Processing the PL construct"
python3 create_network.py PL 845 "$threshold_value_lower"
echo "Finished processing the PL construct"\n

echo "Processing the AP construct"
python3 create_network.py AP 865 "$threshold_value_lower"
echo "Finished processing the AP construct"\n

echo "Processing the APL construct"
python3 create_network.py APL 865 "$threshold_value_lower"
echo "Finished processing the APL construct"\n


# #!/bin/bash

# # A script to run create_network.py for all four constructs

# source header.sh

# # exec &> create_networks.out

# echo "Lower threshold value: $threshold_value_lower"
# echo "Random Walk Steps: $random_walk_steps"
# echo "Random Walk Number Of Walks: $random_walk_number_of_walks"

# echo "Processing the P construct"
# python3 create_network.py P 83 "$threshold_value_lower"
# echo "Finished processing the P construct"\n\n

# echo "Processing the PL construct"
# python3 create_network.py PL 92 "$threshold_value_lower"
# echo "Finished processing the PL construct"\n\n

# echo "Processing the AP construct"
# python3 create_network.py AP 239 "$threshold_value_lower"
# echo "Finished processing the AP construct"\n\n

# echo "Processing the APL construct"
# python3 create_network.py APL 243 "$threshold_value_lower"
# echo "Finished processing the APL construct"\n\n