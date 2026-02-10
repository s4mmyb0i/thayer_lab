#!/bin/bash

# Runs all the absorption probabiltiiy calculations

python3 -m linalg.absorption_probs.absorption_probs \
  /Users/sam/Documents/Wesleyan/Thayer_Lab/2025/summer/bharat/matrices/trans_abs_matrix/P_trans_abs.npz \
  /Users/sam/Documents/Wesleyan/Thayer_Lab/2025/summer/bharat/linalg/absorption_probs

python3 -m linalg.absorption_probs.absorption_probs \
  /Users/sam/Documents/Wesleyan/Thayer_Lab/2025/summer/bharat/matrices/trans_abs_matrix/PL_trans_abs.npz \
  /Users/sam/Documents/Wesleyan/Thayer_Lab/2025/summer/bharat/linalg/absorption_probs

python3 -m linalg.absorption_probs.absorption_probs \
  /Users/sam/Documents/Wesleyan/Thayer_Lab/2025/summer/bharat/matrices/trans_abs_matrix/AP_trans_abs.npz \
  /Users/sam/Documents/Wesleyan/Thayer_Lab/2025/summer/bharat/linalg/absorption_probs

python3 -m linalg.absorption_probs.absorption_probs \
  /Users/sam/Documents/Wesleyan/Thayer_Lab/2025/summer/bharat/matrices/trans_abs_matrix/APL_trans_abs.npz \
  /Users/sam/Documents/Wesleyan/Thayer_Lab/2025/summer/bharat/linalg/absorption_probs#!/bin/bash