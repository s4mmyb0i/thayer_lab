import pandas as pd
from ydata_profiling import ProfileReport

df = pd.read_csv("/Users/sam/Documents/Wesleyan/Thayer_Lab/2025/summer/bharat/random_walk/AP_walks.csv", delimiter=';')

profile = ProfileReport(df, title="Random Walk Profile", explorative=True)
profile.to_file("/Users/sam/Documents/Wesleyan/Thayer_Lab/2025/summer/bharat/random_walk/AP_walks_profile.html")

import os
print("Saving profile to:", os.getcwd())
