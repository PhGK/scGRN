
import pandas as pd
import os

filenames = os.listdir('./results/highest_interactions/')


data = pd.concat([pd.read_csv('./results/highest_interactions/' + file) for file in filenames])
print(data.shape)
data.to_csv('.results/high_values_concat.csv')





