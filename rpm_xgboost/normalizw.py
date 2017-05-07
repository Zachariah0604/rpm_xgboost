
import csv
import pandas as pd

def max_min(X):
    max_x = max(X)
    min_x = min(X)
    return [(x - min_x)/float(max_x - min_x) for x in X]


reader=pd.read_csv('data//data_without_embeding_35.csv')
titles=reader.columns
for title in titles:
    column=reader[str(title)]
    max_x=max(column)
    min_x=min(column)
    for i in range(len(column)):
        x=reader.loc[i,str(title)]
        reader.set_value(i,str(title),(x-min_x)/float(max_x-min_x))
        print('alter ('+ str(i) +',' + str(title) +')')
reader.to_csv('data//nomalize_data_without_embeding_35.csv')