
import csv
import pandas as pd

def max_min(X):
    max_x = max(X)
    min_x = min(X)
    return [(x - min_x)/float(max_x - min_x) for x in X]


reader=pd.read_csv('data//data_embeding.csv')
titles=reader.columns[2:256]
fftile=[]
for ftile  in range(0,60):
    fftile.append(titles[ftile])

for ftile  in range(223,len(titles)):
    fftile.append(titles[ftile])
print(fftile)
for title in fftile:
    column=reader[str(title)]
    max_x=max(column)
    min_x=min(column)
    devalue=float(max_x-min_x)
    reader[str(title)]=(reader[str(title)]-min_x)/devalue
    #for i in range(len(column)):
     #   x=reader.loc[i,str(title)]
     #   reader.set_value(i,str(title),(float(x)-min_x)/float(max_x-min_x))
    print('alter ('+ str(title) +')')
reader.to_csv('data//nomalize_data_embeding.csv',index=False)