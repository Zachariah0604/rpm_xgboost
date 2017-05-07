import pandas as pd
import csv
import numpy as np
reader=pd.read_csv('data//data_without_embeding_35_test.csv')
DataSt=reader.values.astype(float)
mean_line=DataSt.mean(axis=0)
DataSt=DataSt-mean_line
DataSt_cov=np.cov(DataSt,rowvar=0)
i=0