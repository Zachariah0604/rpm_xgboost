import pandas as pd
import numpy as np
import xgboost as xgb
from xgboost.sklearn import XGBRegressor
from sklearn import preprocessing
from sklearn.cross_validation import StratifiedKFold
from sklearn.cross_validation import KFold
from sklearn.externals import joblib
from sklearn.grid_search import GridSearchCV
from math import sqrt
import matplotlib.pyplot as plt
import time
import random
from calculate_libsvm import *
import matplotlib.pyplot as plt
from sklearn.datasets import load_svmlight_file
def get_value(list):
    return float(list[0])
subset = 1
random_seed = 1225
n_rows = 200809
train_rows = int(n_rows * subset)
#random.seed(rand_seed)
#skip = sorted(random.sample(xrange(1,n_rows + 1),n_rows - train_rows))
print('loading data...')
data,label=load_svmlight_file('data//data_embeding.libsvm')
person_list = []
merge_list = []
print('DataShape: ' + str(data.shape))
peptide=data[:,1]
print('get spectrum intensity list...')
merge_list.append(get_merge_list(data,label))
data=data[:,2:].toarray()
params = {}
params['eta']=0.01
params['alpha']=0.005
params['silent'] = 1
params['learning_rate'] = 0.4
params['n_estimators'] = 1000
params['eval_metric']='rmse'
params['max_depth'] = 8
params['min_child_weight'] = 0
params['gamma'] = 0
params['lambda']=500
params['subsample'] = 1
params['colsample_bytree'] = 0.6
params['objective'] = 'reg:linear'
#params['early_stopping_rounds']=100
params['seed'] = -1
#params['updater'] = 'grow_gpu'
plst = list(params.items())
print('K-Folds cross validation...')
cv=KFold(len(merge_list[0]),10,shuffle=True,random_state=1017)
#0.9400
print(cv)
dtrain_predictions = [];idx=[];
print('training model ...')
def get_split_list(array_list):
    list=[]
    for m in array_list:
        for n in range(len(merge_list[0][m])):
            list.append((merge_list[0][m][n][0]-1).astype(int))
    return list
for i,(train_peptide,test_piptide) in enumerate(cv):
    train=get_split_list(train_peptide)
    test=get_split_list(test_piptide) 
    xgtrain = xgb.DMatrix(data[np.array(train)],label=label[np.array(train)])
    xgtest=xgb.DMatrix(data[np.array(test)])
    idx.append([(x+1).astype(int) for x in list(test)])
    tmp = time.time()
    bst = xgb.train(plst,xgtrain,num_boost_round=48)
    boost_time = time.time() - tmp
    rmse = bst.eval(xgb.DMatrix(data[test],label=label[test]))
    
    print('Fold {}:{},Boost Time {}'.format(i+1,rmse,str(boost_time)))
    dtrain_predictions.append(bst.predict(xgtest))
    del bst
print('training complete...')
print('get predict spectrum intensity list...')
pred_result=[];k=0
for i in range(len(idx)):
    for j in range(len(idx[i])):
        k+=1
        pep_id=idx[i][j]
        pred_result.append([pep_id,peptide[pep_id-1,0],dtrain_predictions[i][j]])
pred_result.sort(key=get_value)
pred_label=[pred_result[i][2] for i in range(k)]
merge_list.append(get_merge_list(pred_result,pred_label))      
print('calculate person coefficient..')
sum_person = 0.0
def get_pearson_x(len_of_peptide,m):
    x=[]
    for j in range(len(merge_list[0][len_of_peptide])):
        x.append(merge_list[m][len_of_peptide][j][1])
    return x
for i in range(len(merge_list[0])):
    person_list.append(pearson_r(get_pearson_x(len(merge_list[0][i]),0),get_pearson_x(len(merge_list[1][i]),1)))
for i in range(len(person_list)):
    sum_person+=person_list[i]
person_mean = sum_person / float(len(person_list))
print('r= ' + str(person_mean))


