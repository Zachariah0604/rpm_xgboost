import pandas as pd
import numpy as np
import xgboost as xgb
from xgboost.sklearn import XGBRegressor
from sklearn import preprocessing
from sklearn.cross_validation import StratifiedKFold
from sklearn.cross_validation import KFold
from sklearn.externals import joblib
#from sklearn.model_selection import StratifiedKFold
from sklearn.grid_search import GridSearchCV
from math import sqrt
import matplotlib.pyplot as plt
import time
import random
from calculate import *
import matplotlib.pyplot as plt
def get_value(list):
    return float(list[0])
subset = 1
random_seed = 1225
n_rows = 200809
train_rows = int(n_rows * subset)
#random.seed(rand_seed)
#skip = sorted(random.sample(xrange(1,n_rows + 1),n_rows - train_rows))
print('loading data...')
data = pd.read_csv('data//data_embeding_without0.csv')
person_list = []
merge_list = []
print('DataShape: ' + str(data.shape))
label = data['Intensity'].values
peptide=data['Peptide'].values
print('get spectrum intensity list...')
merge_list.append(get_merge_list(data))
del data['Intensity']
del data['Number']
del data['Peptide']
predictors=[x for x in data]
X = data.values
params = {}

params['eta']=0.01
params['alpha']=0.01
params['silent'] = 0
params['learning_rate'] = 0.4
params['n_estimators'] = 1000
params['eval_metric']='rmse'
params['max_depth'] = 11 
params['min_child_weight'] = 0
params['gamma'] = 0
params['lambda']=500
params['subsample'] = 1
params['colsample_bytree'] = 0.6
params['objective'] = 'reg:linear'
#params['early_stopping_rounds']=100
params['seed'] = -1
params['updater'] = 'grow_gpu'
plst = list(params.items())
print('K-Folds cross validation...')
#cv = StratifiedKFold(label,10)
cv=KFold(len(merge_list[0]),10,shuffle=True,random_state=1017)
print(cv)
dtrain_predictions = [];idx=[];
print('training model ...')
def get_split_list(array_list):
    list=[]
    for m in array_list:
        for n in range(len(merge_list[0][m])):
            list.append(merge_list[0][m][n][0]-1)
    return list
for i,(train_peptide,test_piptide) in enumerate(cv):
    train=get_split_list(train_peptide)
    test=get_split_list(test_piptide) 
    xgtrain = xgb.DMatrix(X[np.array(train)],label=label[np.array(train)])
    xgtest=xgb.DMatrix(X[np.array(test)])
    idx.append([x+1 for x in list(test)])
    tmp = time.time()
    bst = xgb.train(plst,xgtrain,num_boost_round=10)
    #bst.save_model('model2//'+str(i)+'.model')
    ##dump model
    #bst.dump_model('model2//'+str(i)+'.dump.raw.txt')
    ## dump model with feature map
    #bst.dump_model('model2//'+str(i)+'.dump.nice.txt','data//featmap.txt')
   
    boost_time = time.time() - tmp
    rmse = bst.eval(xgb.DMatrix(X[test],label=label[test]))
    
    print('Fold {}:{},Boost Time {}'.format(i+1,rmse,str(boost_time)))
    dtrain_predictions.append(bst.predict(xgtest))
    del bst
print('training complete...')
print('write predict data in file')
pred_result=[];k=0
for i in range(len(idx)):
    for j in range(len(idx[i])):
        k+=1
        pred_result.append([idx[i][j],dtrain_predictions[i][j]])
pred_result.sort(key=get_value)
pred = pd.DataFrame({"Number":[pred_result[i][0] for i in range(k)],"Peptide":peptide,"Intensity":[pred_result[i][1] for i in range(k)]})
pred.to_csv('data//pred2.csv')
print('get predict spectrum intensity list...')
merge_list.append(get_merge_list(pred))      
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


