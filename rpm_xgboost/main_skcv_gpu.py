import pandas as pd
import numpy as np
import xgboost as xgb
from xgboost.sklearn import XGBRegressor
from sklearn import preprocessing
from sklearn.cross_validation import StratifiedKFold
from sklearn.cross_validation import KFold
#from sklearn.model_selection import StratifiedKFold
from sklearn.grid_search import GridSearchCV
from math import sqrt
import matplotlib.pyplot as plt
import time
import random

def get_value(list):
    return float(list[0])

def get_intensity_list(data):
    temp_peptide = ''
    temp_list = []
    intensity_list = []
    for row in data.itertuples():
        peptide = row.Peptide
        if peptide != temp_peptide:
            if temp_peptide == '':
                temp_list.append(row.Intensity)
            else:
                intensity_list.append(temp_list)
                temp_list = []
                temp_list.append(row.Intensity)
            temp_peptide = peptide
        else:
            temp_list.append(row.Intensity)
    intensity_list.append(temp_list)
    return intensity_list
def get_peptide_list(data):
    temp_peptide = ''
    temp_list = []
    peptide_list = []
    for row in data.itertuples():
        peptide = row.Peptide
        if peptide != temp_peptide:
            if temp_peptide == '':
                temp_list.append(row.Peptide)
            else:
                peptide_list.append(temp_list)
                temp_list = []
                temp_list.append(row.Peptide)
            temp_peptide = peptide
       # else:
        #    temp_list.append(row.Peptide)
    peptide_list.append(temp_list)
    return peptide_list
def get_y_Mass_list(data):
    temp_peptide = ''
    temp_list = []
    y_Mass_list = []
    for row in data.itertuples():
        peptide = row.Peptide
        if peptide != temp_peptide:
            if temp_peptide == '':
                temp_list.append(row.y_Mass)
            else:
                y_Mass_list.append(temp_list)
                temp_list = []
                temp_list.append(row.y_Mass)
            temp_peptide = peptide
        else:
            temp_list.append(row.y_Mass)
    y_Mass_list.append(temp_list)
    return y_Mass_list
def multipl(a,b):
    multipl_of_a_b = 0.0
    for i in range(len(a)):
        temp = a[i] * b[i]
        multipl_of_a_b+=temp
    return multipl_of_a_b
# https://wikimedia.org/api/rest_v1/media/math/render/svg/832f0c5c22a0d6f2596c150de811247438a503de
def pearson_r(x,y):
    n = len(x)
    sum_x = sum(x)
    sum_y = sum(y)
    sum_of_xy = multipl(x,y)
    sum_of_x2 = sum([pow(i,2) for i in x])
    sum_of_y2 = sum([pow(j,2) for j in y])
    numerator = sum_of_xy - (float(sum_x) * float(sum_y) / n)
    denominator = sqrt((sum_of_x2 - float(sum_x ** 2) / n) * (sum_of_y2 - float(sum_y ** 2) / n))
    return numerator / denominator

subset = 1
rand_seed = 9
n_rows = 200809
train_rows = int(n_rows * subset)
random.seed(rand_seed)
#skip = sorted(random.sample(xrange(1,n_rows + 1),n_rows - train_rows))
print 'loading data...'
data = pd.read_csv('data/data_embeding.csv')
person_list = []
split_list = []
print 'DataShape: ' + str(data.shape)
label = data['Intensity'].values
peptide=data['Peptide'].values
print 'get spectrum intensity list...'
split_list.append(get_intensity_list(data))
del data['Intensity']
del data['Number']
del data['Peptide']
predictors=[x for x in data]
X = data.values

params = {}
params['silent'] = 0
#param['learning_rate'] = 0.1
params['n_estimators'] = 1000
params['eval_metric']='rmse'
params['max_depth'] = 11
params['min_child_weight'] = 20
params['gamma'] = 0
params['subsample'] = 1
params['colsample_bytree'] = 1
params['objective'] = 'reg:linear'
#params['scale_pos_weight'] = 2
params['seed'] = -1
params['updater'] = 'grow_gpu'
plst = list(params.items())
num_round = 10
print ' K-Folds cross validation...'
#cv = StratifiedKFold(label,10)
cv=KFold(n_rows,10,random_state=0)
print cv
dtrain_predictions = [];idx=[];
print 'training model ...'
i=0
f=open('data//123.txt','w')
f.write('fold\tnum_round\trmse\n')
for train,test in cv:
    xgtrain = xgb.DMatrix(X[train],label=label[train])
    xgtest=xgb.DMatrix(X[test])
    idx.append([x+1 for x in list(test)])
    for round in [x*10 for x in range(1,81)]:
        tmp = time.time()
       # watchlist  = [(xgtest,'test'), (xgtrain,'train')]
        bst = xgb.train(plst,xgtrain,round)
        boost_time = time.time() - tmp
        rmse = bst.eval(xgb.DMatrix(X[test],label=label[test]))
        print 'Fold {}:{},num_round:{},Boost Time {}'.format(i+1,rmse,round,str(boost_time))
        writeStr=str(i)+'\t'+str(round)+'\t:'+str(rmse)+'\n'
        f.write(writeStr)
        del bst
    dtrain_predictions.append(bst.predict(xgtest))
    i+=1
    if i==10:
        break
f.close()
print 'training complete...'
print 'write predict data in file'
pred_result=[];k=0
for i in range(len(idx)):
    for j in range(len(idx[i])):
        k+=1
        pred_result.append([idx[i][j],dtrain_predictions[i][j]])
pred_result.sort(key=get_value)
pred = pd.DataFrame({"idx":[pred_result[i][0] for i in range(k)],"Peptide":peptide,"Intensity":[pred_result[i][1] for i in range(k)]})
pred.to_csv('data//pred2.csv')
print 'get predict spectrum intensity list...'
split_list.append(get_intensity_list(pred))      
print 'calculate person coefficient..'
sum_person = 0.0
cunt2=0
for i in range(len(split_list[0])):
    temp=pearson_r(split_list[0][i],split_list[1][i])
    if np.isnan(temp) or np.isinf(temp):
        cunt2+=1
    person_list.append(pearson_r(split_list[0][i],split_list[1][i]))
for i in range(len(person_list)):
    sum_person+=person_list[i]
person_mean = sum_person / float(len(person_list))
print 'r= ' + str(person_mean)


