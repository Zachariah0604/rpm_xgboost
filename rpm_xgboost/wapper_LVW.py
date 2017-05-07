import pandas as pd
import random
from sklearn.cross_validation import KFold
from calculate import *
from sklearn.tree import DecisionTreeRegressor
import numpy as np
print('get data...')
data=pd.read_csv('data//data_without_embeding_35.csv')
label=data['Intensity'].values
Tscore=0.0
t=0
T=100

merge_list=[]
merge_list.append(get_merge_list(data))
del data['Number']
del data['Peptide']
del data['Intensity']
predictor=[x for x in data.columns if x not in ['Intensity']]
d=len(predictor)
print('score=0.0,T=10,d='+str(d))
def get_split_list(array_list):
    list=[]
    for m in array_list:
        for n in range(len(merge_list[0][m])):
            list.append(merge_list[0][m][n][0]-1)
    return list
cunt=0
while t<T:
    predict_tmp=random.sample(predictor,random.randint(2,34))
    d_tmp=len(predict_tmp)
    print('\n\ntraining model ...#'+str(cunt)+'_'+str(t+1))
    cv=KFold(len(merge_list[0]),10,shuffle=True,random_state=random.randint(1,2000))
    print(cv)
    score_tmp=[]
    for i,(train_peptide,test_peptide) in enumerate(cv):
        train=get_split_list(train_peptide)
        test=get_split_list(test_peptide) 
        model=DecisionTreeRegressor()
        data_tmp=data[predict_tmp]
        DataSt=data_tmp.values
        model.fit(DataSt[train],label[train])
        score_tmp.append(model.score(DataSt[test],label[test]))        
        del model
    print(score_tmp)
    get_score=np.mean(np.array(score_tmp))
    print('score='+str(get_score))
    if get_score>Tscore or (get_score==Tscore and d_tmp<d):
        cunt+=1
        t=0
        Tscore=get_score
        d=d_tmp
        print([x for x in predict_tmp])
        print('d='+str(d)+'\n################################')
    else:
        t+=1