import pandas as pd
import numpy as np
import xgboost as xgb
from xgboost.sklearn import XGBRegressor
from sklearn import preprocessing
from sklearn.cross_validation import train_test_split
from sklearn.grid_search import GridSearchCV
from math import sqrt
import matplotlib.pyplot as plt

def get_value(list):
    return float(list[2])

def get_intensity_list(data):
    temp_peptide=''
    temp_list=[]
    intensity_list=[]
    for row in data.itertuples():
        peptide=row.Peptide
        if peptide != temp_peptide:
            if temp_peptide == '':
                temp_list.append(row.Intensity)
            else:
                intensity_list.append(temp_list)
                temp_list=[]
                temp_list.append(row.Intensity)
            temp_peptide=peptide
        else:
            temp_list.append(row.Intensity)
    intensity_list.append(temp_list)
    return intensity_list
def get_peptide_list(data):
    temp_peptide=''
    temp_list=[]
    peptide_list=[]
    for row in data.itertuples():
        peptide=row.Peptide
        if peptide != temp_peptide:
            if temp_peptide == '':
                temp_list.append(row.Peptide)
            else:
                peptide_list.append(temp_list)
                temp_list=[]
                temp_list.append(row.Peptide)
            temp_peptide=peptide
       # else:
        #    temp_list.append(row.Peptide)
    peptide_list.append(temp_list)
    return peptide_list
def get_y_Mass_list(data):
    temp_peptide=''
    temp_list=[]
    y_Mass_list=[]
    for row in data.itertuples():
        peptide=row.Peptide
        if peptide != temp_peptide:
            if temp_peptide == '':
                temp_list.append(row.y_Mass)
            else:
                y_Mass_list.append(temp_list)
                temp_list=[]
                temp_list.append(row.y_Mass)
            temp_peptide=peptide
        else:
            temp_list.append(row.y_Mass)
    y_Mass_list.append(temp_list)
    return y_Mass_list
def multipl(a,b):
    multipl_of_a_b=0.0
    for i in range(len(a)):
        temp=a[i]*b[i]
        multipl_of_a_b+=temp
    return multipl_of_a_b
# https://wikimedia.org/api/rest_v1/media/math/render/svg/832f0c5c22a0d6f2596c150de811247438a503de
def pearson_r(x,y):
    n=len(x)
    sum_x=sum(x)
    sum_y=sum(y)
    sum_of_xy=multipl(x,y)
    sum_of_x2 = sum([pow(i,2) for i in x])
    sum_of_y2 = sum([pow(j,2) for j in y])
    numerator=sum_of_xy-(float(sum_x)*float(sum_y)/n)
    denominator=sqrt((sum_of_x2-float(sum_x**2)/n)*(sum_of_y2-float(sum_y**2)/n))
    return numerator/denominator

person_list=[];split_list=[]
print('loading data...')
orginal_data=pd.read_csv('data/data_embeding-2.csv')
print('DataShape: '+ str(orginal_data.shape))
train=orginal_data.head(16108)
print('Train dataset: '+str(train.Number.values))
test=orginal_data.loc[16108:orginal_data.shape[0]]
print('Test dataset: '+str(test.Number.values))
print('split test data...')
split_list.append(get_intensity_list(test))
split_list.append(get_peptide_list(test))
split_list.append(get_y_Mass_list(test))
idx=test.Number.values
peptide=test.Peptide.values
train=train.drop(['Number','Peptide'],axis=1)
test=test.drop(['Number','Peptide'],axis=1)

features=list(train.columns[0:])
label=train.Intensity.values.astype(float)

def model_fit(alg,dtrain,dtest,predictors,useTrainCV=True,cv_folds=10,early_stoping_rounds=50):
    if useTrainCV:
        xgb_param=alg.get_xgb_params()
        xgtrain=xgb.DMatrix(dtrain[predictors].values,label=dtrain['Intensity'].values)
        xgtest=xgb.DMatrix(dtest[predictors].values)
        #xgb_param['updater']='grow_gpu'
        cv_result=xgb.cv(xgb_param,xgtrain,num_boost_round=alg.get_params()['n_estimators'],nfold=cv_folds,early_stopping_rounds=early_stoping_rounds, show_stdv=True)
        alg.set_params(n_estimators=cv_result.shape[0])
        #alg.set_params(updater='grow_gpu')
        alg.fit(dtrain[predictors], dtrain['Intensity'],eval_metric='auc')
        dtrain_predictions = alg.predict(dtest[predictors])
        #dtrain_predprob = alg.predict_proba(dtrain[predictors])[:,1]
        print('training complete...')
        print('write predict data in file')

        pred=pd.DataFrame({"id":idx,"Peptide":peptide,"Intensity":dtrain_predictions})
        pred.to_csv('data//pred.csv')
        print('split predict data...')
        split_list.append(get_intensity_list(pred))
        print('calculate person coefficient..')
        for i in range(len(split_list[0])):
            person_list_temp=[]
            person_list_temp.append(split_list[1][i][0])
            person_list_temp.append(split_list[2][i][0])
            person_list_temp.append(pearson_r(split_list[0][i],split_list[3][i]))
            person_list.append(person_list_temp)
        person_list.sort(key=get_value,reverse = True)
        plot_list=[];plost_list_temp=[]
        real_value=[];predictd_value=[];m_dv_z=[];pep=[]
        for i in range(10):
            for j in range(len(split_list[0])):
                if person_list[i][0] == split_list[1][j][0]:
                    real_value.append(split_list[0][j])
                    pep.append(split_list[1][j][0])
                    m_dv_z.append(split_list[2][j])
                    predictd_value.append(split_list[3][j])
        for i in range(10):
            with open('data//plot//'+str(i)+'.txt','w') as f:
                f.write(pep[i]+'\n')
                for j in range(len(real_value[i])):
                    f.write(str(real_value[i][j])+'\t')
                f.write('\n')
                for j in range(len(m_dv_z[i])):
                    f.write(str(m_dv_z[i][j])+'\t')
                f.write('\n')
                for j in range(len(predictd_value[i])):
                    f.write(str(predictd_value[i][j])+'\t')
                f.write('\n')

        #real_value_array=np.array(real_value) 
        #predictd_value_array=np.array(predictd_value)
        #real=[];predict=[]
        #for i in range(10):
        #    real.append(np.mean(real_value_array[i]))
        #    predict.append(np.mean(predictd_value_array[i]))
        #for i in range(10):
        #    plost_list_temp=[]
        #    plost_list_temp.append(real[i])
        #    plost_list_temp.append(pep[i])
        #    plost_list_temp.append(m_dv_z[i])
        #    plost_list_temp.append(predict[i])
        #    plot_list.append(plost_list_temp)
        #plot_list.sort(key=get_value)
        
        sum=0.0
        for i in range(len(person_list)):
            sum+=person_list[i][2]
        person_mean=sum/float(len(person_list))
        print('r= '+str(person_mean))
        #plot(plot_list)
#train=np.array(train)
#test=np.array(test)
predictors=[x for x in train.columns if x not in ['Intensity']]
xgb_model=XGBRegressor(
       
        silent = 0,
        learning_rate=0.1,
        n_estimators=1000,
        max_depth=8,
        min_child_weight=0,
        gamma=0,
        subsample=0.8,
        colsample_bytree=0.8,
        objective='reg:linear',
        scale_pos_weight=1,
        seed=27,
       # updater='grow_gpu'
       )

print('training model ...')
model_fit(xgb_model,train,test,predictors)
