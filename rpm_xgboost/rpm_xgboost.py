import pandas as pd
import numpy as np
import xgboost as xgb
from xgboost.sklearn import XGBRegressor
from sklearn import preprocessing
from sklearn.cross_validation import train_test_split
from sklearn.grid_search import GridSearchCV
from math import sqrt

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
print 'loading data...'
orginal_data=pd.read_csv('data/data_embeding.csv')
print 'DataShape: '+ str(orginal_data.shape)
train=orginal_data.head(105000)
print 'Train dataset: '+str(train.Number.values)
test=orginal_data.loc[105000:orginal_data.shape[0]]
print 'Test dataset: '+str(test.Number.values)
print 'split test data...'
split_list.append(get_intensity_list(test))
idx=test.Number.values
peptide=test.Peptide.values
train=train.drop(['Number','Peptide'],axis=1)
test=test.drop(['Number','Peptide'],axis=1)

features=list(train.columns[0:])
label=train.Intensity.values.astype(float)

def model_fit(alg,dtrain,dtest,predictors,useTrainCV=True,cv_folds=5,early_stoping_rounds=50):
    if useTrainCV:
        xgb_param=alg.get_xgb_params()
        xgtrain=xgb.DMatrix(dtrain[predictors].values,label=dtrain['Intensity'].values)
        xgtest=xgb.DMatrix(dtest[predictors].values)
        cv_result=xgb.cv(xgb_param,xgtrain,num_boost_round=alg.get_params()['n_estimators'],nfold=cv_folds,early_stopping_rounds=early_stoping_rounds, show_stdv=False)
        alg.set_params(n_estimators=cv_result.shape[0])
        alg.fit(dtrain[predictors], dtrain['Intensity'],eval_metric='auc')
        dtrain_predictions = alg.predict(dtest[predictors])
        #dtrain_predprob = alg.predict_proba(dtrain[predictors])[:,1]
        print 'training complete...'
        print 'write predict data in file'
        len1=len(idx) 
        len2=len(peptide)
        len3=len(dtrain_predictions)
        pred=pd.DataFrame({"id":idx,"Peptide":peptide,"Intensity":dtrain_predictions})
        pred.to_csv('data//pred.csv')
        print 'split predict data...'
        split_list.append(get_intensity_list(pred))
        print 'calculate person coefficient..'
        for i in range(len(split_list[0])):
            person_list.append(corrcoef(split_list[0][i],split_list[1][i]))
        person_array=np.array(person_list)
        print 'r= '+str(np.mean(person_array)) 
#train=np.array(train)
#test=np.array(test)
predictors=[x for x in train.columns if x not in ['Intensity']]
xgb_model=XGBRegressor(
        silent = 0,
        learning_rate=0.1,
        n_estimators=1000,
        max_depth=5,
        min_child_weight=1,
        gamma=0,
        subsample=0.8,
        colsample_bytree=0.8,
        objective='reg:linear',
        nthread=4,
        scale_pos_weight=1,
        seed=27
       )

print 'training model ...'
model_fit(xgb_model,train,test,predictors)
