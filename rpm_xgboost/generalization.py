from calculate import *
import pandas as pd
import xgboost as xgb
from sklearn.externals import joblib
test=pd.read_csv('data_embeding_lessthan0.8.csv')
print('DataShape:'+str(test.shape))
merge_list = []
merge_list.append(get_merge_list(test))
peptide=test.Peptide.values
idx=test.Number.values
test=test.drop(['Number','Peptide','Intensity'],axis=1)
predict_list=[];person_list = []
xgbtest=xgb.DMatrix(test)
xgb_model=xgb.Booster(model_file='model//3.model')
predict_list=xgb_model.predict(xgbtest)
pred=pd.DataFrame({'Number':idx,'Peptide':peptide,'Intensity':predict_list})
merge_list.append(get_merge_list(pred))
pred.to_csv('pred.csv')
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
