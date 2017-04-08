from xgboost.sklearn import XGBRegressor
from sklearn.grid_search import GridSearchCV
import pandas as pd
import os
if __name__ =='__main__':
    orginal_data=pd.read_csv('data/data_embeding.csv')
    train=orginal_data
    X=train.drop(['Number','Peptide','Intensity'],axis=1)
    y=train['Intensity']
    xgb_model=XGBRegressor(
        silent=0,
      #  updater='grow_gpu'
        )
    clf=GridSearchCV(xgb_model,
                     {
                     'max_depth':[2,4,5,6,8],
                     'n_estimators':[50,100,200,500,1000],
                     'min_child_weight':1,
                      'gamma':[i/10.0 for i in range(0,5)],
                      'subsample':[0.6,0.8,1],
                      'colsample_bytree':[0.8,1],
                      'reg_alpha':[0, 0.001, 0.005, 0.01, 0.05],
                      'scale_pos_weight':[0.5,0.7,0.8,1],
                     },
                     verbose=1
        )
    clf.fit(X,y)
    print clf.best_score_,clf.best_params_
    f=open('123.txt','w')
    f.write(clf.best_score_+clf.best_params_)
    f.close()
    os.system('pause')