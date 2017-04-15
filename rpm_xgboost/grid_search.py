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
     max_depth=8,
       n_estimators=1000,
      min_child_weight=0,
      gamma=0,
      subsample=0.8,
      colsample_bytree=0.8
        )
    clf=GridSearchCV(xgb_model,
                     {
                     
                    
                    
                      
                      
                      'reg_alpha':[0, 0.001, 0.005, 0.01, 0.05],
                      'scale_pos_weight':[0.5,0.7,0.8,1],
                     },
                     verbose=1
        )
    clf.fit(X,y)
    print clf.best_score_,clf.best_params_
    os.system('pause')
    f=open('123.txt','w')
    f.write(clf.best_score_+clf.best_params_)
    f.close()
    os.system('pause')