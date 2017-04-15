from xgboost import plot_importance
from matplotlib import pyplot
from xgboost import XGBRegressor
import pandas as pd
from numpy import loadtxt

dataset=loadtxt('data//data.csv',delimiter=",")
X=dataset[:,0:34]
y=dataset[:,34]

model =XGBRegressor()
model.fit(X,y)

plot_importance(model)
pyplot.show()