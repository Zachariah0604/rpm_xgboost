from sklearn.datasets import load_svmlight_file
from calculate_libsvm import *
from svmutil import *
dlist=[1, 0.0, 0.22575167]
print(type(dlist))
if type(dlist) is list:
    print(1)
else:
    print(0)
datas,label=load_svmlight_file('D://data_without_embeding_35_test.libsvm')
print(datas[:,1])
#featureList.sort(key=lambda x:int(x.split(':')[0]))
merge_list=[]
merge_list.append(get_merge_list(datas,label))
i=0