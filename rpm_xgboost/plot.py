import numpy as np
import matplotlib.pyplot as plt
from matplotlib import mlab
from xgboost import plot_importance
from matplotlib import rcParams
for i in range(10):
    with open('data//plot//'+str(i)+'.txt','r') as f:
        real_value=[];mz=[];pre_value=[]
        for j in range(4):
            line=f.readline()
            if j==0:
                peptide=line
            if j==1:
                real_value.append([float(x) for x in line.split()])
            if j==2:
                mz.append([float(x) for x in line.split()])
            if j==3:
                pre_value.append([float(x) for x in line.split()])
        plt.figure(1)
        ax1=plt.subplot(211)#在图表2中创建子图1  
        ax2=plt.subplot(212)#在图表2中创建子图2
        plt.figtext(.2,.8, peptide)
        plt.sca(ax1)
        plt.bar(left = mz[0],height = real_value[0] ,width = 0.2)
        plt.ylabel('Experimental Intensity')
        plt.ylim([0,1.0])

        plt.sca(ax2)
        left = mz[0]
        plt.bar(left,height = pre_value[0],width = 0.2)
        plt.ylim([0,1.0])
        ax2.invert_yaxis()
        ax2.xaxis.tick_top()
        plt.ylabel('Predicted Intensity')
        #ax2.set_xticklabels(column_labels, minor=False)
        

        plt.show()
        f.close()
