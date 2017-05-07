import numpy as np
import matplotlib.pyplot as plt
#max_depth
#max_depth=[3,4,5,6,7,8,9,10,11,12,13,14,15]
#pcc=[0.859989007721,0.82234097918,0.89357623581,0.903986025255,0.9281771539,0.937747804166,0.937938862711,0.936537805451,0.951167448437,0.94532226476,0.949008681529,0.944652057635,0.943110920519]
#plt.plot(max_depth,pcc,lw=1.5)
#line1=[(0,0.951167448437),(11,0.951167448437)]
#x1=[0,16]
#y1=[0.951167448437,0.951167448437]
#x2=[11,11]
#y2=[0,0.951167448437]
#x3=[0,16]
#y3=[0.949008681529,0.949008681529]
#x4=[13,13]
#y4=[0,0.949008681529]
#plt.plot(x1,y1,marker='x',linestyle='--')
#plt.plot(x2,y2,marker='x',linestyle='--')
#plt.plot(x3,y3,marker='x',linestyle='--')
#plt.plot(x4,y4,marker='x',linestyle='--')
#plt.annotate('(11,0.9512)',xy=(11,0.951167448437),xytext=(12,0.96),arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2"))
#plt.annotate('(13,0.9490)',xy=(13,0.949008681529),xytext=(13.5,0.92),arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2"))
#plt.xlabel('max_depth')
#plt.ylabel('pcc')
#plt.xlim(1,16)
#plt.ylim(0.80,1.0)
#plt.show()


#subsample
#subsample=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
#pcc=[0.898627876071, 0.91483799622, 0.929968951062, 0.929206005879, 0.937631390935, 0.935071688341, 0.923790747245, 0.926687895459, 0.929761662382, 0.951167448437]
#plt.plot(subsample,pcc,lw=1.5)
#x1=[0,1.0]
#y1=[0.951167448437,0.951167448437]
#x2=[1,1]
#y2=[0,0.951167448437]

#plt.plot(x1,y1,marker='x',linestyle='--')
#plt.plot(x2,y2,marker='x',linestyle='--')

#plt.annotate('(1,0.9512)',xy=(1,0.951167448437),xytext=(0.8,0.96),arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2"))
#plt.xlabel('subsample')
#plt.ylabel('pcc')
#plt.xlim(0,1.1)
#plt.ylim(0.88,0.97)
#plt.show()

#colsample_bytree
#colsample_bytree=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
#pcc=[0.882177506616, 0.911294285035, 0.924326326616, 0.909219763235, 0.947336086283, 0.953925104861, 0.936582004092, 0.909219763235, 0.937539669061, 0.951167448437]
#plt.plot(colsample_bytree,pcc,lw=1.5)
#x1=[0,1.0]
#y1=[0.951167448437,0.951167448437]
#x2=[1,1]
#y2=[0,0.951167448437]
#x3=[0,1.1]
#y3=[0.953925104861,0.953925104861]
#x4=[0.6,0.6]
#y4=[0,0.953925104861]
#plt.plot(x1,y1,marker='x',linestyle='--')
#plt.plot(x2,y2,marker='x',linestyle='--')
#plt.plot(x3,y3,marker='x',linestyle='--')
#plt.plot(x4,y4,marker='x',linestyle='--')
#plt.annotate('(1,0.9512)',xy=(1,0.951167448437),xytext=(0.8,0.96),arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2"))
#plt.annotate('(0.6,0.9539)',xy=(0.6,0.953925104861),xytext=(0.4,0.96),arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2"))
#plt.xlabel('colsample_bytree')
#plt.ylabel('pcc')
#plt.xlim(0,1.1)
#plt.ylim(0.88,0.97)
#plt.show()

#max_depth  num_rounds
with open('data//depth_round.txt','r') as rf:
    num_round=[x for x in range(4,51,2)]
    for x in range(54,103,4):
        num_round.append(x)
    cunt=0
    while cunt<12:
        for i in range(3,15):
            cunt+=1
            pear=[]
            for j in range(37):
                line=rf.readline()
                if not line:
                    break
                list=line.split()
                pear.append(list[2])
            
            plt.plot(num_round,pear,lw=1.5,label="max_depth="+str(i))
    plt.legend(loc='upper left')
    x1=[-5,94]
    y1=[0.962121811227,0.962121811227]
    x2=[94,94]
    y2=[0,0.962121811227]
    plt.plot(x1,y1,marker='x',linestyle='--',color='black')
    plt.plot(x2,y2,marker='x',linestyle='--',color='black')
    plt.annotate('(94,0.9621)',xy=(94,0.9621),xytext=(92,0.97),arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2"))
    plt.xlim(-5,104)
    plt.ylim(0.80,0.98)
    plt.xlabel('num_round')
    plt.ylabel('pcc')
    plt.show()

#learn  lambda
#with open('data//learn_lambda.txt','r') as rf:
#    lambda_=[x for x in range(0,900,100)]
#    cunt=0
#    while cunt<7:
#        for i in [0.001,0.01,0.05,0.1,0.2,0.3,0.4,0.5]:
#            cunt+=1
#            pear=[]
#            for j in range(9):
#                line=rf.readline()
#                if not line:
#                    break
#                list=line.split()
#                pear.append(list[2])
            
#            plt.plot(lambda_,pear,lw=1.5,label="learning_rate="+str(i))
#    plt.legend(loc='lower left')
#    #x1=[-5,52]
#    #y1=[0.959624776977,0.959624776977]
#    #x2=[48,48]
#    #y2=[0,0.959624776977]
#    #plt.plot(x1,y1,marker='x',linestyle='--',color='black')
#    #plt.plot(x2,y2,marker='x',linestyle='--',color='black')
#    #plt.annotate('(48,0.9596)',xy=(48,0.9596),xytext=(46,0.97),arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2"))
#    plt.xlim(-200,900)
#    plt.ylim(0.80,0.98)
#    plt.xlabel('lambda')
#    plt.ylabel('pcc')
#    plt.show()

#eta  alpha
#with open('data//eta_alpha.txt','r') as rf:
#    alpha_=[0.001,0.005,0.01,0.05,0.1,0.2,0.3]
#    cunt=0
#    while cunt<7:
#        for i in [0.001,0.005,0.01,0.05,0.1,0.2,0.3]:
#            cunt+=1
#            pear=[]
#            for j in range(7):
#                line=rf.readline()
#                if not line:
#                    break
#                list=line.split()
#                pear.append(list[2])
            
#            plt.plot(alpha_,pear,lw=1.5,label="eta="+str(i))
#    plt.legend(loc='lower left')
#    #x1=[-5,52]
#    #y1=[0.959624776977,0.959624776977]
#    #x2=[48,48]
#    #y2=[0,0.959624776977]
#    #plt.plot(x1,y1,marker='x',linestyle='--',color='black')
#    #plt.plot(x2,y2,marker='x',linestyle='--',color='black')
#    #plt.annotate('(48,0.9596)',xy=(48,0.9596),xytext=(46,0.97),arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2"))
#    plt.xlim(0.001,0.3)
#    plt.ylim(0.80,0.98)
#    plt.xlabel('alpha')
#    plt.ylabel('pcc')
#    plt.show()
