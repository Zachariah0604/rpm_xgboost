#!usr/bin/env Python
#coding:utf-8

from numpy import *
import numpy
#import cPickle as Pickle
import random
import csv
"""
feature choose,a 25d vector indicate the number of ACDEFGHIKLMNPQRSTVWYcnm(DN)(AN),
23d vector indicate the animo acid,a window of cleavage 4,4*23+2,2 indicate N-term modification,
4d vector indicate the charge of peptide sequence.
"""

dicB={'A':206.4,'C':206.2,'D':208.6,'E':215.6,'F':212.1,\
      'G':202.7,'H':223.7,'I':210.8,'K':221.8,'L':209.6,\
      'M':213.3,'N':212.8,'P':214.4,'Q':214.2,'R':237.0,\
      'S':207.6,'T':211.7,'V':208.7,'W':216.1,'Y':213.1}
dicM={'A':71.03711,'C':103.00919,'D':115.02694,'E':129.04259,'F':147.06841,\
      'G':57.02146,'H':137.05891,'I':113.08406,'K':128.09496,'L':113.08406,\
      'M':131.04048,'N':114.04293,'P':97.05276,'Q':128.05858,'R':156.10111,\
      'S':87.03203,'T':101.04768,'V':99.06841,'W':186.07931,'Y':163.06332}
dicS={'A':0,'C':1,'D':2,'E':3,'F':4,'G':5,'H':6,'I':7,\
      'K':8,'L':9,'M':10,'N':11,'P':12,'Q':13,'R':14,\
      'S':15,'T':16,'V':17,'W':18,'Y':19}
dicHe={'A':1.24,'C':0.79,'D':0.89,'E':0.85,'F':1.26,'G':1.15,'H':0.97,\
       'I':1.29,'K':0.88,'L':1.28,'M':1.22,'N':0.94,'P':0.57,'Q':0.96,\
       'R':0.95,'S':1.00,'T':1.09,'V':1.27,'W':1.07,'Y':1.11}
with open('data\\lx-data_peptide_pif_0.8_y.csv','w') as f1:
    writer=csv.writer(f1)
    CountOfAAInPeptide=''
    IsolationWidth=''
    C_Identity=''
    N_Identity=''
    for i in range(1,21):
        CountOfAAInPeptide+='CountOfAAInPeptide_'+str(i)+','
        C_Identity+='C_Identity_'+str(i)+','
        N_Identity+='N_Identity_'+str(i)+','
    for i in range(6):
        for j in range(1,21):
            IsolationWidth+='IsolationWidth_'+str(i+1)+"_"+str(j)+','
    writeStr='Number,Peptide,'+CountOfAAInPeptide+IsolationWidth+C_Identity+N_Identity+'A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,McRatioOfPeptide,Intensity'
    writer.writerow(writeStr.split(','))
    cunt=0
    f = open('data//ion2_pif_0.8_y_new.txt')
    data = [];labels = [];pep_mz = 0;pep = '1';pep_b = 0;pep_he = 0;k=0
    for line in f:
        attribute = zeros(193)
        f_line = line.split('\t')
        peptide = list(f_line[0])
        ion = list(f_line[2])
        i = 0
        for key in dicS.keys():
            if key in peptide:
                attribute[dicS[key]] = peptide.count(key)/len(peptide)

        ion_format = f_line[5]
        if int(ion_format[1]) > 3 and int(ion_format[1]) < len(peptide) - 3 and ion_format[0] == 'y':
            j = 0
            for i in range(len(peptide) - int(ion_format[1]) - 3,len(peptide) - int(ion_format[1])+3):
                attribute[20 + 20*j+dicS[peptide[i]]] = 1
                j += 1
        elif int(ion_format[1]) > 3 and int(ion_format[1]) < len(peptide) - 3 and ion_format[0] == 'b':
            j = 0
            for i in range(len(ion) - 3,len(ion) + 3):
                attribute[20 + j*20+dicS[peptide[i]]] = 1
                j += 1
        elif int(ion_format[1]) <= 3 and ion_format[0] == 'y' or int(ion_format[1]) >= len(peptide) - 3 and ion_format[0] == 'b':
            j = 0
            if int(ion_format[1]) == 1 or int(ion_format[1]) == len(peptide) - 1:
                for i in range(len(peptide) - 4,len(peptide)):
                    attribute[20 + 20*j + dicS[peptide[i]]] = 1
                    j += 1
            elif int(ion_format[1]) == 2 or int(ion_format[1]) == len(peptide) - 2:
                for i in range(len(peptide) - 5,len(peptide)):
                    attribute[20 + 20*j + dicS[peptide[i]]] = 1
                    j += 1
            else:    
                for i in range(len(peptide) - 6,len(peptide)):
                    attribute[20 + 20*j + dicS[peptide[i]]] = 1
                    j += 1
        elif int(ion_format[1]) <= 3 and ion_format[0] == 'b' or int(ion_format[1]) >= len(peptide) - 3 and ion_format[0] == 'y':
            if int(ion_format[1]) == 1 or int(ion_format[1]) == len(peptide) - 1:
                j = 2
                for i in range(4):
                   attribute[20 + j*20+dicS[peptide[i]]] = 1
                   j += 1
            elif int(ion_format[1]) == 2 or int(ion_format[1]) == len(peptide) - 2:
                j = 1
                for i in range(5):
                   attribute[20 + j*20+dicS[peptide[i]]] = 1
                   j += 1   
            else:
                for i in range(6):
                    attribute[20 + i*20+dicS[peptide[i]]] = 1

        attribute[140 + dicS[peptide[0]]] = 1
        attribute[160 + dicS[peptide[-1]]] = 1

        if ion_format[0] == 'b':
            attribute[180] = float(ion_format[1])/len(peptide)
            attribute[181] = (len(peptide) - float(ion_format[1]))/len(peptide)
            #attribute[182] = abs(len(peptide)/2.0 - float(ion_format[1]))/len(peptide)
        else:
            attribute[181] = float(ion_format[1])/len(peptide)
            attribute[180] = (len(peptide) - float(ion_format[1]))/len(peptide)
            #attribute[182] = abs(len(peptide)/2.0 - float(ion_format[1]))/len(peptide)

        if ion_format[1] == '1':
            attribute[182] = 1
        
        if pep != peptide:
            pep_mz = 0.0;pep_b = 0.0;pep = peptide;pep_he = 0
            for i in range(len(peptide)):
                pep_mz += dicM[peptide[i]]
                pep_b += dicB[peptide[i]]
                pep_he += dicHe[peptide[i]]
        ion_mz = 0.0;ion_b = 0.0
        for i in range(int(ion_format[1])):
            ion_mz += dicM[ion[i]]
            ion_b += dicB[ion[i]]
        attribute[183] = ion_mz/pep_mz
        attribute[184] = (pep_mz - ion_mz)/pep_mz
        if ion_format[0] == 'b':
            attribute[185] = dicB[peptide[int(ion_format[1])]]/pep_b
            attribute[186] = dicB[peptide[int(ion_format[1]) - 1]]/pep_b
            attribute[187] = dicHe[peptide[int(ion_format[1])]]/pep_he
            attribute[188] = dicHe[peptide[int(ion_format[1]) - 1]]/pep_he
        else:
            attribute[185] = dicB[peptide[-int(ion_format[1])]]/pep_b
            attribute[186] = dicB[peptide[-int(ion_format[1]) - 1]]/pep_b
            attribute[187] = dicHe[peptide[-int(ion_format[1])]]/pep_he
            attribute[188] = dicHe[peptide[-int(ion_format[1]) - 1]]/pep_he
        attribute[189] = ion_b/pep_b
        attribute[190] = (pep_b - ion_b)/pep_b
        attribute[191] = (peptide.count('K')+peptide.count("R")+peptide.count('H'))/len(peptide)
        attribute[192]=pep_mz/2.0
        cunt+=1
       # data.append(attribute)
       # labels.append(float(f_line[6].split('\n')[0]))
        data=attribute.tolist()
        data.insert(0,cunt)
        data.insert(1,f_line[0])
        data.append(float(f_line[6].split('\n')[0]))
        writer.writerow(data)
        print('write '+str(cunt)+' line in file')
    f.close()
f1.close()

    #print len(data)
    #print len(labels)

    #Pickle.dump((array(data[:11134]),labels[:11134]),f1,True)
    #Pickle.dump((array(data[11134:22272]),labels[11134:22272]),f1,True)
    #Pickle.dump((array(data[22272:]),labels[22272:]),f1,True)

    #'''Pickle.dump((array(data[:17593]),labels[:17593]),f1,True)
    #Pickle.dump((array(data[17593:35175]),labels[17593:35175]),f1,True)
    #Pickle.dump((array(data[35175:]),labels[35175:]),f1,True)'''

    #f1.close()

    #print "VICTORY!"

