with open('data//ion2_pif_0.8_y.txt')  as f:
    tempStr=''
    list_temp=[];list_less_than_2=[];list_more_than_2=[]
    while True:
        list=[]
        line = f.readline()
        if not line:
            break
        list = line.split()
        if tempStr == '':
            list_temp.append(list)
            tempStr=list[0]
        else:
            if list[0]==tempStr:
                list_temp.append(list)
            elif len(list_temp) > 2:
                list_more_than_2.append(list_temp)
                tempStr=list[0]
                list_temp=[]
                list_temp.append(list)
            else:
                list_less_than_2.append(list_temp)
                tempStr=list[0]
                list_temp=[]
                list_temp.append(list)
f.close()
with open('data//ion2_pif_0.8_y_new.txt','w') as fw:
    for i in range(len(list_more_than_2)):
        for j in range(len(list_more_than_2[i])):
            line_str=''
            k=0
            for a in list_more_than_2[i][j]:
                k+=1
                if k==7:
                    line_str += (str(a)+'\n')
                else:
                    line_str += (str(a)+'\t')
            fw.write(line_str)
fw.close()
