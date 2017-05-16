from math import sqrt
def get_merge_list(data):
    temp_peptide = ''
    temp_list = []
    intensity_list = []
    #print(data)
    for index,row in data.iterrows():
        peptide = row.Peptide
        if peptide != temp_peptide:
            if temp_peptide == '':
                temp_list.append([row.Number,row.Intensity])
            else:
                intensity_list.append(temp_list)
                temp_list = []
                temp_list.append([row.Number,row.Intensity])
            temp_peptide = peptide
        else:
            temp_list.append([row.Number,row.Intensity])
    intensity_list.append(temp_list)
    return intensity_list
def get_peptide_list(data):
    temp_peptide = ''
    temp_list = []
    peptide_list = []
    for row in data.itertuples():
        peptide = row.Peptide
        if peptide != temp_peptide:
            if temp_peptide == '':
                temp_list.append(row.Peptide)
            else:
                peptide_list.append(temp_list)
                temp_list = []
                temp_list.append(row.Peptide)
            temp_peptide = peptide
       # else:
        #    temp_list.append(row.Peptide)
    peptide_list.append(temp_list)
    return peptide_list
def get_y_Mass_list(data):
    temp_peptide = ''
    temp_list = []
    y_Mass_list = []
    for row in data.itertuples():
        peptide = row.Peptide
        if peptide != temp_peptide:
            if temp_peptide == '':
                temp_list.append(row.y_Mass)
            else:
                y_Mass_list.append(temp_list)
                temp_list = []
                temp_list.append(row.y_Mass)
            temp_peptide = peptide
        else:
            temp_list.append(row.y_Mass)
    y_Mass_list.append(temp_list)
    return y_Mass_list
def multipl(a,b):
    multipl_of_a_b = 0.0
    for i in range(len(a)):
        temp = a[i] * b[i]
        multipl_of_a_b+=temp
    return multipl_of_a_b
# https://wikimedia.org/api/rest_v1/media/math/render/svg/832f0c5c22a0d6f2596c150de811247438a503de
def pearson_r(x,y):
    n = len(x)
    sum_x = sum(x)
    sum_y = sum(y)
    sum_of_xy = multipl(x,y)
    sum_of_x2 = sum([pow(i,2) for i in x])
    sum_of_y2 = sum([pow(j,2) for j in y])
    numerator = sum_of_xy - (float(sum_x) * float(sum_y) / n)
    denominator = sqrt((sum_of_x2 - float(sum_x ** 2) / n) * (sum_of_y2 - float(sum_y ** 2) / n))
    return numerator / denominator
