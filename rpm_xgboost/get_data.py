#coding:utf-8
import csv
from chemical_properties import *
class train_feature():
    def __init__(self):
        self.Charge=0
        self.Cleavage_site=0
        self.y_ion=''

        self.Peptide=''
        self.LengthOfPeptide=0
        self.McRatioOfPeptide=0.0
        self.CountOfBasicAA=0
        self.P_Basic=0.0
        self.P_Hydrophobicity=0.0
        self.P_Helicity=0.0
        self.IsolationWidthLeft=''
        self.IsolationWidthRight=''
        self.C_Identity=''
        self.N_Identity=''
        self.L_Basicity=0.0
        self.R_Basicity=0.0
        self.L_Hydrophobicity=0.0
        self.R_Hydrophobicity=0.0
        self.L_Helicity=0.0
        self.R_Helicity=0.0
        self.AvgOfBasic=0.0
        self.AvgOfHydrophobicity=0.0
        self.AvgOfHelicity=0.0
        self.DvalueOfBasic=0.0
        self.DvalueOfHydrophobicity=0.0
        self.DvalueOfHelicity=0.0
        self.b_Distance=0
        self.y_Distance=0
        self.b_Mass=0.0
        self.y_Mass=0.0
        self.DvaleOfMcratioOfyAndPep=0.0
        self.b_Basic=0.0
        self.y_Basic=0.0
        self.b_Hydrophobicity=0.0
        self.y_Hydrophobicity=0.0
        self.b_Helicity=0.0
        self.y_Helicity=0.0
        self.b_CountOfBasicAA=0
        self.y_CountOfBasicAA=0
        self.MassOfPeptide=0.0
        self.RatioOfbAndPeptide=0.0
        self.RatioOfyAndPeptide=0.0
        self.IsTheSamePort=False
        self.Intensity=0.0
global cunt
cunt=0
with open('data//temp_data.csv','wb') as f:
    writer=csv.writer(f)
    writer.writerow(['Number','Peptide','LengthOfPeptide','McRatioOfPeptide','CountOfBasicAA','P_Basic','P_Hydrophobicity','P_Helicity','IsolationWidthLeft','IsolationWidthRight','C_Identity','N_Identity','L_Basicity','R_Basicity','L_Hydrophobicity','R_Hydrophobicity','L_Helicity','R_Helicity','AvgOfBasic','AvgOfHydrophobicity','AvgOfHelicity','DvalueOfBasic','DvalueOfHydrophobicity','DvalueOfHelicity','b_Distance','y_Distance','b_Mass','y_Mass','DvaleOfMcratioOfyAndPep','b_Basic','y_Basic','b_Hydrophobicity','y_Hydrophobicity','b_Helicity','y_Helicity','b_CountOfBasicAA','y_CountOfBasicAA','MassOfPeptide','RatioOfbAndPeptide','RatioOfyAndPeptide','IsTheSamePort','Intensity'])
    orginal_data=open('data/mobile.txt','r')
    while True:
        line=orginal_data.readline()
        if not line:
            break
        flag=False
        try:
            feature=train_feature()
            line_list=line.split()
            feature.Peptide=line_list[0]
            feature.LengthOfPeptide=len(feature.Peptide)
            feature.Charge=line_list[1]
            feature.y_ion=line_list[2]
            feature.Cleavage_site=int(filter(str.isdigit,line_list[5]))
            feature.Intensity=line_list[6]
            feature.IsolationWidthLeft=feature.Peptide[feature.LengthOfPeptide-feature.Cleavage_site-1]
            feature.IsolationWidthRight=feature.Peptide[feature.LengthOfPeptide-feature.Cleavage_site]
            feature.C_Identity=feature.Peptide[0]
            feature.N_Identity=feature.Peptide[feature.LengthOfPeptide-1]
            feature.L_Basicity=dic_basicity[feature.IsolationWidthLeft]
            feature.R_Basicity=dic_basicity[feature.IsolationWidthRight]
            feature.L_Hydrophobicity=dic_hydrophobicity[feature.IsolationWidthLeft]
            feature.R_Hydrophobicity=dic_hydrophobicity[feature.IsolationWidthRight]
            feature.L_Helicity=dic_helicity[feature.IsolationWidthLeft]
            feature.R_Helicity=dic_helicity[feature.IsolationWidthRight]

            feature.AvgOfBasic=(feature.L_Basicity+feature.R_Basicity)/float(2)
            feature.AvgOfHydrophobicity=(feature.L_Hydrophobicity+feature.R_Hydrophobicity)/float(2)
            feature.AvgOfHelicity=(feature.L_Helicity+feature.R_Helicity)/float(2)
            feature.DvalueOfBasic=feature.L_Basicity-feature.R_Basicity
            feature.DvalueOfHydrophobicity=feature.L_Hydrophobicity-feature.R_Hydrophobicity
            feature.DvalueOfHelicity=feature.L_Helicity-feature.R_Helicity

            feature.b_Distance=feature.LengthOfPeptide-feature.Cleavage_site
            feature.y_Distance=feature.Cleavage_site

            for amino_acid_b in feature.Peptide[0:feature.b_Distance]:
                feature.b_Mass+=dic_mass[amino_acid_b]
                feature.b_Helicity+=dic_helicity[amino_acid_b]
                feature.b_Hydrophobicity+=dic_hydrophobicity[amino_acid_b]
                feature.b_Basic+=dic_basicity[amino_acid_b]
                if amino_acid_b == 'R' or amino_acid_b == 'K' or amino_acid_b == 'H':
                    feature.b_CountOfBasicAA+=1
            feature.MassOfPeptide+=feature.b_Mass
            feature.b_Mass += 1
            for amino_acid_y in feature.y_ion:
                feature.y_Mass+=dic_mass[amino_acid_y]
                feature.y_Helicity=dic_helicity[amino_acid_y]
                feature.y_Hydrophobicity=dic_hydrophobicity[amino_acid_y]
                feature.y_Basic=dic_basicity[amino_acid_y]
                if amino_acid_y == 'R' or amino_acid_y == 'K' or amino_acid_y == 'H':
                    feature.y_CountOfBasicAA+=1
            feature.MassOfPeptide+=feature.y_Mass
            feature.y_Mass += 19

            feature.MassOfPeptide+=float(18 + float(feature.Charge))
            feature.McRatioOfPeptide=feature.MassOfPeptide/float(feature.Charge)
            feature.DvaleOfMcratioOfyAndPep=feature.y_Mass-feature.McRatioOfPeptide

            feature.RatioOfbAndPeptide=feature.b_Mass/feature.MassOfPeptide
            feature.RatioOfyAndPeptide=feature.y_Mass/feature.MassOfPeptide

            feature.P_Basic=feature.b_Basic+feature.y_Basic
            feature.P_Helicity=feature.b_Helicity+feature.y_Helicity
            feature.P_Hydrophobicity=feature.b_Hydrophobicity+feature.y_Hydrophobicity
            feature.CountOfBasicAA=feature.b_CountOfBasicAA+feature.y_CountOfBasicAA
            if feature.Cleavage_site == 1 or feature.Cleavage_site == int(feature.LengthOfPeptide - 1):
                feature.IsTheSamePort =True
            else:
                feature.IsTheSamePort = False
            cunt += 1
            flag= True 
        except:
            print('error in ' +str(cunt))
            break
        if flag:
            writer.writerow([cunt,feature.Peptide,feature.LengthOfPeptide,feature.McRatioOfPeptide,feature.CountOfBasicAA,feature.P_Basic,feature.P_Hydrophobicity,feature.P_Helicity,feature.IsolationWidthLeft,feature.IsolationWidthRight,feature.C_Identity,feature.N_Identity,feature.L_Basicity,feature.R_Basicity,feature.L_Hydrophobicity,feature.R_Hydrophobicity,feature.L_Helicity,feature.R_Helicity,feature.AvgOfBasic,feature.AvgOfHydrophobicity,feature.AvgOfHelicity,feature.DvalueOfBasic,feature.DvalueOfHydrophobicity,feature.DvalueOfHelicity,feature.b_Distance,feature.y_Distance,feature.b_Mass,feature.y_Mass,feature.DvaleOfMcratioOfyAndPep,feature.b_Basic,feature.y_Basic,feature.b_Hydrophobicity,feature.y_Hydrophobicity,feature.b_Helicity,feature.y_Helicity,feature.b_CountOfBasicAA,feature.y_CountOfBasicAA,feature.MassOfPeptide,feature.RatioOfbAndPeptide,feature.RatioOfyAndPeptide,feature.IsTheSamePort,feature.Intensity])
            print('write ' + str(cunt) + ' line ')
           
f.close()