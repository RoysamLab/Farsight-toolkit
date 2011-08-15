#!/usr/bin/python
import fnmatch
import os
import shutil
import re
from string import maketrans

rootPath = os.getcwd()

CMK = 'CMakeLists.txt'
TLL = 'TARGET_LINK_LIBRARIES'
ITK = 'ITK'
OPN = '\('
OPR = '( ${ITK_LIBRARIES} '
OOP = '\$ITK_LIBRARIES\}'
OOP1 = '\${\ITK_LIBRARIES\}'
CLC = '\)'
CLS = ')'
CLN = ' ${ITK_LIBRARIES} )'

ATK = 'ITKCommon'
BTK = 'ITKIO'
CTK = 'ITKAlgorithms'
DTK = 'ITKgdcm'
ETK = 'ITKFEM'
FTK = 'ITKNumerics'
GTK = 'ITKStatistics'
HTK = 'ITKsys'
JTK = 'ITKBasicFilters'
EMP = ''

rootPath = os.getcwd()
print rootPath
for root, dirs, files in os.walk(rootPath):
        for filename in fnmatch.filter(files,CMK):
                fileI = os.path.join(root,filename)
                fileO = os.path.join(root,'tmp.txt')
                #print fileI
                #print fileO
                Ifile = open(fileI, 'r')
                Ofile = open(fileO, 'w')
                lin = Ifile.readline()
                while lin:
                        i = lin.find(TLL)
                        j = lin.find(ITK)
                        if i>-1 and j>-1:
                                k = lin.find(CLS)
                                #print lin
                                lin = re.sub(ATK,EMP,lin)
                                lin = re.sub(BTK,EMP,lin)
                                lin = re.sub(CTK,EMP,lin)
                                lin = re.sub(DTK,EMP,lin)
                                lin = re.sub(ETK,EMP,lin)
                                lin = re.sub(FTK,EMP,lin)
                                lin = re.sub(GTK,EMP,lin)
                                lin = re.sub(HTK,EMP,lin)
                                lin = re.sub(JTK,EMP,lin)
                                lin = re.sub(OOP,EMP,lin)
                                lin = re.sub(OOP1,EMP,lin)
                                lin = re.sub(CLC,CLN,lin)

                                print lin
                                Ofile.write(lin)
                                if k==-1:
                                        lin = Ifile.readline()
                                        print lin
                                        lin = re.sub(ATK,EMP,lin)
                                        lin = re.sub(BTK,EMP,lin)
                                        lin = re.sub(CTK,EMP,lin)
                                        lin = re.sub(DTK,EMP,lin)
                                        lin = re.sub(ETK,EMP,lin)
                                        lin = re.sub(FTK,EMP,lin)
                                        lin = re.sub(GTK,EMP,lin)
                                        lin = re.sub(HTK,EMP,lin)
                                        lin = re.sub(JTK,EMP,lin)
                                        lin = re.sub(OOP,EMP,lin)
                                        lin = re.sub(OOP1,EMP,lin)
                                        lin = re.sub(CLC,CLN,lin)

                                        #print lin
                                        Ofile.write(lin)
                                        k = lin.find(CLS)
                        else:
                                Ofile.write(lin)
                        lin  = Ifile.readline()
                Ifile.close()
                Ofile.close()
                os.remove(fileI)
                os.rename(fileO,fileI)
print "Yes he can!"