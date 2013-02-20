#!/usr/bin/python

import fnmatch
import os
import shutil
import re
import sqlite3
from string import maketrans

print "Hello, World!"

rootPath = os.getcwd()
rootPath = rootPath + '\\'
patternN    = '*Dapi.jpg' # Can include any UNIX shell-style wildcards
patternCD34 = '*Qdot 605.jpg'
patternCA9  = '*Qdot 655.jpg'
#patternSMA  = '*Qdot565*'
#patternKI   = 'Image_Cube_C1.tif' --1
#patternRGB  = '*QB_RGB.tif' --2

ppatternN    = 'fullpathDapi.jpg' # Can include any UNIX shell-style wildcards
ppatternCD34 = 'fullpath 605.jpg'
ppatternCA9  = 'fullpath 655.jpg'
#ppatternSMA  = 'fullpath565.jpg'
#ppatternKI   = 'full_path/Image_Cube_C1.tif' --1
#ppatternRGB  = 'full_path/RGB.tif' --2

inpu_xml_img     = 'Histo_Input_Image.xml'
replace_pattern1 = 'full_path'
lab_im           = 'Histo_Input_Image_label.xml'
tab_le           = 'Histo_Input_Image_table.txt'
proc_def_xml     = 'Histo_Project_Definition1.xml'
op_proj          = 'ccrc.xml'
pp = 'C:\\farsight_program\\bin\\projproc'
pppppp = '\"'
for root, dirs, files in os.walk(rootPath):
        count_file = 0
        asd = cmp(os.path.join(root,''),rootPath)
        print asd
        if asd != 0:
                newN = ''
                newCD34 = ''
                newCA9 = ''
                newSMA = ''
                #newKI ='' --1
                #newRGB = '' --1

                intab = "\\"
                outtab = "/"
                trantab = maketrans( intab, outtab )

                for filename in fnmatch.filter(files, patternN):
                        newNN = os.path.join(root,filename)
                        count_file += 1
                        newN = newNN.translate(trantab)

                print newN
                for filename in fnmatch.filter(files, patternCD34):
                        newCD344 = os.path.join(root,filename)
                        count_file += 1
                        newCD34 = newCD344.translate(trantab)

                print newCD34
                for filename in fnmatch.filter(files, patternCA9):
                        newCA99 = os.path.join(root,filename)
                        count_file += 1
                        newCA9 = newCA99.translate(trantab)

                print newCA9
                #for filename in fnmatch.filter(files, patternSMA):
                #        newSMAA = os.path.join(root,filename)
                #        count_file += 1
                #        newSMA = newSMAA.translate(trantab);

                print newSMA
                #for filename in fnmatch.filter(files, patternKI): --1
                #       newKI = os.path.join(root,filename) --1
                #       count_file += 1 --1

                #print newKI --1
                #for filename in fnmatch.filter(files, patternRGB): --2
                #       newRGB = os.path.join(root,filename) --2
                #       count_file += 1 --2

                #print newRGB --2
                        #print replace_pattern
                print os.path.join(root,'')
                print count_file
                if count_file == 3:
                        #Write image file
                        o = open(os.path.join(root,inpu_xml_img),"w")
                        data = open(rootPath+inpu_xml_img).read()
                        print ppatternN
                        print newN
                        print data
                        
                        data = re.sub(ppatternN,newN,data)
                        data = re.sub(ppatternCD34,newCD34,data)
                        data = re.sub(ppatternCA9,newCA9,data)
                        #data = re.sub(ppatternSMA,newSMA,data)
                        #data = re.sub(ppatternKI,newKI,data) --1
                        #data = re.sub(ppatternRGB,newRGB,data) --2
                        o.write( data )
                        o.close()
                        #Write output project file
                        o1 = open(os.path.join(root,op_proj),"w")
                        data1 = open(rootPath+op_proj).read()
                        newPPPP = os.path.join(root,pppppp)
                        newPPP = newPPPP.translate(trantab);
                        data1 = re.sub(replace_pattern1,newPPP,data1)
                        o1.write( data1 )
                        o1.close()
                        o2 = open(os.path.join(root,proc_def_xml),"w")
                        data2 = open(rootPath+proc_def_xml).read()
                        o2.write( data2 )
                        o2.close()
                        #shutil.copy2( rootPath+inpu_xml_img, os.path.join(root,inpu_xml_img))
                        labim1 = os.path.join(root,lab_im)
                        tab1 = os.path.join(root,tab_le)
                        os.chdir( os.path.join(root,'') )
                        com_to_exec = pp + ' \"' + os.path.join(root,inpu_xml_img) + '\" \"' + labim1 + '\" \"' + tab1 + '\" \"' + os.path.join(root,proc_def_xml) + '\"'
                        print com_to_exec
                        #os.system( com_to_exec )
                        os.chdir(rootPath)
                        print com_to_exec
