#!/usr/bin/python

import fnmatch
import os
import shutil
import re
import sqlite3
from string import maketrans

print "Writing Input Tables"

rootPath = os.getcwd()
rootPath = rootPath + '/'
patternN    = '*Histo_Input_Image_table*' # Can include any UNIX shell-style wildcards
ppatternN    = 'full_path/Histo_Input_Image_table.xml' # Can include any UNIX shell-style wildcards
inpu_xml_img     = 'Batch_tables.xml'
replace_pattern1 = 'full_path'


o = open(os.path.join(rootPath,inpu_xml_img),"a+")
o.write('<Table>\n')
o.close()

for root, dirs, files in os.walk(rootPath):
        count_file = 0
        asd = cmp(os.path.join(root,''),rootPath)
        if asd != 0:


                intab = "\\"
                outtab = "/"
                trantab = maketrans( intab, outtab )

                for filename in fnmatch.filter(files, patternN):
                        newNN = os.path.join(root,filename)
                        count_file += 1
                        newN = newNN.translate(trantab)
                        #Write image file
                        o = open(os.path.join(rootPath,inpu_xml_img),"a+")
                        o.write('\t<file>'+newN+'</file>\n')
                        o.close()
o = open(os.path.join(rootPath,inpu_xml_img),"a+")
o.write('</Table>')
o.close()
