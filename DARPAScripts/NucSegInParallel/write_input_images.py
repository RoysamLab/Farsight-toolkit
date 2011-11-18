#!/usr/bin/python

import fnmatch
import os
import shutil
import re
import sqlite3
import subprocess
import multiprocessing
import time
from string import maketrans

print "Hello, World!"

rootPath = os.getcwd()
rootPath = rootPath + '/'
patternN    = 'Chan3*' # Can include any UNIX shell-style wildcards
patternGFP = 'GCHAN3*'

#Things to replace in the input image file
ppatternN    = 'full_pathNuclear.tif' # Can include any UNIX shell-style wildcards
ppatternGFP = 'full_pathGFP.tif'

#Number of charecters to match in the filenames from the end
num_cahr_match = 13

inpu_xml_img     = 'Histo_Input_Image'
inpu_xml_img11   = 'Histo_Input_Imagess'
inpu_xml_img1    = 'Histo_Input_Image.xml'
replace_pattern1 = 'full_path'
lab_im           = 'Histo_Input_Image_label'
tab_le           = 'Histo_Input_Image_table'
proc_def_xml     = 'HistoProjectDef.xml'
op_proj          = 'darpa'
op_proj1         = 'darpa.xml'
exx_ml		 = '.xml'
texxt		 = '.txt'

pp = '/data/kedar/farsight-bin/exe/projproc'
pppppp = '\"'
for root, dirs, files in os.walk(rootPath):
	count_file  = 0
	count_file1 = 0
	asd = cmp(os.path.join(root,''),rootPath)
	#print asd
	if asd != 0:
		newN = []
		newGFP = []

		intab = "\\"
		outtab = "/"
		trantab = maketrans( intab, outtab )

		for filename in fnmatch.filter(files, patternN):
			newNNN = ''
			newNN = os.path.join(root,filename)
			count_file += 1
			newNNN = newNN.translate(trantab)
			newN.append( newNNN )
		for filename in fnmatch.filter(files, patternGFP):
			newGFPPP = ''
			newGFPP = os.path.join(root,filename)
			count_file1 += 1
			newGFPPP = newGFPP.translate(trantab)
			newGFP.append( newGFPPP )
		newN.sort()
		newGFP.sort()
		print newN
		print newGFP

		print os.path.join(root,'')
		#print count_file
		if count_file1 == count_file:
			#Write image file
			com_to_exec = []
			for x in range(len(newGFP)):
				inp_fil_str  = os.path.join(root,inpu_xml_img) + str(x) + exx_ml
				inp_fil_str1 = inpu_xml_img + str(x) + exx_ml
				lab_fil_str  = os.path.join(root,lab_im)       + str(x) + exx_ml
				lab_fil_str1 = lab_im       + str(x) + exx_ml
				tab_fil_str  = os.path.join(root,tab_le)       + str(x) + texxt
				tab_fil_str1 = tab_le       + str(x) + texxt
				prj_fil_str  = os.path.join(root,op_proj)      + str(x) + exx_ml
				prj_fil_str1 = op_proj      + str(x) + exx_ml
				o = open(inp_fil_str,"w")
				data = open(rootPath+inpu_xml_img1).read()
				data = re.sub(ppatternN,newN[x],data)
				data = re.sub(ppatternGFP,newGFP[x],data)
				#print data
				o.write( data )
				o.close()

				#Write output project file
				o1 = open(prj_fil_str,"w")
				data1 = open(rootPath+op_proj1).read()
				newPPPP = os.path.join(root,pppppp)
				newPPP = newPPPP.translate(trantab);
				data1 = re.sub(replace_pattern1, newPPP,       data1)
				data1 = re.sub(lab_im,           lab_fil_str1, data1)
				data1 = re.sub(tab_le,           tab_fil_str1, data1)
				data1 = re.sub(inpu_xml_img11,   inp_fil_str1, data1)
				o1.write( data1 )
				o1.close()

				o2 = open(os.path.join(root,proc_def_xml),"w")
				data2 = open(rootPath+proc_def_xml).read()
				o2.write( data2 )
				o2.close()

				com_to_exec_str = pp + ' \"' + inp_fil_str + '\" \"' + lab_fil_str + '\" \"' + tab_fil_str + '\" \"' + os.path.join(root,proc_def_xml) + '\"'
				com_to_exec.append( com_to_exec_str )
				#print com_to_exec
				#os.system( com_to_exec )
				#os.chdir(rootPath)
			for com_to_exec_str in com_to_exec:
				#os.system( com_to_exec_str )
				subprocess.Popen(com_to_exec_str, shell=True)
				#time in secs between two files in a folder
				time.sleep(2)
				#print com_to_exec_str
			#time in secs between two folders
			time.sleep(600)
