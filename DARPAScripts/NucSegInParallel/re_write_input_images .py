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
rootPath = rootPath + '\\'
patternN    = 'Chan3*' # Can include any UNIX shell-style wildcards
patternGFP  = 'GCHAN3*'
patternLab  = 'Histo_Input_Image_label*_nuc.tif'

#Things to replace in the input image file
ppatternN    = 'full_pathNuclear.tif' # Can include any UNIX shell-style wildcards
ppatternGFP  = 'full_pathGFP.tif'
ppatternLab  = 'full_path_label'

inpu_xml_img     = 'Histo_Input_Image'
inpu_xml_img11   = 'Histo_Input_Imagess'
inpu_xml_img1    = 'Histo_Input_Image.xml'
replace_pattern1 = 'full_path'
lab_im           = 'Histo_Input_Image_label'
lab_im1          = 'Histo_Input_Image_label.xml'
tab_le           = 'Histo_Input_Image_table'
proc_def_xml     = 'HistoProjectDef.xml'
op_proj          = 'darpa'
op_proj1         = 'darpa.xml'
exx_ml		 = '.xml'
texxt		 = '.txt'

pppppp = '\"'
for root, dirs, files in os.walk(rootPath):
	count_file  = 0
	count_file1 = 0
	count_file2 = 0
	asd = cmp(os.path.join(root,''),rootPath)
	#print asd
	if asd != 0:
		newN   = []
		newGFP = []
		newLab = []

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
		for filename in fnmatch.filter(files, patternLab):
			newLABB = ''
			newLABB = os.path.join(root,filename)
			count_file2 += 1
			newLABB = newLABB.translate(trantab)
			newLab.append( newLABB )
		newN.sort()
		newGFP.sort()
		newLab.sort()

		print count_file
		print count_file1
		print count_file2

		print os.path.join(root,'')
		#print count_file
		if count_file1 == count_file:
			print newN
			print newGFP
			print newLab
			#Write image file
			for x in range(len(newGFP)):
				inp_fil_str  = os.path.join(root,inpu_xml_img) + 's' + str(x) + exx_ml
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

				o3 = open(os.path.join(root,lab_fil_str1),"w")
				data3 = open(rootPath+lab_im1).read()
				data3 = re.sub(ppatternLab, newLab[x], data3)
				o3.write( data3 )
				o3.close()
