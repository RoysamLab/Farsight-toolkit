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
patternN    = 'Tile*' # Can include any UNIX shell-style wildcards
patternGFP  = 'Tige*'
patternCY5  = 'Cy5*'

#Things to replace in the input image file
ppatternN    = 'full_pathNuclear.tif' # Can include any UNIX shell-style wildcards
ppatternGFP = 'full_pathGFP.tif'
ppatternCY5 = 'full_pathCY5.tif'

inpu_xml_img     = 'Histo_Input_Image'
inpu_xml_img11   = 'Histo_Input_Imagess'
inpu_xml_img1    = 'Histo_Input_Image.xml'
replace_pattern1 = 'full_path'
lab_im           = 'Histo_Input_Image_label'
lab_im1          = 'Histo_Input_Image_label*'
tab_le           = 'Histo_Input_Image_table'
proc_def_xml     = 'HistoProjectDef.xml'
op_proj          = 'darpa'
op_proj1         = 'darpa.xml'
exx_ml		 = '.xml'
texxt		 = '.txt'

pp = '/data/kedar/farsight-bin/exe/projproc'
pppppp = '\"'
counttt = 0
for root, dirs, files in os.walk(rootPath):
	count_file  = 0
	count_file1 = 0
	count_file2 = 0
	count_file3 = 0
	asd = cmp(os.path.join(root,''),rootPath)
	#print asd
	if asd != 0:
		newN   = []
		newGFP = []
		newCY5 = []
		newDNE = []
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
		for filename in fnmatch.filter(files, patternCY5):
			newCY555 = ''
			newCY55 = os.path.join(root,filename)
			count_file3 += 1
			newCY555 = newCY55.translate(trantab)
			newCY5.append( newCY555 )
		for filename in fnmatch.filter(files, lab_im1):
			newDNEEE = ''
			newDNEE = os.path.join(root,filename)
			count_file2 += 1
			newDNEEE = newDNEE.translate(trantab)
			newDNE.append( newDNEEE )
		newN.sort()
		newGFP.sort()
		newCY5.sort()
		newDNE.sort()
		#print newN
		#print newGFP
		#print newDNE

		print os.path.join(root,'')
		#print count_file
		if count_file1 != count_file3:
			continue
		if count_file1 == count_file:
			#Write image file
			com_to_exec = []
			for x in range(len(newGFP)):
				start_pos = 0
				end_pos   = 0
				curxy     = newGFP[x]
				start_pos = curxy.rfind('/')
				curxy     = curxy[start_pos:]
				#print curxy
				start_pos = curxy.find('_')
				end_pos   = curxy.find('.')
				#print curxy[start_pos:end_pos]
				curxy     = curxy[start_pos:end_pos]
				inp_fil_str  = os.path.join(root,inpu_xml_img) + curxy + exx_ml
				inp_fil_str1 = inpu_xml_img + curxy + exx_ml
				lab_fil_str  = os.path.join(root,lab_im)       + curxy + exx_ml
				lab_fil_str1 = lab_im       + curxy + exx_ml
				tab_fil_str  = os.path.join(root,tab_le)       + curxy + texxt
				tab_fil_str1 = tab_le       + curxy + texxt
				prj_fil_str  = os.path.join(root,op_proj)      + curxy + exx_ml
				prj_fil_str1 = op_proj      + curxy + exx_ml
				o = open(inp_fil_str,"w")
				data = open(rootPath+inpu_xml_img1).read()
				data = re.sub(ppatternN,newN[x],data)
				data = re.sub(ppatternGFP,newGFP[x],data)
				data = re.sub(ppatternCY5,newCY5[x],data)
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

				found_op = 0
				for y in range(len(newDNE)):
					if newDNE[y]==lab_fil_str:
						found_op = 1
						#print lab_fil_str1
				if found_op==0:
					com_to_exec.append( com_to_exec_str )
			#print com_to_exec
				#print com_to_exec
				#os.system( com_to_exec )
				#os.chdir(rootPath)
#			for com_to_exec_str in com_to_exec:
				#os.system( com_to_exec_str )
#				counttt += 1
#				subprocess.Popen(com_to_exec_str, shell=True)
				#Number of instances to run in parallel
#				if counttt == 10:
#					counttt = 0
					#Time to run one instance
#					time.sleep(2700)
				#time in secs between two files in a folder
#				time.sleep(60)
				#print com_to_exec_str

