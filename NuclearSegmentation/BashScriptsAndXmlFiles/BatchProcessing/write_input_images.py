#!/usr/bin/python

import fnmatch
import os
import shutil
import re

print "Hello, World!"

rootPath = os.getcwd()
rootPath = rootPath + '\\'
patternN    = '*Nuclear*' # Can include any UNIX shell-style wildcards
patternCD34 = '*CD34.tif'
patternCA9  = '*CA9.tif'
patternSMA  = '*SMA.tif'
patternKI   = '*pERK.tif'
patternRGB  = '*QB_RGB.tif'

ppatternN    = 'full_path/NuclearH.tif' # Can include any UNIX shell-style wildcards
ppatternCD34 = 'full_path/CD34.tif'
ppatternCA9  = 'full_path/CA9.tif'
ppatternSMA  = 'full_path/SMA.tif'
ppatternKI   = 'full_path/Ki67.tif'
ppatternRGB  = 'full_path/RGB.tif'

inpu_xml_img     = 'Histo_Input_Image.xml'
replace_pattern1 = 'full_path/'
lab_im           = 'Histo_Input_Image_label.xml'
tab_le           = 'Histo_Input_Image_table.txt'
proc_def_xml     = 'HistoProjectDef.xml'
op_proj          = 'ccrc.xml'
pp = 'C:\\farsight_program\\bin\\projproc'

for root, dirs, files in os.walk(rootPath):
	count_file = 0
	asd = cmp(os.path.join(root,''),rootPath)
	print asd
	if asd != 0:
		newN = ''
		newCD34 = ''
		newCA9 = ''
		newSMA = ''
		newKI =''
		newRGB = ''

		for filename in fnmatch.filter(files, patternN):
			newN = os.path.join(root,filename)
			count_file += 1

		print newN
		for filename in fnmatch.filter(files, patternCD34):
			newCD34 = os.path.join(root,filename)
			count_file += 1

		print newCD34
		for filename in fnmatch.filter(files, patternCA9):
			newCA9 = os.path.join(root,filename)
			count_file += 1

		print newCA9
		for filename in fnmatch.filter(files, patternSMA):
			newSMA = os.path.join(root,filename)
			count_file += 1

		print newSMA
		for filename in fnmatch.filter(files, patternKI):
			newKI = os.path.join(root,filename)
			count_file += 1

		print newKI 
		for filename in fnmatch.filter(files, patternRGB):
			newRGB = os.path.join(root,filename)
			count_file += 1

		print newRGB
			#print replace_pattern
		print os.path.join(root,'')
		print count_file
		if count_file == 6:
			#Write image file
			o = open(os.path.join(root,inpu_xml_img),"w")
			data = open(rootPath+inpu_xml_img).read()
			data = re.sub(ppatternN,newN,data)
			data = re.sub(ppatternCD34,newCD34,data)
			data = re.sub(ppatternCA9,newCA9,data)
			data = re.sub(ppatternSMA,newSMA,data)
			data = re.sub(ppatternKI,newKI,data)
			data = re.sub(ppatternRGB,newRGB,data)
			o.write( data )
			o.close()
			#Write output project file
			o1 = open(os.path.join(root,op_proj),"w")
			data1 = open(rootPath+op_proj).read()
			newPPP = os.path.join(root,op_proj)
	 		data1 = re.sub(replace_pattern1,newPPP,data1)
			o1.write( data1 )
			o1.close()
			#shutil.copy2( rootPath+inpu_xml_img, os.path.join(root,inpu_xml_img))
			labim1 = os.path.join(root,lab_im)
			tab1 = os.path.join(root,tab_le)
			os.chdir( os.path.join(root,'') )
			com_to_exec = pp + ' \"' + os.path.join(root,inpu_xml_img) + '\" \"' + labim1 + '\" \"' + tab1 + '\" \"' + rootPath + proc_def_xml + '\"'
			os.system( com_to_exec )
			os.chdir(rootPath)
			print com_to_exec
