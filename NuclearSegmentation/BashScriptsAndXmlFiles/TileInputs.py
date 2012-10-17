import sys
import os
import shutil
import fnmatch
import re
import pdb
import sqlite3
import csv
import pprint
import itk
from string import maketrans
from xml.etree.ElementTree import Element, SubElement, Comment, tostring
import xml.dom.minidom
from xml.dom.minidom import Node
import xml.etree.ElementTree as ET

#Set of Global Vars
Use_recursive_folders = 0
Working_folder = ''
X_limit = 0
Y_limit = 0
Z_limit = 0
Tile_Size = 2048

def usage():
	print \
	"""Incorrect number of arguments"
	Correct usage: *.py Definition.xml (optional 1 to use separate folders for tiles)
	"""
	sys.exit()

def Get_Image_Info( Image_name ):
	pixelType = itk.US
	imageType = itk.Image[pixelType, 3]
	readerType = itk.ImageFileReader[imageType]
	reader = readerType.New()
	print "Reading the image ", Image_name, " to check image sizes"
	reader.SetFileName( Image_name )
	reader.Update()
	global X_limit, Y_limit, Z_limit
	X_limit = reader.GetOutput().GetLargestPossibleRegion().GetSize()[0]
	Y_limit = reader.GetOutput().GetLargestPossibleRegion().GetSize()[1]
	Z_limit = reader.GetOutput().GetLargestPossibleRegion().GetSize()[2]
	print "The image size is ", X_limit, Y_limit, Z_limit, " and the tile size is ",Tile_Size

def Write_Cropped_Image( In_Image_name, Out_Image, Start_X, Start_Y ):
	pixelType  = itk.US
	imageType  = itk.Image[pixelType, 3]
	readerType = itk.ImageFileReader[imageType]
	reader     = readerType.New()
	reader.SetFileName( In_Image_name )
	size  = itk.Size[3](); 
	index = itk.Index[3](); index[0] = Start_X; index[1] = Start_Y; index[2] = 0
	roi   = itk.ImageRegion[3]()
	roi.SetIndex(index)
	if (Start_X+Tile_Size)<X_limit and (Start_Y+Tile_Size)<Y_limit:
		#Crop and write image
		size[0] = Tile_Size; size[1] = Tile_Size; size[2] = Z_limit
		roi.SetSize(size)
		extractType = itk.ExtractImageFilter[imageType,imageType]
		extractfilt = extractType.New()
		extractfilt.SetInput( reader.GetOutput() )
		extractfilt.SetExtractionRegion(roi)
		writerType = itk.ImageFileWriter[imageType]
		writer     = writerType.New()
		writer.SetFileName( Out_Image )
		writer.SetInput( extractfilt.GetOutput() )
		writer.Update()
	else:
		#Crop image
		size[0]=X_limit-Start_X; size[1]=Y_limit-Start_Y; size[2]=Z_limit
		roi.SetSize(size)
		extractType = itk.ExtractImageFilter[imageType,imageType]
		extractfilt = extractType.New()
		extractfilt.SetInput( reader.GetOutput() )
		extractfilt.SetExtractionRegion(roi)
		extractfilt.Update()
		cropIm = itk.Image.US3
		cropIm = extractfilt.GetOutput()

		#Get padded image
		size1 = itk.Size[3](); size1[0] = Tile_Size; size1[1] = Tile_Size; size1[2] = Z_limit
		roi1  = itk.ImageRegion[3]()
		roi1.SetIndex(index)
		roi1.SetSize(size1)
		outIm = itk.Image.US3
		outIm.New(Regions=[Tile_Size, Tile_Size, Z_limit])
		outIm.Allocate()
		outIm.FillBuffer(0)

		#Copy cropped image into padded image and write it ##############
		buf = itk.PyBuffer[itk.Image.US3].GetArrayFromImage(outIm)
		print buf
"""		writerType = itk.ImageFileWriter[imageType]
		writer     = writerType.New()
		writer.SetFileName( Out_Image )
		writer.SetInput( outIm )
		writer.Update() """

def Get_Cropped_Region_Labels(In_Image_name):
	pixelType  = itk.US
	imageType  = itk.Image[pixelType, 3]
	readerType = itk.ImageFileReader[imageType]
	reader     = readerType.New()
	reader.SetFileName( In_Image_name )
	labelStatsType = itk.LabelGeometryImageFilter[imageType,imageType]
	labelfilt = labelStatsType.New()
	labelfilt.SetInput( reader.GetOutput() )
	labelfilt.Update()
	Labels = labelfilt.GetLabels()
	return Labels

def Write_Channel_Tiles( Input_xml_file ):
	#Get channel names, RGB and filenames
	List_Channel_Series = []
	Channel_Names = []; R = []; G = []; B = []; File_Names = [];
	Channels_Xml_Elements = xml.dom.minidom.parse(Working_folder+Input_xml_file)
	for node in Channels_Xml_Elements.getElementsByTagName("file"):
		Channel_Names.append(str(node.getAttribute("chname")))
		R.append(str(node.getAttribute("r")))
		G.append(str(node.getAttribute("g")))
		B.append(str(node.getAttribute("b")))
		for Images in node.childNodes: #There should only be one
			Image_W_Full_Path = str(Images.data)
			Intab = "\\"
			Outtab = "/"
			Trantab = maketrans( Intab, Outtab )
			Image_W_Full_Path = Image_W_Full_Path.translate(Trantab)
			File_Names.append(Image_W_Full_Path)
	X_Start = 0; Y_Start = 0; count=0;
	while Y_Start<Y_limit:
		while X_Start<X_limit:
			#Generate directory string
			Dir_String = Working_folder
			if Use_recursive_folders:
				Dir_String = Dir_String + str(count) + '/'
			if not os.path.exists(Dir_String): os.makedirs(Dir_String)

			#Generate XML elements for each cropped channel
			ImageElement = ET.Element('Image')
			for i in range(len(Channel_Names)):
				FileElement = ET.SubElement( ImageElement,'file' )
				FileElement.set( 'chname', Channel_Names[i] )
				FileElement.set( 'r', R[i] )
				FileElement.set( 'g', G[i] )
				FileElement.set( 'b', B[i] )
				Current_Filename = File_Names[i]
				OutFileName = Dir_String + Current_Filename[Current_Filename.rfind('/')+1:Current_Filename.rfind('.')]
				OutFileName = OutFileName + '_' + str(X_Start) +  '_' + str(Y_Start)
				OutFileName = OutFileName + Current_Filename[Current_Filename.rfind('.'):]
				FileElement.text = OutFileName
				#Write cropped channel
				Write_Cropped_Image(Current_Filename, OutFileName, X_Start, Y_Start)
			#Write XML
			tree = ET.ElementTree(ImageElement)
			OutXmlName = Dir_String + 'ChannelSeries' + '_' + str(X_Start) +  '_' + str(Y_Start) + '.xml'
			tree.write(OutXmlName)
			List_Channel_Series.append(OutXmlName)
			count = count + 1
			X_Start = X_Start+Tile_Size
		X_Start = 0
		Y_Start = Y_Start+Tile_Size
	#Write main XML to link to sub-XMLs
	MainElement = ET.Element('Image')
	for i in range(len(List_Channel_Series)):
		FileElement = ET.SubElement( MainElement,'file' )
		FileElement.text = List_Channel_Series[i]
	tree = ET.ElementTree(MainElement)
	OutXmlName = Working_folder + 'ChannelSeries' + '.xml'
	tree.write(OutXmlName)

def Write_Label_Tiles( Input_file ):
	#Check if it is xml
	IsXML = Input_file.find('xml')
	List_Label_Series = []
	LabelImage = ''
	Labels_List = []
	if IsXML > 0:
		Channels_Xml_Elements = xml.dom.minidom.parse(Working_folder+Input_file)
		for node in Channels_Xml_Elements.getElementsByTagName("file"):
			for Images in node.childNodes: #There should only be one
				Image_W_Full_Path = str(Images.data)
				Intab = "\\"
				Outtab = "/"
				Trantab = maketrans( Intab, Outtab )
				Image_W_Full_Path = Image_W_Full_Path.translate(Trantab)
				LabelImage = Image_W_Full_Path
	else:
		LabelImage = Input_file

	X_Start = 0; Y_Start = 0; count=0;
	while Y_Start<Y_limit:
		while X_Start<X_limit:
			#Generate directory string
			Dir_String = Working_folder
			if Use_recursive_folders:
				Dir_String = Dir_String + str(count) + '/'
			if not os.path.exists(Dir_String): os.makedirs(Dir_String)

			#Generate XML elements for each cropped label
			OutFileName = Dir_String + LabelImage[LabelImage.rfind('/')+1:LabelImage.rfind('.')]
			OutFileName = OutFileName + '_' + str(X_Start) +  '_' + str(Y_Start)
			OutFileName = OutFileName + LabelImage[LabelImage.rfind('.'):]
			ImageElement = ET.Element('Image')
			FileElement = ET.SubElement( ImageElement,'file' )
			FileElement.text = OutFileName
			#Write cropped channel
			Write_Cropped_Image(Working_folder+LabelImage, OutFileName, X_Start, Y_Start)
			#Get_Labels_In_Cropped_Region
			Current_Labels1 = Get_Cropped_Region_Labels(OutFileName)
			Current_Labels = list( Current_Labels1 )
			Current_Labels.sort()
			Labels_List.append(Current_Labels)
			#Write XML
			tree = ET.ElementTree(ImageElement)
			OutXmlName = Dir_String + 'LabelSeries' + '_' + str(X_Start) +  '_' + str(Y_Start) + '.xml'
			tree.write(OutXmlName)
			List_Label_Series.append(OutXmlName)
			count = count + 1
			X_Start = X_Start+Tile_Size
		X_Start = 0
		Y_Start = Y_Start+Tile_Size
	#Write main XML to link to sub-XMLs
	MainElement = ET.Element('Image')
	for i in range(len(List_Label_Series)):
		FileElement = ET.SubElement( MainElement,'file' )
		FileElement.text = List_Label_Series[i]
	tree = ET.ElementTree(MainElement)
	OutXmlName = Working_folder + 'LabelSeries' + '.xml'
	tree.write(OutXmlName)
	return Labels_List

def Write_Table_For_Tiles( Labels_List, Table_name ):
	#Read master table
	In_File_Ptr = open(Working_folder+Table_name, 'rb')
	spamReader = csv.reader(In_File_Ptr, delimiter='\t', quotechar='|')
	index = 1
	tableEntries = []
	Table_List = []
	for row in spamReader:
		row.remove('')
		if index==1:
			columnHeaders = row
			index+=1
		else:
			tableEntries.append(row)
	In_File_Ptr.close()

	X_Start = 0; Y_Start = 0; count=0;
	while Y_Start<Y_limit:
		while X_Start<X_limit:
			#Generate directory string
			Dir_String = Working_folder
			if Use_recursive_folders:
				Dir_String = Dir_String + str(count) + '/'
			if not os.path.exists(Dir_String): os.makedirs(Dir_String)

			OutFileName = Dir_String + Table_name[:Table_name.rfind('.')]
			OutFileName = OutFileName + '_' + str(X_Start) +  '_' + str(Y_Start)
			OutFileName = OutFileName + Table_name[Table_name.rfind('.'):]

			Table_List.append( OutFileName )

			#Open sub-table for writing
			Out_File_Ptr = open(OutFileName, 'wb')
			spamWriter = csv.writer(Out_File_Ptr, delimiter='\t',quotechar='|', quoting=csv.QUOTE_MINIMAL)
			spamWriter.writerow(columnHeaders)

			#Write values in sub-table
			for i in range(len(Labels_List[count])):
				if Labels_List[count][i]==0:
					continue
				for j in range(len(tableEntries)):
					if int(Labels_List[count][i])==int(float(tableEntries[j][0])):
						print 'here'
						spamWriter.writerow(tableEntries[j])

			Out_File_Ptr.close()

			count = count + 1
			X_Start = X_Start+Tile_Size
		X_Start = 0
		Y_Start = Y_Start+Tile_Size
	#Write main XML to link to tables
	MainElement = ET.Element('Table')
	for i in range(len(Table_List)):
		FileElement = ET.SubElement( MainElement,'file' )
		FileElement.text = Table_List[i]
	tree = ET.ElementTree(MainElement)
	OutXmlName = Working_folder + 'TableSeries' + '.xml'
	tree.write(OutXmlName)


def main(argv = None):
	if argv is None:
		argv = sys.argv
	if len(argv) < 2 or len(argv) > 3:
		usage()

	if len(argv) == 3:
		global Use_recursive_folders
		Use_recursive_folders = argv[2]
		print Use_recursive_folders, "-->if > 0, will create new folders for tiles"

	Main_Xml_Elements = xml.dom.minidom.parse(argv[1])

	#Set Working Folder
	global Working_folder
	for node in Main_Xml_Elements.getElementsByTagName("ProjectFiles"):
		Input_folder = node.getAttribute("path")
		#Change windows backslash for folder paths
		Intab = "\\"
		Outtab = "/"
		Trantab = maketrans( Intab, Outtab )
		print "The working folder is ", Input_folder
		Input_folder = str(Input_folder).translate(Trantab)
		Working_folder = Input_folder

	#Get Channels XML filename
	Channel_name = ''
	for node in Main_Xml_Elements.getElementsByTagName("input"):
		Input_channels = node.getAttribute("file")
		Channel_name = str(Input_channels)

	#Get Segmentation Label XML filename
	Label_name = ''
	for node in Main_Xml_Elements.getElementsByTagName("output"):
		Input_label = node.getAttribute("file")
		Label_name = str(Input_label)

	#Get Table filename
	Table_name = ''
	for node in Main_Xml_Elements.getElementsByTagName("table"):
		Input_table = node.getAttribute("file")
		Table_name = str(Input_table)

	#Get the size of one channel
	Channel_Xml_Elements = xml.dom.minidom.parse(Working_folder+Channel_name)
	for node in Channel_Xml_Elements.getElementsByTagName("file"):
		for Images in node.childNodes:
			Image_name = str(Images.data)
			Get_Image_Info( Image_name )
			Out_Image = Image_name[:-5]+'_crop.tif'				#######
			print Out_Image							#######
#			Write_Cropped_Image( Image_name, Out_Image, 10100, 10100 ) 	####### 
			break
		break

	#Write_Channel_Tiles( Channel_name )			##############
	#Labels_List = Write_Label_Tiles  (  Label_name  )	##############
	#Write_Table_For_Tiles( Labels_List, Table_name )	##############

	#Write the main xml file
	MainElement = ET.Element('ProjectFiles')
	MainElement.set( 'path', Working_folder )
	MainElement.set( 'name', 'Segmentation' )
	MainElement.set( 'type', 'Multi' )

	ChannelElement = ET.SubElement( MainElement, 'input' )
	ChannelElement.set( 'file', 'ChannelSeries.xml' )

	LabelElement = ET.SubElement( MainElement, 'output' )
	LabelElement.set( 'file', 'LabelSeries.xml' )

	LogElement = ET.SubElement( MainElement, 'log' )
	Log_File = ''; Validate_Code = '';
	for node in Main_Xml_Elements.getElementsByTagName("log"):
		Log_File = node.getAttribute("file")
		Validate_Code = node.getAttribute("validated")
	LogElement.set( 'file', Log_File )
	LogElement.set( 'validated', Validate_Code )

	DefnElement = ET.SubElement( MainElement, 'definition' )
	Defn_File = ''
	for node in Main_Xml_Elements.getElementsByTagName("definition"):
		Defn_File = node.getAttribute("file")
	DefnElement.set( 'file', Defn_File )

	TableElement = ET.SubElement( MainElement, 'table' )
	TableElement.set( 'file', 'TableSeries.xml' )

	tree = ET.ElementTree(MainElement)
	tree.write(Working_folder+'Tiles.xml')

	sys.exit()

if __name__ == "__main__":
        main()
