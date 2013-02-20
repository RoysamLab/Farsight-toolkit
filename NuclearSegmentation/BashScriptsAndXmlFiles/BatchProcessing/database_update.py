#!/bin/python

import sqlite3
import re
import csv
import math
import fnmatch
import os
databaseFile = "/data/gramak/db_test_pfizer/NE.s3db"
images_path  = "/data/gramak/db_test_pfizer/Alan_Powers_Ki67"
tableFiles   = "*table.txt"
inputImage   = "*_Input_Image.xml"
#Get all the files with pattern tableFiles in images_path
inputFiles  = []
inputImages = []
for root, dirs, files in os.walk(images_path):
	for filename in fnmatch.filter(files, tableFiles):
		inputFiles.append(os.path.join(root,filename))
for root, dirs, files in os.walk(images_path):
	for filename in fnmatch.filter(files, inputImage):
		inputImages.append(os.path.join(root,filename))

for x in range(len(inputFiles)):
	spamReader = csv.reader(open(inputFiles[x], 'rb'), delimiter='\t', quotechar='|')

	#Check column headers
	index = 1
	tableEntries = []
	for row in spamReader:
		row.remove('')
		if index==1:
			columnHeaders = row
			index+=1
		else:
			tableEntries.append(row)
	ColumnHeadersTable = []
	for y in range(len(columnHeaders)):
		ColumnHeadersTable.append(columnHeaders[y].upper())
#	print ColumnHeadersTable
	columnNamesCommand = "PRAGMA table_info(IMAGE_TEST);"
	con1 = sqlite3.connect(databaseFile)
	cur1 = con1.cursor()
	cur1.execute(columnNamesCommand)
	columnHeaders = []
	index = 0
	for row in cur1:
		columnHeaders.append(row[1])
		index+=1
	cur1.close()
	con1.close()
	for y in range(index):
		columnHeaders[y] = str(columnHeaders[y])
#	print columnHeaders
	for y in range(len(ColumnHeadersTable)):
		found = 0
		for z in range(len(columnHeaders)):
			if ColumnHeadersTable[y]==columnHeaders[z]:
				found = 1
		if found==0:
			updateColumn = "ALTER TABLE IMAGE_TEST ADD COLUMN " + ColumnHeadersTable[y] + " FLOAT DEFAULT '0.0'";
			print updateColumn
			con2 = sqlite3.connect(databaseFile)
			cur2 = con2.cursor()
			cur2.execute(updateColumn)
			cur2.close()
			con2.close()

	con22 = sqlite3.connect(databaseFile)
	cur22 = con22.cursor()
	cur22.execute("select * from IMAGE_TEST")
	cur22.close()
	con22.commit()
	con22.close()

	#Check if the image already exists on the database
	ImagePath = inputFiles[x][:inputFiles[x].rfind("/")]
	con3 = sqlite3.connect(databaseFile)
	t  = (ImagePath)
	tt = (inputImages[x])
	print t, tt
	try:
		with con3:
			con3.execute("insert into IMAGE(IMG_NAME, IMG_LOCATION) values (?, ?)", (t, tt))
	except sqlite3.IntegrityError:
		print "ERROR"
	con3.commit()
	con3.close()
	con33 = sqlite3.connect(databaseFile)
	cur3  = con33.cursor()
	ImagePathss = ImagePath + '%'
	cur3.execute("select * from IMAGE where IMG_LOCATION like '%s';"% ImagePathss)
	IDonDB    = []
	for row in cur3:
		IDonDB.append(row[0])
	cur3.close()
	con33.close()
	ttt = []
	for y in range(len(IDonDB)):
		con4 = sqlite3.connect(databaseFile)
		cur4 = con4.cursor()
		ttt = (IDonDB[y],)
		cur4.execute("delete from IMAGE_TEST where IMG_ID=?",ttt)
		cur4.close()
		con4.commit()
		con4.close()
		print "HERE",IDonDB[y]
	con44 = sqlite3.connect(databaseFile)
	cur44  = con44.cursor()
	cur44.execute("select * from IMAGE_TEST where IMG_ID=?",ttt)
	for row in cur44:
		print row

	#Insert path into master table and get image ID if needed
	if len(IDonDB):
		CurImgID = IDonDB[0]
		InsertRowHeaders = "IMG_ID, CELL_ID" 
		for y in range(len(ColumnHeadersTable)):
			if y:
				InsertRowHeaders += ", " + ColumnHeadersTable[y]
		UpdateStr = "insert into IMAGE_TEST("+InsertRowHeaders+") values(?"
		for y in range(len(ColumnHeadersTable)):
			UpdateStr += ", ?"
		UpdateStr += ")"
		print UpdateStr
		con5 = sqlite3.connect(databaseFile)
		cur5 = con5.cursor()
		for y in range(len(tableEntries)):
			CurrentRowFeats = []
			for z in range(len(tableEntries[y])):
				if z==0:
					CellID = tableEntries[y][z]
					CurrentRowFeats.append(CurImgID)
					CurrentRowFeats.append(int(float(CellID)))
				else:
					CurrentRowFeats.append(float(tableEntries[y][z]))
			if y==10:
				print CurrentRowFeats
			cur5.execute(UpdateStr,CurrentRowFeats)
		cur5.close()
		con5.commit()
		con5.close()
		print InsertRowHeaders
		print CurImgID
