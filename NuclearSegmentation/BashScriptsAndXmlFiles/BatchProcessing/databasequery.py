#!/bin/python

import sqlite3
import csv

con = sqlite3.connect("C:/farsight_program/NE.s3db")

images_path = "C:\ccRCC\ccRCC Test Set\Negative"

csv_filename = 'test_op1.csv'

cur = con.cursor()

images_path += '%'

print images_path

print("select * from IMAGE where IMG_LOCATION like '%s%';",images_path)

image_paths = []
image_ids = []
cur.execute("select * from IMAGE where IMG_LOCATION like '%s';"% images_path)
for row in cur:
	image_ids.append( row[0] )
	image_paths.append( row[1] )

for index, object in enumerate(image_paths):
	print image_paths[index]

print ("TOTAL_CELLS")
image_ids1 = []
cell_totals = []
cur.execute("select IMG_ID,count(*) as\"TOTAL_CELLS\" from IMAGE_TEST where IMG_ID in ( select IMG_ID from IMAGE where IMG_LOCATION like '%s' and TRAIN_TUMOR=-1 and TRAIN_ENDOTHELIAL=-1 and TRAIN_ANALYTE=-1 ) group by IMG_ID;"% images_path)
for row in cur:
	image_ids1.append( row[0] )
	cell_totals.append( row[1] )

print ("ENDOTHELIAL_CELLS")
image_ids2 = []
endothelial_totals = []
cur.execute("select IMG_ID,count(*) from (select * from IMAGE_TEST where IMG_ID in ( select IMG_ID from IMAGE where IMG_LOCATION like '%s' ) and PREDICTION_ENDOTHELIAL=1 and TRAIN_TUMOR=-1 and TRAIN_ENDOTHELIAL=-1 and TRAIN_ANALYTE=-1) group by IMG_ID;"% images_path)
for row in cur:
	image_ids2.append( row[0] )
	endothelial_totals.append( row[1] )

print ("ENDOTHELIAL_POSITIVE_CELLS")
image_ids3 = []
endothelial_positive_totals = []
cur.execute("select IMG_ID,count(*) from (select * from IMAGE_TEST where IMG_ID in ( select IMG_ID from IMAGE where IMG_LOCATION like '%s' ) and PREDICTION_ENDOTHELIAL=1 and PREDICTION_ANALYTE=1 and TRAIN_TUMOR=-1 and TRAIN_ENDOTHELIAL=-1 and TRAIN_ANALYTE=-1) group by IMG_ID;"% images_path)
for row in cur:
	image_ids3.append( row[0] )
	endothelial_positive_totals.append( row[1] )

print( "TUMOR_CELLS" )
image_ids4 = []
tumor_totals = []
cur.execute("select IMG_ID,count(*) from (select * from IMAGE_TEST where IMG_ID in ( select IMG_ID from IMAGE where IMG_LOCATION like '%s' ) and PREDICTION_TUMOR=1 and TRAIN_TUMOR=-1 and TRAIN_ENDOTHELIAL=-1 and TRAIN_ANALYTE=-1) group by IMG_ID;"% images_path)
for row in cur:
	image_ids4.append( row[0] )
	tumor_totals.append( row[1] )

print ("TUMOR_POSITIVE_CELLS")
cur.execute("select IMG_ID,count(*) from (select * from IMAGE_TEST where IMG_ID in ( select IMG_ID from IMAGE where IMG_LOCATION like '%s' ) and PREDICTION_TUMOR=1 and PREDICTION_ANALYTE=1  and TRAIN_TUMOR=-1 and TRAIN_ENDOTHELIAL=-1 and TRAIN_ANALYTE=-1) group by IMG_ID;"% images_path)
image_ids5 = []
tumor_positive_totals = []
for row in cur:
	image_ids5.append( row[0] )
	tumor_positive_totals.append( row[1] )

spamWriter = csv.writer(open(csv_filename, 'wb'), delimiter=',',quotechar='"', quoting=csv.QUOTE_MINIMAL)

spamWriter.writerow(['IMAGE_ID']+['IMAGE_PATH']+['TOTAL_CELLS']+['ENDOTHELIAL_CELLS']+['ENDOTHELIAL_POSITIVE_CELLS']+['TUMOR_CELLS']+['TUMOR_POSITIVE_CELLS'])

for k in range(len(image_ids)):
	ids = image_ids[k]
	temp1 = temp2 = temp3 = temp4 = temp5 = 0
	for i in range(len(image_ids1)):
		if image_ids1[i] == ids:
			temp1 = cell_totals[i]
	for i in range(len(image_ids2)):
		if image_ids2[i] == ids:
			temp2 = endothelial_totals[i]
	for i in range(len(image_ids3)):
		if image_ids3[i] == ids:
			temp3 = endothelial_positive_totals[i]
	for i in range(len(image_ids4)):
		if image_ids4[i] == ids:
			temp4 = tumor_totals[i]
	for i in range(len(image_ids5)):
		if image_ids5[i] == ids:
			temp5 = tumor_positive_totals[i]
	print ids, temp1, temp2, temp3, temp4, temp5
	spamWriter.writerow([ids]+[image_paths[k]]+[temp1]+[temp2]+[temp3]+[temp4]+[temp5])

