#!/bin/python

import sqlite3
import csv
import math

con = sqlite3.connect("/home/gramak/sql_test/NE.s3db")

images_path = "C:/STS_Farsight_analysis_CD34/pERK"

csv_filename  = 'summary_output'
csv_extension = 'csv'

sd_constant     = 2
recursion_level = 1 #write code for more levels

cur = con.cursor()

images_path += '%'

#print images_path

#print("select * from IMAGE where IMG_LOCATION like '%s%';",images_path)

image_paths = []
image_ids   = []
cur.execute("select * from IMAGE where IMG_LOCATION like '%s';"% images_path)
#print cur
for row in cur:
	image_ids.append( row[0] )
#	print row[0]
	image_paths.append( row[1] )
#	print row[1]

num_done       = 0
query_str_done = []
for index, object in enumerate(image_paths):
#	print image_paths[index]
	st_ring   = image_paths[index]
	st_r_ng   = st_ring[(len(images_path)+1):len(image_paths[index])]
	end       = st_r_ng.find('/')
	elm_sz    = len(st_r_ng)-end
	#Extend for recursion levels here ***********
	query_str = st_ring[0:(len(st_ring)-elm_sz)]
	exte_str  = st_ring[len(images_path):(len(st_ring)-elm_sz)]
#	print query_str
#	print exte_str
	already_d = 0
	for index1, object in enumerate(query_str_done):
		if query_str_done[index1]==query_str:
			already_d = 1
	if already_d==0:
		query_str_done.append( query_str )
		query_str += '%'
		image_pathas = []
		image_idas   = []
		cur.execute("select * from IMAGE where IMG_LOCATION like '%s';"% query_str)
		for row in cur:
			image_idas.append( row[0] )
			image_pathas.append( row[1] )
		#print ("TOTAL")
		image_ids1 = []
		cell_totals = []
		cur.execute("select IMG_ID,count(*) as\"TOTAL_CELLS\" from IMAGE_TEST where IMG_ID in ( select IMG_ID from IMAGE where IMG_LOCATION like '%s' ) group by IMG_ID;"% query_str)
		for row in cur:
			image_ids1.append( row[0] )
			cell_totals.append( row[1] )

		#print ("EC")
		image_ids2 = []
		endothelial_totals = []
		cur.execute("select IMG_ID,count(*) from (select * from IMAGE_TEST where IMG_ID in ( select IMG_ID from IMAGE where IMG_LOCATION like '%s' ) and PREDICTION_ACTIVE=1) group by IMG_ID;"% query_str)
		for row in cur:
			image_ids2.append( row[0] )
			endothelial_totals.append( row[1] )

		#print ("EC+")
		image_ids3 = []
		endothelial_positive_totals = []
		cur.execute("select IMG_ID,count(*) from (select * from IMAGE_TEST where IMG_ID in ( select IMG_ID from IMAGE where IMG_LOCATION like '%s' ) and PREDICTION_ACTIVE=1 and PREDICTION_ACTIVE_ANALYTE=1) group by IMG_ID;"% query_str)
		for row in cur:
			image_ids3.append( row[0] )
			endothelial_positive_totals.append( row[1] )
		n = 0.0
		n += len(cell_totals)
		mean_cell_totals                 = sum(cell_totals) / n
		sd_cell_totals                   = math.sqrt(sum((x-mean_cell_totals)**2 for x in cell_totals) / n)
		mean_endothelial_totals          = sum(endothelial_totals) / n
		sd_endothelial_totals            = math.sqrt(sum((x-mean_endothelial_totals)**2 for x in endothelial_totals) / n)
		mean_endothelial_positive_totals = sum(endothelial_positive_totals) / n
		sd_endothelial_positive_totals   = math.sqrt(sum((x-mean_endothelial_positive_totals)**2 for x in endothelial_positive_totals) / n)
		print "ID","TOTAL","EC","EC+"
		csv_full   = csv_filename + '_' + exte_str + '.' + csv_extension
		out_fil    = open(csv_full,"wb")
		spamWriter = csv.writer(out_fil, delimiter=',',quotechar='"', quoting=csv.QUOTE_MINIMAL)
		spamWriter.writerow(['ID']+['PATH']+['TOTAL']+['EC']+['EC+'])

		for k in range(len(image_idas)):
			ids = image_idas[k]
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
			print ids, temp1, temp2, temp3
			if temp1<(mean_cell_totals-sd_constant*sd_cell_totals) or \
			   temp1>(mean_cell_totals+sd_constant*sd_cell_totals) or \
			   temp2<(mean_endothelial_totals-sd_constant*sd_endothelial_totals) or \
			   temp2>(mean_endothelial_totals+sd_constant*sd_endothelial_totals) or \
			   temp3<(mean_endothelial_positive_totals-sd_constant*sd_endothelial_positive_totals) or \
			   temp3>(mean_endothelial_positive_totals+sd_constant*sd_endothelial_positive_totals):
				spamWriter.writerow([(k+1)]+[image_pathas[k]]+[temp1]+[temp2]+[temp3]+['*'])
			else:
				spamWriter.writerow([(k+1)]+[image_pathas[k]]+[temp1]+[temp2]+[temp3])
	
		print "sum", sum(cell_totals), sum(endothelial_totals), sum(endothelial_positive_totals)
		print "mean", mean_cell_totals, mean_endothelial_totals, mean_endothelial_positive_totals
		print "std_dev", sd_cell_totals, sd_endothelial_totals, sd_endothelial_positive_totals
		spamWriter.writerow( ['sum']+['']+[sum(cell_totals)]+[sum(endothelial_totals)]+[sum(endothelial_positive_totals)] )
		spamWriter.writerow( ['mean']+['']+[mean_cell_totals]+[mean_endothelial_totals]+[mean_endothelial_positive_totals] )
		spamWriter.writerow( ['std_dev']+['']+[sd_cell_totals]+[sd_endothelial_totals]+[sd_endothelial_positive_totals] )
		spamWriter.writerow( ['Percentage_EC_to_Total'] + [(100.0*sum(endothelial_totals)/sum(cell_totals))] )
		spamWriter.writerow( ['Percentage_EC+_to_EC']   + [(100.0*sum(endothelial_positive_totals)/sum(endothelial_totals))] )
		out_fil.close()
cur.close()

