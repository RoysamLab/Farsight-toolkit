#!/usr/bin/env python
# -*- coding: utf-8 -*-

import shutil
import fnmatch
import os
import subprocess
import sys

MAIN_DATA_FOLDER = '/data/nicolas/dataNew'

SERVER = 'none'
p = os.uname()
if( p[1] == 'Farsight-05.EE.UH.EDU' ):
	SERVER = 'far05'
if( p[1] == 'Farsight-04.EE.UH.EDU' ):
	SERVER = 'far04'

#DATA_FOLDER_ALL = ['/0131_test','/0131_test2'] # For testing dont forget the xTile params
#DATA_FOLDER_ALL = ['/0131_test'] # For testing dont forget the xTile params
#DATA_FOLDER_ALL = ['/0113_NRRD_CROPPED']

#DATA_FOLDER_ALL = ['/0113_NRRD_CROPPED','/0117_NRRD_CROPPED','/0120_NRRD_CROPPED','/0123_NRRD_CROPPED','/0128_NRRD_CROPPED','/0131_NRRD_CROPPED','/0323_NRRD_CROPPED','/0405_NRRD_CROPPED','/0409_NRRD_CROPPED','/0410_NRRD_CROPPED','/0412_NRRD_CROPPED','/1206_NRRD_CROPPED']
#DATA_FOLDER_ALL = ['/0113_NRRD_CROPPED','/0117_NRRD_CROPPED','/0120_NRRD_CROPPED','/0123_NRRD_CROPPED','/0128_NRRD_CROPPED','/0131_NRRD_CROPPED']
#DATA_FOLDER_ALL = ['/0323_NRRD_CROPPED','/0405_NRRD_CROPPED','/0409_NRRD_CROPPED','/0410_NRRD_CROPPED','/0412_NRRD_CROPPED','/1206_NRRD_CROPPED']

#DATA_FOLDER_ALL = ['/0113_NRRD','/0117_NRRD','/0120_NRRD','/0123_NRRD','/0128_NRRD','/0131_NRRD','/0323_NRRD','/0405_NRRD','/0409_NRRD','/0410_NRRD','/0412_NRRD','/1206_NRRD']
TEST_RUN = 0
if len(sys.argv) == 2:
	if(sys.argv[1] =='TEST' ):
		TEST_RUN = 1
	else:
		TEST_RUN = 0

DATA_FOLDER_ALL = ['/0131_test']
if( SERVER == 'far04' ):
	REMOVE_MONTAGES = 1
	#DATA_FOLDER_ALL = ['/0113_NRRD_CROPPED','/0117_NRRD_CROPPED','/0120_NRRD_CROPPED','/0123_NRRD_CROPPED','/0128_NRRD_CROPPED','/0131_NRRD_CROPPED','/0323_NRRD_CROPPED','/0405_NRRD_CROPPED','/0409_NRRD_CROPPED','/0410_NRRD_CROPPED','/0412_NRRD_CROPPED','/1206_NRRD_CROPPED']
	#DATA_FOLDER_ALL = ['/0113_NRRD','/0117_NRRD','/0120_NRRD','/0123_NRRD','/0128_NRRD','/0131_NRRD','/0323_NRRD','/0405_NRRD','/0409_NRRD','/0410_NRRD','/0412_NRRD','/1206_NRRD']
	#DATA_FOLDER_ALL = ['/0323_NRRD','/0405_NRRD','/0409_NRRD','/0410_NRRD','/0412_NRRD','/1206_NRRD']
	DATA_FOLDER_ALL = ['/0113_NRRD','/0117_NRRD','/0120_NRRD','/0123_NRRD','/0128_NRRD','/0131_NRRD']
	#DATA_FOLDER_ALL = ['/0117_NRRD','/0120_NRRD','/0123_NRRD','/0128_NRRD','/0131_NRRD']
	#DATA_FOLDER_ALL = ['/0120_NRRD','/0123_NRRD','/0128_NRRD','/0131_NRRD']
	#DATA_FOLDER_ALL = ['/0113_NRRD']
	#DATA_FOLDER_ALL = ['/0131_test']
	#DATA_FOLDER_ALL = ['/0131_test2']
	#DATA_FOLDER_ALL = ['/0131_test3']
	#DATA_FOLDER_ALL = ['/0131_test4']
	#DATA_FOLDER_ALL = ['/0120_NRRD']
if( SERVER == 'far05' ):
	REMOVE_MONTAGES = 1
	#DATA_FOLDER_ALL = ['/0113_NRRD','/0117_NRRD','/0120_NRRD','/0123_NRRD','/0128_NRRD','/0131_NRRD','/0323_NRRD','/0405_NRRD','/0409_NRRD','/0410_NRRD','/0412_NRRD','/1206_NRRD']
	DATA_FOLDER_ALL = ['/0323_NRRD','/0405_NRRD','/0409_NRRD','/0410_NRRD','/0412_NRRD','/1206_NRRD']
	#DATA_FOLDER_ALL = ['/0405_NRRD','/0409_NRRD','/0410_NRRD','/0412_NRRD','/1206_NRRD']
	#DATA_FOLDER_ALL = ['/0131_test']
	#DATA_FOLDER_ALL = ['/0131_test2']
	#DATA_FOLDER_ALL = ['/0131_test3']
	#DATA_FOLDER_ALL = ['/0131_test4']
	runMake = 0
if( TEST_RUN==1 ):
	DATA_FOLDER_ALL = ['/0131_test4']
	print "THIS IS A TEST"
	print "THIS IS A TEST"


for DATA_FOLDER in DATA_FOLDER_ALL:

	LOCAL_FOLDER = MAIN_DATA_FOLDER+DATA_FOLDER

	if os.path.isdir(LOCAL_FOLDER):
		print 'Folder '+LOCAL_FOLDER+' already exists'
		print 'If you want to remove uncoment this'
		print 'Will be removed'
		shutil.rmtree(LOCAL_FOLDER)
	else:
		print 'Folder '+LOCAL_FOLDER+' will be created'
		os.makedirs(LOCAL_FOLDER)
		LOCAL_FOLDER_PARAMS = LOCAL_FOLDER+'/Parameters'
		os.makedirs(LOCAL_FOLDER_PARAMS)

		BACSUBFILE = LOCAL_FOLDER_PARAMS+'/fijiMacroBackSubst.ijm'
		BACSUBFILE = open(BACSUBFILE, 'w')
		BACSUBFILE.write('file=getArgument();\n')
		BACSUBFILE.write('run("Nrrd ...", "load="+file+".nrrd");\n')
		BACSUBFILE.write('setSlice(5);\n')
		BACSUBFILE.write('run("Subtract Background...", "rolling=50 stack");\n')
		BACSUBFILE.write('run("Nrrd ... ", "nrrd="+file+"_BS.nrrd");\n')
		BACSUBFILE.write('close();\n')
		BACSUBFILE.close()


		MODEL_FILE = LOCAL_FOLDER_PARAMS+'/MG_model.txt'
		MODEL_FILE_NAME = MODEL_FILE
		MODEL_FILE = open(MODEL_FILE, 'w')
		MODEL_FILE.write('bias	volume	mean	maximum	GFP_TOTAL	GFP_AVG\n')
		MODEL_FILE.write('0	31136.5	6.70456	26.5247	224035	8.14896\n')
		MODEL_FILE.write('0	30116.4	20.5461	40.1787	65737.4	8.02105\n')
		MODEL_FILE.write('-0.162498	-0.0615876	1.25E-05	9.09E-06	0.331097	-1.21E-05\n')
		MODEL_FILE.write('0.162498	0.0615875	-1.25E-05	-9.09E-06	-0.331097	1.21E-05\n')
		MODEL_FILE.close()

		CURVE_FILE = LOCAL_FOLDER_PARAMS+'/options_curvelets'
		CURVE_FILE = open(CURVE_FILE, 'w')
		CURVE_FILE.write('-nban\n')
		CURVE_FILE.write('-nsigmas_coarse 2.2\n')
		CURVE_FILE.write('-nsigmas_fine 2.5\n')
		CURVE_FILE.write('-neighb_weight 0.5\n')
		CURVE_FILE.write('-tuning_neighb 0.6\n')
		CURVE_FILE.write('-sigma_ratio 0.01\n')
		CURVE_FILE.write('-num_threads 78\n') # CARAFULLLLLLLLLLLLLLLLLLLLLLLLLLLLL
		CURVE_FILE.write('-tile_size 2048\n')
		CURVE_FILE.write('-border 100\n')
		CURVE_FILE.close()

		PROJE_FILE = LOCAL_FOLDER_PARAMS+'/ProjectDefinition.xml'
		PROJE_FILE = open(PROJE_FILE, 'w')
		PROJE_FILE.write('<ProjectDefinition name="Project DARPA">\n')
		PROJE_FILE.write('    <Inputs>\n')
		PROJE_FILE.write('        <channel number="0" name="dapi" type="NUCLEAR" />\n')
		PROJE_FILE.write('        <channel number="1" name="gfp" type="GFP_MARKER" />\n')
		PROJE_FILE.write('    </Inputs>\n')
		PROJE_FILE.write('    <Pipeline>\n')
		PROJE_FILE.write('        <step name="FEATURE_COMPUTATION" />\n')
		PROJE_FILE.write('        <step name="RAW_ASSOCIATIONS" />\n')
		PROJE_FILE.write('        <step name="CLASSIFY_MCLR" />\n')
		PROJE_FILE.write('    </Pipeline>\n')
		PROJE_FILE.write('    <NuclearSegmentationParameters>\n')
		PROJE_FILE.write('        <parameter name="high_sensitivity" value="0.00" />\n')
		PROJE_FILE.write('        <parameter name="LoG_size" value="30.00" />\n')
		PROJE_FILE.write('        <parameter name="min_scale" value="7.00" />\n')
		PROJE_FILE.write('        <parameter name="max_scale" value="11.00" />\n')
		PROJE_FILE.write('        <parameter name="xy_clustering_res" value="3.00" />\n')
		PROJE_FILE.write('        <parameter name="z_clustering_res" value="2.00" />\n')
		PROJE_FILE.write('        <parameter name="finalize_segmentation" value="0.00" />\n')
		PROJE_FILE.write('        <parameter name="sampling_ratio_XY_to_Z" value="2.00" />\n')
		PROJE_FILE.write('        <parameter name="Use_Distance_Map" value="1.00" />\n')
		PROJE_FILE.write('        <parameter name="refinement_range" value="6.00" />\n')
		PROJE_FILE.write('        <parameter name="min_object_size" value="100.00" />\n')
		PROJE_FILE.write('    </NuclearSegmentationParameters>\n')
		PROJE_FILE.write('    <CytoplasmSegmentationParameters>\n')
		PROJE_FILE.write('        <parameter name="draw_real_boundaries" value="1.00" />\n')
		PROJE_FILE.write('        <parameter name="remove_stromal_cell_boundaries" value="0.00" />\n')
		PROJE_FILE.write('        <parameter name="draw_synthetic_boundaries" value="0.00" />\n')
		PROJE_FILE.write('        <parameter name="radius_of_synthetic_boundaries" value="0.00" />\n')
		PROJE_FILE.write('        <parameter name="number_of_levels" value="1.00" />\n')
		PROJE_FILE.write('        <parameter name="number_of_levels_in_foreground" value="1.00" />\n')
		PROJE_FILE.write('    </CytoplasmSegmentationParameters>\n')
		PROJE_FILE.write('    <AssociationRules>\n')
		PROJE_FILE.write('        <AssociationRule Name="GFP_TOTAL" SegmentationSource="NUCLEAR" Target_Image="gfp" Outside_Distance="3" Inside_Distance="8" Use_Whole_Object="True" Use_Background_Subtraction="True" Use_MultiLevel_Thresholding="False" Number_Of_Thresholds="1" Number_Included_In_Foreground="1" Association_Type="TOTAL" />\n')
		PROJE_FILE.write('        <AssociationRule Name="GFP_AVG" SegmentationSource="NUCLEAR" Target_Image="gfp" Outside_Distance="3" Inside_Distance="8" Use_Whole_Object="True" Use_Background_Subtraction="True" Use_MultiLevel_Thresholding="False" Number_Of_Thresholds="1" Number_Included_In_Foreground="1" Association_Type="AVERAGE" />\n')
		PROJE_FILE.write('    </AssociationRules>\n')
		PROJE_FILE.write('    <Classification_MCLR_Rules>\n')
		PROJE_FILE.write('        <Classification_MCLR_Rule ClassColName="mg" ConfThreshold="0.50" TrainingFileName="'+MODEL_FILE_NAME+'" />\n')
		PROJE_FILE.write('    </Classification_MCLR_Rules>\n')
		PROJE_FILE.write('</ProjectDefinition>\n')
		PROJE_FILE.close()

		SEGME_FILE = LOCAL_FOLDER_PARAMS+'/Seg_Params.ini'
		SEGME_FILE = open(SEGME_FILE, 'w')
		SEGME_FILE.write('high_sensitivity	:	0\n')
		SEGME_FILE.write('LoG_size		:	30\n')
		SEGME_FILE.write('min_scale		:	7\n')
		SEGME_FILE.write('max_scale		:	11\n')
		SEGME_FILE.write('xy_clustering_res	:	3\n')
		SEGME_FILE.write('z_clustering_res	:	2\n')
		SEGME_FILE.write('finalize_segmentation	:	0\n')
		SEGME_FILE.write('sampling_ratio_XY_to_Z	:	2\n')
		SEGME_FILE.write('Use_Distance_Map	:	1\n')
		SEGME_FILE.write('min_object_size		:	100\n')
		SEGME_FILE.close()

		TRACE_FILE = LOCAL_FOLDER_PARAMS+'/options_mnt'
		TRACE_FILE = open(TRACE_FILE, 'w')
		TRACE_FILE.write('-intensity_threshold 0.003\n')
		TRACE_FILE.write('-contrast_threshold 0.003\n')
		TRACE_FILE.write('-cost_threshold 300\n')
		TRACE_FILE.write('-debris_threshold 0.8\n')
		TRACE_FILE.write('-offshoot 15\n')
		TRACE_FILE.write('-device 0\n')
		TRACE_FILE.close()
