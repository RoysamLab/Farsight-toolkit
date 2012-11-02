#!/usr/bin/python
# -*- coding: utf-8 -*-

import getopt,sys, glob, os, shutil, platform, subprocess
import runParallel

FLAGMOSAIC = "0"
if sys.argv[1] == "1":
  print "IT WILL DO THE REGISTRATION AND MOSAIIC ALSO"
  FLAGMOSAIC = "1"
elif sys.argv[1] == "2":
  print "IT WILL DO ONLY THE MOSAIC"
  FLAGMOSAIC = "2"
else:
  print "IT WILL NOT DO THE MOSAIC ONLY THE PREPROCESING"

# Path for the raw data and the destination data
DATASET_RAW = "1005"
DATASET_DEST = DATASET_RAW+"_NRRD"
#GLOBAL_RAW = "/space1/nicolas/FSdata/data/DARPA_RAW/"
#GLOBAL_DEST = "/space1/nicolas/FSdata/data/DARPA_MOSAICS/"

#GLOBAL_RAW = "/data/nicolas/FSdata/DARPA_RAW/"
#GLOBAL_DEST = "/data/nicolas/FSdata/DARPA_MOSAIC/"

# PAth where the data will be copied from (please use the ssh shortcut far-01 to refer to the farsight-01 server, or modified the script accordingly)
GLOBAL_RAW = "/FSdata/data/DARPA_RAW/"
GLOBAL_DEST = "/FSdata/data/DARPA_MOSAIC/"

GLOBAL_RAW_PATH = GLOBAL_RAW+DATASET_RAW+"/"
GLOBAL_DEST_PATH = GLOBAL_DEST+DATASET_DEST+"/"

# Local path
LOCAL_PATH = "/data/prathamesh/data/"


# Path where the scripts and exe are located
SCRIPTS_PATH = "/data/prathamesh/data/montaging_scripts/"
FARSIGHT_EXE = "/data/prathamesh/data/farsight_exes"

# Variables of the image, size, first tile, and most important of all, NUMCOLS
# ORDER: 1: Top left, 2: Bottom right, 3: Top right, 4: Bottom left
ORDER = "3"
NUMCOLS = "21"
XSIZE = "1004"
YSIZE = "1002"
ZSIZE = "350"
ANCHORIMAGE = "0"

#############################################################################
# Remember to keep the standard name of the DARPA project, if you have more colores add them here
COLORS =['NIC','DAPI','Cy5','GFP','TRITC']

#############################################################################
# Create result folder
#if os.path.isdir(GLOBAL_DEST_PATH):
# #print "The path "+GLOBAL_DEST_PATH+" Already exists, will exit"
# #sys.exit(0)
# print "WARNING The path "+GLOBAL_DEST_PATH+" Already exists, will remove it"
# #shutil.rmtree(GLOBAL_DEST_PATH)
#else:
# print "making "+GLOBAL_DEST_PATH
# os.makedirs(GLOBAL_DEST_PATH)

#check if the unused folder exist, this is used as a flag to indicate that the registration has not run already, this folder will be created by the register_automagic script
#if not os.path.isfile(LOCAL_PATH+DATASET_RAW+"/aDAPI.tif") | os.path.isdir(LOCAL_PATH+DATASET_RAW+"/unused"):
if not os.path.isdir(LOCAL_PATH+DATASET_RAW+"/unused"):
  print "The unused folder does not exist, now it will exit !!"

  # Copy images
  print LOCAL_PATH+DATASET_RAW
  if not os.path.isdir(LOCAL_PATH+DATASET_RAW):
    TEMP2 = 'scp -r far-01:'+GLOBAL_RAW_PATH+" "+LOCAL_PATH
    print '\t'+TEMP2
    TEMP8 = subprocess.Popen(TEMP2, shell=True)
    TEMP8.communicate()
  else:
    print "THE DATA ALREADY EXISTS LOCALLY, IT WONT BE COPIED"
  
  # copy scripts
  TEMP2 = "cp "+SCRIPTS_PATH+"register_pairs_parallel.py "+SCRIPTS_PATH+"register_automagic.py "+LOCAL_PATH+DATASET_RAW
  TEMP3 = subprocess.Popen(TEMP2, shell=True)
  TEMP3.communicate()

  ## copy ftkMainDarpa
  #TEMP2 = "cp "+FARSIGHT_EXE+"/ftkMainDarpa "+LOCAL_PATH+DATASET_RAW
  #TEMP3 = subprocess.Popen(TEMP2, shell=True)
  #TEMP3.communicate()

  print "Copied scripts and ftkMainDarpa"

# The file prepScript is used as a flag to tell tha script that the preprocessing is already done.
  if not os.path.isfile(LOCAL_PATH+DATASET_RAW+"/prepScript.ijm"):

    print "Preprocessing....."

    # Preprocess images
    FIJISCRIPT = LOCAL_PATH+DATASET_RAW+"/prepScript.ijm"
    FIJISCRIPT = open(FIJISCRIPT,'w')
    FIJISCRIPT.write('file=getArgument();\n')
    FIJISCRIPT.write('open(file+"DAPIdsu.TIF");\n')
    FIJISCRIPT.write('setSlice(5);\n')
    FIJISCRIPT.write('run("Subtract Background...", "rolling=50 stack");\n')
    FIJISCRIPT.write('run("Median...", "radius=2 stack");\n')
    FIJISCRIPT.write('run("Gaussian Blur...", "sigma=3");\n')
    FIJISCRIPT.write('saveAs("Tiff", file+"NICdsu.TIF");\n')
    FIJISCRIPT.write('close();\n')
    FIJISCRIPT.close()

    DAPITILES=sorted(glob.glob(LOCAL_PATH+DATASET_RAW+"/*DAPIdsu.TIF"))
    #print DAPITILES
    preprocessList = []
    for dapiTile in DAPITILES:
      #print dapiTile
      temp = dapiTile.partition("DAPI")
      #print temp[0]
      TEMP2 = '/data/research/Fiji.app/fiji-linux64 --headless -macro '+LOCAL_PATH+DATASET_RAW+"/prepScript.ijm "+temp[0]+' -batch'
      TEMP3 = subprocess.Popen(TEMP2, shell=True)
      print '\tPrep '+temp[0]+" of "+DAPITILES[-1]
      preprocessList.append(TEMP2)
      TEMP3.communicate()
    #runParallel.main(preprocessList)

    # remove tif in lower case
    NICTILES=sorted(glob.glob(LOCAL_PATH+DATASET_RAW+"/*NIC*dsu.tif"))
    print DAPITILES
    for nicTile in NICTILES:
      temp = nicTile.partition(".tif")
      TEMP2 = "mv "+nicTile+" "+temp[0]+".TIF"
      print TEMP2
      TEMP3 = subprocess.Popen(TEMP2, shell=True)
      TEMP3.communicate()


  #TEMP2 = "for files in /data/nicolas/xyz/"+DATASET_RAW+"/*.tif;do mv \"$files\" `echo $files | sed 's/tif/TIF/g' `;done"
  #TEMP3 = subprocess.Popen(TEMP2, shell=True)
  #TEMP3.communicate()

  ## run testi
  #for color in COLORS:
    #TEMP2 = "./mosaicFake.py /data/nicolas/xyz/"+DATASET_RAW+"/"+" "+ORDER+" "+NUMCOLS+" "+XSIZE+" "+YSIZE+" "+ZSIZE+" "+color
    #print TEMP2
    #TEMP3 = subprocess.Popen(TEMP2, shell=True)
    #TEMP3.communicate()

  ## create the aMosaic
  #for color in COLORS:
    #TEMP2 = LOCAL_PATH+DATASET_RAW+"/ftkMainDarpa TEST_3 "+LOCAL_PATH+DATASET_RAW+"/"+"b"+color+"temp.txt"
    #print TEMP2
    #TEMP3 = subprocess.Popen(TEMP2, shell=True)
    #TEMP3.communicate()

if FLAGMOSAIC == "1":
  # Mosaic
  ACTUAL_PATH = os.getcwd()
  os.chdir(LOCAL_PATH+DATASET_RAW+"/")
  # Real Registration
  TEMP2 = "./register_automagic.py -w "+NUMCOLS+" -r NIC "+" -a "+ANCHORIMAGE
  print TEMP2
  TEMP3 = subprocess.Popen(TEMP2, shell=True)
  TEMP3.communicate()
  os.chdir(ACTUAL_PATH)

#  for color in COLORS:
    #if not "NIC" in color:
      #TEMP2 = "scp "+LOCAL_PATH+DATASET_RAW+"/"+color+"/*.nrrd far-01:"+GLOBAL_DEST_PATH
      #TEMP3 = subprocess.Popen(TEMP2, shell=True)
      #TEMP3.communicate()
      #TEMP2 = "scp "+LOCAL_PATH+DATASET_RAW+"/"+color+"/*2d_proj.png far-01:"+GLOBAL_DEST_PATH
      #TEMP3 = subprocess.Popen(TEMP2, shell=True)
      #TEMP3.communicate()

  #TEMP2 = "scp "+LOCAL_PATH+DATASET_RAW+"/unused/a*.tif far-01:"+GLOBAL_DEST_PATH
  #TEMP3 = subprocess.Popen(TEMP2, shell=True)
  #TEMP3.communicate()

if FLAGMOSAIC == "2":
  # Mosaic
        ACTUAL_PATH = os.getcwd()
        os.chdir(LOCAL_PATH+DATASET_RAW+"/")
        # Real Registration
        TEMP2 = "./register_automagic.py -w "+NUMCOLS+" -r NIC "+" -a "+ANCHORIMAGE+" -m"
        print TEMP2
        TEMP3 = subprocess.Popen(TEMP2, shell=True)
        TEMP3.communicate()
        os.chdir(ACTUAL_PATH)
 
#        for color in COLORS:
#         if not "NIC" in color:
#                        TEMP2 = "scp "+LOCAL_PATH+DATASET_RAW+"/"+color+"/*.nrrd far-01:"+GLOBAL_DEST_PATH
#                        TEMP3 = subprocess.Popen(TEMP2, shell=True)
#                        TEMP3.communicate()
#                        TEMP2 = "scp "+LOCAL_PATH+DATASET_RAW+"/"+color+"/*2d_proj.png far-01:"+GLOBAL_DEST_PATH
#                        TEMP3 = subprocess.Popen(TEMP2, shell=True)
#                        TEMP3.communicate()
 
        #TEMP2 = "scp "+LOCAL_PATH+DATASET_RAW+"/unused/a*.tif far-01:"+GLOBAL_DEST_PATH
        #TEMP3 = subprocess.Popen(TEMP2, shell=True)
        #TEMP3.communicate()











