#!/usr/bin/env python
#############################################################################
# A SCRIPT TO RUN THROUGH THE 5 LABEL TEST IMAGE FROM BEGINNING TO END:
################################################################################
import subprocess
import sys
import os

from farsightutils import *

data = ['1G_EBA_crop1.tif', '1G_EBA_crop2.tif', \
        '1G_GFAP_crop1.tif', '1G_GFAP_crop2.tif', \
        '1G_Iba1_crop1.tif', '1G_Iba1_crop2.tif', \
        '1G_Neuro_crop1.tif', '1G_Neuro_crop2.tif', \
        '1G_Nuclei_crop1.tif', '1G_Nuclei_crop2.tif', \
        'NucleiSegmentationParams.ini', \
        '1G_crop1_assoc_def.xml', '1G_crop2_assoc_def.xml', \
        '1G_GFAP_crop1_TracingParameterValues.xml', '1G_GFAP_crop2_TracingParameterValues.xml',
        '1G_GFAP_crop1.pic','1G_GFAP_crop2.pic']

data_dir = 'C:' + os.sep + 'BADRI_TEST'

###############################################################################

def main():
  print 'You have started the Farsight 5-Label Image Processing Demo'
  yn = raw_input('Would you like to continue with this demonstration (y/n)?:  ')

  if yn!='y' and yn!='Y':
    return;

  print('\nAssuming Test Data is in ' + data_dir)
  yn = raw_input('Is this correct (y/n)?:  ')
  if yn!='y' and yn!='Y':
    print('\nPlease Put the images in ' + data_dir)
    return;
  
  os.chdir(data_dir)
  print('\nChanged Working Directory to ' + data_dir)

  #NOW I WILL SEE IF ALL THE RIGHT IMAGES ARE HERE:
  files = os.listdir(os.getcwd())
  for d in data:
    #See if d is in the list:
    found = False;
    for f in files:
      if f.find(d) != -1: #If I find it break out of loop
        found = True;
        break;
    if not found:   #If I didn't find it print error and exit function
      print('\nCould not find required file: ' + d)
      return;

  print('Found All Required Files')

  #Now BEGIN the Wizard that steps through the modules:

  print('\nLETS BEGIN IMAGE PROCESSING!!!\n')
  image_num = raw_input("Which image would you like to work with (1/2)?")

  if image_num=='1':
    print("GREAT CHOICE")
    nuc_image = data[8]
    ves_image = data[0]
    trace_image = data[15]
    ass_defs = data[11]
    trace_params = data[13]
  elif image_num=='2':
    print("GREAT CHOICE")
    nuc_image = data[9]
    ves_image = data[1]
    trace_image= data[16]
    ass_defs = data[12]
    trace_params = data[14]
  else:
    print("I DON'T KNOW THAT NUMBER, ABORTING")
    return


  #DO THE NUCLEAR SEGMENTATION IF DESIRED:
  (begin,end) = os.path.splitext(nuc_image)
  nuc_result = begin + '_label' + end

  yn = raw_input("\nWould you like to segment the Nuclei in this image (y/n)?")
  if yn=='y' or yn=='Y':
    #RUN SEGMENTATION:
    print("\nSTARTING NUCLEAR SEGMENTATION...")
    subprocess.call(['segment_nuclei.exe', nuc_image, nuc_result, data[10]])
    print("DONE\nCOMPUTING ASSOCIATIVE FEATURES...")
    subprocess.call(['compute_associative_measures', ass_defs])
    print("DONE\nCOMPUTING FEATURES AND CREATING XML RESULT FILE...")
    subprocess.call(['compute_nuclei_features', os.getcwd(), nuc_image, nuc_result])
    print("DONE\nYOU MAY OPEN THE RESULT FROM THE FARSIGHT TOOLBAR")
  elif yn!='n' and yn!='N':
    #UNKNOWN COMMAND ABORT
    print("/n EXITING DEMO")
    return

  
  #DO THE VESSEL SEGMENTATION IF DESIRED:
  (begin,end) = os.path.splitext(ves_image)
  ves_result = begin + '_vessels' + end

  yn = raw_input("\nWould you like to segment the Vessels in this image (y/n)?")
  if yn=='y' or yn=='Y':
    #RUN SEGMENTATION:
    print("\nSTARTING VESSEL SEGMENTATION...")
    subprocess.call(['vessel_segmentation.exe', ves_image, ves_result])
    print("DONE")
  elif yn!='n' and yn!='N':
    #UNKNOWN COMMAND ABORT
    print("/n EXITING DEMO")
    return

  #LOOK FOR RESULT, IF FOUND OFFER TO DO VISUALIZATION
  npts = ves_image + '.npts'
  found = False;
  files = os.listdir(os.getcwd())
  for f in files:
    if f.find(npts) != -1: #If I find it break out of loop
      found = True;
      break;
    
  if found:
    yn=raw_input("\nWould you like to create 3D visualization of vessels (y/n)?")
    if yn=='y' or yn=='Y':
      #RUN VISUALIZATION:
      print("\nSTARTING VISUALIZATION...")
      subprocess.call(['visualize', npts, 'blank'])
      print("DONE")
  

  #DO THE TRACING IF DESIRED:
  yn = raw_input("\nWould you like to Trace the GFAP Channel (y/n)?")
  if yn=='y' or yn=='Y':
    #RUN TRACING:
    print("\nSTARTING TRACE SCRIPT...")
    subprocess.call(["RPITrace3D.exe", trace_params])
    print("\n...DONE")
  elif yn!='n' and yn!='N':
    #UNKNOWN COMMAND ABORT
    print("/n EXITING DEMO")
    return

  #OFFER TO DO PAIRWISE REGISTRATION OF TEST IMAGES:
  yn = raw_input("\nWould you like to Register the Nuclei Channel (y/n)?")
  if yn=='y' or yn=='Y':
    #First create temp file images listed:
    temp = open('temp_list.txt','w+')
    temp.write(data[8] + " " + data[9] + "\n")
    temp.close()
    #RUN registration:
    print("\nSTARTING REGISTRATION SCRIPT...")
    from register_pairs import register
    register(['','temp_list.txt'])
    print("\n...DONE")
    #delete temporary file:
    os.remove('temp_list.txt')
  elif yn!='n' and yn!='N':
    #UNKNOWN COMMAND ABORT
    print("/n EXITING DEMO")
    return

  #LOOK FOR RESULT, IF FOUND OFFER TO DO VISUALIZATION
  (begin,end) = os.path.splitext(data[8])
  m_result = "montage_" + begin + ".xml"
  found = False;
  files = os.listdir(os.getcwd())
  for f in files:
    if f.find(m_result) != -1: #If I find it break out of loop
      found = True;
      break;
    
  if found:
    print("\nYou may VISUALIZE " + m_result + " in Montage Navigator")
    yn=raw_input("Would you like to load this application now (y/n)?")
    if yn=='y' or yn=='Y':
      #RUN VISUALIZATION:
      subprocess.call(['MontageNavigator.exe'])


  #DEMO CONCLUSION
  print("\nTHIS CONCLUDES OUR DEMO. GOODBYE")

################################################################################
#if __name__ == "__main__":
main()
