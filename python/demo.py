#!/usr/bin/env python
#############################################################################
# A SCRIPT TO RUN THROUGH THE 5 LABEL TEST IMAGE FROM BEGINNING TO END:
################################################################################
import subprocess
import sys
import os

from farsightutils import *

data_dir = 'C:' + os.sep + 'BADRI_TEST'

#parts of the original image filenames:
orig_base_name = '100upoint5%25hippo25x1unmixed'
orig_ext = '.lsm'
split_ext = '.pic'

#parts of the cropped Image Filenames:
crop_base_name = 'NM'
crop_ext = '.pic'
crop_id = '_crop';
num_crop = 2;

#channel_names:
nissl_id = '_Nissl'
nuc_id = '_Nuc'
gfap_id = '_GFAP'
iba1_id = '_Iba1'
eba_id = '_EBA'

#Parameter Filenames:
seg_params = 'NucleiSegmentationParams.ini'
ass_def_id = '_AssociationDefs'
trace_params_id = '_TracingParameterValues'

#output file ids:
label_id = '_label'
surf_id = '_surface'
###############################################################################
def find_file(fname):
  found = False;
  files = os.listdir(os.getcwd())
  for f in files:
    if f.find(fname) != -1: #If I find it break out of loop
      found = True;
      break;
    
  if found:
    return True;
  else:
    return False;
###############################################################################
def main_menu():
  print ('\nOPTIONS:')
  print ('  1. RUN SEGMENTATION MODULE')
  print ('  2. VIEW PREVIOUSLY GENERATED RESULTS')
  print ('  3. QUIT')
  choice = raw_input('Please enter selection: ')
  return choice
###############################################################################
def module_menu():
  print("\nMODULES:")
  print("  0. EXECUTE ALL:")
  print("  1. NUCLEAR SEGMENTATION")
  print("  2. VESSEL SEGMENTATION")
  print("  3. TRACING OF ASTROCYTES")
  print("  4. TRACING OF MICROGLIA")
  print("  5. COMPUTE ASSOCIATIVE FEATURES")
  print("  6. COMPUTE INTRINSIC NUCLEAR FEATURES")
  print("  7. CLASSIFICATION OF NUCLEI")
  print("  8. EXIT MENU")
  choice = raw_input('Please enter selection: ')
  return choice
###############################################################################
def result_menu():
  print("\nRESULT VIEWERS:")
  print("  1. VESSELS")
  print("  2. TRACES")
  print("  4. EXIT MENU")
  choice = raw_input("Please enter selection: ")
  return choice      
###############################################################################
def run_wizard():
  print('\nLETS BEGIN IMAGE PROCESSING!!!\n')

  print('Which image would you like to work with:')
  print('  1. ' + orig_base_name + orig_ext) 
  print('  2. ' + crop_base_name + crop_id + '1' + crop_ext)
  print('  3. ' + crop_base_name + crop_id + '2' + crop_ext)
  image_num = raw_input("Which image would you like to work with? ")

  if image_num=='1':
    print("\nThis is a very large image and requires a powerful PC")
    yn = raw_input("Do you want to continue (y/n)? ")
    if yn=='n' or yn=='N':
      return
    nuc_image = orig_base_name + nuc_id + split_ext
    nuc_result = orig_base_name + nuc_id + label_id + split_ext
    eba_image = orig_base_name + eba_id + split_ext
    eba_result = orig_base_name + eba_id + surf_id + split_ext
    trace_astro_xml = orig_base_name + iba1_id + trace_params_id + '.xml'
    trace_micro_xml = orig_base_name + gfap_id + trace_params_id + '.xml'
    ass_defs = orig_base_name + ass_def_id + '.xml'
  elif image_num=='2' or image_num=='3':
    print("GREAT CHOICE")
    crop_num = str(int(image_num)-1)
    nuc_image = crop_base_name + crop_id + crop_num + nuc_id + crop_ext
    nuc_result = crop_base_name + crop_id + crop_num + nuc_id + label_id + crop_ext
    eba_image = crop_base_name + crop_id + crop_num + eba_id + crop_ext
    eba_result = crop_base_name + crop_id + crop_num + eba_id + surf_id + crop_ext
    trace_astro_xml = crop_base_name + crop_id + crop_num + iba1_id + trace_params_id + '.xml'
    trace_micro_xml = crop_base_name + crop_id + crop_num + gfap_id + trace_params_id + '.xml'
    ass_defs = crop_base_name + crop_id + crop_num + ass_def_id + '.xml'
  else:
    print("I DON'T KNOW THAT NUMBER, ABORTING")
    return
  
  while (1):
    choice = module_menu()
    
    if choice == '0':
      
      print("\nSTARTING NUCLEAR SEGMENTATION...")
      if find_file(nuc_image) and find_file(seg_params):
        subprocess.call(['segment_nuclei.exe', nuc_image, nuc_result, seg_params])
        print("\n...DONE")
      else:
        print("COULD NOT FIND INPUT FILES")

      print("\nSTARTING VESSEL SEGMENTATION...")  
      if find_file(eba_image):
        subprocess.call(['vessel_segmentation.exe', eba_image, eba_result])
        print("\n...DONE")
      else:
        print("COULD NOT FIND INPUT FILE")

      print("\nSTARTING TRACING OF ASTROCYTES...")  
      if find_file(trace_astro_xml):
        subprocess.call(["RPITrace3D.exe", trace_astro_xml])
        print("\n...DONE")
      else:
        print("COULD NOT FIND INPUT FILE")

      print("\nSTARTING TRACING OF MICROGLIA...")
      if find_file(trace_micro_xml):
        subprocess.call(["RPITrace3D.exe", trace_micro_xml])
        print("\n...DONE")
      else:
        print("COULD NOT FIND INPUT FILE")

      print("\nCOMPUTING ASSOCIATIVE FEATURES...")
      if find_file(ass_defs):
        subprocess.call(['compute_associative_measures.exe', ass_defs])
        print("\n...DONE")
      else:
        print("COULD NOT FIND INPUT FILE")

      print("\nCOMPUTING INTRINSIC FEATURES AND CREATING XML RESULT FILE...")  
      if find_file(nuc_image) and find_file(nuc_result):
        subprocess.call(['compute_nuclei_features.exe', os.getcwd(), nuc_image, nuc_result])
        print("\n...DONE\n  YOU MAY OPEN THE RESULT FROM THE FARSIGHT TOOLBAR")
      else:
        print("COULD NOT FIND INPUT FILES")
        
    elif choice == '1':
      print("\nSTARTING NUCLEAR SEGMENTATION...")
      if find_file(nuc_image) and find_file(seg_params):
        subprocess.call(['segment_nuclei.exe', nuc_image, nuc_result, seg_params])
        print("\n...DONE")
      else:
        print("COULD NOT FIND INPUT FILES")
        
    elif choice == '2':
      print("\nSTARTING VESSEL SEGMENTATION...")  
      if find_file(eba_image):
        subprocess.call(['vessel_segmentation.exe', eba_image, eba_result])
        print("\n...DONE")
      else:
        print("COULD NOT FIND INPUT FILE")
        
    elif choice == '3':
      print("\nSTARTING TRACING OF ASTROCYTES...")  
      if find_file(trace_astro_xml):
        subprocess.call(["RPITrace3D.exe", trace_astro_xml])
        print("\n...DONE")
      else:
        print("COULD NOT FIND INPUT FILE")
        
    elif choice == '4':
      print("\nSTARTING TRACING OF MICROGLIA...")
      if find_file(trace_micro_xml):
        subprocess.call(["RPITrace3D.exe", trace_micro_xml])
        print("\n...DONE")
      else:
        print("COULD NOT FIND INPUT FILE")
        
    elif choice == '5':
      print("\nCOMPUTING ASSOCIATIVE FEATURES...")
      if find_file(ass_defs) and find_file(nuc_result):
        subprocess.call(['compute_associative_measures.exe', ass_defs])
        print("\n...DONE")
      else:
        print("COULD NOT FIND INPUT FILES")
        
    elif choice == '6':
      print("\nCOMPUTING INTRINSIC FEATURES AND CREATING XML RESULT FILE...")  
      if find_file(nuc_image) and find_file(nuc_result):
        subprocess.call(['compute_nuclei_features.exe', os.getcwd(), nuc_image, nuc_result])
        print("\n...DONE\n  YOU MAY OPEN THE RESULT FROM THE FARSIGHT TOOLBAR")
      else:
        print("COULD NOT FIND INPUT FILES")
        
    elif choice == '7':
      print("\nTHIS FUNCTION NOT YET IMPLEMENTED")
      
    elif choice == '8':
      return
    else:
      print("\nUNRECOGNIZED OPTION")
###############################################################################
def view_results():
  from farsightutils import GetFilename
  while (1):
    choice = result_menu()
    if choice == '1':
      fname = GetFilename('NPTS','.npts')
      if(fname != ''):
        print("\nSTARTING VISUALIZATION...")
        subprocess.call(['visualize.exe', fname, 'blank'])
        print("\n...DONE")
        
    elif choice == '2':
      fname = GetFilename('XML','.xml')
      if(fname != ''):
        base = os.path.basename(fname)
        loc = base.find('TracedPoints.xml')
        dname = base[0:loc-1]
        print("SEARCHING FOR DATA FILE(" + dname + ")...")
        found = False;
        files = os.listdir(os.getcwd())
        for f in files:
          if f.find(dname) != -1: #If I find it break out of loop
            dname = f;
            found = True;
            break;
        
        if found:
          subprocess.call(['vtkTraceViewer.exe', fname, dname])
        else:
          subprocess.call(['vtkTraceViewer.exe', fname])

        print("\n...DONE")
        
    elif choice == '3':
      return
    else:
      print("\nUNRECOGNIZED OPTION")
###############################################################################  
def main():
  
  print ('\nYou have started the FARSIGHT 5-Label Image Processing DEMO')
  print ('\nAssuming Test Data is in ' + data_dir)
  os.chdir(data_dir)

  while (1):
    choice = main_menu()
    if choice == '1':
      run_wizard()
    elif choice == '2':
      view_results()
    elif choice == '3':
      print("\nGOODBYE")
      return

################################################################################
#if __name__ == "__main__":
main()
