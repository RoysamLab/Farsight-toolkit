#!/usr/bin/env python
#############################################################################
# A SCRIPT TO RUN THROUGH THE 5 LABEL TEST IMAGE FROM BEGINNING TO END:
# THESE IMAGES CAN BE FOUND ON THE CENSSIS SERVER AT: /data/SHARE/Yousef/NN_Images
################################################################################
import subprocess
import sys
import os

from farsightutils import *

data_dir = 'C:' + os.sep + 'BADRI_TEST'
full_image_dir = data_dir + os.sep + 'ORIGINAL'
crop_image_dir = data_dir + os.sep + 'CropDemo'
tissuenets_dir = data_dir + os.sep + 'TissueNetsDemo'

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
seg_params_id = '_SegParams'
ass_def_id = '_AssociationDefs'
trace_params_id = '_TracingParameterValues'
rend_params_id = '_RenderParameters'

#output file ids:
label_id = '_label'
surf_id = '_surface'
ass_feat_id = '_AssocFeatures'
traced_id = 'TracedPoints'

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
  print ('\nDEMO OPTIONS:')
  print ('  1. COMPLETE IMAGE VISUALIZATION')
  print ('  2. 5-LABEL IMAGE PROCESSING')
  print ('  3. TRACE EDITOR')
  print ('  4. REGISTRATION DEMO')
  print ('  5. TISSUE NETS DEMO')
  print ('  6. QUIT')
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
  print("  8. VIEW RESULT RENDERING")
  print("  9. EXIT MENU")
  choice = raw_input('Please enter selection: ')
  return choice          
###############################################################################
def classify(test_data):
  print("\nSCALING TRAINING SET...")
  training_set = "NM_Training_Set.txt"
  training_scale = "NM_Training_Set_Scaled.txt"
  training_range = "NM_Training_Set_Ranges.txt"
  cout = file(training_scale,'w')
  arg = "svmscale -s " + training_range + " " + training_set
  subprocess.Popen(arg,stdout=cout).wait()
  cout.close()
  print("...DONE")

  print("\nSCALING TESTING SET...")
  if not find_file(test_data):
    print("COULD NOT FIND INPUT")
    return

  base = os.path.basename(test_data)
  loc = base.find('.')
  dname = base[0:loc]
  testing_scale = dname + "_Scaled.txt"
  cout = file(testing_scale,'w')
  arg = "svmscale -r " + training_range + " " + test_data
  subprocess.Popen(arg,stdout=cout).wait()
  cout.close()
  print("...DONE")

  print("\nTRAINING BASED ON SCALED TRAINING SET...")
  cout = file('train.log','w')
  arg = "svmtrain -t 2 " + training_scale
  subprocess.Popen(arg,stdout=cout).wait()
  cout.close()
  print("...DONE")

  print("\nCLASSIFYING TEST DATA...")
  testing_predict = dname + "_Classes.txt"
  arg = "svmpredict " + testing_scale + " " + training_scale + ".model" + " " + testing_predict
  subprocess.Popen(arg).wait()
  print("...DONE")

  print("\nIMPORTING CLASS INFO INTO XML...")
  loc = dname.find('_libSVM')
  xname = dname[0:loc] + '.xml'
  arg = "classify_nuclei " + os.getcwd() + " " + xname + " " + testing_predict
  subprocess.Popen(arg).wait();
  print("...DONE")

###############################################################################
def run_wizard():
  print('\nLETS BEGIN IMAGE PROCESSING!!!\n')

  print('Which image would you like to work with:') 
  print('  1. ' + crop_base_name + crop_id + '1')
  print('  2. ' + crop_base_name + crop_id + '2')
  image_num = raw_input("Which image would you like to work with? ")

  if image_num=='1' or image_num=='2':
    print("GREAT CHOICE")
    #crop_num = str(int(image_num)-1)
    crop_num = image_num
    nuc_image = crop_base_name + crop_id + crop_num + nuc_id + crop_ext
    nuc_result = crop_base_name + crop_id + crop_num + nuc_id + label_id + ".tiff"
    seg_params = crop_base_name + crop_id + crop_num + nuc_id + seg_params_id + ".ini"
    eba_image = crop_base_name + crop_id + crop_num + eba_id + crop_ext
    eba_result = crop_base_name + crop_id + crop_num + eba_id + surf_id + crop_ext
    trace_astro_xml = crop_base_name + crop_id + crop_num + iba1_id + trace_params_id + '.xml'
    trace_micro_xml = crop_base_name + crop_id + crop_num + gfap_id + trace_params_id + '.xml'
    trace_astro_out = crop_base_name + crop_id + crop_num + iba1_id + traced_id + '.xml'
    trace_micro_out = crop_base_name + crop_id + crop_num + gfap_id + traced_id + '.xml'
    ass_defs = crop_base_name + crop_id + crop_num + ass_def_id + '.xml'
    ass_feats = crop_base_name + crop_id + crop_num + ass_def_id + ass_feat_id + '.XML'
    rend_params = crop_base_name + crop_id + crop_num + rend_params_id + '.txt'
    svm_file = crop_base_name + crop_id + crop_num + nuc_id + "_libSVM.txt"
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
        if find_file(ass_feats):
          subprocess.call(['compute_nuclei_features.exe', os.getcwd(), nuc_image, nuc_result, ass_feats])
        else:
          subprocess.call(['compute_nuclei_features.exe', os.getcwd(), nuc_image, nuc_result])
        print("\n...DONE\n  YOU MAY OPEN THE RESULT FROM THE FARSIGHT TOOLBAR")
      else:
        print("COULD NOT FIND INPUT FILES")

      #DO CLASSIFICATION:
      classify(svm_file)

      print("\nSTARTING VISUALIZATION...")
      if find_file(nuc_result) and find_file(eba_result) and find_file(trace_astro_out) and find_file(trace_micro_out) and find_file(rend_params):
        subprocess.call(['render.exe', rend_params])
        print("\n...DONE")
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
        if find_file(ass_feats):
          subprocess.call(['compute_nuclei_features.exe', os.getcwd(), nuc_image, nuc_result, ass_feats])
        else:
          subprocess.call(['compute_nuclei_features.exe', os.getcwd(), nuc_image, nuc_result])
        print("\n...DONE\n  YOU MAY OPEN THE RESULT FROM THE FARSIGHT TOOLBAR")
      else:
        print("COULD NOT FIND INPUT FILES")
        
    elif choice == '7':
      classify(svm_file)

    elif choice == '8':
      if find_file(nuc_result) and find_file(eba_result) and find_file(trace_astro_out) and find_file(trace_micro_out) and find_file(rend_params):
        print("\nSTARTING VISUALIZATION...")
        subprocess.call(['render.exe', rend_params])
        print("\n...DONE")
      else:
        print("COULD NOT FIND INPUT FILES")
        
    elif choice == '9':
      return
    else:
      print("\nUNRECOGNIZED OPTION")
###############################################################################  
def main():
  
  print ('\nYou have started the FARSIGHT DEMO')
  print ('\nAssuming Test Data is in ' + data_dir)
  os.chdir(data_dir)

  while (1):
    choice = main_menu()
    if choice == '1':
      
      os.chdir(full_image_dir)
      
      nuc_result = orig_base_name + nuc_id + label_id + ".tiff"
      eba_result = orig_base_name + eba_id + surf_id + split_ext
      trace_astro_out = orig_base_name + iba1_id + traced_id + '.xml'
      trace_micro_out = orig_base_name + gfap_id + traced_id + '.xml'
      rend_params = orig_base_name + rend_params_id + '.txt'
      
      if find_file(nuc_result) and find_file(eba_result) and find_file(trace_astro_out) and find_file(trace_micro_out) and find_file(rend_params):
        print("\nSTARTING VISUALIZATION...")
        subprocess.call(['render.exe', rend_params])
        print("\n...DONE")
      else:
        print("COULD NOT FIND INPUT FILES")

      os.chdir(data_dir)
      
    elif choice == '2':
      os.chdir(crop_image_dir)
      run_wizard()
      os.chdir(data_dir)
      
    elif choice == '3':
      os.chdir(crop_image_dir)
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
          subprocess.call(['trace_editor.exe', fname, dname])
        else:
          subprocess.call(['trace_editor.exe', fname])
        print("\n...DONE")
        os.chdir(data_dir)
        
    elif choice == '4':
      os.chdir(crop_image_dir)
      from register_pairs import register
      print("\nSTARTING REGISTRATION OF PAIRS")
      register([os.getcwd()+os.sep,'NM_RegistrationPairs.txt'])
      print("\n...DONE")
      print("\nOPENING MONTAGE BROWSER")
      subprocess.call(["MontageNavigator.exe"])
      print("\nMONTAGE BROWSER CLOSED")
      os.chdir(data_dir)
      
    elif choice == '5':
      os.chdir(tissuenets_dir)
      print("\nIMPORTING tissuenets_demo.py")
      #os.environ['CLASSPATH'] = "C:\\Program Files (x86)\\Farsight 0.1.1\\bin\\saxon9.jar"
      import tissuenets_demo
      os.chdir(data_dir)
      
    elif choice == '6':
      print("\nGOODBYE")
      return

################################################################################
#if __name__ == "__main__":
main()
