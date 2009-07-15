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
trace_demo_dir = data_dir + os.sep + 'TraceDemo'
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
def classify(test_data):
  print("\nSCALING TRAINING SET...")
  training_set = "NM_Training_Set.txt"
  training_scale = "NM_Training_Set_Scaled.txt"
  training_range = "NM_Training_Set_Ranges.txt"
  cout = file(training_scale,'w')
  arg = "svm-scale -s " + training_range + " " + training_set
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
  arg = "svm-scale -r " + training_range + " " + test_data
  subprocess.Popen(arg,stdout=cout).wait()
  cout.close()
  print("...DONE")

  print("\nTRAINING BASED ON SCALED TRAINING SET...")
  cout = file('train.log','w')
  arg = "svm-train -t 2 " + training_scale
  subprocess.Popen(arg,stdout=cout).wait()
  cout.close()
  print("...DONE")

  print("\nCLASSIFYING TEST DATA...")
  testing_predict = dname + "_Classes.txt"
  arg = "svm-predict " + testing_scale + " " + training_scale + ".model" + " " + testing_predict
  subprocess.Popen(arg).wait()
  print("...DONE")

  print("\nIMPORTING CLASS INFO INTO XML...")
  loc = dname.find('_libSVM')
  xname = dname[0:loc] + '.xml'
  arg = "classify_nuclei " + os.getcwd() + " " + xname + " " + testing_predict
  subprocess.Popen(arg).wait();
  print("...DONE")
###############################################################################
def bioformats_menu():
  print("\nBIO-FORMATS Command Line Tools: ")
  print("  1. showinf - print image info and display image")
  print("  2. bfview - launch bio-format image viewer")
  print("  3. ijview - dispay image in ImageJ using plugin")
  print("  4. bfconvert - convert image from one format to another")
  print("  5. EXIT MENU")
  choice = raw_input('Please enter selection: ')
  return choice
###############################################################################
def run_bioformats():
  print("\nBio-Formats is a standalone Java library for reading \
and writing life sciences image file formats. It is capable of \
parsing both pixels and metadata for a large number of formats, \
as well as writing to several formats")

  while (1):
    choice = bioformats_menu()

    #showinf - allow user to browse for file
    if choice == '1':
      ftypes = [('all files', '.*')]
      fname = GetFilename(ftypes)
      if(fname != ''):
        args = "showinf.bat " + fname
        subprocess.Popen(args).wait()

    #open bfview    
    elif choice == '2':
      subprocess.Popen("bfview.bat").wait()

    #open ijview
    elif choice == '3':
      ftypes = [('all files', '.*')]
      fname = GetFilename(ftypes)
      if(fname != ''):
        args = "ijview.bat " + fname
        subprocess.Popen(args).wait()

    #bfconvert - all selects file to convert, then types in new name    
    elif choice == '4':
      ftypes = [('all files', '.*')]
      fname = GetFilename(ftypes)
      if(fname != ''):
        print("Base File: " + fname)
        fdir1 = os.path.abspath( os.path.dirname(fname) ) #absolute directory of the file
        fnam1 = os.path.basename(fname)			                #name of file
        os.chdir(fdir1)
        new_name = raw_input('Enter New Name: ')
        fdir2 = os.path.abspath( os.path.dirname(new_name) ) #absolute directory of the file
        fnam2 = os.path.basename(new_name)			                #name of file
        
        args = "bfconvert.bat " + fnam1 + " " + fnam2
        subprocess.Popen(args).wait()
        os.chdir(data_dir)
        
    elif choice == '5':
      return
    
    else:
      print("\nUNRECOGNIZED OPTION")  
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
    trace_astro_xml = crop_base_name + crop_id + crop_num + gfap_id + trace_params_id + '.xml'
    trace_micro_xml = crop_base_name + crop_id + crop_num + iba1_id + trace_params_id + '.xml'
    trace_astro_out = crop_base_name + crop_id + crop_num + gfap_id + traced_id + '.xml'
    trace_micro_out = crop_base_name + crop_id + crop_num + iba1_id + traced_id + '.xml'
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
        #now make sure a cache folder exists
        if not os.path.exists('cache'):
          os.mkdir('cache')
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
        #now make sure a cache folder exists
        if not os.path.exists('cache'):
          os.mkdir('cache')
        subprocess.call(['render.exe', rend_params])
        print("\n...DONE")
      else:
        print("COULD NOT FIND INPUT FILES")
        
    elif choice == '9':
      return
    else:
      print("\nUNRECOGNIZED OPTION")
      
###############################################################################
def main_menu():
  print ('\nDEMO OPTIONS:')
  print ('  1. RENDER 5-LABEL IMAGE PROCESSING RESULTS')
  print ('  2. 5-LABEL IMAGE PROCESSING DEMO')
  print ('  3. 3D TRACING DEMO') 
  print ('  4. 3D NUCLEAR SEGMENTATION DEMO')
  print ('  5. 3D VESSEL SEGMENTATION DEMO')
  print ('  6. REGISTRATION DEMO')
  print ('  7. TISSUE NETS DEMO')
  print ('  8. RENDER GRAYSCALE IMAGE')
  print ('  9. RENDER SEGMENTATION RESULT FILE')
  print (' 10. OPEN TRACE EDITOR')
  print (' 11. BIO-FORMATS COMMAND-LINE TOOLS')
  print (' 12. QUIT')
  choice = raw_input('Please enter selection: ')
  return choice
###############################################################################
def main():
  
  print ('\nYou have started the FARSIGHT DEMO')
  print ('\nAssuming Test Data is in ' + data_dir)
  os.chdir(data_dir)

  while (1):
    choice = main_menu()

    #DO 3D RENDERING OF ALL 4 CHANNELS, AND COLOR NUCLEI BY CLASS:
    if choice == '1':
      os.chdir(full_image_dir)
      #now make sure a cache folder exists
      if not os.path.exists('cache'):
        os.mkdir('cache')
      args = "render.exe 100upoint5%25hippo25x1unmixed_RenderParameters.txt";
      cout = file("cout.log",'w')
      print("\nSTARTING VISUALIZATION (IN NEW THREAD - PLEASE BE PATIENT)...")
      subprocess.Popen(args,stdout=cout);
      os.chdir(data_dir)

    #RUN THE PROCESSING DEMO THAT ALLOWS FOR EACH STEP IN PROCESSING TO BE PERFORMED ON TWO IMAGES
    elif choice == '2':
      os.chdir(crop_image_dir)
      run_wizard()
      os.chdir(data_dir)

    #TRACING DEMO OF SINGLE NEURON.  TRACE, then OPEN Trace Editor
    elif choice == '3':
      os.chdir(trace_demo_dir)
      print("\nTRACING...")
      cout = file("trace.log",'w')
      args = "RPITrace3D.exe ETR052Z1_TracingParameterValues.xml"
      subprocess.Popen(args,stdout=cout).wait()
      cout.close()
      print("...DONE")
      print("OPENING TRACE EDITOR...")
      args = "trace_editor.exe ETR052Z1TracedPoints.xml ETR052Z1.pic"
      subprocess.Popen(args,stdout=subprocess.PIPE)
      os.chdir(data_dir)

    #SHOW 3D RENDING OF ORIGINAL DATA OF SMALL IMAGE WITH NUCLEI.
    #SEGMENT THE NUCLEI, THEN SHOW A 3D RENDERING OF THE RESULT
    elif choice == '4':
      os.chdir(crop_image_dir)
      crop_num = '1'
      nuc_image = crop_base_name + crop_id + crop_num + nuc_id + crop_ext
      nuc_result = crop_base_name + crop_id + crop_num + nuc_id + label_id + ".tiff"
      seg_params = crop_base_name + crop_id + crop_num + nuc_id + seg_params_id + ".ini"
      
      if find_file(nuc_image) and find_file(seg_params):
        print("\nRENDERING NUCLEI CHANNEL...")
        args = "trace_editor.exe " + nuc_image
        subprocess.Popen(args,stdout=subprocess.PIPE)
        print("STARTING NUCLEAR SEGMENTATION...")
        args = 'segment_nuclei.exe ' + nuc_image + " " + nuc_result + " " + seg_params
        cout = file("seg_nuclei.log",'w')
        subprocess.Popen(args,stdout=cout).wait()
        cout.close()
        print("\n...DONE")
        print("RENDING SEGMENTATION RESULT...")
        #now make sure a cache folder exists
        if not os.path.exists('cache'):
          os.mkdir('cache')
        args = "render.exe NM_crop1_Nuc_RenderParameters.txt"
        subprocess.Popen(args,stdout=subprocess.PIPE)
        
      os.chdir(data_dir)

    #SHOW 3D RENDING OF VESSEL DATA, SEGMENT THE VESSELS, THEN SHOW 3D RENDING OF RESULT
    elif choice == '5':
      os.chdir(crop_image_dir)
      crop_num = '1'
      eba_image = crop_base_name + crop_id + crop_num + eba_id + crop_ext
      eba_result = crop_base_name + crop_id + crop_num + eba_id + surf_id + crop_ext

      if find_file(eba_image):
        print("\nRENDERING VESSEL CHANNEL...")
        args = "trace_editor.exe " + eba_image
        subprocess.Popen(args,stdout=subprocess.PIPE)
        print("STARTING VESSEL SEGMENTATION...")
        args = 'vessel_segmentation.exe ' + eba_image + " " + eba_result
        cout=file("seg_vessels.log",'w')
        subprocess.Popen(args,stdout=cout).wait()
        cout.close()
        print("...DONE")
        print("RENDERING SEGMENTATION RESULT...")
        #now make sure a cache folder exists
        if not os.path.exists('cache'):
          os.mkdir('cache')
        args = "render.exe NM_crop1_EBA_RenderParameters.txt"
        subprocess.Popen(args,stdout=subprocess.PIPE)
        
      os.chdir(data_dir)

    #RUN THE REGISTRATION DEMO - CREATES MONTAGE OF 2 CROPPED REGIONS  
    elif choice == '6':
      os.chdir(crop_image_dir)
      from register_pairs import register
      print("\nSTARTING REGISTRATION OF PAIRS")
      register([os.getcwd()+os.sep,'NM_RegistrationPairs.txt'])
      print("\n...MONTAGE CREATED")
      print("LOAD THE MONTAGE IMAGE (CropDemo/pairwise_montage_NM_crop2_Nuc.tiff)")
      print("IN FARSIGHT TO VIEW RESULTS")
      #print("\nOPENING MONTAGE BROWSER")
      #subprocess.call(["MontageNavigator.exe"])
      #subprocess.call(['trace_editor.exe', "montage_NM_crop2_Nuc.tiff"])
      #print("\nMONTAGE BROWSER CLOSED")
      os.chdir(data_dir)

    #RUN THE TISSUENETS DEMO  
    elif choice == '7':
      os.chdir(tissuenets_dir)
      print("\nIMPORTING tissuenets_demo.py")
      import tissuenets_demo
      os.chdir(data_dir)

    #ALLOW USER TO BROWSE FOR A DATA FILE, THEN RENDER IT USING trace_edit program
    #FILE MUST BE A GRAYSCALE IMAGE
    elif choice == '8':
      ftypes = [('PIC', '.pic'),('TIFF', '.tiff .tif'),('all files', '.*')]
      fname = GetFilename(ftypes)
      if(fname != ''):
        args = 'trace_editor.exe ' + fname
        subprocess.Popen(args,stdout=subprocess.PIPE)

    #ALLOW USER TO BROWSE FOR A RESULT FILE, THEN RENDER THE OBJECTS using render program
    #FILE MY BE label image, binary image, or ..TracedPoints.xml.
    elif choice == '9':
      os.chdir(data_dir)
      ftypes = [('XML','.xml'),('PIC', '.pic'),('TIFF', '.tiff .tif'),('all files', '.*')]
      fname = GetFilename(ftypes)
      if(fname != ''):
        fdir = os.path.abspath( os.path.dirname(fname) ) #absolute directory of the file
        fnam = os.path.basename(fname)			                #name of file
        os.chdir(fdir)
        rfile = file('r.tmp','w')
        rfile.write("1 1 1\n")
        rfile.write(fnam + " 0 1 1")
        rfile.close()
        #now make sure a cache folder exists
        if not os.path.exists('cache'):
          os.mkdir('cache')
        args = 'render.exe r.tmp'
        subprocess.Popen(args,stdout=subprocess.PIPE)
        os.chdir(data_dir)


    #ALLOW USER TO BROWSE FOR A FILE AND THEN OPEN THE TRACE EDITOR FOR THIS FILE:
    elif choice == '10':
      os.chdir(data_dir)
      ftypes = [('XML','.xml'),('all files', '.*')]
      fname = GetFilename(ftypes)
      if(fname != ''):
        fdir = os.path.abspath( os.path.dirname(fname) ) #absolute directory of the file
        os.chdir(fdir)
        fnam = os.path.basename(fname)			                #name of file
        loc = fnam.find('TracedPoints.xml')
        dname = fnam[0:loc-1]
        print("SEARCHING FOR DATA FILE(" + dname + ")...")
        found = False;
        files = os.listdir(os.getcwd())
        for f in files:
          if f.find(dname) != -1: #If I find it break out of loop
            dname = f;
            found = True;
            break;
        
        if found:
          args = 'trace_editor.exe ' + fname + ' ' + dname
        else:
          args = 'trace_editor.exe ' + fname
          
        subprocess.Popen(args,stdout=subprocess.PIPE)
        os.chdir(data_dir)

    elif choice == '11':
      run_bioformats()
      os.chdir(data_dir)
      
    #EXIT  
    elif choice == '12':
      print("\nGOODBYE")
      return

################################################################################
#if __name__ == "__main__":
main()
