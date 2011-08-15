# Python scrip which register all possible pairs and can be
# initialized with prior transformations if given in an xml file
# containing initial transformations

import os, sys, re, platform

#initialize executables for platform
regp = ''
regj = ''
moim = ''
moimp = ''
much = ''

if platform.system() == 'Windows':
    regp = 'register_pair.exe'
    regj = 'register_joint.exe'
    moim = 'mosaic_images.exe'
    moimp = 'mosaic_image_pair.exe'
    much = 'multi_channels_2D.exe'
else:
    regp = 'register_pair'
    regj = 'register_joint'
    moim = 'mosaic_images'
    moimp = 'mosaic_image_pair'
    much = 'multi_channels_2D'

# A function to convert a string to a number of base 26. The input is
# assumed to be a string of character
def char_to_num(index):
    index_upper=index.upper()
    power = len(index_upper)
    total = 0
    pos = 0
    for char in index:
        total = total + (ord(char)-ord('A')) * pow(26,power-pos-1)
        pos = pos+1
    return total

# The function to create image pairs for registration
def create_pairs(argv):
    f=open(argv[0],'r');
    names=[]
    rows=[]
    cols=[]
    pairs=[]

    #pairwise registration
    for line in f:
        s_line = line.rstrip().rstrip('\n')
        names.append(s_line);
        
    # Now do the pairing for from_name in names:
    from_ind = -1
    for from_name in names:
        from_ind = from_ind + 1
        to_ind = from_ind
        for to_name in names[from_ind+1:]:
            to_ind = to_ind + 1
            pairs.append(from_name+' '+to_name)
            print 'pair = '+from_name+' '+to_name
    return pairs

# The function to perform pairwise registration, joint registration,
# and montaging by taking the first image
def register(pair_list, argv):
    image_dir = argv[0]
    #f=open(sys.argv[1],'r');
    f_o = open(argv[1]+'.failed_pairs','w')
    f_xforms = open('xxx_123.txt','w')
    names=[]
    xforms=[]

    #pairwise registration
    for line in pair_list:
        s_line = line.rstrip().rstrip('\n')
        pos = s_line.find(' ')
        from_image = s_line[:pos]
        to_image = s_line[pos+1:]
            
        # perform registration
        if (len(argv)>2):
             success = os.system(regp+' '+image_dir+from_image+' '+ image_dir+to_image +' -remove_2d');
        else:
            success = os.system(regp+' '+image_dir+from_image+' '+ image_dir+to_image +' -remove_2d');
        if success == 0:
            # Add the names to the list if not already there. The list
            # keeps potential images which can be the anchor images for
            # montaging
            found = False 
            for name in names:
                if name == from_image:
                    found = True
                    break
            if not found: names.append(from_image)
        
            found = False 
            for name in names:
                if name == to_image:
                    found = True
                    break
            if not found: names.append(to_image)
        
            # Add the transformed.xml to the list
            from_dot = from_image.find('.')
            to_dot = to_image.find('.');
            f_xforms.write(from_image[:from_dot]+"_to_"+to_image[:to_dot]+"_transform.xml\n")
        else:
            f_o.write("From "+from_image+" to "+to_image+" failed\n")
            print "From "+from_image+" to "+to_image+" failed"
            
    # Perform joint registration write the xform list to a
    # temporary file and remove it after joint registration
    f_xforms.close()
    f_o.close();
    print("\nSTART register_joint.exe...")
    os.system(regj+' xxx_123.txt -multiplier 4')
    
    print("DONE")

    # perform montaging using the first image as the anchor
    print("\nSTART mosaic_images ...")
    cmd_executed = False;
    numPairs = len(pair_list);
    if (numPairs > 1):
        if (len(argv) > 2): #multiple channels
            fc = open(argv[2],'r');
            fc_o = open(argv[2]+'_123.txt','w')
            channel_count = 0;
            for line in fc:
                dot_pos = names[0].find('.');
                name_no_ext = names[0][:dot_pos]
                os.system(moim + ' joint_transforms.xml '+names[0]+ ' -3d -path ' + argv[0]+' -channel '+str(channel_count)+' -output montage_'+name_no_ext+'_Ch'+str(channel_count));
                fc_o.write('montage_'+name_no_ext+'_Ch'+str(channel_count)+'_2d_proj.png '+ line)
                channel_count += 1
            fc_o.close()
            os.system(much + ' '+argv[2]+'_123.txt '+'montage_'+name_no_ext+'_color_2d_proj.png')
            os.remove(argv[2]+'_123.txt')
            cmd_executed = True;
        else :
            cmd = moim + " joint_transforms.xml -3d " + names[0] + " -path " + argv[0]
    else : # just for one image pair, so call mosaic_image_pair
        cmd = moimp + " joint_transforms.xml " + names[0] +" "+ names[1] + " -path " + argv[0]
    if ( not cmd_executed):
        os.system(cmd)
        
    print("DONE")
    
    # TEMP FILE CLEANUP:
    print("\nCLEANING TEMP FILES...")
    os.remove('xxx_123.txt')
    print("DONE")

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print '\nUsage: '+sys.argv[0]+' image_dir image_list_file [channel_color_file]\n '
        sys.exit(1)
    pairs = create_pairs(sys.argv[2:])
    register(pairs,sys.argv[1:])


