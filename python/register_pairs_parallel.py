#!/usr/bin/python

#python script to register pairs specified in a file, which contains a pair of
#image names per line.

import os, sys, platform

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

# the main function
def register(argv):
    if len(argv) < 2:
        print 'Arguments: image_dir pair_list'
        sys.exit(1)

    f=open(argv[1],'r');
    f_o = open(argv[1]+'.failed_pairs','w')
    f_xforms = open('xxx_123.txt','w')
    names=[]
    xforms=[]
    from_image_list=[]
    to_image_list=[]
    subprocess_command_list=[]
    subprocess_list=[]
    file_handle_list=[]

    #pairwise registration
    for line in pair_list:
        s_line = line.rstrip().rstrip('\n')
        pos = s_line.find(' ')
        from_image = s_line[:pos]
        to_image = s_line[pos+1:]
        
        from_image_list.append(from_image)
        to_image_list.append(to_image)
        
        subprocess_command_list.append(regp+' '+image_dir+from_image+' '+ image_dir+to_image +' -remove_2d')
    
    #spawn and add subprocesses (each representing an execution of register_pair.exe) into a list
    idx = 0
    for subprocess_command in subprocess_command_list:
        fh = open('debug_from_' + from_image_list[idx] + '_to_' + to_image_list[idx] + '.txt', 'w')
        subprocess_list.append(subprocess.Popen(subprocess_command, stdout = fh, stderr = fh))
        file_handle_list.append(fh)
        idx = idx + 1
        
    #loop through entire list to check if process are done
    processes_done = False
    while not processes_done:
        processes_done = True
        time.sleep(1) #sleep 1 second before polling all processes
        for subp in subprocess_list:
            if subp.poll() == None:
                processes_done = False
    
    #close all the file handles since we are done writing to them
    for fh in file_handle_list:
        fh.close()
    
    #all processes are done, so read return values and take appropriate action
    idx = 0
    for subp in subprocess_list:
        if subp.poll() == 0: #success
            found = False
            for name in names:
                if name == from_image_list[idx]:
                    found = True
                    break
            if not found:
                names.append(from_image_list[idx])
            
            found = False
            for name in names:
                if name == to_image_list[idx]:
                    found = True
                    break
            if not found:
                names.append(to_image)
                
            #add the transformed.xml to the list
            from_dot = from_image_list[idx].find('.')
            to_dot = to_image_list[idx].find('.')
            f_xforms.write(from_image_list[idx][:from_dot]+"_to_"+to_image_list[idx][:to_dot]+"_transform.xml\n")
            print "From "+from_image+" to "+to_image+" succeeded"
        else:
            f_o.write("From "+from_image_list[idx]+" to "+to_image_list[idx]+" failed\n")
            print "From "+from_image+" to "+to_image+" failed"
            
        idx = idx + 1
    
    
    # Perform joint registration write the xform list to a
    # temporary file and remove it after joint registration
    f_xforms.close()
    f_o.close();
    print("\nSTART register_joint...")
    os.system(regj + ' xxx_123.txt -multiplier 4')
    print("DONE")
    #os.system('rm xxx_123.txt')    #not cross platform!!

    # perform montaging using the first image as the anchor
    print("\nSTART...")
    cmd_executed = False;
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
        print 'Usage: '+sys.argv[0]+' image_dir pair_list_file [channel_color_file]\n'
        sys.exit(1)
    register(sys.argv[1:])
