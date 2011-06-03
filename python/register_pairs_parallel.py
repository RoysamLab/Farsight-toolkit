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
    launched_processes_list=[]
    launched_subp_list=[]
    still_launched_subp_list=[]

     #pairwise registration
    for line in pair_list:
        s_line = line.rstrip().rstrip('\n')
        pos = s_line.find(' ')
        from_image = s_line[:pos]
        to_image = s_line[pos+1:]
        
        from_image_list.append(from_image)
        to_image_list.append(to_image)
        
        subprocess_command_list.append(regp+' '+image_dir+from_image+' '+ image_dir+to_image +' -remove_2d')

    if (not (os.path.exists(os.getcwd() + '\\debug\\'))):
        os.makedirs(os.getcwd() + '\\debug\\')
    
    #spawn and add subprocesses (each representing an execution of register_pair.exe) into a list
    idx = 0
    threads_launched = 0
    for subprocess_command in subprocess_command_list:
        #while loop to check active processes to see if we can launch more threads
        while threads_launched >= multiprocessing.cpu_count():
            threads_launched = 0
            still_launched_subp_list = []
            for launched_subp in launched_processes_list:                  
                if launched_subp.poll() == None:
                    threads_launched = threads_launched + 1
                    still_launched_subp_list.append(launched_subp)
            launched_processes_list = still_launched_subp_list[:]
            #print "Threads launched: " + str(threads_launched)
            time.sleep(1)
            
                        
        fh = open(os.getcwd() + '\debug\debug_from_' + from_image_list[idx] + '_to_' + to_image_list[idx] + '.txt', 'w')
        print "Launching " + subprocess_command
        subp = subprocess.Popen(subprocess_command, stdout = fh, stderr = fh)
        
        subprocess_list.append(subp)
        launched_processes_list.append(subp)
        file_handle_list.append(fh)
        threads_launched = threads_launched + 1;
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
    start_time = time.clock()
    register(sys.argv[1:])
    print 'Registration and Montaging took: ' + time.clock() - start_time + ' seconds'
