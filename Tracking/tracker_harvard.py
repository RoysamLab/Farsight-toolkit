import os, re, operator, pdb, subprocess


def populate_directories(path, dataset_id):
    if os.path.lexists(path) == False:
        os.mkdir(path)
    if os.path.lexists(os.path.join(path, 'cache')) == False:
        os.mkdir(os.path.join(path, 'cache'))
    if os.path.lexists(os.path.join(path, 'cache', dataset_id)) == False:
        os.mkdir(os.path.join(path, 'cache', dataset_id))

if __name__ == '__main__':

    # default values for all variables defined here
   # data_directory = 'C:\\Users\\Arun\\Research\\piexoto_data\\TSeries-02102009-1455-624-3Dmovies'
    data_directory = 'C:\\Users\\Arun\\Research\\Tracking\\data\\TSeries-02102009-1455-624-3Dmovies'
    #cwd = 'C:\\Users\\Arun\\Research\\Tracking\\harvard'
    cwd = 'L:\\Tracking\\'
    exe_dir = 'C:\\Users\\Arun\\Research\\Farsight\\exe\\bin'
    dataset_id = 'test'
    number_of_cores = 7;
    dirList = os.listdir(data_directory)
    list_of_filenames =[]
    channels = []
    time_points = []
    #pdb.set_trace()
    for fname in dirList:
        if fname.endswith('.tif') and os.path.isfile(os.path.join(data_directory,fname)):
            list_of_filenames.append(fname)
    #pdf.set_trace()
    filenames = {(0,0):'Null'}
    #TSeries-02102009-1455-624_Cycle001_CurrentSettings_Ch2
    for fname in list_of_filenames:
        m = re.search('^(.*)_Cycle([0-9]+).*_Ch([0-9]+)\.tif',fname) 
        dataset_temp = m.group(1) # set dataset_id
        channels.append(int(m.group(3))) # collect all channel values
        time_points.append(int(m.group(2))) # collect all time_points
        filenames[int(m.group(3)), int(m.group(2))] = fname # we establish a reverse dictionary lookup of the filenames from channel number and time point
        
    channels = sorted(list(set(channels)))
    time_points = sorted(list(set(time_points)))
    #pdb.set_trace()
    dataset_id += dataset_temp
    num_channels = len(channels)
    num_time_points = len(time_points)
    print (num_channels, num_time_points)

    populate_directories(cwd,dataset_id)
    os.chdir(cwd)

    cache_prefix = os.path.join(cwd,'cache',dataset_id)
    
    # it should be easy to just call the exe's with the filenames from now on
    # use filenames dictionary lookup to get the filenames
    # prefix the filanames with necessary things like
    # unmixed_, labeled_, labeled_tracks_, vessel_binarized_, etc..

    time_points = time_points[0:30] # DEBUG
    #channels = [2,3,4]
    #pdb.set_trace()
    ######################### Delete slices #############################
##    for w in channels:
##        for t in time_points:
##            temp_fname = [];
##            temp_fname.append(os.path.join(exe_dir,'extract_slices'))
##            temp_fname.append(os.path.join(data_directory,filenames[(w,t)]))
##            temp_fname.append('2')
##            temp_fname.append('44')
##            temp_fname.append(os.path.join(cache_prefix, 'extracted_' + filenames[(w,t)]))
##            subprocess.call(temp_fname)


    ######################### Smooth ####################################

    popen_objs = [];
    for t in time_points:
        for w in channels:
            temp_fname = []
            output_filename =  os.path.join(cache_prefix, 'smoothed_' + filenames[(w,t)])
            if os.path.exists(output_filename):
                continue
            # for anisodiff
            # temp_fname.append(os.path.join(exe_dir,'anisodiff'))
            # for medianfilter
            temp_fname.append(os.path.join(exe_dir,'medianfilter'))
            temp_fname.append(os.path.join(data_directory,filenames[(w,t)]))
            #for anisodiff
            #temp_fname.append('10,0.0625,2')
            #for medianfilter
            temp_fname.append('2,2,1')
            temp_fname.append(os.path.join(cache_prefix, 'smoothed_' + filenames[(w,t)]))
            popen_objs.append(subprocess.Popen(temp_fname))
##            subprocess.call(temp_fname)
            if len(popen_objs) == number_of_cores:
                for obj in popen_objs:
                    out = obj.communicate()[0]    
                popen_objs = []
    for obj in popen_objs:
        out = obj.communicate()[0] 
    #####################################################################
    ######################### Background Subraction #####################
##    for w in channels:
##        temp_fname = [];
##        temp_fname.append(os.path.join(exe_dir,'background_subtraction'))
##        for t in time_points:
##            temp_fname.append(os.path.join(data_directory,filenames[(w,t)]))
##        for t in time_points:
##            temp_fname.append(os.path.join(cache_prefix,'bg_sub_' + filenames[(w,t)]))
##        subprocess.call(temp_fname)

    ######################### Unmixing ##################################
##    for t in time_points:
##        temp_fname = [];
##        temp_fname.append(os.path.join(exe_dir,'unmix'))
##        
##        #first add input filenames then add output filenames
##        for w in channels:
##            temp_fname.append(os.path.join(cache_prefix,'extracted_' + filenames[(w,t)]))
##        for w in channels:
##            temp_fname.append(os.path.join(cache_prefix,'unmixed_' + filenames[(w,t)]))
##        
##        # call the unmix executable
##        subprocess.call(temp_fname)

    ######################## Segmentation ###############################

    nuclei_segmentation_cfg = 'CellSegmentation.ini'
    params = [];
    params.append('70,30,5');
    params.append('2,50,5');

    
    channels_to_segment = [2,3,4]
    popen_objs = [];
##    tpr = int(number_of_cores/len(channels_to_segment));
    # segment all the channels first. Then trace the vessels.
    for t  in time_points:
        for w in channels_to_segment:
            temp_fname = [];
            output_filename = os.path.join(cache_prefix, 'labeled_' + filenames[(w,t)])
            if os.path.exists(output_filename): # dont do anything if the output file exists
                continue;
            if w == 3 or w == 4:
                # call yousef segmentation
                temp_fname.append(os.path.join(exe_dir,'segment_nuclei_harvard'))
                temp_fname.append(os.path.join(cache_prefix, 'smoothed_' + filenames[(w,t)]))
                temp_fname.append(output_filename)
                temp_fname.append(nuclei_segmentation_cfg)# parameters filename
            else:
                # call simple binary morphological operators based segmentation
                temp_fname.append(os.path.join(exe_dir,'segmentation'))
                temp_fname.append(os.path.join(cache_prefix, 'smoothed_' + filenames[(w,t)]))
                temp_fname.append('2') # type of segmentation
                temp_fname.append(params[0]) # paramerters
                temp_fname.append(output_filename)
            popen_objs.append(subprocess.Popen(temp_fname));
            if len(popen_objs) == number_of_cores:
                for obj in popen_objs:
                        out = obj.communicate()[0];
                popen_objs = []
    for obj in popen_objs:
        out = obj.communicate()[0];
    #pdb.set_trace()
    # vessel tracing
##    vessel_w = 4; #vessel channel
##    temp_fname = [];
##    temp_fname.append(os.path.join(exe_dir,'vessel_thinning'))
##    for t in time_points:
##        temp_fname.append(os.path.join(cache_prefix, 'labeled_' + filenames[(vessel_w,t)]))
##    temp_fname.append(os.path.join(cache_prefix, 'vessel_binarized_' + dataset_id + '_w' + str(vessel_w) + '.tif'))
##    temp_fname.append(os.path.join(cache_prefix, 'vessel_trace_'+ dataset_id + '_w' + str(vessel_w) + '.tif'))
##    subprocess.call(temp_fname)
      ############################# Clean up segmentation ################

    for t in time_points:
        for w in channels_to_segment:
            temp_fname = []
            output_filename = os.path.join(cache_prefix, 'clabeled_' + filenames[(w,t)])
            if os.path.exists(output_filename):
                continue;
            temp_fname.append(os.path.join(exe_dir,'remove_small_components'))
            temp_fname.append(os.path.join(cache_prefix, 'labeled_' + filenames[(w,t)]))
            temp_fname.append('100')
            temp_fname.append(os.path.join(cache_prefix, 'clabeled_' + filenames[(w,t)]))
            #pdb.set_trace()
            popen_objs.append(subprocess.Popen(temp_fname))
            if len(popen_objs) == number_of_cores :
                for obj in popen_objs:
                    out = obj.communicate()[0]
                popen_objs = []
    for obj in popen_objs:
        out = obj.communicate()[0]        

##    ############################# Tracking ##############################
##
##    # track channel 1 and 2 only
    channels_to_track = [4]
    for w in channels_to_track:
        temp_fname = [];
        temp_fname.append(cache_prefix)
        for t in time_points:
            temp_fname.append('smoothed_' + filenames[(w,t)])
        temp_fname.append(cache_prefix)
        for t in time_points:
            temp_fname.append('clabeled_' + filenames[(w,t)])
        temp_fname.append(cache_prefix)
        for t in time_points:
            temp_fname.append('labeled_tracks_' + filenames[(w,t)])
        print temp_fname
        f = open(os.path.join(cache_prefix,'tracking_harvard_filenames.txt'),'w')
        for x in temp_fname:
            f.write(x)
            f.write('\r\n')
        f.close()
        temp_fname1 = [];
        temp_fname1.append(os.path.join(exe_dir,'tracking_multiframe'))
        temp_fname1.append(os.path.join(cache_prefix,'tracking_harvard_filenames.txt'))
        temp_fname1.append(os.path.join(cwd,'Tracking_harvard_parameters.txt'))
        #subprocess.call(temp_fname1);
        subprocess.call(temp_fname1);
##    
##    ####################### Feature computation #########################
##
##    # compute features for channel 1,2 against channel 3,4
    for w in channels_to_track:
        temp_fname = [];
        temp_fname.append('0.96 0.96 4.0')
        temp_fname.append(str(len(time_points)))
        temp_fname.append('1'); # number of associated channels to compute features with
        for t in time_points:
            temp_fname.append(os.path.join(cache_prefix, 'smoothed_' + filenames[(w,t)]))
        for t in time_points:
            temp_fname.append(os.path.join(cache_prefix, 'labeled_tracks_' + filenames[(w,t)]))
        temp_fname.append('DC') # type of channel
        for t in time_points:
            temp_fname.append(os.path.join(cache_prefix, 'clabeled_' + filenames[(2,t)])) # add DC segmented files too
##        temp_fname.append('Vessel') # type of channel
##        temp_fname.append(os.path.join(cache_prefix, 'vessel_trace_' + dataset_id + '_w' + str(vessel_w) + '.tif'))
        temp_fname.append(os.path.join(cache_prefix, 'track_summary_' + dataset_id + '_w' + str(w) + '.txt'))
        temp_fname.append(os.path.join(cache_prefix, 'track_points_summary_' + dataset_id + '_w' + str(w) + '.txt'))
        f = open(os.path.join(cache_prefix,'summary_harvard_filenames.txt'),'w')
        for x in temp_fname:
            f.write(x)
            f.write('\r\n')
        f.close()
        temp_fname1 = [];
        temp_fname1.append(os.path.join(exe_dir,'summary'))
        temp_fname1.append(os.path.join(cache_prefix,'summary_harvard_filenames.txt'))
        #subprocess.call(temp_fname1);

  ####################### Rendering #########################
##  temp_fname = [];
##  colors = [];
##  colors.append('1,0,0');
##  colors.append('0,1,0');
##  colors.append('0,0,1');
##  temp_fname.append(os.path.join(exe_dir,'render_tracks'))
##  temp_fname.append(str(len(time_points)))
##  temp_fname.append(str(len(channels)))
##  for co in range(1,len(channels),1):
##      temp_fname.append(colors(co))
##  for w in channels:
##      for t in time_points:

          
