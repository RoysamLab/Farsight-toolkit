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
    data_directory = 'D:\\ucb dataset\\data\\051309B6CFPLATGFPinB3rd matrixlinearUnmix_median4x4movie4 track_w1'
    cwd = 'D:\\ucb dataset\\output\\ena\\run'
    exe_dir = 'D:\\farsight\\bin'
    dataset_id = 'second_'
    dirList = os.listdir(data_directory)
    list_of_filenames =[]
    channels = []
    time_points = []
    
    for fname in dirList:
        if fname.endswith('.tif') and os.path.isfile(os.path.join(data_directory,fname)):
            list_of_filenames.append(fname)

    filenames = {(0,0):'Null'}
    for fname in list_of_filenames:
        m = re.search('^(.*)_w([0-9]+).*_t([0-9]+)\.tif',fname) 
        dataset_temp = m.group(1) # set dataset_id
        channels.append(int(m.group(2))) # collect all channel values
        time_points.append(int(m.group(3))) # collect all time_points
        filenames[int(m.group(2)), int(m.group(3))] = fname # we establish a reverse dictionary lookup of the filenames from channel number and time point
        
    channels = sorted(list(set(channels)))
    time_points = sorted(list(set(time_points)))

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

    time_points = time_points[0:5] # DEBUG
    channels = [1,2]
##    ######################### Delete slices #############################
##    for w in channels:
##        for t in time_points:
##            temp_fname = [];
##            temp_fname.append(os.path.join(exe_dir,'extract_slices'))
##            temp_fname.append(os.path.join(data_directory,filenames[(w,t)]))
##            temp_fname.append('2')
##            temp_fname.append('44')
##            temp_fname.append(os.path.join(cache_prefix, 'extracted_' + filenames[(w,t)]))
##            subprocess.call(temp_fname)
##
##
##    ######################### Smooth ####################################
##    for w in channels:
##        for t in time_points:
##            temp_fname = [];
##            temp_fname.append(os.path.join(exe_dir,'anisodiff'))
##            temp_fname.append(os.path.join(cache_prefix, 'extracted_' + filenames[(w,t)]))
##            temp_fname.append('10,0.0625,2')
##            temp_fname.append(os.path.join(cache_prefix, 'smoothed_' + filenames[(w,t)]))
##            subprocess.call(temp_fname)
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
##
##    nuclei_segmentation_cfg = 'CellSegmentation.ini'
##    params = [];
##    params.append('15,50,10');
##    params.append('2,50,5');
##
##    channels_to_segment = [1]
##    # segment all the channels first. Then trace the vessels.
##    for t  in time_points:
##        for w in channels_to_segment:
##            temp_fname = [];
##            temp_fname.append(os.path.join(exe_dir,'segmentation'))
##            temp_fname.append(os.path.join(cache_prefix, 'smoothed_' + filenames[(w,t)]))
##            if w == 1 or w == 2:
##                # call yousef segmentation
##                temp_fname.append('1') # type of segmentation
##                temp_fname.append(nuclei_segmentation_cfg)# parameters filename
##            else:
##                # call simple binary morphological operators based segmentation
##                temp_fname.append('2') # type of segmentation
##                temp_fname.append(params[int(w)-3]) # paramerters
##            temp_fname.append(os.path.join(cache_prefix, 'labeled_' + filenames[(w,t)]))
##            subprocess.call(temp_fname)

    # vessel tracing
##    vessel_w = 4; #vessel channel
##    temp_fname = [];
##    temp_fname.append(os.path.join(exe_dir,'vessel_thinning'))
##    for t in time_points:
##        temp_fname.append(os.path.join(cache_prefix, 'labeled_' + filenames[(vessel_w,t)]))
##    temp_fname.append(os.path.join(cache_prefix, 'vessel_binarized_' + dataset_id + '_w' + str(vessel_w) + '.tif'))
##    temp_fname.append(os.path.join(cache_prefix, 'vessel_trace_'+ dataset_id + '_w' + str(vessel_w) + '.tif'))
##    subprocess.call(temp_fname) 

    ############################# Tracking ##############################

    # track channel 1 and 2 only
    channels_to_track = [1,2]
    for w in channels_to_track:
        temp_fname = [];
        temp_fname.append(os.path.join(exe_dir,'tracking'))
        for t in time_points:
            temp_fname.append(os.path.join(cache_prefix, 'smoothed_' + filenames[(w,t)]))
        for t in time_points:
            temp_fname.append(os.path.join(cache_prefix, 'labeled_' + filenames[(w,t)]))
        for t in time_points:
            temp_fname.append(os.path.join(cache_prefix, 'labeled_tracks_' + filenames[(w,t)]))
        subprocess.call(temp_fname);
    
    ####################### Feature computation #########################

    # compute features for channel 1,2 against channel 3,4
    for w in channels_to_track:
        temp_fname = [];
        temp_fname.append(os.path.join(exe_dir,'summary'))
        temp_fname.append(str(len(time_points)))
        temp_fname.append('0'); # number of associated channels to compute features with
        for t in time_points:
            temp_fname.append(os.path.join(cache_prefix, 'smoothed_' + filenames[(w,t)]))
        for t in time_points:
            temp_fname.append(os.path.join(cache_prefix, 'labeled_tracks_' + filenames[(w,t)]))
##        temp_fname.append('DC') # type of channel
##        for t in time_points:
##            temp_fname.append(os.path.join(cache_prefix, 'labeled_' + filenames[(3,t)])) # add DC segmented files too
##        temp_fname.append('Vessel') # type of channel
##        temp_fname.append(os.path.join(cache_prefix, 'vessel_trace_' + dataset_id + '_w' + str(vessel_w) + '.tif'))
        temp_fname.append(os.path.join(cache_prefix, 'track_summary_' + dataset_id + '_w' + str(w) + '.txt'))
        temp_fname.append(os.path.join(cache_prefix, 'track_points_summary_' + dataset_id + '_w' + str(w) + '.txt'))
        subprocess.call(temp_fname);
