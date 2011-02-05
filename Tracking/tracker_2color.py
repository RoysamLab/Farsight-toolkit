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
    data_directory = 'L:\\\Tracking\\\data\\p14neg\\movie5c\\'
    cwd = 'L:\\Tracking'
    exe_dir = 'C:\\Users\\Arun\\Research\\Farsight\\exe\\bin'
    dataset_id = 'p14neg'
    number_of_cores = 7;
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

    time_points = time_points # DEBUG
    channels = channels;
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
    popen_objs = [];
    for t in time_points:
        temp_fname = [];
        temp_fname.append(os.path.join(exe_dir,'unmix_satfix_smoothassign'))
        
        #first add input filenames then add output filenames
        temp_fname.append(str(len(channels)))
        temp_fname.append(str(len(channels)))
        for w in channels:
           temp_fname.append(os.path.join(data_directory,filenames[(w,t)]))
        for w in channels:
            temp_fname.append(os.path.join(cache_prefix,'unmixed_' + filenames[(w,t)]))
        popen_objs.append(subprocess.Popen(temp_fname))
        if len(popen_objs) == number_of_cores:
                for obj in popen_objs:
                    out = obj.communicate()[0]    
                popen_objs = []
    for obj in popen_objs:
        out = obj.communicate()[0]       
        
    ######################## Segmentation ###############################

    nuclei_segmentation_cfg = '2color_CellSegmentation.ini'
    params = [];
    params.append('15,50,10');
    params.append('2,50,5');

    vessel_w = 4;
    dc_w = 1;
    channels_to_segment = channels
    popen_objs = [];
    # segment all the channels first. Then trace the vessels.
    for t  in time_points:
        for w in channels_to_segment:
            temp_fname = [];
            output_filename = os.path.join(cache_prefix, 'labeled_' + filenames[(w,t)])
            if os.path.exists(output_filename): # dont do anything if the output file exists
                continue;
            if w == 2:
                # call yousef segmentation
                temp_fname.append(os.path.join(exe_dir,'segment_nuclei_harvard'))
                temp_fname.append(os.path.join(cache_prefix, 'unmixed_' + filenames[(w,t)]))
                temp_fname.append(output_filename)
                temp_fname.append(nuclei_segmentation_cfg)# parameters filename
            else:
                # call simple binary morphological operators based segmentation
                temp_fname.append(os.path.join(exe_dir,'segmentation'))
                temp_fname.append(os.path.join(cache_prefix, 'unmixed_' + filenames[(w,t)]))
                temp_fname.append('2') # type of segmentation
                temp_fname.append(params[int(w)]) # paramerters
                temp_fname.append(output_filename)
            popen_objs.append(subprocess.Popen(temp_fname));
            if len(popen_objs) == number_of_cores:
                for obj in popen_objs:
                        out = obj.communicate()[0];
                popen_objs = []
    for obj in popen_objs:
        out = obj.communicate()[0];

    # vessel tracing
     #vessel channel
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
            temp_fname.append('40')
            temp_fname.append(os.path.join(cache_prefix, 'clabeled_' + filenames[(w,t)]))
            #pdb.set_trace()
            popen_objs.append(subprocess.Popen(temp_fname))
            if len(popen_objs) == number_of_cores :
                for obj in popen_objs:
                    out = obj.communicate()[0]
                popen_objs = []
    for obj in popen_objs:
        out = obj.communicate()[0]   
    ############################# Tracking ##############################

    # track channel 1 and 2 only
    channels_to_track = [2]
    for w in channels_to_track:
        temp_fname = [];
        temp_fname.append(cache_prefix)
        for t in time_points:
            temp_fname.append('unmixed_' + filenames[(w,t)])
        temp_fname.append(cache_prefix)
        for t in time_points:
            temp_fname.append('clabeled_' + filenames[(w,t)])
        temp_fname.append(cache_prefix)
        for t in time_points:
            temp_fname.append('labeled_tracks_' + filenames[(w,t)])
        print temp_fname
        f = open(os.path.join(cache_prefix,'tracking_berkeley_filenames.txt'),'w')
        for x in temp_fname:
            f.write(x)
            f.write('\r\n')
        f.close()
        temp_fname1 = [];
        temp_fname1.append(os.path.join(exe_dir,'tracking_multiframe'))
        temp_fname1.append(os.path.join(cache_prefix,'tracking_berkeley_filenames.txt'))
        temp_fname1.append(os.path.join(cwd,'tracking_2color_parameters.txt'))
        #subprocess.call(temp_fname1);
        subprocess.call(temp_fname1);
    
    ####################### Feature computation #########################

    for w in channels_to_track:
        temp_fname = [];
        #temp_fname.append(os.path.join(exe_dir,'summary'))
        temp_fname.append('0.64 0.64 2.0')
        temp_fname.append(str(len(time_points)))
        temp_fname.append('1'); # number of associated channels to compute features with
        for t in time_points:
            temp_fname.append(os.path.join(cache_prefix, 'unmixed_' + filenames[(w,t)]))
        for t in time_points:
            temp_fname.append(os.path.join(cache_prefix, 'labeled_tracks_' + filenames[(w,t)]))
        temp_fname.append('DC') # type of channel
        for t in time_points:
            temp_fname.append(os.path.join(cache_prefix, 'clabeled_' + filenames[(dc_w,t)])) # add DC segmented files too
        #temp_fname.append('Vessel') # type of channel
        #temp_fname.append(os.path.join(cache_prefix, 'vessel_trace_' + dataset_id + '_w' + str(vessel_w) + '.tif'))
        temp_fname.append(os.path.join(cache_prefix, 'track_summary_' + dataset_id + '_w' + str(w) + '.txt'))
        temp_fname.append(os.path.join(cache_prefix, 'track_points_summary_' + dataset_id + '_w' + str(w) + '.txt'))
        f = open(os.path.join(cache_prefix,'summary_berkeley_filenames.txt'),'w')
        for x in temp_fname:
            f.write(x)
            f.write('\r\n')
        f.close()
        temp_fname1 = [];
        temp_fname1.append(os.path.join(exe_dir,'summary'))
        temp_fname1.append(os.path.join(cache_prefix,'summary_berkeley_filenames.txt'))
        subprocess.call(temp_fname1);

    ######################## Segmentation feature computation ##############


    for w in channels_to_track:
        temp_fname = [];
        #temp_fname.append(os.path.join(exe_dir,'summary'))
        temp_fname.append('0.64 0.64 2.0')
        temp_fname.append(str(len(time_points)))
        temp_fname.append('1'); # number of associated channels to compute features with
        for t in time_points:
            temp_fname.append(os.path.join(cache_prefix, 'unmixed_' + filenames[(w,t)]))
        for t in time_points:
            temp_fname.append(os.path.join(cache_prefix, 'clabeled_' + filenames[(w,t)]))
        temp_fname.append('DC') # type of channel
        for t in time_points:
            temp_fname.append(os.path.join(cache_prefix, 'clabeled_' + filenames[(dc_w,t)])) # add DC segmented files too
        #temp_fname.append('Vessel') # type of channel
        #temp_fname.append(os.path.join(cache_prefix, 'vessel_trace_' + dataset_id + '_w' + str(vessel_w) + '.tif'))
        temp_fname.append(os.path.join(cache_prefix, 'segmentation_summary_' + dataset_id + '_w' + str(w) + '.txt'))
        temp_fname.append(os.path.join(cache_prefix, 'segmentation_points_summary_' + dataset_id + '_w' + str(w) + '.txt'))
        f = open(os.path.join(cache_prefix,'summary_berkeley_filenames.txt'),'w')
        for x in temp_fname:
            f.write(x)
            f.write('\r\n')
        f.close()
        temp_fname1 = [];
        temp_fname1.append(os.path.join(exe_dir,'summary'))
        temp_fname1.append(os.path.join(cache_prefix,'summary_berkeley_filenames.txt'))
        subprocess.call(temp_fname1);
        
