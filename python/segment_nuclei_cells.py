import os, sys, platform, subprocess, multiprocessing, time

def seg_nuclei(basename_list):
    command_list=[]
    list=[]
    launched_processes_list=[]
    file_handle_list=[]
    
    idx = 0

    for line in basename_list:
        command_list.append("segment_nuclei " + basename_list[idx].rstrip('\n') + ".mhd " + basename_list[idx].rstrip('\n') + "_label.mhd Seg_Params.ini")
        idx = idx + 1

    #make "debug" directory to store the console output of each subprocess
    if (not (os.path.exists(os.getcwd() + '/debug/'))):
        os.makedirs(os.getcwd() + '/debug/')

    #spawn and add subprocesses (each representing an execution of register_pair.exe) into a list
    idx = 0
    threads_launched = 0
    for command in command_list:
        #while loop to check active processes to see if we can launch more threads
        while threads_launched >= multiprocessing.cpu_count():
        #while threads_launched >= 4:
            threads_launched = 0
            still_launched_subp_list = []
            for launched_subp in launched_processes_list:                  
                if launched_subp.poll() == None:
                    threads_launched = threads_launched + 1
                    still_launched_subp_list.append(launched_subp)
            launched_processes_list = still_launched_subp_list[:]
            #print "Threads launched: " + str(threads_launched)
            time.sleep(1)
                          
        #fh = open(os.getcwd() + '/debug/debug_' + basename_list[idx].rstrip('\n') + ".txt", 'w')
        print "Launching " + command
        subp = subprocess.Popen(command + " > " + os.getcwd() + '/debug/debug_' + basename_list[idx].rstrip('\n') + ".txt 2>&1" , shell=True)
        #subp = subprocess.Popen(command)
        
        list.append(subp)
        launched_processes_list.append(subp)
        #file_handle_list.append(fh)
        threads_launched = threads_launched + 1;
        idx = idx + 1

    #loop through entire list to check if process are done
    processes_done = False
    while not processes_done:
        processes_done = True
        time.sleep(1) #sleep 1 second before polling all processes
        for subp in list:
            if subp.poll() == None:
                processes_done = False

    #close all the debug file handles since we are done writing to them
    #for fh in file_handle_list:
        #fh.close()
    
    

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print 'Usage: '+sys.argv[0]+' <cell_list.txt>\n'
        sys.exit(1)
    start_time = time.clock()
    basename_list = open(sys.argv[1], 'r').readlines()
    seg_nuclei(basename_list)
    print 'seg_nuclei took: ' + str(time.clock() - start_time) + ' seconds'
