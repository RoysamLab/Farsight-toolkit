#!/usr/bin/env python
# -*- coding: utf-8 -*-

import shutil, fnmatch, os, subprocess, os.path, time, glob, sys, inspect, filecmp, datetime, getpass, multiprocessing
#from subprocess 


def main(subprocess_command_list):
  #print subprocess_command_list
  

  #spawn and add subprocesses (each representing an execution of register_pair.exe) into a list
  subprocess_list=[]
  file_handle_list=[]
  launched_processes_list=[]
  launched_subp_list=[]
  still_launched_subp_list=[]

  idx = 0
  threads_launched = 0
  contador = 0
  flagAllStart = 0
  for subprocess_command in subprocess_command_list:
    #while loop to check active processes to see if we can launch more threads
    process = subprocess.Popen('free -k|grep Mem:', shell=True, stdout=subprocess.PIPE)
    lista = process.communicate()[0].split()

    memoryAvilable = 1
    if int(lista[3]) < 150000000:
      memoryAvilable = 0
      memoryAvilable = 1
    print lista[3]+str(memoryAvilable)	

    while (threads_launched >= max(1,multiprocessing.cpu_count() - 70) ) | (memoryAvilable == 0): # at least one thread, but no more than 3/4 of the machine
    #while (threads_launched >= max(1,1) ) | (memoryAvilable == 0): # at least one thread, but no more than 3/4 of the machine
      if flagAllStart == 0:
        flagAllStart = 1
        print "\tALL THE PROCESS ARE STARTED"
      threads_launched = 0
      still_launched_subp_list = []
      for launched_subp in launched_processes_list:                  
        if launched_subp.poll() == None:
          threads_launched = threads_launched + 1
          still_launched_subp_list.append(launched_subp)
      launched_processes_list = still_launched_subp_list[:]
      #print "Threads launched: " + str(threads_launched)
      time.sleep(0.1)
      process = subprocess.Popen('free -k|grep Mem:', shell=True, stdout=subprocess.PIPE)
      lista = process.communicate()[0].split()
      memoryAvilable = 1
      if int(lista[3]) < 150000000:
        memoryAvilable = 0
        memoryAvilable = 1
    #fh = open(os.getcwd() + '/debug/debug_from_' + from_image_list[idx] + '_to_' + to_image_list[idx] + '.txt', 'w')
    contador = contador + 1
    print "Launching " + subprocess_command + ", "+str(contador)+" of "+str(len(subprocess_command_list))
    subp = subprocess.Popen(subprocess_command, shell=True)
    #subp = subprocess.Popen(subprocess_command)
    
    subprocess_list.append(subp)
    launched_processes_list.append(subp)
    #file_handle_list.append(fh)
    threads_launched = threads_launched + 1;
    idx = idx + 1
    time.sleep(0.1)
    
  #loop through entire list to check if process are done
  processes_done = False
  while not processes_done:
    processes_done = True
    time.sleep(1) #sleep 1 second before polling all processes
    for subp in subprocess_list:
      if subp.poll() == None:
        processes_done = False

   

if __name__ == "__main__":
  main()
