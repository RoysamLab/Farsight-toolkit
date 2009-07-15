########################################################
## FILE: FARSIGHTUTILS
##
## A MODULE WITH SOME USEFULL FUNCTIONS FOR FARSIGHT
##
########################################################

##############
## A FUNCTION TO BROWSE FOR A FILENAME
##############
def GetFilename(ftypes):

    import os
    import Tkinter,Tkconstants,tkFileDialog
    
    origdir = os.getcwd()

    root = Tkinter.Tk()
    root.withdraw()
  
    options = {}
    options['defaultextension'] = '' # couldn't figure out how this works
    options['filetypes'] = ftypes
    options['initialdir'] = '/'
    options['initialfile'] = ''
    options['parent'] = root
    options['title'] = 'Select File'

    fname = tkFileDialog.askopenfilename(**options)

    return fname
    
##############


##############
## A FUNCTION TO ALLOW ONE TO BROWSE FOR A DIRECTORY and make it
## the current working directory
##############
def SetWorkingDirectory():

    import os,sys
    import Tkinter,Tkconstants,tkFileDialog
    
    origdir = os.getcwd()
    sys.argv = ['']

    root = Tkinter.Tk()
    root.withdraw()
  
    options = {}
    options['initialdir'] = origdir
    options['mustexist'] = False
    ##options['parent'] = root
    options['title'] = 'Choose Working Directory'
    wrkdir = tkFileDialog.askdirectory(**options)
    os.chdir(wrkdir)

    print '\nWorking directory changed to: ' + wrkdir + '\n'
    
###############

##############
## A FUNCTION TO ALLOW A NEW DIRECTORY TO BE ADDED TO THE PYTHON PATH
##############
def AppendPythonPath():

    import sys
    import Tkinter,Tkconstants,tkFileDialog
    
    root = Tkinter.Tk()
    root.withdraw()
  
    options = {}
    options['initialdir'] = '/'
    options['mustexist'] = False
    options['parent'] = root
    options['title'] = 'Choose directory to add to PATH'
    newdir = tkFileDialog.askdirectory(**options)
    sys.append(newdir)

    print '\nDirectory add to PATH: ' + newdir + '\n'
    
###############

##############
## A FUNCTION TO PRINT ALL FILES IN THE CURRENT WORKING DIRECTORY
##############
def ls():
    import glob
    print glob.glob('*')
    
###############
