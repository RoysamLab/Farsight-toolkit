#!/usr/bin/python

# Run this script from a folder with the raw data coming out of Bill's
# lab. If you do not know who Bill is, you'll probably have to walk through and 
# reorganize the code to suit your new needs.

import getopt,sys, glob, os, shutil, platform, subprocess
import register_pairs_parallel

width=7
colors=None
regisColor=None
anchorImage=0

#scriptsfolder='/FSdata/bbusse/ftk/src/DARPAScripts/montaging' 
scriptsfolder='/data/nicolas/src-bin/bin/farsight-rel/exe/' #The folder with all the registration scripts and stuff
gdbicpFolder='/data/nicolas/src-bin/src/farsight/trunk/DARPAScripts/montaging/'
#execfolder='/FSdata/bbusse/ftk/release/Farsight/exe' #The folder with all the executables
imgpairname='imgpairs.txt' #keeping windows names so as not to confuse people
index_length=5 #number of digits in image num, used for leading zeros

#initialize executables for platform
if platform.system() == 'Windows':
    regpair = 'register_pair.exe'
    regjoint = 'register_joint.exe'
    mosaicim = 'mosaic_images.exe'
    gdbicp = 'gdbicp.exe'
    ftkMainDarpa = 'ftkMainDarpa.exe'
else:
    regpair = 'register_pair_16'
    regjoint = 'register_joint_16'
    mosaicim = 'mosaic_images_16'
    gdbicp = 'gdbicp'
    ftkMainDarpa = 'ftkMainDarpa'

def setupVars():
    #estimate colors, if not specified already
    global colors
    if not colors: #this bit is mostly just a lot of globs and filtering
        alltiles=glob.glob('*dsu.tif')
        alltiles.extend(glob.glob('*dsu.TIF'))
        if alltiles==[]: #no images in folder, grab all subfolders as colors (probably already split up)
            colors = os.listdir('.')
	    print "PP"
	    #print colors
	    newColor = []
            for color in colors:
	    	#print color
                if 'unused' == color: continue
                if not os.path.isdir(color): continue
		newColor.append(color)
		#print newColor
	    colors = newColor	
        else:
            if not 'kt' in alltiles[0] or not '_w' in alltiles[0]:
                print "all folder images not in Shain lab format (ktNUM_wCHANNELdsu.TIF)"
                sys.exit() 
            firstnum= alltiles[0].lstrip('kt').partition('_w')[0]
            if not firstnum.isdigit():
                print "all folder images not in Shain lab format (ktNUM_wCHANNELdsu.TIF)"
                sys.exit() 
            firsttile=glob.glob('kt'+firstnum+'*dsu.tif')
            firsttile.extend(glob.glob('kt'+firstnum+'*dsu.TIF'))
             
            colors=[]#now we can do the colors themselves
            for t in firsttile:
                colors.append(t.partition('dsu')[0].partition('_w')[2].lstrip('0123456789'))

    print colors
    global regisColor
    if (not regisColor) and colors: #if not regiscolor, grab the first color        
        regisColor = colors[0]

def prepFolders():    
    #move tiles, make imgpairlist
    for c in colors: #move tiles
        try: 
            os.mkdir(c) 
        except: 
            pass
        tiles=glob.glob('*'+c+'dsu.tif')
        tiles.extend(glob.glob('*'+c+'dsu.TIF'))
        for t in tiles:
            shutil.move(t,c)
 
    try:  #move everything else
        os.mkdir('unused')   
    except:
        pass     
    tiles=glob.glob('*.*')+glob.glob('*ftk*')
    print "todo lo demas"
    print tiles
    for t in tiles:
        if not "register" in t:
            shutil.move(t,'unused')

    os.chdir(regisColor)

    prefix='kt'
    registiles=glob.glob('*dsu.tif')
    registiles.extend(glob.glob('*dsu.TIF'))
    suffix='_w'+registiles[0].partition('_w')[2]
    
    registiles.sort()

    first=registiles[0].partition('_w')[0].lstrip('kt')
    last=registiles[-1].partition('_w')[0].lstrip('kt')

    tst='0'+str(index_length)+'d'
    f=open(imgpairname,'w') #make image pair list    

    print first
    print last
    for i in range((int(last))-int(first)):
	#print i
    	if ((i+1)%width != 0) :
		#print i
	        f.write((prefix+'%(first)'+tst+suffix+' '+prefix+'%(second)'+tst+suffix+'\n') % {'first': int(first)+i+1, 'second': int(first)+i})
        if i+1 >= width:
            f.write((prefix+'%(first)'+tst+suffix+' '+prefix+'%(second)'+tst+suffix+'\n') % {'first': int(first)+i+1, 'second': int(first)+i+1-width})
    f.close()

    os.chdir('..')


def register():
    #copy over any necessary registration stuff
    shutil.copy(os.path.join(gdbicpFolder,gdbicp), regisColor)
    shutil.copy(os.path.join(scriptsfolder,regpair), regisColor)
    shutil.copy(os.path.join(scriptsfolder,regjoint), regisColor)

    os.chdir(regisColor)
    register_pairs_parallel.register(['./',imgpairname])
    os.chdir('..')
    
def mosaic():
    #copy joint transform files from regiscolor folder to the rest, modify them
    os.chdir(regisColor)
    tiles=glob.glob('*dsu.tif')
    tiles.extend(glob.glob('*dsu.TIF'))
    tiles.sort()
    regisSuffix='_w'+tiles[0].partition('_w')[2]
    os.chdir('..')
    for color in colors:
        if color != regisColor:
            os.chdir(color)
            tiles=glob.glob('*dsu.tif')
            tiles.extend(glob.glob('*dsu.TIF'))
            tiles.sort()
            suffix='_w'+tiles[0].partition('_w')[2]    
            print regisSuffix, suffix        
            os.chdir('..')
            f=open(os.path.join(regisColor,'joint_transforms.xml'),'r')
            o=open(os.path.join(color,'joint_transforms.xml'),'w')
            for line in f:
                o.write(line.replace(regisSuffix,suffix))
            f.close()
            o.close()

    for color in colors: #run mosaicing
	if not "NIC" in color:
            print "\tNOW IT WILL MOSAIC"
            os.chdir(color)
            tiles=glob.glob('*dsu.tif')
            tiles.extend(glob.glob('*dsu.TIF')) #yes I'm globbing three times, deal with it
            tiles.sort()
            shutil.copy(os.path.join(scriptsfolder,mosaicim), '.')
      	    print "\tThe anchor IMAGE WILL BE: "+tiles[anchorImage]
            #os.system('' if platform.system() == 'Windows' else './' + mosaicim + ' joint_transforms.xml '+tiles[anchorImage]+ ' -3d'+' -blending 0')
	    TEMP2 = './' + mosaicim + ' joint_transforms.xml '+tiles[anchorImage]+ ' -3d'+' -blending 3'
            print TEMP2
            TEMP3 = subprocess.Popen(TEMP2, shell=True)
            TEMP3.communicate()
	    #print './' + mosaicim + ' joint_transforms.xml '+tiles[anchorImage]+ ' -3d'+' -blending 0'
            os.chdir('..')

def usage():
    print '''Usage:  register_automagic.py [Options]
    -w, --width W  How many horizontal tiles exist (default 7)
    -c, --colors C1,C2,etc  Channel names (otherwise attempt to figure them out)
    -r, --regisColor C1  Channel used for registration computation
    -m, --mosaicOnly  Assume the tiles are already registered, and just stitch the full montage
    -h, --help  Print this message and exit
    -anchor, is the anchor image

    This registration process is designed to cater to the Shain Lab naming convention.  If "kt01234_w311DAPIdsu.TIF" does not make sense to you, you will probably need to examine this code by hand, because there will certainly be naming conventions that need changing at the very least, along with other reorganization.'''

def main(argv):  
    reg=True      
    prep=True
    mos=True                   
    try:                                
        opts, args = getopt.getopt(argv, "a:w:c:hmr:d", ['anchorImage=',"width=", "colors=",'help','mosaicOnly','regisColor='])
    except getopt.GetoptError:           
        usage()                          
        sys.exit(2)  
    print opts
    for opt,arg in opts:               
        if opt in ("-h", "--help"):      
            usage()                     
            sys.exit()                  
        elif opt in ("-w", "--width"):   
            global width                           
            width = int(arg)
	elif opt in ("-a", "--anchor"):
	    global anchorImage
	    anchorImage = int(arg)
        elif opt in ("-r", "--regisColor"):   
            global regisColor           
            regisColor = arg
        elif opt in ("-c", "--colors"): 
            global colors
            colors = arg.split(',')
        elif opt in ('-m','--mosaicOnly'):
            reg=False
            prep=False
    print "\t\t ANCHOR: "+str(anchorImage) 
    print "\t\t regisColor: "+str(regisColor)
    setupVars()
    if prep: prepFolders()         
    if reg: register()
    if mos: mosaic()

if __name__ == "__main__":
    head= os.path.split(sys.argv[0])[0]
    #scriptsfolder = head if head else '.'
    main(sys.argv[1:])

