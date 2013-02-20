import ij
import os
import math
 
global foldername
foldername=None
tilesize=1024

##############################
''' Start of actual script '''
##############################

#fd = OpenDialog("Save Tiles In",None)
prefixval = 'Chan3_'
fd= DirectoryChooser("Save Tiles In")
foldername= fd.getDirectory()
if foldername:

	#turn off IJ gui stuff, makes it a bit faster but less pretty
	#Interp = ij.macro.Interpreter()
	#Interp.batchMode = True

	im=ij.IJ.getImage() #big montage image we want to work on
	im.show() #make sure the thing's on top

	#make a dir to put the tiles into
	try:
		os.mkdir(foldername+"/"+im.getTitle()+"_Tiles/")
	except:
		pass #it probably already exists

	tx= math.ceil(1.0*im.getWidth()/tilesize) #number of x-tiles
	ty= math.ceil(1.0*im.getHeight()/tilesize) #number of y-tiles

	maxz= im.getStackSize() #number of slices

	numtiles=tx*ty #how many tiles we'll have in total
	current_tilenum = 1 #current tile we're on
	IJ.showProgress(0.0)

	IJ.setBackgroundColor(0,0,0)
	
	#document everything in a farsight xml format file
	xml=open(foldername+"/"+im.getTitle()+"_Tiles/LoadNuclei.xml", "w")
	xml.write('<?xml version="1.0" ?>\n')
	xml.write('<Source>\n')

	for y in range(ty):
		for x in range(tx):
			#limit our x and/or y if the tile is an edge tile
			dx=min(tilesize,im.getWidth()-(tilesize*x+1))
			dy=min(tilesize,im.getHeight()-(tilesize*y+1))
			print (tilesize*x+1,tilesize*y+1,dx,dy)
			r=ij.gui.Roi(tilesize*x+1,tilesize*y+1,dx,dy) #define our tile ROI: x,y,w,h
					
			#set the roi
			im.setRoi(r)
			
			#extract the tile
			tilexy = ij.plugin.Duplicator().run(im,1,maxz)
			tilexy.getProcessor().setColor(0)
			
			#maybe help keep some memory clear, who knows if this does anything
			im.trimProcessor()
			
			#pad the tile if need be
			if tilexy.getWidth() < tilesize or tilexy.getHeight() < tilesize:
				print "padding tile"
				resizer=ij.plugin.CanvasResizer()
				newstack = resizer.expandStack(tilexy.getImageStack(),tilesize,tilesize,0,0)
				tilexy.setStack(newstack)
			
			#save the tile
			IJ.saveAs(tilexy,"Tiff",foldername+"/"+im.getTitle()+"_Tiles/"+str(prefixval)+str(current_tilenum)+"_"+str(int(tilesize*x))+"_"+str(int(tilesize*y)))
			
			#append filenames to xml document
			xml.write('\t<File FileName="'+foldername+"/"+im.getTitle()+'_Tiles/'+prefixval+str(current_tilenum)+"_"+str(int(tilesize*x))+"_"+str(int(tilesize*y))+'.tif" Type="Image" tX="'+str(int(tilesize*x))+'" tY="'+str(int(tilesize*y))+'" tZ="0"/>\n')
			xml.write('\t<File FileName="'+foldername+"/"+im.getTitle()+'_Tiles/'+prefixval+str(current_tilenum)+"_"+str(int(tilesize*x))+"_"+str(int(tilesize*y))+'_table.txt" Type="Nuclei_Table" tX="'+str(int(tilesize*x))+'" tY="'+str(int(tilesize*y))+'" tZ="0"/>\n')
			
			#delete tile to free up memory
			tilexy.close()
			
			current_tilenum += 1
			IJ.showProgress(1.0*current_tilenum/numtiles)

	#finish up xml
	xml.write('</Source>\n')
	xml.close()
	
	#Interp.batchMode = False
	#print "Done"


