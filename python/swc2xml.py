import sys

if (len(sys.argv) > 1):
    filename=sys.argv[1]
else:
    print "Usage: swc2xml.py cell_list.txt"
    print "Ask Ho about formatting of cell_list and intermediate files"
    print "If you don't know Ho, tough luck"
    quit()

#document everything in a farsight xml format file
xml=open("LoadNuclei.xml", "w")
xml.write('<?xml version="1.0" ?>\n')
xml.write('<Source>\n')

cellsfile = open(filename, 'r')
cells = cellsfile.readlines()
for cell in cells:
    locfile= open(cell.strip()+".txt",'r')
    glbcoord=locfile.readline().split(' ')
    loccoord=locfile.readline().split(' ')
    offset=[int(glbcoord[0])-int(loccoord[0]),int(glbcoord[1])-int(loccoord[1]),int(glbcoord[2])-int(loccoord[2])]
	
    xml.write('\t<File FileName="'+cell.strip()+'_ANT.swc" Type="Trace" tX="'+str(offset[0])+'" tY="'+str(offset[1])+'" tZ="'+str(offset[2])+'"/>\n')
	
#finish up xml
xml.write('</Source>\n')
xml.close()