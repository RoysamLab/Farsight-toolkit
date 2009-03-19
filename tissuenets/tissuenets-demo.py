import os
import sys

from Tkinter import *
from tkFileDialog   import asksaveasfilename

from tkFileDialog   import askopenfilename
from tkMessageBox import *

#Adds n empty rows starting from line rowno in self which might be a frame
def emptyRows(n, rowno,self):
      for i in range(n):
         Label(self,text = "").grid(row=rowno+i, column=0,columnspan = 2, sticky =W)

fileName1 = ""
fileName2 = ""
fileName3 = ""
fileName4 = ""
flag = 0 #

class dlgCreateNode(Toplevel):
   def __init__(self, master):

      master.title('Create Nodes')
      Toplevel.__init__(self, master)
      #This window and its parent are required to be related
      self.transient(master)
      self.master=master
      dlgNodeFrame = Frame(self,width=50, height=60, bd=1)
      dlgNodeFrame.pack()
      #call grabset to make dialog modal= Only 1 dialog can appear on the screen
      #We don't want to confuse users
      self.grab_set()

      self.geometry("+%d+%d" % (master.winfo_rootx()+350,
                        master.winfo_rooty()+120))

      Label(dlgNodeFrame,text = "").grid(row=0, column=0,columnspan = 2, sticky =W)
              
      Label(dlgNodeFrame,
            text = "Input1: Filename for reading all objects :"
            ).grid(row=1, column=0,columnspan = 2, sticky =W)
      self.txtNodeInFileName=dlgNodeFrame.txtNodeInFileName = Entry(dlgNodeFrame)
      dlgNodeFrame.txtNodeInFileName.grid(row=1, column= 3, sticky = W)

      Button(dlgNodeFrame,
             text= "Browse...",
             command = self.cmdBrowse2
             ).grid(row=1,column=4,sticky=W)

      #Add 1 empty row starting from row 2
      emptyRows(1,2,dlgNodeFrame)

      Label(dlgNodeFrame,
            text = "Output: Enter the Node File name to be created :"
            ).grid(row=3, column=0,columnspan = 2, sticky =W)
      self.txtNodeOutFileName=dlgNodeFrame.txtNodeOutFileName = Entry(dlgNodeFrame)
      dlgNodeFrame.txtNodeOutFileName.grid(row=3, column= 3, sticky = W)

      Button(dlgNodeFrame,
             text= "Browse...",
             command = self.cmdBrowse1
             ).grid(row=3,column=4,sticky=W)  

      #Add 1 empty row starting from row 4
      emptyRows(1,4,dlgNodeFrame)
      
      # Ask for Object type
      Label(dlgNodeFrame,
            text = " Input2:Object Type :"
            ).grid(row=5, column=0, sticky =W)
      self.txtObjectType=dlgNodeFrame.txtObjectType = StringVar()
      
      allObjects = ["Astrocytes", "Microglia", "Neurons","Endothelials","Electrode"]
      dlgNodeFrame.txtObjectType.set(2)
      column = 2
      row=6
      for fobject in allObjects:
         self.objectTypes=Radiobutton(dlgNodeFrame,
                     text=fobject,
                     variable=dlgNodeFrame.txtObjectType,
                     value=fobject
                     ).grid(row=row,column=0,sticky=W)
         row+=1


      #DEMO
      self.txtNodeInFileName.delete(0,END)  
      self.txtNodeInFileName.insert(0,"SupplementD.xml")

      self.txtNodeOutFileName.delete(0,END)  
      self.txtNodeOutFileName.insert(0,"Neurons.xml")
      dlgNodeFrame.txtObjectType.set("Neurons")
      #DEMO
      

      #Add 1 empty row starting from row 2
      #emptyRows(1,2,dlgNodeFrame)
      #Label(dlgNodeFrame,text = "").grid(row=4, column=0,columnspan = 2, sticky =W)

      #Add 4 empty rows starting from row 5
      emptyRows(4,5,dlgNodeFrame)
      #photo1 = PhotoImage(file="b3.gif")
         
      self.buttonCreate= Button(dlgNodeFrame,
             text= "  Create Nodes...  ",
             #image=photo1,                   
             command = self.cmdCreateNodes
             ).grid(row=11,column=1,sticky=W)
      # Do this for garbage collection
      #self.buttonCreate.image=photo1
      #button1.pack(side=tk.LEFT, padx=2, pady=2)
 
      self.buttonCancel=Button(dlgNodeFrame,
             text= "  Cancel  ",
             command = self.cmdQuit
             ).grid(row=11,column=3,sticky=W)

      dlgNodeFrame.pack(side=TOP,fill=X)

                   
   def cmdBrowse1(self):
      fileName1=asksaveasfilename()
      if fileName1:
            self.txtNodeOutFileName.delete(0,END)
            self.txtNodeOutFileName.insert(0,fileName1)
   
   def cmdBrowse2(self):
      fileName2 = askopenfilename()
      if fileName2:
            self.txtNodeInFileName.delete(0,END)
            self.txtNodeInFileName.insert(0,fileName2)
      

   def cmdCreateNodes(self):

      # Get the object types from the radio button group
      objType = ""
      objType = self.txtObjectType.get()

      fileName1=self.txtNodeOutFileName.get()
      fileName2=self.txtNodeInFileName.get()

      if (objType != "") and (fileName1 != "") and (fileName2 != ""):
            print "Creating "+ objType + " Nodes ..."
            javaCommand = "java -Xmx1024m net.sf.saxon.Query -o:"+fileName1+" "
            javaCommand = javaCommand + "defineNodes.xquery allnodes=" + fileName2+ " "
            javaCommand = javaCommand + "objecttype=" + objType
            os.system(javaCommand)
            showinfo('Warning!', 'Nodes are created. Output is written to: '+fileName1)
      else:
            showerror('Error!', 'Inputs are Incomplete')
                         
   def cmdQuit(self):
      self.master.focus_set()
      self.destroy()

##########################################################################################
##########################################################################################
# Usage: python defineLinks output nodes1 nodes2 numofneighbors

class dlgCreateLinks(Toplevel):
   def __init__(self, master):

      master.title('Create Links')
      Toplevel.__init__(self, master)
      #This window and its parent are required to be related
      self.transient(master)
      self.master=master
      dlgNodeFrame = Frame(self,width=50, height=60, bd=1)
      #dlgNodeFrame.title('Create Links')
      dlgNodeFrame.pack()

      #call grabset to make dialog modal= Only 1 dialog can appear on the screen
      #We don't want to confuse users
      self.grab_set()

      self.geometry("+%d+%d" % (master.winfo_rootx()+350,
                        master.winfo_rooty()+120))

      Label(dlgNodeFrame,text = "").grid(row=0, column=0,columnspan = 2, sticky =W)

      fileName1=""
      fileName2=""
      fileName3=""
      Label(dlgNodeFrame,
            text = "Input1: Filename for reading nodes of specific object type :"
            ).grid(row=1, column=0,columnspan = 2, sticky =W)
      self.txtNodeInFileName=dlgNodeFrame.txtNodeInFileName = Entry(dlgNodeFrame)
      dlgNodeFrame.txtNodeInFileName.grid(row=1, column= 3, sticky = W)

      Button(dlgNodeFrame,
             text= "Browse...",
             command = self.cmdBrowse2
             ).grid(row=1,column=4,sticky=W)

      #Add 1 empty row starting from row 2
      emptyRows(1,2,dlgNodeFrame)
 
      Label(dlgNodeFrame,
            text = " Input2: Filename for reading nodes of specific object type :"
            ).grid(row=3, column=0,columnspan = 2, sticky =W)
      self.txtNodeInFileName2=dlgNodeFrame.txtNodeInFileName2 = Entry(dlgNodeFrame)
      dlgNodeFrame.txtNodeInFileName2.grid(row=3, column= 3, sticky = W)

      Button(dlgNodeFrame,
             text= "Browse...",
             command = self.cmdBrowse3
             ).grid(row=3,column=4,sticky=W)

      #Add 1 empty row starting from row 4
      emptyRows(1,4,dlgNodeFrame)
      
      Label(dlgNodeFrame,
            text = "Output: Enter the Links File name to be created :"
            ).grid(row=5, column=0,columnspan = 2, sticky =W)
      self.txtNodeOutFileName=dlgNodeFrame.txtNodeOutFileName = Entry(dlgNodeFrame)
      dlgNodeFrame.txtNodeOutFileName.grid(row=5, column= 3, sticky = W)

      Button(dlgNodeFrame,
             text= "Browse...",
             command = self.cmdBrowse1
             ).grid(row=5,column=4,sticky=W)
      
      #Add 1 empty row starting from row 6
      emptyRows(1,6,dlgNodeFrame)
      
      #Label(dlgNodeFrame,text = "").grid(row=4, column=0,columnspan = 2, sticky =W)
      
      # Ask for Object type
      k = -1 
      Label(dlgNodeFrame,
            text = " Input3:Number of nearest Neighbors for k-Nearest Neighbors :"
            ).grid(row=7, column=0, sticky =W)
      self.NumofNeighbors=dlgNodeFrame.NumofNeighbors = Entry(dlgNodeFrame)
      dlgNodeFrame.NumofNeighbors.grid(row=7, column= 3, columnspan=1,sticky = W)

      #DEMO
      self.txtNodeInFileName.delete(0,END)  
      self.txtNodeInFileName.insert(0,"Neurons.xml")

      self.txtNodeInFileName2.delete(0,END)  
      self.txtNodeInFileName2.insert(0,"Microglia.xml")

      self.txtNodeOutFileName.delete(0,END)  
      self.txtNodeOutFileName.insert(0,"links-Neurons-Microglia-5NN.xml")

      self.NumofNeighbors.delete(0,END);
      self.NumofNeighbors.insert(0,"5")
      
      #DEMO

      #Add 4 empty rows starting from row 5
      emptyRows(4,5,dlgNodeFrame)
         
      self.buttonCreate= Button(dlgNodeFrame,
             text= "  Create Links...  ",
             command = self.cmdCreateLinks
             ).grid(row=11,column=1,sticky=W)

      self.buttonCancel=Button(dlgNodeFrame,
             text= "  Cancel  ",
             command = self.cmdQuit
             ).grid(row=11,column=3,sticky=W)

      dlgNodeFrame.pack(side=TOP,fill=X)

                   
   def cmdBrowse1(self):
      fileName1=asksaveasfilename()
      if fileName1:
            self.txtNodeOutFileName.delete(0,END)  
            self.txtNodeOutFileName.insert(0,fileName1)
   
   def cmdBrowse2(self):
      fileName2 = askopenfilename()
      if fileName2:
            self.txtNodeInFileName.delete(0,END)
            self.txtNodeInFileName.insert(0,fileName2)
      
   def cmdBrowse3(self):
      fileName3 = askopenfilename()
      if fileName3:
            self.txtNodeInFileName2.delete(0,END)
            self.txtNodeInFileName2.insert(0,fileName3)
      
   def cmdCreateLinks(self):

      # k for kNN
      k = self.NumofNeighbors.get()

      fileName1=self.txtNodeOutFileName.get()
      fileName2=self.txtNodeInFileName.get()
      fileName3=self.txtNodeInFileName2.get()
      
      if (k != -1) and (fileName1 != "") and (fileName2 != "") and (fileName3 != "") :
            #java -Xmx1024m net.sf.saxon.Query -o:"+sys.argv[1]+" kneighbor.xquery nodes1=" + sys.argv[2]+ " nodes2="+sys.argv[3]+" n=" + sys.argv[4])
            javaCommand = "java -Xmx1024m net.sf.saxon.Query -o:"+fileName1+" "
            javaCommand = javaCommand + "kneighbor.xquery nodes1=" + fileName2+ " "
            javaCommand = javaCommand + "nodes2=" + fileName3 + " "
            javaCommand = javaCommand + "n=" + k
            os.system(javaCommand)
            showinfo('Warning!', 'Links are created. Output is written to: '+fileName1)
      else:
            showerror('Error!', 'Inputs are Incomplete')
                         
   def cmdQuit(self):
      self.master.focus_set()
      self.destroy()

   def defineLinks(self):
     crtNodeObj=dlgCreateLinks(root)

##########################################################################################
##########################################################################################
# Usage: python defineNetwork output nodes1 nodes2 links

class dlgCreateNetwork(Toplevel):
   def __init__(self, master):

      master.title('Create Network')
      Toplevel.__init__(self, master)
      #This window and its parent are required to be related
      self.transient(master)
      self.master=master
      dlgNodeFrame = Frame(self,width=50, height=60, bd=1)
      dlgNodeFrame.pack()
      #call grabset to make dialog modal= Only 1 dialog can appear on the screen
      #We don't want to confuse users
      self.grab_set()

      self.geometry("+%d+%d" % (master.winfo_rootx()+350,
                        master.winfo_rooty()+120))

      Label(dlgNodeFrame,text = "").grid(row=0, column=0,columnspan = 2, sticky =W)

      fileName1=""
      fileName2=""
      fileName3=""
      fileName4=""
 
      Label(dlgNodeFrame,
            text = "Input1: Filename for nodes of specific object type :"
            ).grid(row=1, column=0,columnspan = 2, sticky =W)
      self.txtNodeInFileName=dlgNodeFrame.txtNodeInFileName = Entry(dlgNodeFrame)
      dlgNodeFrame.txtNodeInFileName.grid(row=1, column= 3, sticky = W)

      Button(dlgNodeFrame,
             text= "Browse...",
             command = self.cmdBrowse2
             ).grid(row=1,column=4,sticky=W)

      #Add 1 empty row starting from row 2
      emptyRows(1,2,dlgNodeFrame)
      
      Label(dlgNodeFrame,
            text = " Input2: Filename for nodes of specific object type :"
            ).grid(row=3, column=0,columnspan = 2, sticky =W)
      self.txtNodeInFileName2=dlgNodeFrame.txtNodeInFileName2 = Entry(dlgNodeFrame)
      dlgNodeFrame.txtNodeInFileName2.grid(row=3, column= 3, sticky = W)

      Button(dlgNodeFrame,
             text= "Browse...",
             command = self.cmdBrowse3
             ).grid(row=3,column=4,sticky=W)

      #Add 1 empty row starting from row 4
      emptyRows(1,4,dlgNodeFrame)
      
      #Label(dlgNodeFrame,text = "").grid(row=4, column=0,columnspan = 2, sticky =W)
      
      # Ask for Links
      Label(dlgNodeFrame,
            text = " Input3: Filename for links defined between two object types: "
            ).grid(row=5, column=0, sticky =W)
      self.txtNodeInFileName3=dlgNodeFrame.txtNodeInFileName3 = Entry(dlgNodeFrame)
      dlgNodeFrame.txtNodeInFileName3.grid(row=5, column= 3, sticky = W)

      Button(dlgNodeFrame,
             text= "Browse...",
             command = self.cmdBrowse4
             ).grid(row=5,column=4,sticky=W)

      #Add 1 empty row starting from row 2
      emptyRows(1,6,dlgNodeFrame)
      
      Label(dlgNodeFrame,
            text = "Output: Enter the Network File name to be created :"
            ).grid(row=7, column=0,columnspan = 2, sticky =W)
      self.txtNodeOutFileName=dlgNodeFrame.txtNodeOutFileName = Entry(dlgNodeFrame)
      dlgNodeFrame.txtNodeOutFileName.grid(row=7, column= 3, sticky = W)

      Button(dlgNodeFrame,
             text= "Browse...",
             command = self.cmdBrowse1
             ).grid(row=7,column=4,sticky=W)

      #DEMO
      self.txtNodeInFileName.delete(0,END)  
      self.txtNodeInFileName.insert(0,"Neurons.xml")

      self.txtNodeInFileName2.delete(0,END)  
      self.txtNodeInFileName2.insert(0,"Microglia.xml")

      self.txtNodeInFileName3.delete(0,END)  
      self.txtNodeInFileName3.insert(0,"links-Neurons-Microglia-5NN.xml")

      self.txtNodeOutFileName.delete(0,END)  
      self.txtNodeOutFileName.insert(0,"network-Neurons-Microglia-5NN.xml")
      
      #DEMO

      #Add 4 empty rows starting from row 5
      emptyRows(4,5,dlgNodeFrame)
         
      self.buttonCreate= Button(dlgNodeFrame,
             text= "  Create Network...  ",
             command = self.cmdCreateNetwork
             ).grid(row=11,column=1,sticky=W)

      self.buttonCancel=Button(dlgNodeFrame,
             text= "  Cancel  ",
             command = self.cmdQuit
             ).grid(row=11,column=3,sticky=W)

      dlgNodeFrame.pack(side=TOP,fill=X)

   def cmdBrowse1(self):
      fileName1=asksaveasfilename()
      if fileName1:
            self.txtNodeOutFileName.delete(0,END)
            self.txtNodeOutFileName.insert(0,fileName1)
   
   def cmdBrowse2(self):
      fileName2 = askopenfilename()
      if fileName2:
            self.txtNodeInFileName.delete(0,END)
            self.txtNodeInFileName.insert(0,fileName2)
      
   def cmdBrowse3(self):
      fileName3 = askopenfilename()
      if fileName3:
            self.txtNodeInFileName2.delete(0,END)
            self.txtNodeInFileName2.insert(0,fileName3)

   def cmdBrowse4(self):
      fileName4 = askopenfilename()
      if fileName4:
            self.txtNodeInFileName3.delete(0,END)
            self.txtNodeInFileName3.insert(0,fileName4)        
      
   def cmdCreateNetwork(self):

      fileName1=self.txtNodeOutFileName.get()
      fileName2=self.txtNodeInFileName.get()
      fileName3=self.txtNodeInFileName2.get()
      fileName4=self.txtNodeInFileName3.get()
      
      if (fileName1 != "") and (fileName2 != "") and (fileName3 != "") and (fileName4 != "") :
            #java -Xmx1024m net.sf.saxon.Query -o:"+sys.argv[1]+" kneighbor.xquery nodes1=" + sys.argv[2]+ " nodes2="+sys.argv[3]+" n=" + sys.argv[4])
            #java -Xmx1024m net.sf.saxon.Query -o:"+sys.argv[1]+" defineNetwork.xquery nodes1=" + sys.argv[2]+" nodes2="+sys.argv[3]+" links=" + sys.argv[4])
            
            javaCommand = "java -Xmx1024m net.sf.saxon.Query -o:"+fileName1+" "
            javaCommand = javaCommand + "defineNetwork.xquery nodes1=" + fileName2+ " "
            javaCommand = javaCommand + "nodes2=" + fileName3 + " "
            javaCommand = javaCommand + "links=" + fileName4
            os.system(javaCommand)
            showinfo('Warning!', 'Network has been created. Output is written to: '+fileName1)
      else:
            showerror('Error!', 'Inputs are Incomplete')
            
   def cmdQuit(self):
      self.master.focus_set()
      self.destroy()

##########################################################################################
##########################################################################################

class dlgDisplayNetwork(Toplevel):
   def __init__(self, master):

      master.title('Display Network')
      Toplevel.__init__(self, master)
      #This window and its parent are required to be related
      self.transient(master)
      self.master=master
      dlgNodeFrame = Frame(self,width=50, height=60, bd=1)
      dlgNodeFrame.pack()
      #call grabset to make dialog modal= Only 1 dialog can appear on the screen
      #We don't want to confuse users
      self.grab_set()

      self.geometry("+%d+%d" % (master.winfo_rootx()+350,
                        master.winfo_rooty()+120))

      Label(dlgNodeFrame,text = "").grid(row=0, column=0,columnspan = 2, sticky =W)

      fileName1=""
      fileName2=""
      fileName3=""
      fileName4=""

      Label(dlgNodeFrame,
            text = "Input: Filename that defines the network :"
            ).grid(row=1, column=0,columnspan = 2, sticky =W)
      self.txtNetworkInFileName=dlgNodeFrame.txtNetworkInFileName = Entry(dlgNodeFrame)
      dlgNodeFrame.txtNetworkInFileName.grid(row=1, column= 3, sticky = W)

      Button(dlgNodeFrame,
             text= "Browse...",
             command = self.cmdBrowse
             ).grid(row=1,column=4,sticky=W)

      #Add 2 empty rows starting from row 2
      #emptyRows(2,2,dlgNodeFrame)

      # Ask for Object type
      Label(dlgNodeFrame,
            text = " Select Network Type :"
            ).grid(row=2, column=0, sticky =W)
      self.txtObjectType=dlgNodeFrame.txtObjectType = StringVar()
      
      #graphTypes = ["Standard Network", "Minimum Spanning Tree"]
      #self.txtObjectType.set("Standard Network")
      column = 2
      row=3

      # RadioButtons for selecting raw Network or MST
      self.objectTypes=Radiobutton(dlgNodeFrame,
                     text="Standard Network",                     
                     variable=self.txtObjectType,
                     value="stdNet",
                     command=self.setNestedRadioBtns
                  ).grid(row=row,column=0,sticky=W)

      #The standard Network is selected by default   
      self.txtObjectType.set("stdNet")

      #emptyRows(2,4,dlgNodeFrame)

      self.stdNetType= dlgNodeFrame.stdNetParam=StringVar()
      #self.stdNetParam.set(2)
      #dlgNodeFrame.stdNetParam.set(2)
      self.stdNetParam=StringVar()
      self.stdNetParam.set(1)
      # Nested RadioButtons for selecting raw Network Type
      self.stdNetType=Radiobutton(dlgNodeFrame,
                     text="Pyramidial",
                     variable=self.stdNetParam,
                     value="pyramidialNet",
                     command=self.setStandard
                     ).grid(row=4,column=1,sticky=W)

      #For pyramidial, a cutoff point has to be entered!
      Label(dlgNodeFrame,
            text = "cutoff point :"
            ).grid(row=4, column=2, sticky =W)
      self.cutOff=Entry(dlgNodeFrame)
      self.cutOff.grid(row=4, column= 3, sticky = W)
      
      #dlgNodeFrame.txtNetworkInFileName.grid(row=1, column= 3, sticky = W)
      
      #dlgNodeFrame.txtNetworkInFileName = Entry(dlgNodeFrame)
      
      
      self.stdNetType=Radiobutton(dlgNodeFrame,
                                  text="Distance from cell-2",
                                  variable=self.stdNetParam,
                                  value="distanceNet",
                                  command=self.setStandard
                                  ).grid(row=5,column=1,sticky=W)

      #For every cell-1, draw links to cell2's within (<) a given distance
      Label(dlgNodeFrame,
            text = "distance :"
            ).grid(row=5, column=2, sticky =W)
      self.distance=Entry(dlgNodeFrame)
      self.distance.grid(row=5, column= 3, sticky = W)
      
      
      self.objectTypes=Radiobutton(dlgNodeFrame,
                                  text="Minimum Spanning Tree",
                                  variable=self.txtObjectType,
                                  value="mstNet",
                                  command=self.setNestedRadioBtns
                                  ).grid(row=6,column=0,sticky=W)

      #DEMO
      self.txtNetworkInFileName.delete(0,END)  
      self.txtNetworkInFileName.insert(0,"network-Neurons-Microglia-5NN.xml")

      self.stdNetParam.set("pyramidialNet")

      self.cutOff.delete(0,END)  
      self.cutOff.insert(0,"50")
      #DEMO

      self.buttonCreate= Button(dlgNodeFrame,
             text= "  Display Network...  ",
             command = self.cmdDisplayNetwork
             ).grid(row=7,column=1,sticky=W)

      self.buttonCancel=Button(dlgNodeFrame,
             text= "  Cancel  ",
             command = self.cmdQuit
             ).grid(row=7,column=3,sticky=W)

      dlgNodeFrame.pack(side=TOP,fill=X)

   def setStandard(self):
         self.txtObjectType.set("stdNet")
         if (self.stdNetParam.get()=="pyramidialNet"):
               self.distance.delete(0,END)
         if (self.stdNetParam.get()=="distanceNet"):
               self.cutOff.delete(0,END)
                  
   def setNestedRadioBtns(self):
         # If mst is selected we need to clear nested radio buttons
         # and associated texboxes that are for only standard networks
         if (self.txtObjectType.get() == "mstNet") :
               self.stdNetParam.set(1)               
               self.cutOff.delete(0,END)
               self.distance.delete(0,END)
               
   def cmdDisplayNetwork(self):

      fileName2=self.txtNetworkInFileName.get()
      stdMst=IntVar() #whether standard or an mst network is asked for!
      #osCommand = StringVar()
      osCommand = " "
      distance= "-1"  # No distance by default
      cutoff  = "-1"  # No cutoff by default
      stdMst  = "-1"  # 0->standard  1->Minimum Spanning Tree
                      # 2->standard and Pyramidial 3->standard and distance    
      stdNet = self.txtObjectType.get()

      if (fileName2 != "") :
            if (stdNet == "mstNet") :
                  stdMst = "1" #Display the mst of the given network
                  osCommand = "bionet "+fileName2+" "+distance+" "+stdMst+" "+cutoff
            else:   #The user wants to display the standard network
                  if (self.stdNetParam.get()=="pyramidialNet"):                        
                        if self.cutOff.get():
                              stdMst = "2"
                              cutoff = self.cutOff.get()
                              osCommand = "tissuenets "+fileName2+" "+distance+" "+stdMst+" "+cutoff
                        else: 
                             showerror('Error!', 'Enter a CutOff point first')       
                  else:
                        if (self.stdNetParam.get()=="distanceNet"):                              
                              if self.distance.get():
                                    stdMst = "3"
                                    distance = self.distance.get()
                                    osCommand = "tissuenets "+fileName2+" "+distance+" "+stdMst+" "+cutoff
                              else: 
                                   showerror('Error!', 'Enter a distance first')     

                        else: #No parameter selections for standard network selection
                              stdMst="0"
                              osCommand = "tissuenets "+fileName2+" "+distance+" "+stdMst+" "+cutoff
      else:
            showerror('Error!', 'Enter a file Name')

      if (stdMst != "-1"):
            os.system(osCommand)

   def cmdBrowse(self):
      fileName2 = askopenfilename()
      if fileName2:
            self.txtNetworkInFileName.delete(0,END)
            self.txtNetworkInFileName.insert(0,fileName2)
        
                         
   def cmdQuit(self):
      self.master.focus_set()
      self.destroy()

   def defineLinks(self):
     crtNodeObj=dlgCreateLinks(root)

   def getOptions(self):
     row=row+1

   #self.objectTypes.set(1)
     
     #crtNodeObj=dlgCreateLinks(root)     

##########################################################################################
##########################################################################################

class dlgCreateHistograms(Toplevel):
   def __init__(self, master):

      master.title('Histogram Creation')
      Toplevel.__init__(self, master)
      #This window and its parent are required to be related
      self.transient(master)
      self.master=master
      dlgNodeFrame = Frame(self,width=50, height=60, bd=1)
      dlgNodeFrame.pack()
      #call grabset to make dialog modal= Only 1 dialog can appear on the screen
      #We don't want to confuse users
      self.grab_set()

      self.geometry("+%d+%d" % (master.winfo_rootx()+350,
                        master.winfo_rooty()+120))

      Label(dlgNodeFrame,text = "").grid(row=0, column=0,columnspan = 2, sticky =W)

      fileName1=""
      fileName2=""
      fileName3=""
      fileName4=""

      #Add 2 empty rows starting from row 2
      #emptyRows(2,2,dlgNodeFrame)

      # Ask for Object type
      Label(dlgNodeFrame,
            text = " Select Feature :"
            ).grid(row=2, column=0, sticky =W)
      self.txtObjectType=dlgNodeFrame.txtObjectType = StringVar()
            
      self.features = ["Volume",
                  "Intensity",
                  "Intensity_variation",
                  "Skew_of_intensity_distribution",
                  "Energy_of_intensity_distribution",
                  "Texture",
                  "Convexity",
                  "Surface_area",
                  "Shape",
                  "Surface_curvature",
                  "Surface_gradient",
                  "Boundary_sharing"
                  "Eccentricity",
                  "Radius_variation",
                  "Orientation",
                  "Interior_gradient",
                  "Intensity_ratio",
                  "Depth",
                  "Number_of_poles_on_surface",
                  "Iba1_Signal",
                  "GFAP_Signal",
                  "Nissl_sig",
                  "Laminin_Signal",
                  "Distance_to_blood_vessel",
                  "Distance_to_convergence_point"]
    
      column = 0
      row=3
      row1=3

      self.feature=StringVar()  
      for feature in self.features:
            # RadioButtons for selecting Features
            self.objectTypes=Radiobutton(dlgNodeFrame,
                     text=feature,                     
                     variable=self.feature,
                     value=feature,
                     command=self.setDefaultFileName
                  ).grid(row=row,column=column,sticky=W)
            row+=1
      self.feature.set("Volume")

      # Ask for Object type
      Label(dlgNodeFrame,
            text = " Select Cell Type :"
            ).grid(row=row1-1, column=column+2, sticky =W)
      self.txtObjectType=dlgNodeFrame.txtObjectType = StringVar()

      allObjects = ["Astrocytes", "Microglia", "Neurons","Endothelials"]
      
      self.cellType=StringVar()
      for fobject in allObjects:
         self.objectTypes=Radiobutton(dlgNodeFrame,
                     text=fobject,
                     variable=self.cellType,
                     value=fobject
                     ).grid(row=row1,column=column+2,sticky=W)
         row1+=1
      self.cellType.set("Astrocytes")
         
      self.stdNetType= dlgNodeFrame.stdNetParam=StringVar()
      self.stdNetParam=StringVar()
      self.stdNetParam.set(1)

      #Add 4 empty rows starting from row 5
      emptyRows(4,row,dlgNodeFrame)

      
      Label(dlgNodeFrame,
            text = " Input: Filename for Histogram Data :"
            ).grid(row=row+1, column=column,columnspan = 2, sticky =W)
      self.txtNetworkInFileName=dlgNodeFrame.txtNetworkInFileName = Entry(dlgNodeFrame)
      dlgNodeFrame.txtNetworkInFileName.grid(row=row+1, column= column+2, sticky = W)

      Button(dlgNodeFrame,
             text= "Browse...",
             command = self.cmdBrowse
             ).grid(row=row+1,column=column+3,sticky=W)

      row+=1

      self.numofBins=IntVar()  
      Label(dlgNodeFrame,
            text = " Input: Number of bins ( 2 <= bins <= 10 ) :"
            ).grid(row=row+1, column=column,columnspan = 2, sticky =W)
      self.numofBins= dlgNodeFrame.numofBins=Entry(dlgNodeFrame)
      dlgNodeFrame.numofBins.grid(row=row+1, column= column+2, sticky = W)

      row+=1

      Label(dlgNodeFrame,
            text = "Output: File name to store histogram data :"
            ).grid(row=row+1, column=column,columnspan = 2, sticky =W)
      self.txtNodeOutFileName=dlgNodeFrame.txtNodeOutFileName = Entry(dlgNodeFrame)
      dlgNodeFrame.txtNodeOutFileName.grid(row=row+1, column= column+2, sticky = W)

      #DEMO
      self.txtNetworkInFileName.delete(0,END)  
      self.txtNetworkInFileName.insert(0,"SupplementD.xml")
      
      self.numofBins.delete(0,END)  
      self.numofBins.insert(0,"6")

      self.txtNodeOutFileName.delete(0,END)  
      self.txtNodeOutFileName.insert(0,"Microglia_Intensity.xml")

      self.feature.set("Intensity")
      self.cellType.set("Microglia")
      
      #DEMO
      
      Button(dlgNodeFrame,
             text= "Browse...",
             command = self.cmdBrowse1
             ).grid(row=row+1,column=column+3,sticky=W)
      
      #Add 1 empty row starting from row
      emptyRows(1,row,dlgNodeFrame)
      
      self.buttonCreate= Button(dlgNodeFrame,
             text= "  Create Histogram Data...",
             command = self.cmdCreateHistogramData
             ).grid(row=row+2,column=column+1,sticky=W)

      self.buttonCancel=Button(dlgNodeFrame,
             text= "  Cancel  ",
             command = self.cmdQuit
             ).grid(row=row+2,column=column+2,sticky=W)

      dlgNodeFrame.pack(side=TOP,fill=X)
               
   def cmdCreateHistogramData(self):

      fileName1=self.txtNetworkInFileName.get()                
      fileName2=self.txtNodeOutFileName.get()
      numofBins = self.numofBins.get()

      if (fileName1 != "") and (fileName2 != "") and (numofBins !="")  :
            if (int(numofBins) < 2) or (int(numofBins) > 10):
                  showerror('Error!', 'Number of bins can be a number from only 2 to 10')
            else:                  
                  feature = self.feature.get()
                  cellType = self.cellType.get()
                  javaCommand = "java -Xmx1024m net.sf.saxon.Query -o:"+fileName2+" "
                  javaCommand = javaCommand + "defineHistogram.xquery allnodes=" + fileName1+ " "
                  javaCommand = javaCommand + "objecttype=" + cellType + " "
                  javaCommand = javaCommand + "featurename=" + feature
                  os.system(javaCommand)            
                  osCommand = "histogram " + fileName2 + " " + numofBins + " " + fileName2
                  os.system(osCommand) 
      else:   #The user wants to display the standard network
            showerror('Error!', 'Empty Fields! Enter inputs...')

   def setDefaultFileName(self):
      self.txtNodeOutFileName.delete(0,END)
      self.txtNodeOutFileName.insert(0, self.cellType.get()+"_"+self.feature.get()+".xml")

   def cmdBrowse(self):
      fileName1 = askopenfilename()
      if fileName1:
            self.txtNetworkInFileName.delete(0,END)
            self.txtNetworkInFileName.insert(0,fileName1)

   def cmdBrowse1(self):
      fileName2=asksaveasfilename()
      if fileName2:
            self.txtNodeOutFileName.delete(0,END)
            self.txtNodeOutFileName.insert(0,"")
            self.txtNodeOutFileName.insert(0,fileName2)
          
   def cmdQuit(self):
      self.master.focus_set()
      self.destroy()

   def getOptions(self):
     row=row+1

   #self.objectTypes.set(1)
     
     #crtNodeObj=dlgCreateLinks(root)
     
##########################################################################################
##########################################################################################
      
class MainMenu:

   def notdone(self):
         doNothing=0     
     
   def defineNode(self):
     crtNodeObj=dlgCreateNode(root)

   def defineLinks(self):
     crtNodeObj=dlgCreateLinks(root)
     
   def defineNetwork(self):
     crtNodeObj=dlgCreateNetwork(root)

   def displayNetwork(self):
     crtNodeObj=dlgDisplayNetwork(root)
     
   def createHistograms(self):
     crtNodeObj=dlgCreateHistograms(root)

   def defineAll(self):
     #flag =1    
     crtNodeObj=dlgCreateNode(root)
     crtNodeObj=dlgCreateLinks(root)
     crtNodeObj=dlgCreateNetwork(root)     
     
   def __init__(self,win):
      top = Menu(win)
      win.config(menu=top)
    
      menuCreateNetwork = Menu(top,tearoff=0)
      menuCreateNetwork.add_command(label='Define Node ...', command=self.defineNode, underline=0)
      menuCreateNetwork.add_command(label='Define Links ...', command=self.defineLinks, underline=0)
      menuCreateNetwork.add_command(label='Create Network ...', command=self.defineNetwork, underline=0)
      #menuCreateNetwork.add_command(label='Display Network ...', command=self.displayNetwork, underline=0)
      #menuCreateNetwork.add_command(label='Create Everything ...', command=self.defineAll, underline=0)

      menuDisplayNetwork = Menu(top, tearoff=0)
      menuDisplayNetwork.add_command(label='Tissue Network with k-NN ...', command=self.displayNetwork, underline=0)
      #menuDisplayNetwork.add_command(label='Bio Network with k-NN and "less than some distance" ...', command=self.notdone, underline=0)
      #menuDisplayNetwork.add_command(label='Bio Network with "less than some distance" ...', command=self.notdone, underline=0)
      #menuDisplayNetwork.add_command(label='Minimum Spanning Tree of a Bio Network ...', command=self.notdone, underline=0)

      #menuDisplayNetwork.add_command(label='Bio Network with k-NN (NO LABELS) ...', command=self.displayNetwork, underline=0)
      #menuDisplayNetwork.add_command(label='Bio Network with k-NN and "less than some distance" (NO LABELS) ...', command=self.notdone, underline=0)
      #menuDisplayNetwork.add_command(label='Bio Network with "less than some distance" (NO LABELS) ...', command=self.notdone, underline=0)
      #menuDisplayNetwork.add_command(label='Minimum Spanning Tree of a Bio Network (NO LABELS) ...', command=self.notdone, underline=0)
      
      menuQueries = Menu(top, tearoff=0)
      menuQueries.add_command(label='Define Histograms ...', command=self.createHistograms, underline=0)
      menuQueries.add_command(label='Enter an XQUERY Query...(NOT WORKING YET!)', command=self.notdone, underline=0)
      menuQueries.add_command(label='Built-in Query-1 ...(NOT WORKING YET!)', command=self.notdone, underline=0)
      menuQueries.add_command(label='Built-in Query-1 ...(NOT WORKING YET!)', command=self.notdone, underline=0)

      top.add_cascade(label='Create Network', menu=menuCreateNetwork, underline=0)
      top.add_cascade(label='Display Network', menu=menuDisplayNetwork, underline=0)
      top.add_cascade(label='Queries and Histograms', menu=menuQueries, underline=0)
      
      
root =Tk()
root.title('F A R S I G H T')
root.geometry("700x700")

farsightMenu = MainMenu(root)
msg = Label(root, text='F A R S I G H T')
msg.pack(expand=YES, fill=BOTH)
msg.config(relief=SUNKEN, width=40, height=7, bg='light green')
root.mainloop()
