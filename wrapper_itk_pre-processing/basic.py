# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# basic.py
# This module contains important function definitions that are used frequently in the various python scripts written.

# Author	: Adarsh K. Ramasubramonian
# Date 		: 21 May, 2009
#
# Last modified : 21 May, 2009
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

DEBUGMODE = False

# Debug print statement
def dprint(message):
	if DEBUGMODE:
		print message
 		
# Checks if the file name passed has the desired extension.                     
def CheckFileExtension(name, extension):			 
	lengthOfName = len(name)
	return (name[lengthOfName-4:] == extension)	 
