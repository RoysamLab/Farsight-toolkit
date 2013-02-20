# This script imports the Hyperlink module for the Python-based interface for ITK pre-processing algorithms. One of the reasons for specifying this separately is that some of the older versions of wxPyhon had the hyperlink module stored in wx.lib, and not in wx.lib.agw (in the latest version).

import sys
try:
    import wx.lib.agw.hyperlink as hyperlink
except:
    try:
        import wx.lib.hyperlink as hyperlink
    except:
        print "Hyperlink module not found. Please make sure that you have the wxPython version with the Hyperlink module."
