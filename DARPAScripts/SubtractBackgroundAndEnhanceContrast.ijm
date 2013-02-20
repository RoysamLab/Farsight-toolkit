macro "Batch Preprocess" 
{
	convert();
} 

function convert() 
{
	requires("1.33s");
	dir1 = getDirectory("Choose Source Directory "); 
	dir2 = getDirectory("Choose Destination Directory "); 
	list = getFileList(dir1);
	//setBatchMode(true);
	for (i=0; i<list.length; i++) 
	{ 
		showProgress(i+1, list.length); 
		open(dir1+list[i]);
		run("Enhance Contrast", "saturated=0.00 use");
		run("Median...", "radius=1");
		run("Subtract Background...", "rolling=50 stack");
		run("Enhance Contrast", "saturated=0.00 use");
		//run("16-bit");
		run("8-bit");
		
		saveAs("tiff", dir2+list[i]);
		close();
	}
}

