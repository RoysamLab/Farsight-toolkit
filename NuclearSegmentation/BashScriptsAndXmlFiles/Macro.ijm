dir = getDirectory("Choose a Directory ");
print("asd");
listFiles(dir);
function listFiles(dir) {
	print("dsa");
	list = getFileList(dir);
	for (i=0; i<list.length; i++) {
		print(": " + dir + list[i]);
		open( dir + list[i] );
		vinayisaloser=File.nameWithoutExtension;
		run("8-bit");
		selectWindow(list[i]);
		saveAs("Tiff",dir + "M" + vinayisaloser + ".tif");
		close();
	}
}