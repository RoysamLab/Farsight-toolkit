dir = getDirectory("Choose a Directory ");
print("asd");
listFiles(dir);
function listFiles(dir) {
	print("dsa");
	dir1 = "C:\\Lab\\ArunFiles\\Data\\Peixoto\\TSeries-09012011-A-009\\Unzipped\\Data\\asd\\";
	for (i=1; i<181; i++) {
		for (j=2; j<5; j++) {
			string1 = "";
			for (k=1; k<12; k++) {
				string = "TSeries-09012011-A-009_Cycle";
				if( i<10 ){
					string = string + "00"+ i;
				}
				else if( i<100){
					string = string + "0" + i;
				}
				else{
					string = string + i;
				}
				string = string + "_CurrentSettings_Ch" + j;
				string1 = string;
				string = string + "_0000";
				if( k<10 ){
					string = string + "0" + k;
				}
				else{
					string = string + k;
				}
				string = string + ".tif";
				print( dir + string );
				open( dir + string );
			}
			run("Images to Stack", "name=[] title=TSer use");
			string1 = string1+".tif";
			print( dir1 + string1 );
			saveAs("Tiff", dir1 + string1 );
			run("Close");
		}
	}
}