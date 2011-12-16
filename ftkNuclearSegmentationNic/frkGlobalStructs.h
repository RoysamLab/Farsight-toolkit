// ############################################################################################################################################################################
#ifndef _frkGlobalStructs_h_
#define _frkGlobalStructs_h_
// ############################################################################################################################################################################


namespace ftk{
	/**
	* namespace corresponding to the segmentation algorithm 
	* implemented by nicolas
	*/
	namespace nucSecNic{

		struct ConnComp
		{
			int x1;
			int x2;
			int y1;
			int y2;
			int z1;
			int z2;
		};

		class Seed
		{
		public:
			Seed(){
				xVal = 0;
				yVal = 0;
				zVal = 0;
				idVal = 0;
				ccVal = 0;
			}
			Seed( int x, int y, int z, int ID, int cc) {
				xVal = x;
				yVal = y;
				zVal = z;
				idVal = ID;
				ccVal = cc;
			}
			
			void setX(int x ){ xVal = x; }
			void setY(int y ){ yVal = y; }
			void setZ(int z ){ zVal = z; }
			void setID(int i){ idVal = i;}
			void setCC(int cc) { ccVal = cc;}

			int x() const { return xVal; }
			int y() const { return yVal; }
			int z() const { return zVal; }
			int ID() const {return idVal;}
			int CC() const {return ccVal;}

		private:
			int xVal;
			int yVal;
			int zVal;
			int idVal;
			int ccVal;
		};
	};
};


#endif