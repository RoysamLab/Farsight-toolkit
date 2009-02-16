
#ifndef Connected_h
#define Connected_h

#define MaxNumOfSomas 50
#define MaxNumOfSegments 100
#define MaxNumOfColors 512

class EquivalentColors 
{
public:
	EquivalentColors();

	inline void SetLabel(int color) 
	{
		label = color;
	}
	
	void Merge(int color, int recursiveFlag);

	int label;
	int NumOfColors;
	int myColors[MaxNumOfColors];
	int validColor;
};

class CTree;

// a soma has a center, a color, and a list of other somas and segments it 
// is directly or indirectly connected to
class CSoma
{
public:
	
	CSoma();

	int AddSegment(int id);
	void AddIntersectionPoint(int id);
	
	void AddConnectedSomaID(int id);
	
	void Print(FILE *);
	
	void PrintTrees(FILE *);

	// Draw the trees hanging from all somas in the projection
	// images 
	void DrawTrees();

	void WriteID(CImage &, unsigned char color = LetterColor);
	void WriteIDXZ(CImage &, unsigned char color = LetterColor);
	void WriteIDYZ(CImage &, unsigned char color = LetterColor);
	void WriteLabel(CImage &, int, int, int, int, unsigned char color = LetterColor);

	void SetLabel(int i);
	
	void ConstructTrees();
	
	// write the structures in Neurolucida format
	void PrintNeurolucidaFormat(ofstream &outFile);


	// data
	CTree *m_aTrees;

	CPoint m_Center;
	int m_iLabel;
	int m_iVolume;
	int m_aiConnectedSomaIDs[MaxNumOfSomas];
	int m_aiIntersectionPoints[MaxNumOfSegments];
	int m_aiSegments[MaxNumOfSegments];
	int m_iNumOfIntersectionPoints;
	int m_iNumOfSegments;
	int m_iNumOfSomas;
	int m_iInitialized;
	int m_iHitsImageBoundary;  // does the soma hit the image boundary?
	int m_iSumOfAllSegments;
	char m_achID[3];
};	


// Collection class of somas
// the class contains some extra functionality as well
class CSomas
{
public:
	// CTOR and DTOR
	CSomas(int = 0);
	~CSomas();

	// Write the soma IDs to an image
	void WriteIDs(CImage &, unsigned char color = LetterColor);
	void WriteIDsXZ(CImage &, unsigned char color = LetterColor);
	void WriteIDsYZ(CImage &, unsigned char color = LetterColor);


	// for each of the somas, construct its trees
	void ConstructTrees();

	void Print(char *fName);
	void PrintTrees(char *fName = NULL);

	// inform the two given somas that they are connected
	void ConnectSomas(int id1, int id2);

	// inform the some with the given id that it hits the image boundary
	void SomaHitsImageBoundary(int somaID);

	// add the segment id to the soma's list of ids
	int AddSegment(int somaID, int segmentID);

	// Draw the trees hanging from all somas in the projection
	// images 
	void DrawTrees();

	// write the structures in Neurolucida format
	void PrintNeurolucidaFormat(char *fName);

	// data members

	// an array of somas
	CSoma *m_aData;

	// the number of somas in my collection
	int    m_iNumOfSomas;
};

#endif
