/////////////////////////////////////////////////////////////////////
// File: Cvector.h
//
// Contains the declaration for the class vector.h. It represents a 3D vector 
//

#ifndef CVector_h
#define CVector_h

class CVector
{
public:
	CVector(int slices, int col, int row, float , float , int len);
	~CVector();

	// draw the template in the given image
	void DisplayVector(C3DImage &aImage, int color = 255);

	// set the base of the template to that position
	void Position(int x, int y, int z);
	
	// data	
	// the directions of the vector (in degrees)
	// the direction of avector is defined by two angles,
	// the one with the +ve x-axis (m_fHDir), and the one with
	// the +ve z-axiz (m_fVDir);
	float m_fHDir;
	float m_fVDir;
	// the length of the vector
	int m_iLength;
	// the number of rows and columns, in the corresponding image
	int m_iRows;
	int m_iCols;
	// the vector.
	// notice that each vector cells represents an abosolute shift from
	// the vector's base to the cell.
	int *m_piData;
	int m_iBase;
	// each vector also contains the (x, y, z) coordinates of shifts.
	CPoint *m_pIndices;
};

#endif
