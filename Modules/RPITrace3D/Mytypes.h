// all colors below this color are reserved for the tracking process

#define HighestReservedColor 5

typedef enum TopEndMiddleNeither {  Neither = 0, OnTop = 1, OnEnd = 2, Middle = 3};

// a node in a tree can have any of these types
typedef enum TreeNodeType { ROOT = 1, NORMAL = 2, LEAF = 3 };

typedef enum BranchingOrIntersection { Branching = 0, Intersection = 1};


// these color values are used to indicate the type for all pixels
// > 250; 

typedef enum PixelColors { SomaColor = 1, CenterlineColor = 2, ValidSeedColor = 2,
								   InvalidSeedColor = 2, VesselColor = 3, 
									InvalidSeedGrayLevel = 5, FillingColor = 2,
									BoundaryColor = 4, IDColor = 0,						
									HighSomaColor = 255, StartingSomaColor = 230,
									StartingIntersectionColor = 150,
									LetterColor = 4, IntersectionPointColor = 3};


