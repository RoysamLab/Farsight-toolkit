/*=========================================================================
Copyright 2009 Rensselaer Polytechnic Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License. 
=========================================================================*/

// all colors below this color are reserved for the tracking process

#define HighestReservedColor 5

enum TopEndMiddleNeither {  Neither = 0, OnTop = 1, OnEnd = 2, Middle = 3};

// a node in a tree can have any of these types
enum TreeNodeType { ROOT = 1, NORMAL = 2, LEAF = 3 };

enum BranchingOrIntersection { Branching = 0, Intersection = 1};


// these color values are used to indicate the type for all pixels
// > 250; 
enum PixelColors { SomaColor = 1, CenterlineColor = 2, ValidSeedColor = 2,
								   InvalidSeedColor = 2, VesselColor = 3, 
									InvalidSeedGrayLevel = 5, FillingColor = 2,
									BoundaryColor = 4, IDColor = 0,						
									HighSomaColor = 255, StartingSomaColor = 230,
									StartingIntersectionColor = 150,
									LetterColor = 4, IntersectionPointColor = 3};

