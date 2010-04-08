
#include "TraceObject.cxx"
#include "ftkCommon/ftkProjectManager.h"
#include <string>
#include <set>


typedef std::pair<int,int> PairT;
typedef struct { int xMin; int xMax; int yMin; int yMax; } RegionType;
std::vector<RegionType> FindRegions( ftk::ProjectManager * project );
void MergeInRegion(TraceObject *tObj, RegionType region);
int IsParallelEndBegin(TraceLine * tLine1, TraceLine * tLine2);
int IsParallelEndEnd(TraceLine * tLine1, TraceLine * tLine2);
int IsParallelBeginBegin(TraceLine * tLine1, TraceLine * tLine2);
int IsParallelBeginEnd(TraceLine * tLine1, TraceLine * tLine2);
void MergePair(TraceObject *tObj, TraceLine * tLine1, TraceLine * tLine2);
void MergeBeginBegin(TraceObject *tObj, TraceLine * tLine1, TraceLine * tLine2);
void MergeEndBegin(TraceObject *tObj, TraceLine * tLine1, TraceLine * tLine2);
void MergeBeginEnd(TraceObject *tObj, TraceLine * tLine1, TraceLine * tLine2);
void MergeEndEnd(TraceObject *tObj, TraceLine * tLine1, TraceLine * tLine2);
double EuclideanDist(TraceBit b1, TraceBit b2, double z_penalty_factor=1.0);

const double maxDist = 10; //Worked with 8,10 here
const int minCount = 5;    //And 5 here
const double z_penalty_factor = 1.0;
const int imageX = 1024;
const int imageY = 1024;

int main(int argc, char *argv[])
{
	if(argc < 3)
	{
		std::cout << " Usage: " << argv[0] << " InputProject(xml) OutputFile(swc) <AlternateTranslations(txt)>" << std::endl;
	}

	char * projFilename = argv[1];
	char * outFilename = argv[2];
	char * altTrans = NULL;

	if(argc == 4)
	{
		altTrans = argv[3];
	}

	int len = strlen(projFilename);
	if( strcmp(projFilename+len-3, "xml" ) != 0)
	{
		std::cout << "Project file must be xml" << std::endl;
	}

	TraceObject allTraces;
	ftk::ProjectManager project(projFilename);

	if(altTrans)
	{
		project.ReplaceTranslations(altTrans);
		//project.writeProject("newTproject.xml");
	}

	for (unsigned int i = 0; i < project.size(); i++)
	{ 
		std::string fname = project.GetFileName(i);
		std::string type = project.GetFileType(i);
		if ( type.find("Trace") != std::string::npos )
		{
			allTraces.SetTraceOffset(project.GetTranslationX(i), project.GetTranslationY(i), project.GetTranslationZ(i));

			len = fname.size();
			if(fname.find(".swc") != std::string::npos )
			{
				allTraces.ReadFromSWCFile((char*)fname.c_str());
			}
			else if(fname.find(".xml") != std::string::npos )
			{
				allTraces.ReadFromRPIXMLFile((char*)fname.c_str());
			}
			else if (fname.find(".vtk") != std::string::npos )
			{
				allTraces.ReadFromVTKFile((char*)fname.c_str());
			}
		}
	}

	allTraces.GetVTKPolyData();	//This is necessary for markers (and merge)

	std::vector<RegionType> regions = FindRegions( &project );

	for(unsigned int i=0; i < regions.size(); ++i)
	{
		std::cout << "Region " << i+1 << " of " << regions.size() << "...\n";
		MergeInRegion(&allTraces, regions.at(i) );
		std::cout << "...done\n";
	}

	allTraces.WriteToSWCFile(outFilename);

	//std::cerr << "PRESS ENTER TO EXIT\n";
	//getchar();

	return EXIT_SUCCESS;
}

std::vector<RegionType> FindRegions( ftk::ProjectManager * project )
{
	std::vector<RegionType> regions;

	for (unsigned int i = 0; i < project->size(); i++)
	{ 
		std::string type = project->GetFileType(i);
		if ( type.find("Trace") != std::string::npos )
		{
			int tX1 = project->GetTranslationX(i);
			int tY1 = project->GetTranslationY(i);

			for (unsigned int j = i+1; j < project->size(); j++)
			{ 
				type = project->GetFileType(j);
				if ( type.find("Trace") != std::string::npos )
				{
					RegionType r;

					int tX2 = project->GetTranslationX(j);
					int tY2 = project->GetTranslationY(j);
					
					r.xMin = tX2 > tX1 ? tX2 : tX1;
					r.yMin = tY2 > tY1 ? tY2 : tY1;
					
					int xBeg = tX2 < tX1 ? tX2 : tX1;
					int yBeg = tY2 < tY1 ? tY2 : tY1;
					r.xMax = xBeg + imageX;
					r.yMax = yBeg + imageY;

					if(r.xMin < r.xMax && r.yMin < r.yMax)
						regions.push_back(r);
				}
			}
		}
	}

	return regions;
}

void MergeInRegion(TraceObject *tObj, RegionType region)
{
	int xMin = region.xMin;
	int xMax = region.xMax;
	int yMin = region.yMin;
	int yMax = region.yMax;

	std::vector<TraceLine*> lines = tObj->GetTraceLines();

	std::vector<PairT> closePairs;
	std::vector<int> count;

	//Iterate through lines and compare Trace Bits to see if the lines are running very close to each other within the region
	for(int i=0; i<(int)lines.size(); ++i)
	{
		int bestCount = 0;
		int bestJ = -1;

		for(int j=i+1; j<(int)lines.size(); ++j)
		{
			TraceLine * tLine1 = lines.at(i);
			TraceLine * tLine2 = lines.at(j); 

			//If Root start at Begin
			//If Leaf start at End
			//If free could be either
			TraceLine::TraceBitsType::iterator tb_it;
		
			TraceBit beg1 = *tLine1->GetTraceBitIteratorBegin();
			tb_it = tLine1->GetTraceBitIteratorEnd();
			tb_it--;
			TraceBit end1 = *tb_it;
			
			TraceBit beg2 = *tLine2->GetTraceBitIteratorBegin();
			tb_it = tLine2->GetTraceBitIteratorEnd();
			tb_it--;
			TraceBit end2 = *tb_it;

			//Find end of Line1:
			bool beg1_ok = false;
			bool end1_ok = false;
			bool beg2_ok = false;
			bool end2_ok = false;
			if(beg1.x >= xMin-1 && beg1.x <= xMax && beg1.y >= yMin-1 && beg1.y <= yMax)
				beg1_ok = true;
			if(end1.x >= xMin-1 && end1.x <= xMax && end1.y >= yMin-1 && end1.y <= yMax)
				end1_ok = true;
			if(beg2.x >= xMin-1 && beg2.x <= xMax && beg2.y >= yMin-1 && beg2.y <= yMax)
				beg2_ok = true;
			if(end2.x >= xMin-1 && end2.x <= xMax && end2.y >= yMin-1 && end2.y <= yMax)
				end2_ok = true;

			if(tLine1->isFree() && beg1_ok && end1_ok)
			{
				//Fragment so delete??
				continue;
			}
			if(tLine2->isFree() && beg2_ok && end2_ok)
			{
				//Fragment so delete??
				continue;
			}
			
			// 0 = can't do it
			// 1 = beginning
			// 2 = end
			int start1 = 0;
			int start2 = 0;

			//Deal with line 1:
			if( tLine1->isFree() && beg1_ok )
				start1 = 1;	//start at beginning
			else if( tLine1->isFree() && end1_ok )
				start1 = 2; //start at end
			else if( tLine1->isLeaf() && end1_ok )
				start1 = 2;
			else if( tLine1->isRoot() && beg1_ok )
				start1 = 1;

			//Line 2:
			if( tLine2->isFree() && beg2_ok )
				start2 = 1;	//start at beginning
			else if( tLine2->isFree() && end2_ok )
				start2 = 2; //start at end
			else if( tLine2->isLeaf() && end2_ok )
				start2 = 2;
			else if( tLine2->isRoot() && beg2_ok )
				start2 = 1;

			//Now run appropriate search function:
			int c;
			if(start1 == 0 || start2 == 0)
				continue;
			else if(start1 == 1 && start2 == 1)
				c = IsParallelBeginBegin(tLine1, tLine2);
			else if(start1 == 1 && start2 == 2)
				c = IsParallelBeginEnd(tLine1, tLine2);
			else if(start1 == 2 && start2 == 1)
				c = IsParallelEndBegin(tLine1, tLine2);
			else if(start1 == 2 && start2 == 2)
				c = IsParallelEndEnd(tLine1, tLine2);

			if(c > bestCount)
			{
				bestCount = c;
				bestJ = j;
			}

			//std::cerr << tLine1->GetId() << ", " << tLine2->GetId() << "   c = " << c << std::endl;

		} // end for j

		if(bestJ != -1 && bestCount >= minCount)
		{
			PairT p(i,bestJ);
			//PairT p(lines.at(i)->GetId(), lines.at(bestJ)->GetId());
			closePairs.push_back(p);
			count.push_back(bestCount);
			std::cerr << "    Found Pair = " << p.first << ", " << p.second << " at " << bestCount << std::endl;
		}

	} // end for i

	//Merge the line pairs:
	std::set<int> usedIds;

	for(int p=0; p<(int)closePairs.size(); ++p)
	{
		int id1 = closePairs.at(p).first;
		int id2 = closePairs.at(p).second;

		if( usedIds.find(id1) != usedIds.end() || usedIds.find(id2) != usedIds.end() )
			continue;

		TraceLine * tLine1 = lines.at(id1);
		TraceLine * tLine2 = lines.at(id2);
		
		MergePair(tObj, tLine1, tLine2);

		usedIds.insert(id1);
		usedIds.insert(id2);
	}

}

//Start at end of TraceLine 1 and at the Beginning of Trace Line 2
int IsParallelEndBegin(TraceLine * tLine1, TraceLine * tLine2)
{
	int count = 0;

	//I want to find points that are really close between these lines
	TraceLine::TraceBitsType::iterator it1 = tLine1->GetTraceBitIteratorEnd();
	it1--;
	for( ; it1 != tLine1->GetTraceBitIteratorBegin(); it1--)
	{
		TraceBit tBit1 = *it1;
			
		TraceBit closeBit;
		closeBit.id = -1;
		double closeDist = maxDist;

		TraceLine::TraceBitsType::iterator it2 = tLine2->GetTraceBitIteratorBegin();
		for( ; it2 != tLine2->GetTraceBitIteratorEnd(); it2++)
		{
			TraceBit tBit2 = *it2;

			//I am looking for the closest bit
			double dist = EuclideanDist(tBit1, tBit2, z_penalty_factor);
			if(dist < closeDist)
			{
				closeDist = dist;
				closeBit = tBit2;
			}
		} // end for it2

		if(closeBit.id != -1)
			count++;

	} // end for it1

	return count;
}

//Start at end of TraceLine 1 and at the end of Trace Line 2
int IsParallelEndEnd(TraceLine * tLine1, TraceLine * tLine2)
{
	int count = 0;

	//I want to find points that are really close between these lines
	TraceLine::TraceBitsType::iterator it1 = tLine1->GetTraceBitIteratorEnd();
	it1--;
	for( ; it1 != tLine1->GetTraceBitIteratorBegin(); it1--)
	{
		TraceBit tBit1 = *it1;
			
		TraceBit closeBit;
		closeBit.id = -1;
		double closeDist = maxDist;

		TraceLine::TraceBitsType::iterator it2 = tLine2->GetTraceBitIteratorEnd();
		it2--;
		for( ; it2 != tLine2->GetTraceBitIteratorBegin(); it2--)
		{
			TraceBit tBit2 = *it2;

			//I am looking for the closest bit
			double dist = EuclideanDist(tBit1, tBit2, z_penalty_factor);
			if(dist < closeDist)
			{
				closeDist = dist;
				closeBit = tBit2;
			}
		} // end for it2

		if(closeBit.id != -1)
			count++;

	} // end for it1

	return count;
}

//Start at Beginning of TraceLine 1 and at the end of Trace Line 2
int IsParallelBeginEnd(TraceLine * tLine1, TraceLine * tLine2)
{
	int count = 0;

	//I want to find points that are really close between these lines
	TraceLine::TraceBitsType::iterator it1 = tLine1->GetTraceBitIteratorBegin();
	for( ; it1 != tLine1->GetTraceBitIteratorEnd(); it1++)
	{
		TraceBit tBit1 = *it1;
			
		TraceBit closeBit;
		closeBit.id = -1;
		double closeDist = maxDist;

		TraceLine::TraceBitsType::iterator it2 = tLine2->GetTraceBitIteratorEnd();
		it2--;
		for( ; it2 != tLine2->GetTraceBitIteratorBegin(); it2--)
		{
			TraceBit tBit2 = *it2;

			//I am looking for the closest bit
			double dist = EuclideanDist(tBit1, tBit2, z_penalty_factor);
			if(dist < closeDist)
			{
				closeDist = dist;
				closeBit = tBit2;
			}
		} // end for it2

		if(closeBit.id != -1)
			count++;

	} // end for it1

	return count;
}

//Start at Beginning of TraceLine 1 and at the Beginning of Trace Line 2
int IsParallelBeginBegin(TraceLine * tLine1, TraceLine * tLine2)
{
	int count = 0;

	//I want to find points that are really close between these lines
	TraceLine::TraceBitsType::iterator it1 = tLine1->GetTraceBitIteratorBegin();
	for( ; it1 != tLine1->GetTraceBitIteratorEnd(); it1++)
	{
		TraceBit tBit1 = *it1;
			
		TraceBit closeBit;
		closeBit.id = -1;
		double closeDist = maxDist;

		TraceLine::TraceBitsType::iterator it2 = tLine2->GetTraceBitIteratorBegin();
		for( ; it2 != tLine2->GetTraceBitIteratorEnd(); it2++)
		{
			TraceBit tBit2 = *it2;

			//I am looking for the closest bit
			double dist = EuclideanDist(tBit1, tBit2, z_penalty_factor);
			if(dist < closeDist)
			{
				closeDist = dist;
				closeBit = tBit2;
			}
		} // end for it2

		if(closeBit.id != -1)
			count++;

	} // end for it1

	return count;
}


void MergePair(TraceObject *tObj, TraceLine * tLine1, TraceLine * tLine2)
{
	//Find closest end points:
	TraceLine::TraceBitsType::iterator tb_it;
	TraceBit beg1 = *tLine1->GetTraceBitIteratorBegin();
	tb_it = tLine1->GetTraceBitIteratorEnd();
	tb_it--;
	TraceBit end1 = *tb_it;
	TraceBit beg2 = *tLine2->GetTraceBitIteratorBegin();
	tb_it = tLine2->GetTraceBitIteratorEnd();
	tb_it--;
	TraceBit end2 = *tb_it;

	double d1 = EuclideanDist(beg1, beg2);
	double d2 = EuclideanDist(beg1, end2);
	double d3 = EuclideanDist(end1, beg2);
	double d4 = EuclideanDist(end1, end2);

	if(d1<d2 && d1<d3 && d1<d4)
		MergeBeginBegin(tObj, tLine1, tLine2);
	else if(d2<d1 && d2<d3 && d2<d4)
		MergeBeginEnd(tObj, tLine1, tLine2);
	else if(d3<d1 && d3<d2 && d3<d4)
		MergeEndBegin(tObj, tLine1, tLine2);
	else if(d4<d1 && d4<d2 && d4<d3)
		MergeEndEnd(tObj, tLine1, tLine2);
}

void MergeBeginBegin(TraceObject *tObj, TraceLine * tLine1, TraceLine * tLine2)
{
	TraceBit closeBit1;
	TraceBit closeBit2;
	double closeDist = maxDist;

	TraceLine::TraceBitsType::iterator it1;
	TraceLine::TraceBitsType::iterator it2;

	//I want to find points that are really close between these lines
	it1 = tLine1->GetTraceBitIteratorBegin();
	for( ; it1 != tLine1->GetTraceBitIteratorEnd(); it1++)
	{
		TraceBit tBit1 = *it1;

		it2 = tLine2->GetTraceBitIteratorBegin();
		for( ; it2 != tLine2->GetTraceBitIteratorEnd(); it2++)
		{
			TraceBit tBit2 = *it2;

			//I am looking for the closest bit
			double dist = EuclideanDist(tBit1, tBit2);
			if(dist < closeDist)
			{
				closeDist = dist;
				closeBit1 = tBit1;
				closeBit2 = tBit2;
			}
		} // end for it2
	} // end for it1

	//Now shrink Trace 1 back to closeBit1:
	int count = 0;
	it1 = tLine1->GetTraceBitIteratorBegin();
	for( ; it1 != tLine1->GetTraceBitIteratorEnd(); it1++)
	{
		TraceBit tBit1 = *it1;
		if( tBit1.id == closeBit1.id )
			break;

		count++;
	}
	for(int c=0; c<=count; ++c)
	{
		tLine1->removeFirstBit();
	}


	//Now shrink Trace 2 back to closeBit2:
	count = 0;
	it2 = tLine2->GetTraceBitIteratorBegin();
	for( ; it2 != tLine2->GetTraceBitIteratorEnd(); it2++)
	{
		TraceBit tBit2 = *it2;
		if( tBit2.id == closeBit2.id )
			break;

		count++;
	}
	for(int c=0; c<=count+1; ++c)
	{
		tLine2->removeFirstBit();
	}

	it1 = tLine1->GetTraceBitIteratorBegin();
	it2 = tLine2->GetTraceBitIteratorBegin();
	tObj->mergeTraces((*it1).marker,(*it2).marker);

}

void MergeEndBegin(TraceObject *tObj, TraceLine * tLine1, TraceLine * tLine2)
{
	TraceBit closeBit1;
	TraceBit closeBit2;
	double closeDist = maxDist;

	TraceLine::TraceBitsType::iterator it1;
	TraceLine::TraceBitsType::iterator it2;

	//I want to find points that are really close between these lines
	it1 = tLine1->GetTraceBitIteratorEnd();
	it1--;
	for( ; it1 != tLine1->GetTraceBitIteratorBegin(); it1--)
	{
		TraceBit tBit1 = *it1;

		it2 = tLine2->GetTraceBitIteratorBegin();
		for( ; it2 != tLine2->GetTraceBitIteratorEnd(); it2++)
		{
			TraceBit tBit2 = *it2;

			//I am looking for the closest bit
			double dist = EuclideanDist(tBit1, tBit2);
			if(dist < closeDist)
			{
				closeDist = dist;
				closeBit1 = tBit1;
				closeBit2 = tBit2;
			}
		} // end for it2
	} // end for it1

	

	//Now shrink Trace 1 back to closeBit1:
	int count = 0;
	it1 = tLine1->GetTraceBitIteratorEnd();
	it1--;
	for( ; it1 != tLine1->GetTraceBitIteratorBegin(); it1--)
	{
		TraceBit tBit1 = *it1;
		if( tBit1.id == closeBit1.id )
			break;

		count++;
	}
	for(int c=0; c<=count; ++c)
	{
		tLine1->removeLastBit();
	}


	//Now shrink Trace 2 back to closeBit2:
	count = 0;
	it2 = tLine2->GetTraceBitIteratorBegin();
	for( ; it2 != tLine2->GetTraceBitIteratorEnd(); it2++)
	{
		TraceBit tBit2 = *it2;
		if( tBit2.id == closeBit2.id )
			break;

		count++;
	}
	for(int c=0; c<=count+1; ++c)
	{
		tLine2->removeFirstBit();
	}

	it1 = tLine1->GetTraceBitIteratorEnd();
	it1--;
	it2 = tLine2->GetTraceBitIteratorBegin();
	tObj->mergeTraces((*it1).marker,(*it2).marker);

}

void MergeBeginEnd(TraceObject *tObj, TraceLine * tLine1, TraceLine * tLine2)
{
	TraceBit closeBit1;
	TraceBit closeBit2;
	double closeDist = maxDist;

	TraceLine::TraceBitsType::iterator it1;
	TraceLine::TraceBitsType::iterator it2;

	//I want to find points that are really close between these lines
	it1 = tLine1->GetTraceBitIteratorBegin();
	for( ; it1 != tLine1->GetTraceBitIteratorEnd(); it1++)
	{
		TraceBit tBit1 = *it1;

		it2 = tLine2->GetTraceBitIteratorEnd();
		it2--;
		for( ; it2 != tLine2->GetTraceBitIteratorBegin(); it2--)
		{
			TraceBit tBit2 = *it2;

			//I am looking for the closest bit
			double dist = EuclideanDist(tBit1, tBit2);
			if(dist < closeDist)
			{
				closeDist = dist;
				closeBit1 = tBit1;
				closeBit2 = tBit2;
			}
		} // end for it2
	} // end for it1

	//Now shrink Trace 1 back to closeBit1:
	int count = 0;
	it1 = tLine1->GetTraceBitIteratorBegin();
	for( ; it1 != tLine1->GetTraceBitIteratorEnd(); it1++)
	{
		TraceBit tBit1 = *it1;
		if( tBit1.id == closeBit1.id )
			break;

		count++;
	}
	for(int c=0; c<=count; ++c)
	{
		tLine1->removeFirstBit();
	}


	//Now shrink Trace 2 back to closeBit2:
	count = 0;
	it2 = tLine2->GetTraceBitIteratorEnd();
	it2--;
	for( ; it2 != tLine2->GetTraceBitIteratorBegin(); it2--)
	{
		TraceBit tBit2 = *it2;
		if( tBit2.id == closeBit2.id )
			break;

		count++;
	}
	for(int c=0; c<=count+1; ++c)
	{
		tLine2->removeLastBit();
	}

	it1 = tLine1->GetTraceBitIteratorBegin();
	it2 = tLine2->GetTraceBitIteratorEnd();
	it2--;
	tObj->mergeTraces((*it1).marker,(*it2).marker);

}

void MergeEndEnd(TraceObject *tObj, TraceLine * tLine1, TraceLine * tLine2)
{
	TraceBit closeBit1;
	TraceBit closeBit2;
	double closeDist = maxDist;

	TraceLine::TraceBitsType::iterator it1;
	TraceLine::TraceBitsType::iterator it2;

	//I want to find points that are really close between these lines
	it1 = tLine1->GetTraceBitIteratorEnd();
	it1--;
	for( ; it1 != tLine1->GetTraceBitIteratorBegin(); it1--)
	{
		TraceBit tBit1 = *it1;

		it2 = tLine2->GetTraceBitIteratorEnd();
		it2--;
		for( ; it2 != tLine2->GetTraceBitIteratorBegin(); it2--)
		{
			TraceBit tBit2 = *it2;

			//I am looking for the closest bit
			double dist = EuclideanDist(tBit1, tBit2);
			if(dist < closeDist)
			{
				closeDist = dist;
				closeBit1 = tBit1;
				closeBit2 = tBit2;
			}
		} // end for it2
	} // end for it1

	//Now shrink Trace 1 back to closeBit1:
	int count = 0;
	it1 = tLine1->GetTraceBitIteratorEnd();
	it1--;
	for( ; it1 != tLine1->GetTraceBitIteratorBegin(); it1--)
	{
		TraceBit tBit1 = *it1;
		if( tBit1.id == closeBit1.id )
			break;

		count++;
	}
	for(int c=0; c<=count; ++c)
	{
		tLine1->removeLastBit();
	}


	//Now shrink Trace 2 back to closeBit2:
	count = 0;
	it2 = tLine2->GetTraceBitIteratorEnd();
	it2--;
	for( ; it2 != tLine2->GetTraceBitIteratorBegin(); it2--)
	{
		TraceBit tBit2 = *it2;
		if( tBit2.id == closeBit2.id )
			break;

		count++;
	}
	for(int c=0; c<=count+1; ++c)
	{
		tLine2->removeLastBit();
	}

	it1 = tLine1->GetTraceBitIteratorEnd();
	it1--;
	it2 = tLine2->GetTraceBitIteratorEnd();
	it2--;
	tObj->mergeTraces((*it1).marker,(*it2).marker);

}

double EuclideanDist(TraceBit b1, TraceBit b2, double z_penalty_factor)
{
	double dx = b1.x - b2.x;
	double dy = b1.y - b2.y;
	double dz = z_penalty_factor*(b1.z - b2.z);
	double dist = sqrt(dx*dx + dy*dy + dz*dz);
	return dist;
}