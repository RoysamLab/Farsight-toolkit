#include "TraceObject.h"
#include "CellTraceModel.h"
#include "CellTrace.h"
#include <string>
#include <limits>

int main(int argc, char* argv[])
{
	TraceObject* tobj = new TraceObject();

	if (argc < 2)
	{
		std::cout << "Must run program with correct number of arguments" << std::endl;
		return -1;
	}
	
	char* xml_file = argv[1];
	std::cout << "Reading XML file: " << xml_file << std::endl;
	
	tobj->ReadFromSWCFile(xml_file);
	CellTraceModel* cells = new CellTraceModel(tobj->CalculateCellFeatures());

	double maxX = 0;
	double minX = std::numeric_limits<double>::max();
	double maxY = 0;
	double minY = std::numeric_limits<double>::max();
	double maxZ = 0;
	double minZ = std::numeric_limits<double>::max();


	for (int k = 0; k < cells->getCellCount(); k++)
	{
		CellTrace* cell = cells->GetCellAt(k);
		
		double bounds[6];
		cell->getCellBounds(bounds);

		std::cout << "Min X: " << bounds[0] << " Max X: " << bounds[1] << " Min Y: " << bounds[2] << " Max Y: " << bounds[3] << " Min Z: " << bounds[4] << " Max Z: " << bounds[5] << std::endl;
	
		if (bounds[1] > maxX)
			maxX = bounds[1];
		if (minX > bounds[0])
			minX = bounds[0];

		if (bounds[3] > maxY)
			maxY = bounds[3];
		if (minY > bounds[2])
			minY = bounds[2];

		if (bounds[5] > maxZ)
			maxZ = bounds[5];
		if (minZ > bounds[4])
			minZ = bounds[4];
	}

	std::cout << "Image bounds: (" << minX << ", " << minY << ", " << minZ << ") to (" << maxX << ", " << maxY << ", " << maxZ << ")" << std::endl; 

	return 0;
}