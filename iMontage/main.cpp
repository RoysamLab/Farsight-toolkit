#include "TraceObject.h"
#include "CellTraceModel.h"
#include "CellTrace.h"
#include <string>

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

	for (int k = 0; k < cells->getCellCount(); k++)
	{
		CellTrace* cell = cells->GetCellAt(k);
		
		double bounds[6];
		cell->getCellBounds(bounds);

		//std::cout << "Min X: " << bounds[0] << " Max X: " << bounds[1] << " Min Y: " << bounds[2] << " Max Y: " << bounds[3] << " Min Z: " << bounds[4] << " Max Z: " << bounds[5] << std::endl;
	}

	return 0;
}