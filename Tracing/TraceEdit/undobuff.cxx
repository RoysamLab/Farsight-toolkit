#include "undobuff.h"

/*
int main(int argc, char *argv[])
{
	undoBuffer<std::string> buffer;
	std::string testString = "ABCDEFGHIJ";
	for(unsigned int i = 0; i < 10; i++)
	{
		buffer.Add(testString.substr(i, 1));
		buffer.print();
	}
	std::cout << std::endl;
	
	for(unsigned int i = 0; i < 11; i++)
	{
		if(buffer.UndoOrRedo(0))
		{
			std::cout << "Current Pos: " << buffer.getCurrPos() << std::endl;
			std::cout << "Current state: " << buffer.getState() << std::endl;
			buffer.print();
			std::cout << std::endl;
		}
	}
		
	for(unsigned int i = 0 ; i < 12; i++)
	{
		if(buffer.UndoOrRedo(1))
		{
			std::cout << "Current Pos: " << buffer.getCurrPos() << std::endl;
			std::cout << "Current State: " << buffer.getState() << std::endl;
			buffer.print();
			std::cout << std::endl;
		}
	}
		
	std::string test3 = "1234567890";
	for(unsigned int i = 0; i < 10; i++)
	{
		buffer.Add(test3.substr(i,1));
		buffer.print();
	}
	std::cout << std::endl;
	return 0;
}
*/
