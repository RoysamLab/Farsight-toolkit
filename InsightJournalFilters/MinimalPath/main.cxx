#if defined(_MSC_VER)
//Warning about: identifier was truncated to '255' characters in the debug information (MVC6.0 Debug)
#pragma warning( disable : 4786 )
#endif

// General includes
#include <iostream>

// ITK includes
#include "Testing/MinimalPathTests.cxx"

int main(int argc, char* argv[])
{
    return Test_SpeedToPath_RegularStepGradientDescent_2D( argc, argv );
}
