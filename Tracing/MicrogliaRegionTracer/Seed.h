#include <cstddef>

class Seed
{
public:
	Seed(size_t x, size_t y, size_t z);
	
	size_t getX() const;
	size_t getY() const;
	size_t getZ() const;

private:
	size_t x;
	size_t y;
	size_t z;
};
