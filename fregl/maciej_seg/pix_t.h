/** @file pix_t.h
*   @brief class for pix
*   This is the class that represents the pix in the image
*
*   @author Maciej Wotjon
*/

#ifndef __PIX_T_H_
#define __PIX_T_H_

/** @brief pix class
*/
class pix_t
{
public:
	//coordinates
	int x_, y_, z_;

	pix_t( int x= -1, int y= -1, int z= -1)
		:x_(x)
		, y_(y)
		, z_(z)
	{
	};

	pix_t(const pix_t& p )
		:x_(p.x_)
		, y_(p.y_)
		, z_(p.z_)
	{
	};

	~pix_t(void)
	{
	};
};

#endif
