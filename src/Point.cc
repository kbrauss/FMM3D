/** Point.cc
 *  Created on: June 20, 2016
 *      Author: Keith D Brauss
 * */

#include <cmath>
#include <sstream>
#include <string>
#include <iostream>
#include <limits>  // DBL_EPSILON
#include <vector>

#include "Point.h"
#include "Util.h"

/**
 * Header Interface for Class Point
 *
class Point
{
  public:
	std::vector<double> coord; // coordinates for point in three-dimensional space

	// constructors
	Point();
    //Point(double x, double y, double z);
    Point(std::vector<double> coord);
    //Point(const Point &p);

    std::vector<double> getCoord() { return coord; };
    void setCoord(std::vector<double> coord) { this->coord = coord; };
    void   setX(double x_coord) { this->coord[0]=x_coord; };
    double getX() { return coord[0]; };
    void   setY(double y_coord) { this->coord[1]=y_coord; };
    double getY() { return coord[1]; };
    void   setZ(double z_coord) { this->coord[2]=z_coord; };
    double getZ() { return coord[2]; };

    std::string coordToString();
    bool equals(Point &p);
    int getBoxIndex(unsigned int level);
};
 */

Point::Point()
    :
    coord(3)
{
  this->coord[0] = 0.0;
  this->coord[1] = 0.0;
  this->coord[2] = 0.0;
}

Point::Point(double x, double y, double z)
    :
    coord(3)
{
  coord[0] = x;
  coord[1] = y;
  coord[2] = z;
}


Point::Point(std::vector<double> coord)
    :
    coord(coord)
{}

//inline Point::Point(double x, double y)
//    : coord(x,y)
//{}

//inline Point::Point(std::complex<double> c)
//    : coord(c.real(),c.imag())
//{}

//inline Point::Point(const Point &p)
//    : coord(p.coord.real(),p.coord.imag())
//{}

std::string Point::coordToString()
{
  std::string result;     // string which will contain the result
  std::ostringstream convert;  // stream used for the conversion
  convert  << "Point("<< this->coord[0] << "," << this->coord[1] << "," << this->coord[2] << ")";
  result = convert.str(); // set 'result' to the contents of the stream
  return result;
}

bool Point::equals(Point &p)
{
  std::vector<double> p_Coords = p.getCoord();
  double diff = 0.0;
  for (unsigned int i=0; i<3; ++i)
	  diff += std::abs(this->coord[i] - p_Coords[i]);
  if (diff < std::numeric_limits<double>::epsilon() )
    return true;
  else
    return false;
}

/**
 * Explanation of getBoxIndex member function
 *
 * The function uses the Util class member function interleave to
 * get the index number (n) of the cell (box) that this Point point is located in
 * for the refinement level that has been specified.
 *
 * getBoxIndex consists of just one call to the interleave member function of
 * the class Util.  The member function interleave of class Util is therefore
 * explained here using the example(s) below.  Our examples are based on
 * a level l = 2 refinement.
 * The indexing for level 1 and 2 is shown below as well as the location
 * of the two example points x and o.  We show how getBoxIndex determines
 * the level l=2 cells where the points are located.
 *
 * Before going to the examples below, it may be helpful to look at the
 * explanation for the setbit and getbit member functions of the class Util
 * in the Util.cc file.  There you can find a review of the conversion from
 * binary to decimal.  The setting of the 8 bits (assuming 8 bit integer
 * - probably really 32 bit integer) that define the decimal value of an
 * integer (corresponds to the index number of the cell here) is the main tool used in
 * locating the cell where the point lies.
 *
 *
 *  This is level l=2 refinement
 *
 *  * At level l=2, we have five vertical layers (layers 0,1,2,3,4) that slice through the domain
 *  * The top layer (layer 0) that would be at the "front" of the domain and located
 *    at the bottom left below is missing.
 *  * This is so that we can "peer into the domain" (and save space on the screen) to see the ordering
 *    that is inside.
 *  * The numbering that is shown lies between vertical layers and "on top" of the vertical layer
 *    where it is located.
 *  * At level l=1 there are a total of 8 cells.  The cells are labeled below as well.
 *  * At level l=1 there are only three layers (layers 0,2, and 4 - again layer 0 is removed)
 *    We would not see layers 1 and 3 at l=1.  And so you see repeats of the numbers ordering
 *    these cells, as these cells pass through layers 1 and 3.
 *                                                                                                                   y = 1.00
 *                                                                                                         1.00 ______ ______ ______ ______
 *                                                                                                             |      |      |      |      |
 *                                                                                 y = 0.75                    |  27  |  31  |  59  |  63  |
 *                                                                                                         0.75|______3______|______7______|
 *                                                                       1.00 ______ ______ ______ ______      |      |      |      |      |
 *       |                                                                   |      |      |      |      |     |  26  |  30  |  58  |  62  |
 *                                               y = 0.50                    |  25  |  29  |o 57  |  61  | 0.50|______|______|______|______|
 *                                                                       0.75|______3______|______7______|     |      |      |      |      |
 *                                     1.00 ______ ______ ______ ______      |      |      |      |      |     |  19  |  23  |  51  |  55  |
 *      z                                  |      |      |      |      |     |  24  |  28  |  56  |  60  | 0.25|______2______|______6______|
 *       |      y = 0.25                   |  11  |  15  |  43  |  47  | 0.50|______|______|______|______|     |      |      |      |      |
 *       |                             0.75|______1______|______5______|     |      |      |      |      |     |  18  |  22  |  50  |  54  |
 *  1.00 |______ ______ ______ ______      |      |      |      |      |     |  17  |  21  |  49  |  53  |     |______|______|______|______|
 *       |      |      |      |      |     |  10  |  14  |  42  |  46  | 0.25|______2______|______6______|    0.0    0.25   0.5    0.75   1.0       x
 *       |  9   |  13  |  41  |  45  | 0.50|______|______|______|______|     |      |      |      |      |
 *  0.75 |______1______|______5______|     |    x |      |      |      |     |  16  |  20  |  48  |  52  |
 *       |      |      |      |      |     |  3   |   7  |  35  |  39  |     |______|______|______|______|
 *       |  8   |  12  |  40  |  44  | 0.25|______0______|______4______|    0.0    0.25   0.5    0.75   1.0       x
 *  0.50 |______|______|______|______|     |      |      |      |      |
 *       |      |      |      |      |     |  2   |   6  |  34  |  38  |
 *       |  1   |  5   |  33  |  37  |     |______|______|______|______|
 *  0.25 |______0______|______4______|    0.0    0.25   0.5    0.75   1.0       x
 *       |      |      |      |      |
 *       |  0   |  4   |  32  |  36  |
 *       |______|______|______|______|______
 *      0.0    0.25   0.5    0.75   1.0       x
 *
 *
 *     Example 1:
 *
 *     Let
 *        x - has coordinates (0.1875,0.375,0.3125)
 *        x is in cell (n,l) = (3,2)
 *        where cell number n = 3 and refinement level l = 2
 *
 *     Note that in terms of bits of an integer
 *     00000001 = 2^0 = 1
 *     00000010 = 2^1 = 2
 *     00000100 = 2^2 = 4
 *
 *     We review the bit interleaving that is being performed here.
 *
 *     On the first pass through the for loop
 *        floor(coord[0]*std::pow(2,level)) = floor(0.1875*2^2)
 *                                              = floor(0.1875*4)
 *                                              = floor(0.75)
 *                                              = 0
 *        floor(coord[1]*std::pow(2,level)) = floor(0.375*2^2)
 *                                              = floor(0.375*4)
 *                                              = floor(1.5)
 *                                              = 1
 *        floor(coord[2]*std::pow(2,level)) = floor(0.3125*2^2)
 *                                              = floor(0.3125*4)
 *                                              = floor(1.25)
 *                                              = 1
 *     If we start at the (0,0,0) corner at the front bottom left of the
 *     figure above and count cell lengths in each coordinate direction
 *     corresponding to the values returned from the floor function calls,
 *     we end up at the front-lower-left corner of the cell containing the
 *     point of interest.
 *     Starting the front-lower-left corner of cell 0 (where (0,0,0) is located)
 *     and counting 0 cell lengths in the x-direction, 1 cell length in the y-direction
 *     and 1 cell length in the z-direction we arrive at the front-lower-left corner of
 *     cell 3 where the point x with coordinates (0.1875,0.375,0.3125) is located.
 *
 *     Applying getBoxIndex to this point results in
 *     the call to the Util class's member function interleave
 *        Util::interleave(int x, int y, int z, int level)
 *        Util.interleave(0,1,1,2)
 *     The interleave function passes in x=0, y=1, z=1, level=2
 *     Then calls setbit of the Util class in the for loop below
 *
 *     Note that we will only loop two passes since we are at level l=2
 *
 *          int ans = 0;
 *          for (unsigned int i=0; i<level=2; ++i)
 *          {
 *            ans = setbit(ans, (level-i)*3-1, getbit(x, level-i-1));
 *            ans = setbit(ans, (level-i)*3-2, getbit(y, level-i-1));
 *            ans = setbit(ans, (level-i)*3-3, getbit(z, level-i-1));
 *          }
 *          return ans;
 *
 *     On the first pass (i=0) through the for loop we have
 *        ans = setbit(int n, int pos, int setto)
 *            = setbit(ans, (level-i)*3-1, getbit(x, level-i-1))
 *            = setbit(0, (2-0)*3-1, getbit(0, 2-0-1))
 *            = setbit(0, 2*3-1, getbit(0,1))     // 0 = 00000000
 *            = setbit(0, 5, getbit(000000*0, 1)) // getting bit 1 of int 0 (starred bit)
 *            = setbit(0, 5, 0)  // set int 0 bit 5 to 0 (but already 0)
 *            = 0
 *
 *        ans = setbit(int n, int pos, int setto)
 *            = setbit(ans, (level-i)*3-2, getbit(y, level-i-1))
 *            = setbit(0, (2-0)*3-2, getbit(1,2-0-1))
 *            = setbit(0, 4, getbit(1,1))         // 1 = 00000001
 *            = setbit(0, 4, getbit(000000*1, 1)) // getting bit 1 of int 1 (starred)
 *            = setbit(0, 4, 0)  // set int 0 bit 4 to 0 (but already 0)
 *            = setbit(0, 4, 0)
 *            = 0
 *
 *        ans = setbit(int n, int pos, int setto)
 *            = setbit(ans, (level-i)*3-3, getbit(z, level-i-1))
 *            = setbit(0, (2-0)*3-3, getbit(1,2-0-1))
 *            = setbit(0, 3, getbit(1,1))         // 1 = 00000001
 *            = setbit(0, 3, getbit(000000*1, 1)) // getting bit 1 of int 1 (starred)
 *            = setbit(0, 3, 0)  // set int 0 bit 3 to 0 (but already 0)
 *            = setbit(0, 3, 0)
 *            = 0
 *
 *     On the second pass (i=1) through the for loop we have
 *        ans = setbit(int n, int pos, int setto)
 *            = setbit(ans, (level-i)*3-1, getbit(x, level-i-1))
 *            = setbit(0, (2-1)*3-1, getbit(0, 2-1-1))
 *            = setbit(0, 1*3-1, getbit(0,0))     // 0 = 00000000
 *            = setbit(0, 2, getbit(00000000, 0)) // getting bit 0 of int 0 (this returns 0)
 *            = setbit(0, 2, 0)  // set int 0 bit 2 to 0 (results in 00000000 set to 00000000)
 *            = 0  (00000100 = 1*2^2 = 1*4 = 4)
 *
 *        ans = setbit(int n, int pos, int setto)
 *            = setbit(ans, (level-i)*3-2, getbit(y, level-i-1))
 *            = setbit(0, (2-1)*3-2, getbit(1,2-1-1))
 *            = setbit(0, 1, getbit(1,0))         // 1 = 00000001
 *            = setbit(0, 1, getbit(00000001,0))  // getting bit 0 of int 1 (this returns 1)
 *            = setbit(0, 1, 1)  // set int 0 bit 1 to 1 (results in 00000000 set to 00000010
 *            = 2  (00000010 = 1*2^1 = 2)
 *
 *        ans = setbit(int n, int pos, int setto)
 *            = setbit(ans, (level-i)*3-3, getbit(z, level-i-1))
 *            = setbit(2, (2-1)*3-3, getbit(1,2-1-1))
 *            = setbit(2, 0, getbit(1,0))         // 1 = 00000001
 *            = setbit(2, 0, getbit(00000001,0))  // getting bit 0 of int 1 (this returns 1)
 *            = setbit(2, 0, 1)  // set int 2 bit 0 to 1 (results in 00000010 set to 00000011
 *            = 3  (00000011 = 1*2^1 + 1*2^0 = 2 + 1 = 3)
 *
 *     And the position of point x at level two is in cell (n,l)
 *     with box index n = 3 and level l = 2 (see figure above for location of x)
 *
 *
 *     Example 2:
 *
 *     Let
 *        o - has coordinates (0.5625,0.625,0.875)
 *        x is in cell (n,l) = (57,2)
 *        where cell number n = 57 and refinement level l = 2
 *     Note that
 *     00111001 = 1*2^5 + 1*2^4 + 1*2^3 + 0*2^2 + 0*2^1 + 1*2^0 = 32 + 16 + 8 + 1 = 57
 *
 *     That is, we expect the interleave function to return 57.  There is only one way
 *     that it can do this in terms of bits - shown above.
 *
 *     From examples we deduce that the first loop locates the cell that the
 *     point is in at level 1. (At level zero, there is only one cell.
 *     So there is not much to determine in terms of what cell).
 *     Another way to view the first pass is statement one (x) determines whether
 *     the point is on the left half or right half of the domain.  The second
 *     command (y) determines whether the  point is on the closer or further half
 *     (into the screen - y) of the domain.  The third command of the first pass (z) determines
 *     if the point is in the upper half or lower half of the domain.  For level 2, this correlates
 *     to 0 or 1 at bit 5 (x), 0 or 1 at bit 4 (y), 0 or 1 at bit 3 (z).
 *
 *     Once we have reached this location, then we
 *
 *
 *     Then
 *        floor(coord[0]*std::pow(2,level)) = floor(0.5625*2^2)
 *                                              = floor(0.5625*4)
 *                                              = floor(2.25)
 *                                              = 2
 *        floor(coord[1]*std::pow(2,level)) = floor(0.625*2^2)
 *                                              = floor(0.625*4)
 *                                              = floor(2.5)
 *                                              = 2
 *        floor(coord[2]*std::pow(2,level)) = floor(0.875*2^2)
 *                                              = floor(0.875*4)
 *                                              = floor(3.5)
 *                                              = 3
 *     Looking at the illustration and the last two examples, it looks
 *     We can see the product and floor functions return values
 *     (x,y,z) = (2,2,3) are giving the increments necessary
 *     to reach the cell from the lower left corner of the domain (0,0,0).
 *
 *     Therefore applying getBoxIndex to this point results in
 *     the call to the Util class's member function interleave
 *        Util::interleave(int x, int y, int z, int level)
 *        Util.interleave(2,2,3,2)
 *     The interleave function passes in x=2, y=2, z=3, level=2
 *     Then calls setbit of the Util class in the for loop below
 *
 *          int ans = 0;
 *          for (unsigned int i=0; i<level=2; ++i)
 *          {
 *            ans = setbit(ans, (level-i)*3-1, getbit(x, level-i-1));
 *            ans = setbit(ans, (level-i)*3-2, getbit(y, level-i-1));
 *            ans = setbit(ans, (level-i)*3-3, getbit(y, level-i-1));
 *          }
 *          return ans;
 *
 *     On the first pass (i=0) through the for loop we have
 *        ans = setbit(int n, int pos, int setto)
 *            = setbit(ans, (level-i)*3-1, getbit(x, level-i-1))
 *            = setbit(0, (2-0)*3-1, getbit(2, 2-0-1))
 *            = setbit(0, 2*3-1, getbit(2,1))     // 2 = 00000010
 *            = setbit(0, 5, getbit(00000010, 1)) // getbit 1 of int 2 (which is second bit from right - 1 )
 *            = setbit(0, 5, 1)  // set int 0 bit 5 to 1 -> setting 00000000 to 00100000 = 2^5 = 32)
 *            = 32
 *
 *        ans = setbit(int n, int pos, int setto)
 *            = setbit(ans, (level-i)*3-2, getbit(y, level-i-1))
 *            = setbit(32, (2-0)*3-2, getbit(2,2-0-1))
 *            = setbit(32, 4, getbit(2,1))
 *            = setbit(32, 4, getbit(000000010, 1)) // getbit 1 of int 2 (which is second bit from right - 1)
 *            = setbit(32, 4, 1)  // set int 32 bit 4 to 1 -> setting 00100000 to 00110000
 *            = 48  (1*2^5 + 1*2^4 = 32 + 16 = 48)
 *
 *        ans = setbit(int n, int pos, int setto)
 *            = setbit(ans, (level-i)*3-2, getbit(z, level-i-1))
 *            = setbit(48, (2-0)*3-3, getbit(3,2-0-1))
 *            = setbit(48, 3, getbit(3,1))
 *            = setbit(48, 3, getbit(000000011, 1)) // getbit 1 of int 3 (which is second bit from right - 1)
 *            = setbit(48, 3, 1)  // set int 48 bit 3 to 1 -> setting 00110000 to 00111000
 *            = 56  (1*2^5 + 1*2^4 + 1*2^3 = 32 + 16 + 8 = 56)
 *
 *     On the second pass (i=1) through the for loop we have
 *        ans = setbit(int n, int pos, int setto)
 *            = setbit(ans, (level-i)*3-1, getbit(2, level-i-1))
 *            = setbit(56, (2-1)*3-1, getbit(2, 2-1-1))
 *            = setbit(56, 1*3-1, getbit(2,0))     // 2 = 00000010
 *            = setbit(56, 2, getbit(00000010, 0)) // getbit 0 of int 2 (which is first bit from right - 0)
 *            = setbit(56, 2, 0)  // set int 56 bit 2 to 0 (00111000 to 00111000 - no change)
 *            = 56
 *
 *        ans = setbit(int n, int pos, int setto)
 *            = setbit(ans, (level-i)*3-2, getbit(y, level-i-1))
 *            = setbit(56, (2-1)*3-2, getbit(2,2-1-1))
 *            = setbit(56, 1, getbit(2,0))         // 2 = 00000010
 *            = setbit(56, 1, getbit(00000010, 0)) // getbit 0 of int 2 (which is first bit from right - 0)
 *            = setbit(56, 1, 0)  // set int 56 bit 1 to 0 (00111000 to 00111000 - no change)
 *            = 56
 *
 *        ans = setbit(int n, int pos, int setto)
 *            = setbit(ans, (level-i)*3-3, getbit(z, level-i-1))
 *            = setbit(56, (2-1)*3-3, getbit(3,2-1-1))
 *            = setbit(56, 0, getbit(3,0))         // 2 = 00000010
 *            = setbit(56, 0, getbit(00000011, 0)) // getbit 0 of int 3 (which is first bit from right - 1)
 *            = setbit(56, 0, 1)  // set int 56 bit 0 to 1 (00111000 to 00111001 - no change)
 *            = 57 (1*2^5 + 1*2^4 + 1*2^3 + 1*2^0 = 32 + 16 + 8 +1 = 57)
 *
 *     And the position of point o at level two is in cell (n,l)
 *     with box index n = 57 and level l = 2 (see figure above for location of o)
 *
 *
 *     Remarks:
 *
 *     From the two examples above can see that for level l = 2
 *     the for loop goes from i = 0 up to i = 1 < 2
 *
 *     Within the three statements in the for loop
 *            ans = setbit(ans, (level-i)*3-1, getbit(x, level-i-1));
 *            ans = setbit(ans, (level-i)*3-2, getbit(y, level-i-1));
 *            ans = setbit(ans, (level-i)*3-3, getbit(z, level-i-1));
 *     are the arguments
 *       (level-i)*3-1
 *       (level-i)*3-2
 *       (level-i)*3-3
 *
 *     We can see from the examples that these two arguments locate the bit that is
 *     to be set to zero or one
 *     For the level l = 2
 *       For i = 0
 *         we have (level - i)*3 - 1 = (2 - 0)*3 - 1 = 6 - 1 = 5
 *         and     (level - i)*3 - 2 = (2 - 0)*3 - 2 = 6 - 2 = 4
 *         and     (level - i)*3 - 2 = (2 - 0)*3 - 3 = 6 - 3 = 3
 *       For i = 1
 *         we have (level - i)*3 - 1 = (2 - 1)*3 - 1 = 3 - 1 = 2
 *         and     (level - i)*3 - 2 = (2 - 1)*3 - 2 = 3 - 2 = 1
 *         and     (level - i)*3 - 3 = (2 - 1)*3 - 3 = 3 - 3 = 0
 *
 *     Further, the loops count down the bit locations as it progresses.
 *     Here for level l=2, the count is from bit location 5 to bit location 0
 *     In terms of bit positions of an integer written in binary, this is
 *     00*00000 (bit position 5), 000*0000 (bit position 4), 0000*000 (bit position 3)
 *     00000*00 (bit position 2), 000000*0 (bit position 1), 0000000* (bit position 0)
 *
 *     Depending on the coordinates of the point of interest, each of these bits will
 *     be set to 0 or 1.  We can see that the largest number possible for level l=2
 *     is when all 6 of the bits are set to 1 (from position 5 to position 0) resulting in
 *     the number
 *                00111111 = 0*2^7 + 0*2^6 + 1*2^5 + 1*2^4 + 1*2^3 + 1*2^2 + 1*2^1 + 1*2^0
 *                         = 0     + 0     + 32    + 16    + 8     + 4     + 2     + 1
 *                         = 40 + 20 + 3 = 63
 *     And this is exactly the number of cells (64 cells counted from 0 to 63) in the refinement level l=2
 *
 *     We also see if the data type int holds 32 bits (used to be amount for unsigned long - can find out with
 *     std::numeric_limits<unsigned>::digits, but need #include <limits>, then the largest amount of cells
 *     possible occurs at level
 *     l=7 and this
 *     is
 *     11111111111111111111111111111111 = 1*2^31     + 1*2^30     + 1*2^29    + ... + 1*2^4 + 1*2^3 + 1*2^2 + 1*2^1 + 1*2^0
 *                                      = 2147483648 + 1073741824 + 536870912 + ... + 16    + 8     + 4     + 2     + 1
 *                                      = 4294967295 < 2^32 = 4294967296
 *     So at level 32 (l = 32) we have 4294967296 cells counted from 0 to 4294967295
 *
 *     Notes
 *
 *     The first pass through the for loop
 *      - sets the 5th bit or does not
 *        - if the bit is set to 0 (0*2^5 = 0)
 *          - the location of the point will lie in the left half of the domain
 *            looking at it above (left half of x-values)
 *            (set bit to 0 implies cell number n = 0*2^5 = 0 or higher (up to 31))
 *        - if the bit is set to 1 (1*2^5 = 32)
 *          - the location of the point will be in right half of the domain
 *            (set bit to 1 implies cell number n = 1*2^5=32 or higher (up to 63))
 *      - sets the 4th bit or does not
 *        - if the bit is set to 1 ( 1*2^4 = 16) the location of the point will be in the second half (back half)
 *          of the y-values - the further away half of the domain w.r.t the y-values (into the screen).
 *          - (i) the location of the point will now be in the left half of the domain (x) and furthest away (y)
 *            (cell number n = 0*2^5 + 1*2^4 = 0 + 16 = 16 or higher (up to n = 16 + 1*2^3 + 1*2^2 + 1*2^1 + 1*2^0 = 31) )
 *          - (ii) the location of the point will be in the right half of the domain (x) and furthest away (y)
 *            (cell number n = 1*2^5 + 1*2^4 = 32 + 16 = 48 or higher (up to n = 48 + 1*2^3 + 1*2^2 + 1*2^1 + 1*2^0 = 63))
 *        - if the bit is set to 0 (0*2^4 = 0) the location of the point will be in the closer half w.r.t. y (front half)
 *          - (i) the location of the point will now be in the left half and front half
 *            (cell number n = 0*2^5 + 0*2^4 = 0 + 0 = 0 or higher (up to n = 0 + 1*2^3 + 1*2^2 + 1*2^1 + 1*2^0 = 15) )
 *          - (ii) the location of the point will be in right half (x) and front half (y)
 *            (cell number n = 1*2^5 + 0*2^4 = 32 + 0 = 32 or higher (up to 32 + 1*2^3 + 1*2^2 + 1*2^1 + 1*2^0 = 47))
 *      - sets the 3rd bit or does not
 *        - if the bit is set to 1 ( 1*2^3 = 8) the location of the point will be in the upper half
 *          of the z-values - vertical upper half according to figure.
 *          - (i) the location of the point will now be in the left half (x) and front half (y) and upper half (z)
 *            (cell number n = 0*2^5 + 0*2^4 + 1*2^3 = 0 + 0 + 8 = 8 or higher (up to n = 16 + 1*2^2 + 1*2^1 + 1*2^0 = 23) )
 *          - (ii) the location of the point will be in the left half (x) and back half (y) and upper half (z)
 *            (cell number n = 0*2^5 + 1*2^4 + 1*2^3 = 0 + 16 + 8 = 24 or higher (up to n = 24 + 1*2^2 + 1*2^1 + 1*2^0 = 31))
 *          - (iii) the location of the point will now be in the right half (x) and front half (y) and upper half (z)
 *            (cell number n = 1*2^5 + 0*2^4 + 1*2^3 = 32 + 0 + 8 = 40 or higher (up to n = 40 + 1*2^2 + 1*2^1 + 1*2^0 = 47) )
 *          - (iv) the location of the point will be in the right half (x) and back half (y) and upper half (z)
 *            (cell number n = 1*2^5 + 1*2^4 + 1*2^3 = 32 + 16 + 8 = 56 or higher (up to n = 56 + 1*2^2 + 1*2^1 + 1*2^0 = 63))
 *        - if the bit is set to 0 (0*2^3 = 0) the location of the point will be in the lower half w.r.t. z (vertical lower half)
 *          - (i) the location of the point will now be in the left half (x) and front half (y) and lower half (z)
 *            (cell number n = 0*2^5 + 0*2^4 + 0*2^3 = 0 + 0 + 0 = 0 or higher (up to n = 0 + 1*2^2 + 1*2^1 + 1*2^0 = 7) )
 *          - (ii) the location of the point will be in right half (x) and front half (y) and lower half (z)
 *            (cell number n = 1*2^5 + 0*2^4 + 0*2^3 = 32 + 0 + 0 = 32 or higher (up to 32 + 1*2^2 + 1*2^1 + 1*2^0 = 39))
 *          - (iii) the location of the point will now be in the right half (x) and back half (y) and lower half (z)
 *            (cell number n = 1*2^5 + 1*2^4 + 0*2^3 = 32 + 16 + 0 = 48 or higher (up to n = 48 + 1*2^2 + 1*2^1 + 1*2^0 = 55) )
 *          - (iv) the location of the point will be in left half (x) and back half (y) and lower half (z)
 *            (cell number n = 0*2^5 + 1*2^4 + 0*2^3 = 0 + 16 + 0 = 16 or higher (up to 16 + 1*2^2 + 1*2^1 + 1*2^0 = 23))
 *
 *     The second pass through the loop
 *      - first statment sets the 2rd bit or does not (1*2^2 = 4 )
 *        - from examples, we can see that the second pass does exactly the same thing as the first pass;
 *          however, at the next refinement level.
 *          - setting the second bit to 1 will put the point on the right half w.r.t. the x-axis of
 *            the cell that we found in the previous pass (starting cell value to ending cell value) at the
 *            previous level
 *          - setting to zero will place the point on the left half of that cell
 *      - second statement sets the 1st bit or does not (1*2^1 = 2)
 *        - from examples, we can see that the second pass does exactly the same thing as the first pass;
 *          however, at the next refinement level.
 *          - setting the first bit to 1 will put the point on the back half w.r.t. the y-axis of
 *            the cell that we found in the previous pass (starting cell value to ending cell value) at the
 *            previous level
 *          - setting to zero will place the point on the front half of that cell
 *      - third statment sets the 0th bit or does not (1*2^0 = 1 )
 *        - from examples, we can see that the second pass does exactly the same thing as the first pass;
 *          however, at the next refinement level.
 *          - setting the zeroth bit to 1 will put the point on the upper half w.r.t. the z-axis of
 *            the cell that we found in the previous pass (starting cell value to ending cell value) at the
 *            previous level (adding 1 to any of base cells in that set will give index of cell above)
 *          - setting to zero will place the point on the lower half of that parent cell (does nothing still base cell)
 *
 */


int Point::getBoxIndex(unsigned int level)
{
  // here refinement level count starts on l = 0
  Util util;
  return util.interleave(std::floor(coord[0]*std::pow(2,level)),
		                 std::floor(coord[1]*std::pow(2,level)),
		                 std::floor(coord[2]*std::pow(2,level)),
		                 level);
}
