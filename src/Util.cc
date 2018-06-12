/**  Util.cc
 *
 *   Created on: Jun 20, 2016
 *       Author: Keith D Brauss
 */

#include <cmath>
#include <vector>

#include "Util.h"

/**
 * Header Interface for Class Util
 *
 *
class Util
{
  public:
    Util() {};
    int interleave(int x, int y, int z, int level);
    std::vector<double> uninterleave(int n, int L);
    int setbit(int n, int pos, int setto);
    int getbit(int n, int pos);

};
*
*/

/** The interleave function operation is explained in Point.cc */
// Interleave is used by class Point to return the cell
// (for a given refinement level) where the point is located
// In the getBoxIndex function of class Point, int x, int y, and int z
// are integers that are the floor of the product of the
// x, y, and z coordinates of the point and the base 2^{level}
// These values appear to indicate the increments (cell lengths)
// necessary to reach the lower left hand corner of the cell
// of interest (containing the point with desired coordinates)
// Specifically, int x, int y, and int z appear to correspond to the
// horizontal increments (x) to the right of the lower left hand corner
// of the domain, forward increments (y) into the screen from that
// location, and vertical increments (z)
// upward to reach the lower left hand corner of the cell of interest
int Util::interleave(int x, int y, int z, int level)
{
  if (x == 0 && y == 0 && z == 0)
    return 0;
  else
  {
    int ans = 0;
    for (int i=0; i<level; ++i)
    {
      ans = setbit(ans, (level-i)*3-1, getbit(x, level-i-1));
      ans = setbit(ans, (level-i)*3-2, getbit(y, level-i-1));
      ans = setbit(ans, (level-i)*3-3, getbit(z, level-i-1));
    }
    return ans;
  }
}

/**
 * Explanation of Uninterleave
 *
 * This function returns the location to the lower left corner of the cell whose index n
 * and level L are given.  The location is returned with respect to the lower-left-front
 * corner of the domain (located at coordinates (0,0,0), see the figure below)
 * in units of cell lengths.  For example, at L = 2 and cell index n = 3, the vector
 * [0, 1, 1] is returned.  The vector indicates that zero units (horizontal steps) in the
 * x-direction are required to reach the lower-front-left corner of the cell with index n = 3,
 * 1 cell-length from the lower-front-left corner of the domain in the y-direction
 * (into the screen - "to the back") is required to reach the lower-front-left corner
 * of the cell with index n = 3, and 1 cell-length in the z-direction (upward)
 * is also required to reach the lower-front-left corner of the cell with index n = 3.
 *
 *
 *
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
 * Our examples below match up with the interleave examples shown in Point.cc
 *
 *
 * Example:
 *
 * Consider the call to uninterleave(index n, level L) = uninterleave(3, 2)
 * The refinement level L = 2 and index n = 3 are passed into the function.
 * In the refinement level L = 2 above, we can see the location of the cell with
 * index n = 3.
 *
 * If the index of the cell is n == 0, then nothing needs to be done to reach the
 * location of the lower-left-front corner of this cell.  We know that the
 * vector returned will be [0,0,0].
 * If the index of the cell is not n == 0, then the first step initializes
 * xInt, yInt, and zInt to 0.
 * We then pass through the for loop shown.  Here (n,L) = (3,2)
 *
 *   for (int i=0; i<L; ++i)
 *   {
 *     xInt = setbit(xInt, L-i-1, getbit(n, (L-i-1)*3+2));
 *     xInt = setbit(yInt, L-i-1, getbit(n, (L-i-1)*3+1));
 *     yInt = setbit(zInt, L-i-1, getbit(n, (L-i-1)*3));
 *   }
 *
 * Passing through the for loop - first pass
 * i = 0:
 *   xInt = setbit(xInt, L-i-1, getbit(n, (L-i-1)*3+2));
 *   yInt = setbit(yInt, L-i-1, getbit(n, (L-i-1)*3+1));
 *   zInt = setbit(zInt, L-i-1, getbit(n, (L-i-1)*3));
 *   becomes
 *   xInt = setbit(0, 2-0-1, getbit(3, (2-0-1)*3+2));
 *   yInt = setbit(0, 2-0-1, getbit(3, (2-0-1)*3+1));
 *   zInt = setbit(0, 2-0-1, getbit(3, (2-0-1)*3));
 *   implies
 *   xInt = setbit(0, 1, getbit(3, (1)*3+2)); //
 *   yInt = setbit(0, 1, getbit(3, (1)*3+1)); //
 *   zInt = setbit(0, 1, getbit(3, (1)*3));   // 3 = 2^1 + 2^0 = 00000011
 *   implies
 *   xInt = setbit(0, 1, getbit(3, 5));     // getbit(000*0011,5) = 0 (get bit 6th from right)
 *   yInt = setbit(0, 1, getbit(3, 4));     // getbit(000*0011,4) = 0 (get bit 5th from right)
 *   zInt = setbit(0, 1, getbit(3, 3));     // getbit(0000*011,3) = 0 (get 4th bit from right)
 *   implies								// getbit(n,m) - getting number n's bit m (count of bits starts with 0)
 *   xInt = setbit(0, 1, 0);                // setbit(0,2,0) = 00000000 = 0
 *   yInt = setbit(0, 1, 0);                // setbit(0,2,0) = 00000000 = 0
 *   zInt = setbit(0, 1, 0);                // setbit(0,2,0) = 00000000 = 0
 *                                          // setbit(n,m,l) - setting number n's bit m to value l
 *
 * i = 1:
 *   xInt = setbit(xInt, L-i-1, getbit(n, (L-i-1)*3+2));
 *   yInt = setbit(yInt, L-i-1, getbit(n, (L-i-1)*3+1));
 *   zInt = setbit(zInt, L-i-1, getbit(n, (L-i-1)*3));
 *   becomes
 *   xInt = setbit(0, 2-1-1, getbit(3, (2-1-1)*3+2));
 *   yInt = setbit(0, 2-1-1, getbit(3, (2-1-1)*3+1));
 *   zInt = setbit(0, 2-1-1, getbit(3, (2-1-1)*3));
 *   implies
 *   xInt = setbit(0, 0, getbit(3, (0)*3+2));
 *   yInt = setbit(0, 0, getbit(3, (0)*3+1));
 *   zInt = setbit(0, 0, getbit(3, (0)*3)); // 3 = 2^1 + 2^0 = 00000011
 *   implies
 *   xInt = setbit(0, 0, getbit(3, 2));     // getbit(3,1) = 0
 *   yInt = setbit(0, 0, getbit(3, 1));     // getbit(3,1) = 1
 *   zInt = setbit(0, 0, getbit(3, 0));     // getbit(3,0) = 1
 *   implies
 *   xInt = setbit(0, 0, 0);                // setbit(0,0,0) = 00000000 = 0
 *   yInt = setbit(0, 0, 1);                // setbit(0,0,1) = 00000001 = 1
 *   zInt = setbit(0, 0, 1);                // setbit(0,0,1) = 00000001 = 1
 *
 *   the for loop ends
 *   and therefore
 *   xInt = 1
 *   yInt = 1
 *   Therefore, the vector [0,1,1] is returned by the function.
 *   This matches what we expected from our discussion above.
 *
 * Example:
 *
 * Consider the call to uninterleave(index n, level L) = uninterleave(57, 2)
 * The refinement level L = 2 and index n = 57 are passed into the function.
 * In the refinement level L = 2 shown above, we can see the location of the cell with
 * index n = 57.
 *
 * We then pass through the for loop shown.  Here (n,L) = (3,2)
 *
 *   for (int i=0; i<L; ++i)
 *   {
 *     xInt = setbit(xInt, L-i-1, getbit(n, (L-i-1)*3+2));
 *     yInt = setbit(yInt, L-i-1, getbit(n, (L-i-1)*3+1));
 *     zInt = setbit(zInt, L-i-1, getbit(n, (L-i-1)*3));
 *   }
 *
 * Passing through the for loop - first pass
 * i = 0:
 *   xInt = setbit(xInt, L-i-1, getbit(n, (L-i-1)*3+2));
 *   yInt = setbit(yInt, L-i-1, getbit(n, (L-i-1)*3+1));
 *   zInt = setbit(zInt, L-i-1, getbit(n, (L-i-1)*3));
 *   becomes
 *   xInt = setbit(0, 2-0-1, getbit(57, (2-0-1)*3+2));
 *   yInt = setbit(0, 2-0-1, getbit(57, (2-0-1)*3+1));
 *   zInt = setbit(0, 2-0-1, getbit(57, (2-0-1)*3));
 *   implies
 *   xInt = setbit(0, 1, getbit(57, (1)*3+2)); //
 *   yInt = setbit(0, 1, getbit(57, (1)*3+1)); //
 *   zInt = setbit(0, 1, getbit(57, (1)*3));   // 57 = 2^5 + 2^4 + 2^3 + 2^0 = 00111001
 *   implies
 *   xInt = setbit(0, 1, getbit(57, 5));     // getbit(00111001,5) = 1 (get bit 6th from right)
 *   yInt = setbit(0, 1, getbit(57, 4));     // getbit(00111001,4) = 1 (get bit 5th from right)
 *   zInt = setbit(0, 1, getbit(57, 3));     // getbit(00111001,3) = 1 (get 4th bit from right)
 *   implies								// getbit(n,m) - getting number n's bit m (count of bits starts with 0)
 *   xInt = setbit(0, 1, 1);                // setbit(0,1,1) = 00000010 = 2
 *   yInt = setbit(0, 1, 1);                // setbit(0,1,1) = 00000010 = 2
 *   zInt = setbit(0, 1, 1);                // setbit(0,1,1) = 00000010 = 2
 *                                          // setbit(n,m,l) - setting number n's bit m to value l
 *
 * i = 1:
 *   xInt = setbit(xInt, L-i-1, getbit(n, (L-i-1)*3+2));
 *   yInt = setbit(yInt, L-i-1, getbit(n, (L-i-1)*3+1));
 *   zInt = setbit(zInt, L-i-1, getbit(n, (L-i-1)*3));
 *   becomes
 *   xInt = setbit(2, 2-1-1, getbit(57, (2-1-1)*3+2));
 *   yInt = setbit(2, 2-1-1, getbit(57, (2-1-1)*3+1));
 *   zInt = setbit(2, 2-1-1, getbit(57, (2-1-1)*3));
 *   implies
 *   xInt = setbit(2, 0, getbit(57, (0)*3+2));
 *   yInt = setbit(2, 0, getbit(57, (0)*3+1));
 *   zInt = setbit(2, 0, getbit(57, (0)*3)); // 57 = 2^5 + 2^4 + 2^3 + 2^0 = 00111001
 *   implies
 *   xInt = setbit(2, 0, getbit(57, 2));     // getbit(00111001,2) = 0
 *   yInt = setbit(2, 0, getbit(57, 1));     // getbit(00111001,1) = 0
 *   zInt = setbit(2, 0, getbit(57, 0));     // getbit(00111001,0) = 1
 *   implies
 *   xInt = setbit(2, 0, 0);                // setbit(2,0,0) = 00000010 = 2
 *   yInt = setbit(2, 0, 0);                // setbit(2,0,0) = 00000010 = 2
 *   zInt = setbit(2, 0, 1);                // setbit(2,0,1) = 00000011 = 3
 *
 *
 *   the for loop ends
 *   and we have the result
 *   xInt = 2
 *   yInt = 2
 *   zInt = 3
 *   Therefore, the vector [2,2,3] is returned by the function.
 *   This matches what we expected from our discussion above.
 *   If from the lower-left-front corner (0,0,0) we move two cell
 *   lengths in the x-direction, 2 cell-lengths in the y-direction,
 *   and 3 cell-lengths in the z-direction, then we reach the
 *   lower-left-front corner of the cell having index n=57 at level
 *   l = 2.
 *
 *
 * In the for loop below, setbit turns on/off (sets to 1 or 0)
 * the L-i-1 bit in the integer
 * - xInt if getbit(n, (L-i-1)*2+1) is 1/0 (on or off)
 * - yInt if getbit(n, (L-i-1)*2) is 1/0 (on or off)
 *
 * The largest bit positions that can be set is L-i-1 = 2 - 0 - 1 = 1
 * The largest that xInt, yInt, and zInt can be after passing through the
 * for loop occurs if all the bits in positions 1 and 0 are turned on (set to 1)
 * This results in the value 2^1 + 2^0 = 2 + 1 = 3
 * The vector returned from the function would be [xInt,yInt,zInt] = [3,3,3]
 * This is the location in cell-lengths from the lower-front-left corner of the domain
 * to the lower-front-left corner of the cell in the upper-back-right corner
 * of the domain - the cell with the last index of the refinement level L = 2.
 * We cane see from this that the method is able to locate the lower-front-left corners
 * of all the cells at this refinement level
 *
 * To obtain a particular point in the cell we can multiply the values of xInt, yInt,
 * and zInt by the length of the cells in the refinement.
 * For the example with L = 2 and (xInt,yInt,zInt = (2,2,3)
 * - The cell lengths are 0.25 (see figure above).
 * - The lower-left corner of the cell having index n = 57 of the refinement
 *   L = 2 has coordinates (x,y) = (0.50,0.50,0.75)
 * - This can be obtained by multiplying (xInt,yInt,zInt) = (2,2,3) by cell-length = 0.25
 *   i.e., (2*0.25,2*0.25,3*0.25) = (0.50,0.50,0.75)
 * - Using this location, we can map to any position in cell n = 57
 *
 * Note: would make sense that return vector is std::vector<int>
 */

std::vector<double> Util::uninterleave(int n, int L)
{
  std::vector<double> llf_box_corner(3);  // llf - lower left front
  if (n==0)
  {
    llf_box_corner[0] = 0.0; llf_box_corner[1] = 0.0;
    llf_box_corner[2] = 0.0;
  }
  else
  {
    int xInt = 0;
    int yInt = 0;
    int zInt = 0;
    for (int i=0; i<L; ++i)
    {
      xInt = setbit(xInt, L-i-1, getbit(n, (L-i-1)*3+2));
      yInt = setbit(yInt, L-i-1, getbit(n, (L-i-1)*3+1));
      zInt = setbit(zInt, L-i-1, getbit(n, (L-i-1)*3));
    }
    llf_box_corner[0] = xInt; llf_box_corner[1] = yInt;
    llf_box_corner[2] = zInt;
  }
  return llf_box_corner;
}




/** Explanation of setbit and getbit member functions
 *
 *  Recall the formula for going from binary to decimal
 *  Ex: byte 10010110   < - - - - - - - bit zero  : 0 x 2^0
 *           |   |  |                   bit one   : 1 x 2^1
 *           |   |  bit 0               bit two   : 1 x 2^2
 *           |   bit 3                  bit three : 0 x 2^3
 *           bit 7                      bit four  : 1 x 2^4
 *                                      bit five  : 0 x 2^5
 *                                      bit six   : 0 x 2^6
 *                                  +   bit seven : 1 x 2^7
 *                                -----------------------------
 *                                      2^1 + 2^2 + 2^4 + 2^7 = 2 + 4 + 16 + 128
 *                                                            = 22 + 128
 *                                                            = 150
 *
 *  1u is an unsigned value (int) with the single bit 0 set (00000001 - 8 bits w/ bit 0 set = 1)
 *  here we are using 1 instead of 1u and ints instead of unsigned ints
 *  now 1u << 0 means shift 1=00000001 to the left zero units
 *      1u << 1 means shift 1=00000001 to the left one unit = 00000010 = 1 x 2^1 = 2
 *      1u << 2 means shift 1=00000001 to the left two units = 00000100 = 1 x 2^2 = 4
 *
 *  Therefore, for the setbit function below
 *  In the if statement
 *    n | (1u<<pos) means "n union (00000001 shifted to the left pos units)"
 *    and results in a union of n with 1u<<pos (a union of all the bits that are on (ones))
 *    ex:         10100001 | 1u<<3 = 10100001 | 00001000 = 10101001
 *    in decimal:      161 | 8     =                     = 128 + 32 + 8 + 1 = 169
 *  And the code below in the else statement
 *    n & ~(1u<<pos) = n and (00000001 shifted to the left pos units)
 *    results in an intersection of n with 1u<<pos (an intersection of all bits that are on (ones))
 *    ex:         10101001 & ~(1u<<3) = 10101001 & ~(00001000) = 10101001 & 11110111 = 10100001
 */

/** turns on (set to 1) bit pos of int n if setto == 1 else turns off (set to 0) bit pos of int n */
int Util::setbit(int n, int pos, int setto)
{
  if (setto == 1)
    return n | (1<<pos);  // turning on bit 'pos' in number n (setting bit pos equal to 1)
  else
    return n & ~(1<<pos); // turning off bit 'pos' in number n (setting bit pos equal to 0)
}

/** checks whether a bit at position 'pos' of an int 'n' is on or off and returns 1 if on and 0 if off */
int Util::getbit(int n, int pos)
{
  if ((n & (1 << pos)) != 0)  // checking whether bit pos of int n is on or off (1 or 0)
    return 1;                 // returning 1 if bit pos of int n is 1 (on)
  else
    return 0;                 // returning 0 if bit pos of int n is 0 (off)
}
