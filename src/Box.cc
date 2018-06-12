/** Box.cc
 *
 *  Created on: Jun 21, 2016
 *      Author: Keith D Brauss
 */

#include <complex>
#include <vector>
#include <iostream>
#include <string>

#include "Box.h"
#include "Util.h"

using namespace std;


/**
 * The Problem Type
 *
 * We want to apply the FMM to a problem that has the following properties
 * - The domain is the unit square [0,1]x[0,1]
 * - The refinement is quadtree where one cell is a square and is split into
 *   4 squares of equal size
 * - The source points x and the target points y will be the same
 *   (the source points x will also be the target points y).
 *   The cells at each refinement level are squares
 * - The number of source (or target) particles in each cell at the lowest refinement level (l = L) are known
 *   - We will set this number of particles per cell at l = L to nParticlesPerCell = 4
 * - The particles' locations are determined by formulas that use the coordinates of the corners of their cells
 *   - Example: - Suppose that the corners are ll=(0.0,0.0), lr=(1.0,0.0), ul=(0.0,1.0), and ur=(1.0,1.0),
 *                where ll stands for lower-left, lr stands for lower-right, ul stands for upper left,
 *                ur stands for upper right
 *              - We wish the four particles to be inside the cell and 1/4 away from each corner w.r.t. the x and y axis
 *              - Therefore, assuming that we are working with squares for cells, we can measure the length
 *                of one side of the cell \f$s = \sqrt{ (ll[0] - lr[0])^2 + (ll[1] - lr[1])^2 }\f$.
 *              - We want 1/4 of this distance.  Let d = 0.25s
 *              - Lets refer to the particles as points in the source vector x, say x1, x2, x3, x4
 *              - Then let x1[0] = ll[0] + d;  x1[1] = ll[1] + d;
 *                         x2[0] = lr[0] - d;  x2[1] = lr[1] + d;
 *                         x3[0] = ul[0] + d;  x3[1] = ul[1] - d;
 *                         x4[0] = ur[0] - d;  xr[1] = ur[1] - d;
 *
 *                ul(0,1)  *------------------------* ur(1,1)
 *                         |                        |
 *                         |                        |
 *                         |      *          *      |
 *                         |      x3         x4     |
 *                         |                        |
 *                         |                        |
 *                         |                        |
 *                         |      *          *      |
 *                         |      x1         x2     |
 *                         |                        |
 *                ll(0,0)  *------------------------* lr(1,0)
 *
 *                         |------------s-----------|
 *                         |--d---|          |---d--|
 *
 *
 * The series representing the mother function is truncated at the index value p.
 * - Therefore, the coefficient arrays: c, dtilde, and d will have a size of p+1
 *   (since the series index starts counting at zero).
 * For each refinement level l, there are 4^l cells
 * - for l = 1, there are 4^1 = 4 cells
 * - for l = 2, there are 4^2 = 16 cells
 * - for l = 3, there are 4^3 = 64 cells
 * - ...
 * - for l = L, there are 4^L cells
 * Therefore, since the number of particles per cell at the lowest refinement level is set,
 * the total number of points in the domain is known and equals nParticlesPerCell*pow(4,L)
 * - Example: If the lowest refinement level is L = 3 and nParticlesPerCell = 4,
 *            then the total number of particles is 4*4^3 = 4*64 = 256 particles in the domain
 *
 */

/**
 * Header Interface for Class Point
 *
class Box
{
  public:
    int DEFAULT_P = 16;
    int DEFAULT_LEVEL=3;
    int DEFAULT_INDEX=0;

    int level;                             // refinement level of box (start count on 0)
    int index;                             // cell index of box (or cell)
    int p;                                 // p is the index at which the series are truncated
    int order_of_approximation             // the size of the coefficient matrices and the
                                           // highest order derivatives used in series approximation
                                           // also known as abs_alpha

//    unsigned int nParticlesPerCell = 4;    // points in each cell at lowest refinement level (l = L)
//    int totalParticles;                     // total particles in domain (unit square)


    bool empty;                            // is box empty


    std::vector<std::vector<std::complex<double> > > c;
    std::vector<double> dtilde;
    std::vector<double> d;

    //
    // The number of source points (and target points) in a box will depend
    // on the refinement level.  The lowest level l = L will have 4 particles
    // per box.  As we work up through the levels to level l = 0, the number
    // of particles will increase by a multiple of 4 when going from one level
    // up to the next (four children for each parent)

    std::vector<Point> x;  // source points
    std::vector<Point> y;  // target points

    // Creates a new instance of Node
    Box();
    Box(int level, int index, int p, int order_of_approximation);

    int       getLevel() { return level; };
    void      printLevel() { std::cout << "Box level is " << level << '\n'; };
    void      setLevel(int i) { this->level = i; };
    int       getIndex() { return index; };
    void      printIndex() { std::cout << "Box index is " << index << '\n'; };
    void      setIndex(int i) { this->index = i; };
    Point     getCenter();
    double    getSize() { return std::pow(2.0, -level); };

    void      setP(int p) { this->p = p; this->c.resize(p); this->d.resize(p); this->dtilde.resize(p);};
    bool      isEmpty() { return empty; };

    std::vector<double>                getC() {return c; };
    void                               addToC(std::vector<std::vector<std::complex<double> > > &increment);
    std::string                        cToString();
    void                               printC();

    std::vector<double>                getD() {return d; };
    void                               addToD(std::vector<double> &increment);
    std::string                        dToString();
    void                               printD();

    void                               addToDtilde(const std::vector<double> &increment);
    std::string                        dtildeToString();
    void                               printDtilde();


    void               addX(Point &p);
    int                getSizeX() { return this->x.size(); };
    std::vector<Point> getX() { return this->x; };
    void               printSizeX() { std::cout << "Box sizeX is " << x.size() << "\n"; };

    void               addY(Point &p);
    int                getSizeY() { return this->y.size(); };
    std::vector<Point> getY() { return this->y; };
    void               printSizeY() { std::cout << "Box sizeY is " << y.size() << "\n"; };

    std::string        toString();
    int                getParentIndex();
    void               getNeighborsIndex        (std::vector<int> &neighbor_indexes);
    void               getParentsNeighborsIndex (std::vector<int> &parents_neighbor_indexes);

    void               getNeighborsE4Index      (std::vector<int> &neighborE4_indexes);

    void               getChildrenIndex(std::vector<int> &children_indexes);

  private:

    void               getChildrenIndexOfBox(int levelOfBox, int indexOfBox, std::vector<int> &children_indexes_of_box);

};
*/


Box::Box()
   :
   level(DEFAULT_LEVEL),
   index(DEFAULT_INDEX),
   p(DEFAULT_P),
   order_of_approximation(DEFAULT_ORDER_OF_APPROXIMATION),
   empty(true),
   c(order_of_approximation+1,std::vector<std::complex<double> > (order_of_approximation+1)),
   dtilde(order_of_approximation+1,std::vector<std::complex<double> > (order_of_approximation+1)),
   d(order_of_approximation+1,std::vector<std::complex<double> > (order_of_approximation+1))
//   dtilde_grad(order_of_approximation,std::vector<std::complex<double> >(order_of_approximation * 3)),
//   d_grad(order_of_approximation,std::vector<std::complex<double> >(order_of_approximation * 3))
{

  // initializing coefficients
  for (int j=0; j<order_of_approximation+1; ++j)
    for(int i=0; i<order_of_approximation+1; ++i)
    {
      c[i][j] = 0.0;
      dtilde[i][j] = 0.0;
      d[i][j] = 0.0;
//      if (i < order_of_approximation)
//        if (j < order_of_approximation)
//          for (int k=0; k<3; ++k)
//          {
//            dtilde_grad[i][j][k] = 0.0;
//            d_grad[i][j][k] = 0.0;
//          }
    }
}


Box::Box(int level, int index, int p, int order_of_approximation)
   :
   level(level),
   index(index),
   p(p),
   order_of_approximation(order_of_approximation),
   empty(true),
   c(order_of_approximation+1,std::vector<std::complex<double> > (order_of_approximation+1)),
//   c(p+1),
   dtilde(order_of_approximation+1,std::vector<std::complex<double> > (order_of_approximation+1)),
   d(order_of_approximation+1,std::vector<std::complex<double> > (order_of_approximation+1))
//   dtilde_grad(order_of_approximation,std::vector<std::complex<double> >(order_of_approximation * 3)),
//   d_grad(order_of_approximation,std::vector<std::complex<double> >(order_of_approximation * 3))
{

  // initializing coefficients c, dtilde and d
  for (int j=0; j<order_of_approximation+1; ++j)
    for(int i=0; i<order_of_approximation+1; ++i)
    {
      c[i][j] = 0.0;
      dtilde[i][j] = 0.0;
      d[i][j] = 0.0;
//      if (i < order_of_approximation)
//        if (j < order_of_approximation)
//          for (int k=0; k<3; ++k)
//          {
//            dtilde_grad[i][j][k] = 0.0;
//            d_grad[i][j][k] = 0.0;
//          }
    }

}

Point Box::getCenter()
{
  Util util;
  std::vector<double> llf_corner = util.uninterleave(this->index, this->level);
  std::vector<double> middle = {0.5,0.5,0.5};
  llf_corner[0] += middle[0];  llf_corner[1] += middle[1];  llf_corner[2] += middle[2];
  llf_corner[0] *= getSize();  llf_corner[1] *= getSize();  llf_corner[2] *= getSize();
  Point center(llf_corner);

  return center;
}


void Box::addToC(std::vector<std::vector<std::complex<double> > > &increment)
{
//std::cout << "size of c is " << c[0].size() << std::endl;
  for (unsigned int i=0; i<(unsigned int)order_of_approximation+1; ++i)
    for (unsigned int j=0; j<(unsigned int)order_of_approximation+1; ++j)
    c[j][i] = c[j][i] + increment[j][i];
}


/*
std::string Box::cToString()
{
  std::string ansc = "        C: ";
  for (unsigned int i=0; i<c.size(); ++i)
  {
    ansc+= "(" + std::to_string(c[i].real()) + " " + std::to_string(c[i].imag()) + ") ";
  }
  return ansc + "\n";
}

void Box::printC()
{
  std::cout << cToString() << "\n";
}
*/

void Box::addToD(std::vector<std::vector<std::complex<double> > > &increment)
{
  for (unsigned int i=0; i<(unsigned int)order_of_approximation+1; ++i)
    for (unsigned int j=0; j<(unsigned int)order_of_approximation+1; ++j)
      d[j][i] = d[j][i] + increment[j][i];
}

/*
std::string Box::dToString()
{
  std::string ansd = "        D: ";
  for (unsigned int i=0; i<d.size(); ++i)
  {
    ansd+= "(" + std::to_string(d[i].real()) + " " + std::to_string(d[i].imag()) + ") ";
  }
  return ansd + "\n";
}

void Box::printD()
{
  std::cout << dToString() << "\n";
}
*/

void Box::addToDtilde(const std::vector<std::vector<std::complex<double> > > &increment)
{
  for (unsigned int i=0; i<(unsigned int)order_of_approximation+1; ++i)
    for (unsigned int j=0; j<(unsigned int)order_of_approximation+1; ++j)
      dtilde[j][i] = dtilde[j][i] + increment[j][i];
}

/*
std::string Box::dtildeToString()
{
  std::string ansdt = "        Dtilde: ";
  for (unsigned int i=0; i<dtilde.size(); ++i)
  {
    ansdt+= "(" + std::to_string(dtilde[i].real()) + " " + std::to_string(dtilde[i].imag()) + ") ";
  }
  return ansdt + "\n";
}

void Box::printDtilde()
{
  std::cout << dtildeToString() << "\n";
}
*/

void Box::addX(Point &p)
{
  this->x.push_back(p);
}

void Box::addY(Point &p)
{
  this->y.push_back(p);
}

/*
std::string Box::toString()
{
  std::string ans = "box (l = " + std::to_string(level) + ", n = " + std::to_string(index) + ") \n";
  std::string ansc = "        C: ";
  std::string ansdt = "   Dtilde: ";
  std::string ansd = "        D: ";
  for (unsigned int i=0; i<c.size(); ++i)
  {
    ansc+= "(" + std::to_string(c[i].real()) + " " + std::to_string(c[i].imag()) + ") ";
    ansdt+="(" + std::to_string(dtilde[i].real()) + " " + std::to_string(dtilde[i].imag()) + ") ";
    ansd+= "(" + std::to_string(d[i].real()) + " " + std::to_string(d[i].imag()) + ") ";
  }
  ans += ansc + "\n" + ansdt + "\n" + ansd + "\n";
  return ans;
}
*/

/**
 * Explanation of getParentIndex
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
 *  1.00 |______ ______ ______ ______      |      |      |      |      |     |  17  |* 21  |  49  |  53  |     |______|______|______|______|
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
 *
 *     This is level l=2 refinement
 *
 *
 * going to parent level (level - 1)
 * shifting the index of this box to right 2 bits (syntax is index >> 2)
 *
 * Example: index n = 8 and level l = 2
 *
 *          At level l = 2 have 2^(l*d) = 2^(2*3) = 2^6 = 64 cells
 *          At parent level l = 2-1 = 1 have 2^(1*3) = 2^3 = 8 cells
 *
 *          In binary n = 8 is  00001000
 *          Shifting 00001000 to right 3 bits gives 0000001 = 2^0 = 1
 *          We can see above that cell n=8 at l=2 has parent cell n=1 at l=1
 *
 * Example: index n = 29 and level l = 2
 *
 *          At level l = 2 have 2^(l*d) = 2^(2*3) = 2^6 = 64 cells
 *          At parent level l = 2-1 = 1 have 2^(1*3) = 2^3 = 8 cells
 *
 *          In binary n = 29 is  00011101
 *          Shifting 00011101 to right 3 bits gives 00000011 = 2^1 + 2^0 = 2 + 1 = 3
 *          We can see above that cell n=29 at l=2 has parent cell n=3 at l=1
 *
 */
int Box::getParentIndex()
{
  return index >> 3;
}

/**
 * Explanation of getNeighborsIndex
 *
 * The nested for loop
 *  for (int i=-1; i<=1; ++i)
      for (int j=-1; j<=1; ++j)
        for (int k=-1; k<=1; ++k)
          if (  (i!=0||j!=0) && x+i>=0 && x+i<Math.pow(2, level) && y+j>=0 && y+j<Math.pow(2, level))
            ans.addElement( t.getBox(level, Util.interleave(x+i, y+j, level)) );
 *
 * Looking at the nested for loop, we see that there are 3*3*3 = 27 cases and that the case
 * (i,j,k) = (0,0,0) is removed as a possibility by the if statement (the box we are working with is removed
 * from consideration as a possible neighbor).
 * Therefore, there are 26 cases and these match up with the 8 possible neighbors.
 * Again, the case (0,0,0) matches up with the center cell for which we are trying to find the nearest neighbors.
 * We can see that -1, 0, and 1 for i corresponds to left, center (zero) and right one cell length (one unit)
 * in the x direction
 * and -1, 0, and 1 for j corresponds to backward, center (zero) and forward one unit in the y direction from the cell of interest.
 *
 * Example: Let tmp =  (0.3125, 0.6875, 0.375)
 *          Then uninterleave function of class Util returns the location of the lower-left-front corner
 *          of the cell where the point tmp is located.  The location is given in terms of cell lengths from the
 *          lower-left-front corner of the domain (0,0,0).  For this example the point (0.3125, 0.6875, 0.375) is located by
 *          the lower-left-front corner of cell n = 21. The uninterleave function returns x = 1 and y = 2 and z = 1
 *          implying right two cell lengths from the lower-left-front corner of the domain, forward three cell lengths
 *          (parallel to the y-axis), and up one cell length (parallel to z-axis) to reach the lower left corner of cell n = 21.
 *
 *          Let level = 2
 *          For the cell n = 21 containing the point (0.3125, 0.6875, 0.375) and having x = 1 and y = 2 and z = 1
 *          The neighbors are n = 17, (<--x-->) 49, 7, (<--y-->) 23, 20, (<--z-->) 28,		6 neighbors in axial directions
 *          leftlower: 2, 16, 18,   leftmiddle: 3, 19,       leftupper: 10, 24, 26,         8 on left (into y axis)
 *          middlelower: 6, 22,                              middleupper: 14, 30            4 in middle (into y axis)
 *          rightlower: 34, 48, 50, rightmiddle: 35, 51,     rightupper: 42, 56, 58         8 on right (into y axis)
 *																					   sum: 26 neighbors
 * i = -1
 *   j = -1
 *     k = -1
 *       ((i=)-1!=0 || (j=)-1!=0 || (k=)-1!=0) && ((x+i=1+-1=)0>=0) && ((x+i=1+-1=)0<4(=2^3))
 *                T                             T                    T
 *                                             && ((y+j=2+-1=)1>=0) && ((y+j=2+-1=)1<4(=2^3))
 *                                              T                    T         = False
 *                                             && ((z+k=1+-1=)0>=0) && ((z+k=1+-1=)0<4(=2^3))
 *                                              T                    T                 = True
 * (True) There is a neighbor diagonally off the lower-left-front corner of the cell n = 21
 *  This neighbor is cell n = 2
 *  Note: At level l = 2, four cell-lengths from the lower-left-front corner in either of the 3 axial directions
 *        will take you outside the domain.  Therefore, the lower-left-front corner of any cell will be located
 *        at most 3 cell lengths away from the lower-left-front corner of the domain (origin). Hence, if x+i,
 *        y+j, or z+k is greater than or equal to 4 ( not < 2^2 = 4 ) we are at a location on the edge of the
 *        domain (on the boundary).  And there is no neighboring cell at that location.
 *
 * i = -1
 *   j = -1
 *     k = 0
 *       ((i=)-1!=0 || (j=)-1!=0 || (k=)0!=0)  && ((x+i=1+-1=)0>=0) && ((x+i=1+-1=)0<8(=2^3))
 *                T                             T                    T
 *                                             && ((y+j=2+-1=)1>=0) && ((y+j=2+-1=)1<8(=2^3))
 *                                              T                    T
 *                                             && ((z+k=1+0=)1>=0) && ((z+k=1+0=)1<8(=2^3))
 *                                              T                    T                 = True
 * (True) There is a neighbor diagonally off the front-left edge of the cell n = 21
 *  This neighbor is cell n = 3
 * Recall that ans is a vector of boxes
 *    ans.addElement( t.getBox(level, Util.interleave(x+i, y+j, z+k, level)) );
 *    ans.addElement( t.getBox(3, Util.interleave(1+-1, 2+-1, 1+0, 2)) );
 *    ans.addElement( t.getBox(3, Util.interleave(0, 1, 1, 2)) );
 *
 * The interleave function returns the index n = 17.
 * Indeed, Util.interleave(0,1,1,2) returns the box that is 0 increments right horizontally (x-direction) and
 * 1 increments forward (y-direction), and 1 cell-length (increment) upward vertically (z-direction) from
 * the lower-left-front corner of the domain.  The index for this cell is n = 17.
 *
 * i = -1
 *   j = -1
 *     k = 1
 *       ((i=)-1!=0 || (j=)-1!=0 || (k=)1!=0)  && ((x+i=1+-1=)0>=0) && ((x+i=1+-1=)0<8(=2^3))
 *                T                             T                    T
 *                                             && ((y+j=2+-1=)1>=0) && ((y+j=2+-1=)1<8(=2^3))
 *                                              T                    T
 *                                             && ((z+k=1+1=)2>=0) && ((z+k=1+1=)2<8(=2^3))
 *                                              T                    T                 = True
 * (True) There is a neighbor diagonally off the upper-left-front corner of the cell n = 21
 *  This neighbor is cell n = 10
 * Recall that ans is a vector of boxes
 *    ans.addElement( t.getBox(level, Util.interleave(x+i, y+j, z+k, level)) );
 *    ans.addElement( t.getBox(3, Util.interleave(1+-1, 2+-1, 1+1, 2)) );
 *    ans.addElement( t.getBox(3, Util.interleave(0, 1, 2, 2)) );
 *
 * The interleave function returns the index n = 10.
 * Indeed, Util.interleave(0,1,2,2) returns the box that is zero increments (cell-lengths) horizontally to the right
 * (in x-direction), and one increment forward (in y-direction), and two increments vertically upward (in z-direction)
 * from the lower-left-front corner of the domain (0,0,0).  The index for this cell is n = 1.
 *
 * Next, i = -1, j = 0 and we iterate over k=-1,0,1
 *   This covers lower to upper cells on the left-middle side of the cell of interest (at the center of the 27 cell neighbors)
 * Next, i = -1, j = 1 and we iterate over k=-1,0,1
 *   This covers lower to upper (vertical z-direction for k=-1,0,1) cells on the left-back edge of the cell of interest
 * Next, i = 0, and we iterate over j =-1,0,1 and  k=-1,0,1
 *   This covers the neighbors located at the middle-front, middle-middle and middle-back of the cell of interest
 * Last, i = 1, and we iterate over j =-1,0,1 and  k=-1,0,1
 *   This covers the neighbors located at the right-front, right-middle and right-back of the cell of interest
 *
 */


void Box::getNeighborsIndex(std::vector<int> &neighbor_indexes)
{
  // Finding x,y increments (cell lengths) from the lower left corner of
  // the domain to the lower left corner of this Box (cell) using the uninterleave
  // function from class Util and placing these increments in tmp in terms of
  // real term = horizontal units from lower left corner of domain and
  // imag term = vertical units from lower left corner of the domain
  neighbor_indexes.resize(0);
  Util util;

  std::vector<double> tmp = util.uninterleave(index,level);
  int x = (int)tmp[0];
  int y = (int)tmp[1];
  int z = (int)tmp[2];
  for (int i=-1; i<=1; ++i)
    for (int j=-1; j<=1; ++j)
      for (int k=-1; k<=1; ++k)
        if (  (i!=0||j!=0||k!=0) && x+i>=0 && x+i<std::pow(2, level)
                                 && y+j>=0 && y+j<std::pow(2, level)
                                 && z+k>=0 && z+k<std::pow(2, level)
           )
        neighbor_indexes.push_back(util.interleave(x+i, y+j, z+k, level));
}

 // see getNeighborsIndex above for explanation
 void Box::getParentsNeighborsIndex(std::vector<int> &parents_neighbor_indexes)
 {
   // Finding x,y, and z increments (cell lengths) from the lower-left-front corner of
   // the domain to the lower-left-front corner of this Box (cell) using the uninterleave
   // function from class Util and placing these increments in tmp in terms of
   // tmp[0] term = horizontal units to right from lower-left-front corner of domain and
   // tmp[1] term = forward units (into y-direction) from lower-left-front corner of the domain
   // tmp[2] term = vertically upward (z-direction)  from lower-left-front corner of the domain
   parents_neighbor_indexes.resize(0);
   Util util;

   int parent_index;
   parent_index = this->getParentIndex();

   std::vector<double> tmp = util.uninterleave(parent_index, level-1);
   int x = (int)tmp[0];
   int y = (int)tmp[1];
   int z = (int)tmp[2];
   for (int i=-1; i<=1; ++i)
     for (int j=-1; j<=1; ++j)
       for (int k=-1; k<=1; ++k)
         if (  (i!=0||j!=0||k!=0) && x+i>=0 && x+i<std::pow(2, level-1)
                                  && y+j>=0 && y+j<std::pow(2, level-1)
                                  && z+k>=0 && z+k<std::pow(2, level-1))
           parents_neighbor_indexes.push_back(util.interleave(x+i, y+j, z+k, level-1));
 }


/**
 * Explanation of getNeighborsE4Index
 *
 * The member function gets the interaction list for this box
 * The interaction list or set E_4 is made up of this box's parent's
 * nearest neighbor's children (peers of this box).
 * However, the list does not include the nearest neighbors of this box.
 * We obtain the interaction list E_4 by
 * (1) getting the indices of the near neighbors of this box
 *     - we will check to make sure that these cell indexes are
 *     - not in the interaction list set E_4 that will be returned
 * (2) getting the indexes of this box's parent's nearest neighbors
 *     - we will then get the children of each of these parent boxes
 *     - the children then form the interaction list - minus the near neighbors of this box
 * (3) getting the set of indices for the children of this box's parent's nearest neighbors
 * (4) comparing and removing the indices of this box's near neighbors from the list (set)
 *     of step (3)
 */


void Box::getNeighborsE4Index(std::vector<int> &neighborE4_indexes)
{
  neighborE4_indexes.resize(0);

  // getting the indexes of the near neighbors of this box
  std::vector<int> neighbor_indexes;
  neighbor_indexes.resize(0);
  this->getNeighborsIndex(neighbor_indexes);

  // getting the (parent) indexes of the near neighbors of the parent of this box
  std::vector<int> parents_neighbor_indexes;
  parents_neighbor_indexes.resize(0);
  this->getParentsNeighborsIndex(parents_neighbor_indexes);

//for (unsigned int i=0; i<neighbor_indexes.size(); ++i)
//    std::cout << "neighbor_indexes[" << i << "] = " << neighbor_indexes[i] << "\n";
//for (unsigned int j=0; j<parents_neighbor_indexes.size(); ++j)
//    std::cout << "parents_neighbor_indexes[" << j << "] = " << parents_neighbor_indexes[j] << "\n";


  // getting the (peer) indexes of the children of the neighbors of the parent
  std::vector<int> parent_neighbor_children_indexes;
  std::vector<int> box_children_indexes;
  parent_neighbor_children_indexes.resize(0);
  for (unsigned int i=0; i<parents_neighbor_indexes.size(); ++i)
  {
	box_children_indexes.resize(0);
    getChildrenIndexOfBox(level-1, parents_neighbor_indexes[i],
    		              box_children_indexes);
    for (unsigned int m=0; m<box_children_indexes.size(); ++m)
    parent_neighbor_children_indexes.push_back(box_children_indexes[m]);
  }

  // comparing and removing the indexes from the parent_neighbor_children that
  // are also near neighbors of this box
  bool flag_not_in_list = false;  // true that parent_neighbor_children_indexes[i] is not
                                  // in the interaction list (is a near neighbor) of this box
  for (unsigned int i=0; i<parent_neighbor_children_indexes.size(); ++i)
  {
	flag_not_in_list =  false;

	for (unsigned int j=0; j<neighbor_indexes.size(); ++j)
      if (parent_neighbor_children_indexes[i] == neighbor_indexes[j])
        flag_not_in_list = true;

    if (flag_not_in_list == true)
    {
      // not adding the index 'parent_neighbor_children_indexes[i]'
      // to the interaction list neighborE4_indexes
    }
    else
    {
      neighborE4_indexes.push_back(parent_neighbor_children_indexes[i]);
    }

  }

}


/** Explanation of getChildrenIndex
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
 *  1.00 |______ ______ ______ ______      |      |      |      |      |     |  17  |* 21  |  49  |  53  |     |______|______|______|______|
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
 *
 *     This is level l=2 refinement
 *
 *  For 3D octree structure each box has 8 children and therefore member
 *  function getChildrenIndex returns 8 indices
 *
 *  In the for loop
 *    for (int i=0; i<8; ++i)
 *      children_indexes.push_back((index<<3)+i);
 *
 *  the command
 *      (index<<3)+i;
 *
 *  shifts the index to the left three bits
 *
 *  Example: l = 1 and index n = 7
 *
 *           n = 7 = 2^2 + 2^1 + 2^0 = 00000111
 *
 *           Shifting n = 9 to the left two bits (index << 3)
 *           results in 00111000 = 2^3 + 2^4 + 2^5 = 8 + 16 + 32 = 56
 *
 *           From the diagram above, this is the index of the lower-left-front
 *           child cell at level l = 2 with index n = 56 of the parent cell l = 1 and n = 7
 *
 *           The other 7 children are obtained by adding 1, 2, 3, 4, 5, 6, and 7 to the
 *           index of this child cell (56 + 1 = 57 (ufl - upper front left), 56 + 2 = 58 (lbl),
 *           56 + 3 = 59 (ubl), 56 + 4 = 60 (lfr), 56 + 5 = 61 (ufr), 56 + 6 = 62 (lbr - lower back right),
 *           56 + 7 = 63 (ubr)
 *
 *
 *
 */

void Box::getChildrenIndex(std::vector<int> &children_indexes)
{
  children_indexes.resize(0);
  for (int i=0; i<8; ++i)
    children_indexes.push_back((this->index<<3)+i);
}

void Box::getChildrenIndexOfBox(int levelOfBox, int indexOfBox, std::vector<int> &children_indexes_of_box)
{
  // need to assert that levelOfBox is between 2 and 8
  // if refinement level is greater than 8, then bitwise operations will be incorrect (8 bits)
  // if refinement level is less than 2, then FMM does not apply (series approximation convergence not guaranteed)
  children_indexes_of_box.resize(0);
  for (int i=0; i<8; ++i)
    children_indexes_of_box.push_back((indexOfBox<<3)+i);
}

