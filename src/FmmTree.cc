/*
 * FmmTree.cc
 *
 *  Created on: Jun 21, 2016
 *      Author: Keith D Brauss
 */

#include <complex>
#include <vector>
#include <iostream>
#include <limits>
#include <cmath>
#include <cassert>

#include "FmmTree.h"
#include "Box.h"
#include "Point.h"


using namespace std;

/**
 * Header Interface for Class Point
 *
class FmmTree
{
  public:
    int MAX_NUM_LEVEL=8;
    int DEFAULT_NUM_LEVEL=3;
    int DEFAULT_ORDER_OF_APPROXIMATION = 3;

    int dimension = 2;

    int numOfLevels;                         // number of levels is l - 1
    int currLevel;                           // current level is l

    std::vector<Point> x;
    std::vector<Point> y;

    Potential potential;
    int order_of_approximation;

    std::vector<std::vector<Box> > tree_structure;              // an array of structs

    long numOpsIndirect;
    long numOpsDirect;

    FmmTree();                                // Constructor
    FmmTree(int level, Potential potential);  // Constructor
    FmmTree(int level, std::vector<Point> source, std::vector<Point> target, Potential potential);

    void initStruct();

    int getNumOfLevels() { return this->numOfLevels; };
    int getClusterThreshold();
    int getIndex(std::vector<Point> &z, Point &p);
    Box getBox(int level, int index) { return tree_structure[level][index]; };

    void printX ();
    void printY ();
    void printBoxInformation();
    void printTreeStructure();
    std::vector<double> solve(std::vector<double> &u);
    std::vector<double> solveDirect(std::vector<double> &u);

  private:
    void upwardPass(std::vector<double> &u);
    void downwardPass1();
    void downwardPass2();

};
*/

// construct the most basic tree (number of refinement levels is 4)
// that can still use FMM
FmmTree::FmmTree()
       :
       numOfLevels(DEFAULT_NUM_LEVEL),
       //currLevel(numOfLevels-1),
       tree_structure(numOfLevels),
       numOpsIndirect(0),
       numOpsDirect(0)
{}


// Explanation of Constructor FmmTree:
//
// the first argument 'level' of the constructor passes in the largest
// refinement level L.  However, this value is only true when counting of
// levels starts at the index l = 1
// Indexing in C++ starts with zero, therefore with respect to
// starting the count at l = 0, the largest refinement level
// (the first argument 'level' of the constructor below is)
// is L = level - 1
// We see this taken into account with currLevel and the
// numOfLevels-1 in the initializer list
FmmTree::FmmTree(int nlevels, std::vector<Point> &sources, std::vector<Point> &targets, Potential &potential)
       :
       numOfLevels(nlevels),
       //currLevel(numOfLevels-1),
       potential(potential),
       order_of_approximation(potential.getOrderApprox()),
       tree_structure(numOfLevels),
       numOpsIndirect(0),
       numOpsDirect(0)
{
  // need to assert that levelOfBox is between 0 and 8
  // if refinement level is greater than 8, then bitwise operations will be incorrect (8 bits)
  // a refinement level less than 0 does not make sense (not defined)
  // (series approximation convergence not guaranteed)
  // Note: level gives refinement level based on count beginning with 1
  //       levels l = 1, 2, 3, 4, 5, ....
  assert(nlevels>0 && "FmmTree level < 1");
  assert(nlevels<=8 && "FmmTree level > 8");

  unsigned int length = sources.size();
  x.resize(length);
  for (unsigned int i=0; i<length; ++i)
    x[i] = sources[i];

  length = targets.size();
  y.resize(length);
  for (unsigned int i=0; i<length; ++i)
	y[i] = targets[i];

  initStruct();
}

/**
 * Explanation of initStruct()
 *
 * First For Loop:
 *
 * Incrementing on i (i is row number of the struct matrix as well as refinement level)
 *
 * Each element of struct (let's say row number i) is initialized as a vector (column corresponding to the row)
 * of Box objects with its size being the number of cells for the refinement level i (4^i)
 *
 * Inside the first for loop is a second loop incrementing on j
 * (j is the column number of the struct matrix as well as the cell index n)
 *
 * Thinking of struct as a matrix.  Each element of the (matrix) struct is also
 * a Box object - Box(i,j,potential.getP()) where
 * - i is the row number
 *   - with respect to the Box constructor this the refinement level
 * - j is the column number
 *   - with respect to the Box constructor this is the cell index
 * - p is an integer from the potential object
 *   - with respect to the Box constructor this is the truncation index for the series approximation
 *
 * Second and Third For Loop:
 *
 * The second loop increments over i through each the source particle x[i].
 *  - For the refinement level numOfLevels-1, getBoxIndex determines the cell index n for each source particle x[i]
 *  - The addX member function of Box is called to add the source point x[i] to that cell (or box)
 *
 * The third loop increments over i through each the target particle y[i].
 *  - For the refinement level numOfLevels-1, getBoxIndex determines the cell index n for each source particle y[i]
 *  - The addY member function of Box is called to add the target point y[i] to that cell (or box)
 *
 */

void FmmTree::initStruct()
{
  int cells_in_tree_structure_i;

  for (unsigned int i=0; i<tree_structure.size(); ++i)
  {
    // creating all 8^i cells (boxes) for refinement level i
	cells_in_tree_structure_i = std::pow(8,i);
std::cout << "tree_structure[" << i << "] = " << tree_structure[i].size() << std::endl;
    tree_structure[i].resize(cells_in_tree_structure_i);
std::cout << "tree_structure[" << i << "] = " << tree_structure[i].size() << std::endl;
    // setting index, level and last term p in series sum for cell
    // in tree_structure at level i
    for (unsigned int j=0; j<tree_structure[i].size(); ++j)
    {
      tree_structure[i][j].setLevel(i);              // refinement level of box (cell) is i
      tree_structure[i][j].setIndex(j);              // index of box (cell) is j
      tree_structure[i][j].setP(potential.getP());
      tree_structure[i][j].setOrderApprox(potential.getOrderApprox());
//      std::cout << "I am here?" << std::endl;

    }
  }
  // using getBoxIndex to perform sorting of source and target particles
  // into boxes (cells) for currLevel (numOfLevel-1)
//  std::cout << "x.size() = " << x.size() << "\n";
  for (unsigned int i=0; i<x.size(); ++i)
  {
    tree_structure[numOfLevels-1][x[i].getBoxIndex(numOfLevels-1)].addX(x[i]);
  }
  for (unsigned int i=0; i<y.size(); ++i)
    tree_structure[numOfLevels-1][y[i].getBoxIndex(numOfLevels-1)].addY(y[i]);

}

// getting the largest number of points (source x or target y) in a cell
// for level numOfLevels-1 (highest refinement level when counting of l
// starts with l = 0
int FmmTree::getClusterThreshold()
{
  int ans = 0;
  for (unsigned int i=0; i<tree_structure[numOfLevels-1].size(); ++i)
  {
    int xlength = tree_structure[numOfLevels-1][i].getSizeX();
    int ylength = tree_structure[numOfLevels-1][i].getSizeY();
    if (xlength>ans)
      ans = xlength;
    if (ylength>ans)
      ans = ylength;
  }
  return ans;
}

// returns index of p in vector z
// point p must be a point in vector z
// else index returned is ans = -1
int FmmTree::getIndex(std::vector<Point> &z, Point &p)
{
  int ans = -1;
  for (unsigned int i=0; i<z.size(); ++i)
    if (z[i].equals(p))
    {
      ans = i;
      i = z.size();
    }
  return ans;
}

void FmmTree::printX()
{
  for (unsigned int i=0; i<x.size(); ++i)
  {
    std::string x_coords = this->x[i].coordToString();
	std::cout << x_coords << '\n';
  }

}

void FmmTree::printY()
{
  for (unsigned int i=0; i<x.size(); ++i)
  {
    std::string y_coords = this->y[i].coordToString();
	std::cout << y_coords << '\n';
  }

}

void FmmTree::printBoxInformation()
{
  for (unsigned int i=0; i<tree_structure[numOfLevels-1].size(); ++i)
  {
    std::cout << "For the Box at tree_structure[" << numOfLevels-1 << "][" << i << "]" << "\n";
    tree_structure[numOfLevels-1][i].printLevel();
    tree_structure[numOfLevels-1][i].printIndex();
    tree_structure[numOfLevels-1][i].printSizeX();
    tree_structure[numOfLevels-1][i].printSizeY();
  }
}


//void FmmTree::printTreeStructure()
//{
//  for (int i=0; i<numOfLevels; ++i)
//  {
//    std::cout << 'tree structure level' << i
//    		  << 'has ' << tree_structure[i].size() << 'cells' << '\n';
    //for (unsigned int j=0; j<tree_structure[i].size(); ++j)
    //  std::cout << 'cell[' << i << '][' << j << '] has '
    //            << tree_structure[i][j].getSizeX() << 'points.'
    //            << '\n';
//  }
//}


//std::vector<double> FmmTree::solve(std::vector<double> &u)
std::vector<std::vector<std::complex<double> > > FmmTree::solve(std::vector<double> &u)
{
  // v is the answer to the potential calculation using FMM
  // for each target y[i] (element of y), we will have calculated the potential
  // v[i] due to all the sources x using FMM
//  std::vector<double> v;
//  v.resize(y.size());
  std::vector<std::vector<std::complex<double> > > v(4, std::vector<std::complex<double> >());
  v[0].resize(y.size());
  v[1].resize(y.size());
  v[2].resize(y.size());
  v[3].resize(y.size());
  //std::vector<double> answer;

  std::cout << "Starting updward pass..." << "\n";
  upwardPass(u);
  std::cout << "Completed updward pass..." << "\n";

  std::cout << "Starting downward pass 1..." << "\n";
  downwardPass1();
  std::cout << "Completed downward pass 1..." << "\n";

  std::cout << "Starting downward pass 2..." << "\n";
  downwardPass2();
  std::cout << "Completed downward pass 2..." << "\n";

  std::cout << "Completing FMM..." << "\n";

  // Explanation of Nested For Loops in Code Below - could be separate member function of FmmTree
  //
  // [0] - for loop through the cells at the highest refinement level(lowest level) numOfLevels
  //       (index starts on zero, so numOfLevels-1) determining and summing the regular and singular part
  //       of the potential function.  The regular part of the potential function is the approximated
  //       series due to transformed far-field approximations of sources in the interaction list and the rest
  //       of the domain (beyond the list - not including near neighbors and made up of the interaction lists
  //       of the parental hierarchy).  The singular part of the potential function calculation consists of
  //       direct calculations of the potential function for sources in the near neighbors acting on the
  //       target of the cell.
  //   [1] - getting a reference thisBox for the box to be worked on
  //   [2] - getting the target points yPoints of this box
  //   [3] - if there are target points in this box
  //     [4] - for each target point yPoints[j]
  //       [5] - getting a reference thisY for the target point
  //       [6] - declaring the regular part of the potential calculation
  //             where the source points x[i] are far enough away from thisBox
  //             that the potential calculation can be approximated by a series
  //       [7] - retrieving series coefficients D (could be done outside this loop)
  //       [8] - retrieving powers of R-expansion for thisY
  //       [9] - for each of p terms of R-expansion series
  //        [10] - calculating and adding the first p terms of the series - a
  //               truncated approximation to the infinite series
  //      [11] - initializing the singular part of the potential calculation
  //             where the source points x[i] are too close to approximate
  //             the potential calculation with a series and the calculation
  //             must be done directly
  //      [12] - declaring and collecting the neighbors and this box's indices
  //             in this line and the two lines above
  //      [13] - for each near neighbor (including this box)
  //        [14-15] - obtaining a reference thisNeighborsBox to neighor's box
  //                  and x-values (sources) thisNeighborsX in neighbor's box
  //        [16] - if there are source terms thisNeighborsX in neighbor's box
  //          [17] - for each of the source terms thisNeighborsX[q]
  //            [18-19] - getting reference thisX and thisU for thisNeighborsX[q]
  //                      and the charge u for that source particle thisX
  //            [20-22] - declaring and initializing coordinates for thisX and thisY
  //            [23-25] - if thisX and thisY are not the same
  //                      using relative and absolute comparison for the cases
  //                      where thisXCoord and thisYCoord may be large or small
  //                      if thisXCoord and thisYCoord are small then maxXYOne is 1
  //                      and we are comparing absolutely and not relatively
  //                      (relative comparison - percentages - could be smaller than
  //                      machine epsilon
  //              [26-28] - calculating the potential directly (singular part)
  //                        for thisX on thisY (and taking only real part?)
  //      [29] - once completing direct calculations on thisY for all sources thisX
  //             in the near neighbors, then adding the result (singular part) to
  //             the result from the series approximations to the potential calculation
  //             for sources far enough away (regular part)
  //             Making sure to put this final result in the same location (have same index value)
  //             as the corresponding location (index value) of yPoints[j] = thisY in the vector
  //             of target points y
  //
  for (unsigned int i=0; i<tree_structure[numOfLevels-1].size(); ++i)                       // 0
  {
    Box& thisBox = tree_structure[numOfLevels-1][i];                                        // 1
    std::vector<Point> yPoints = thisBox.getY();                                            // 2
    if (yPoints.size() > 0)                                                                 // 3
    {
      for (unsigned int j=0; j<yPoints.size(); ++j)                                         // 4
      {
        Point& thisY = yPoints[j];                                                          // 5
        std::complex<double> regPart = 0.0;                                                 // 6

        std::vector<std::complex<double> > regPartGrad(3, 0.0);                                             // 6

        std::vector<std::vector<std::complex<double> > >
                            d = thisBox.getD();                                             // 7

        std::vector<std::vector<std::complex<double> > >
          rVec = potential.getRVector(thisY.getCoord(), thisBox.getCenter().getCoord());    // 8

        std::vector<std::vector<std::vector<std::complex<double> > > >
          rVecGrad = potential.getRVectorGrad(rVec);    // 8

        numOpsIndirect+=potential.getP();
        for (unsigned int k=0; k<((unsigned int)order_of_approximation + 1); ++k)           // 9
          for (unsigned int i=0; i<((unsigned int)order_of_approximation + 1); ++i)
          {
            regPart += (d[i][k] * rVec[i][k]);                                              // 10
            regPartGrad[0] += (d[i][k] * rVecGrad[0][i][k]);                                              // 10
            regPartGrad[1] += (d[i][k] * rVecGrad[1][i][k]); // * std::complex<double>{0.0,-1.0});                                              // 10
            regPartGrad[2] += (d[i][k] * rVecGrad[2][i][k]); // * std::complex<double>{-1.0, 0.0});                                              // 10
            numOpsIndirect++;

            std::cout << "d[" << i << "][" << k << "] = " << d[i][k] << std::endl;
            std::cout << "rVecGrad[1][" << i << "][" << k << "] = " << rVecGrad[1][i][k] << std::endl;
            std::cout << "regPartGrad[1] = " << regPartGrad[1] << std::endl;

          }
        std::cout << "regPart = " << regPart << std::endl;
        std::cout << "regPartGrad[0] = " << regPartGrad[0] << std::endl;
        std::cout << "regPartGrad[1] = " << regPartGrad[1] << std::endl;
        std::cout << "regPartGrad[2] = " << regPartGrad[2] << std::endl;

        double sinPart = 0.0;                                                               // 11

        std::vector<double> sinPartGrad(3, 0.0);                                                               // 11

        std::vector<int> neighbors_indexes_and_me;
        thisBox.getNeighborsIndex(neighbors_indexes_and_me);
        neighbors_indexes_and_me.push_back(thisBox.getIndex());                             // 12

        for (unsigned int m=0; m<neighbors_indexes_and_me.size(); ++m)                      // 13
        {
          std::cout << "We are at neighbors" << std::endl;
          Box& thisNeighborsBox
              = tree_structure[numOfLevels-1][neighbors_indexes_and_me[m]];                 // 14
          std::vector<Point> thisNeighborsX = thisNeighborsBox.getX();                      // 15
 //         if (thisNeighborsX.size() > 0)                                                    // 16
 //         {
        	std::cout << "We are in singular part" << std::endl;
            for (unsigned int q=0; q<thisNeighborsX.size(); ++q)                            // 17
            {
              Point& thisX = thisNeighborsX[q];                                             // 18
              double thisU = u[getIndex(x, thisX)];                                         // 19
              std::vector<double> thisYCoord, thisXCoord, diffXY;                           // 20
              double inc, length_thisYCoord, length_thisXCoord, length_diffXY;

              std::vector<double> incGrad(3, 0.0);

              thisYCoord = thisY.getCoord();                                                // 21
              thisXCoord = thisX.getCoord();                                                // 22
              diffXY.resize(3);
              diffXY[0] = thisYCoord[0] - thisXCoord[0];
              diffXY[1] = thisYCoord[1] - thisXCoord[1];
              diffXY[2] = thisYCoord[2] - thisXCoord[2];
              length_thisYCoord = pow(thisYCoord[0],2.0)
            		              + pow(thisYCoord[1],2.0) + pow(thisYCoord[2],2.0);
              length_thisYCoord = pow(length_thisYCoord,0.5);
              length_thisXCoord = pow(thisXCoord[0],2.0)
            		              + pow(thisXCoord[1],2.0) + pow(thisXCoord[2],2.0);
              length_thisXCoord = pow(length_thisXCoord,0.5);
              length_diffXY = pow(diffXY[0],2.0) + pow(diffXY[1],2.0) + pow(diffXY[2],2.0);
              length_diffXY = pow(length_diffXY,0.5);

              double maxXY = std::max(length_thisXCoord,
            	                         length_thisYCoord);
              double maxXYOne = std::max(1.0,maxXY);                                        // 23
              if (length_diffXY
                      <= std::numeric_limits<double>::epsilon()*maxXYOne)                   // 24
              {
                // Do nothing - thisYCoord and thisXCoord are the same point
            	// (up to machine epsilon)
            	// This can happen when the target and source points are the same
            	// Specifically, this happens when thisNeighborsBox is thisBox
            	// (last q-value), and thisX and thisY are a target point and
            	// a source point in the same box (and the target and source
            	// points are the same).
              }
              else // target and source points thisY and thisX are not the same             // 25
              {
                  inc = potential.direct(thisY.getCoord(), thisX.getCoord());               // 26

                  incGrad = potential.directGrad(thisY.getCoord(), thisX.getCoord());               // 26

                  inc = inc * thisU;                                                        // 27
                  sinPart += inc;                                                           // 28

                  incGrad[0] = incGrad[0] * thisU;                                                        // 27
                  incGrad[1] = incGrad[1] * thisU;                                                        // 27
                  incGrad[2] = incGrad[2] * thisU;                                                        // 27
                  sinPartGrad[0] += incGrad[0];                                                           // 28
                  sinPartGrad[1] += incGrad[1];                                                           // 28
                  sinPartGrad[2] += incGrad[2];                                                           // 28

                  numOpsIndirect++;
              }
            }
            std::cout << "sinPart = " << sinPart << std::endl;
            std::cout << "sinPartGrad[0] = " << sinPartGrad[0] << std::endl;
            std::cout << "sinPartGrad[1] = " << sinPartGrad[1] << std::endl;
            std::cout << "sinPartGrad[2] = " << sinPartGrad[2] << std::endl;

 //         }
        }

        v[0][getIndex(y, thisY)] = sinPart + regPart.real();                                   // 29
        v[1][getIndex(y, thisY)] = sinPartGrad[0] + regPartGrad[0];                                   // 29
        v[2][getIndex(y, thisY)] = sinPartGrad[1] + regPartGrad[1];                                   // 29
        v[3][getIndex(y, thisY)] = sinPartGrad[2] + regPartGrad[2];                                   // 29
      }
    }
  }

  std::cout << "Completed FMM..." << "\n";
//  std::cout << "numOpsIndirect = " << numOpsIndirect << std::endl;

  return v;
}


/* Explanation of upwardPass
 *
 * The 1st part of the upward pass starts at the highest level of refinement (lowest level l = L)
 * l = numOfLevels-1.  The first step of the upward pass is to obtain the far-field
 * expansions (really just the coefficients of the expansions) for all the sources
 * w.r.t. the center of the cell in which they are located at this refinement level.
 * A for loop over the cells at this refinement level starts this member function.
 * While looping over a cell the potential functions corresponding to all source points
 * of the cell are approximated by far-field expansions whose center is the center of
 * the cell.  Corresponding coefficients of these approximations are added together
 * to form one series for all the sources of the cell (since the powers of each source's series
 * were the same).  Once this loop is complete we have far-field s-expansion for each cell
 * (all the sources of each cell) at the lowest level.
 *
 * The 2nd part of the upward pass is a nested loop.  The outer loop is over each refinement
 * level, working upward.   The inner loop is over the cells at that particular refinement level.
 * The purpose of the second part of the upward pass is to shift the center of the source expansion
 * for a cell at the current refinement level (child level) up to the center at the parent level.
 * Once this nested loop is complete, we have a far-field expansion for each box (cell) at any
 * level l >= 2 that incorporates all of the sources located within that cell (box).
 *
 * The purpose of the 2nd part of the upward pass is to minimize the work in obtaining approximations to
 * the potential function between a target and a source that are far from each other.  For example,
 * suppose that at level l = 2, a source and target are one cell away.  Then at this level, the source cell
 * is in the interaction list of the target cell and we would approximate the potential function with
 * a far field expansion.  Further, we can approximate all the sources of the source cell with an s-expansion
 * Therefore, to save operations we take the s-expansion for the cell (all source incorporated) instead of going
 * down to the lowest refinement level taking each source point expansion individually.  The purpose of the downward
 * pass is to determine when distances between cells are far enough to perform these operations.
 *
 *
 */


void FmmTree::upwardPass(std::vector<double> &u)
{
  std::cout << "Entering Upward Pass " << std::endl;

  // 1st Part of the Upward Pass
  // recall that tree structure is an array.  Each row (row associated with first index)
  // contains all the boxes (cells) for the particular level.
  // The for loop below iterates through all of the cells at
  // the highest level l = numOfLevels-1 to obtain the power series expansion
  // for the potential function associated with each point in the particular cell.
  for (unsigned int i=0; i<tree_structure[numOfLevels-1].size(); ++i)
  {
	// The first two steps set up an alias for the ith box of the refinement
	// to allow for the Box to be worked on, and pull the points of the box.
    Box& thisBox = tree_structure[numOfLevels-1][i];
    std::vector<Point> xPoints = thisBox.getX();

    // Here we make sure that there are actually points in the box
    // and then loop through the points
    // We then make an alias for the point to allow work on its variables.
    // Our first goal in this loop is to get the coefficients of the individual
    // expansions that approximate the potential functions corresponding to the
    // source points of this cell.
    // The call to potential.getSCoeff
    if (xPoints.size() > 0)
      for (unsigned int j=0; j<xPoints.size(); ++j)
      {
    	Point& thisX = xPoints[j];
        std::cout << "xPoints[" << j << "] = " << xPoints[j].coordToString() << std::endl;
        std::vector<std::vector<std::complex<double> > >
           B = potential.getSCoeff(thisX.getCoord(), thisBox.getCenter().getCoord());

//std::cout << "Got back B in Part One upward pass" << std::endl;
//std::cout << "Size of B is " << B.size() << std::endl;
//std::cout << "Size of C is " << C.size() << std::endl;
        numOpsIndirect++;
        numOpsIndirect+=potential.getP();        // O(p) flops in getSCoeff
                                                 // p*(1 subtraction, 1 pow, 1 division, 1 mult by -1)

        // determining charge associated with point thisX
        double thisU = u[getIndex(x,thisX)];
        std::cout << "thisU = " << thisU << std::endl;

        std::vector<std::vector<std::complex<double> > > getBackCV1 = thisBox.getC();
        for (unsigned int i=0; i<getBackCV1.size(); ++i)
          for (unsigned int j=0; j<getBackCV1.size(); ++j)
            std::cout << "CofBoxV1[" << i << "][" << j << "] = " << getBackCV1[i][j] << std::endl;

        for (unsigned int k=0; k<B.size(); ++k)
          for (unsigned int j=0; j<B.size(); ++j)
        {
          B[j][k] = B[j][k]*thisU;
          //std::cout << "B[" << j << "][" << k << "] = " << B[j][k] << std::endl;
          numOpsIndirect++;
        }
        thisBox.addToC(B);

        for (unsigned int i=0; i<B.size(); ++i)
          for (unsigned int j=0; j<B.size(); ++j)
              std::cout << "B[" << i << "][" << j << "] = " << B[i][j] << std::endl;

        std::vector<std::vector<std::complex<double> > > getBackCV2 = thisBox.getC();
        for (unsigned int i=0; i<getBackCV2.size(); ++i)
          for (unsigned int j=0; j<getBackCV2.size(); ++j)
            std::cout << "CofBoxV2[" << i << "][" << j << "] = " << getBackCV2[i][j] << std::endl;

        numOpsIndirect+=potential.getP();
//std::cout << "Added B in Part One upward pass" << std::endl;
      }

  }
  // 2nd Part of the Upward Pass
  // Translating the far-field expansion from the child centers
  // to parent centers upward the parent hierarchy
  // Here el stands for refinement level and
  //      k  stands for box (cell) index
  for (int el = numOfLevels-1; el>=2; --el)
  {
    std::cout << "Upward pass level " << +el << "\n";
    for (unsigned int k=0; k<tree_structure[el].size(); ++k)
    {
      Box& thisBox = tree_structure[el][k];
      int parentBoxLevel = el-1;
      Box& parentBox = tree_structure[parentBoxLevel][thisBox.getParentIndex()];
      numOpsIndirect++;
      std::vector<double> from = thisBox.getCenter().getCoord();
      numOpsIndirect++;
      std::vector<double> to = parentBox.getCenter().getCoord();
      numOpsIndirect++;

      // translating the thisBox's series that has coeffs thisBoxC
      // from its center at location 'from' = thisBox.getCenter().getCoord()
      // to its parent's center at location 'to' = parentBox.getCenter().getCoord()
      // The new series with parent center can be added to the parent's
      // C series since the powers for each term of the two series are now the same
      std::vector<std::vector<std::complex<double> > >
        newCoeffs = potential.getSS(from, to, thisBox.getC());
      parentBox.addToC(newCoeffs);

      if (newCoeffs[0][0] != 0.0 && thisBox.getIndex() == 511 && el == numOfLevels-1)
      {
    	  std::cout << "level = " << el << " and " << "thisBox.getIndex = " << thisBox.getIndex() << std::endl;
          std::cout << "from = [" << from[0] << ", " << from[1] << ", " << from[2] << "]" << std::endl;
          std::cout << "to = [" << to[0] << ", " << to[1] << ", " << to[2] << "]" << std::endl;

          std::vector<std::vector<std::complex<double> > > CofBox = thisBox.getC();
          for (unsigned int i=0; i<CofBox.size(); ++i)
            for (unsigned int j=0; j<CofBox.size(); ++j)
              std::cout << "CofBox[" << i << "][" << j << "] = " << CofBox[i][j] << std::endl;

        for (unsigned int i=0; i<newCoeffs.size(); ++i)
          for (unsigned int j=0; j<newCoeffs.size(); ++j)
            std::cout << "parent newCoeffs[" << i << "][" << j << "] = " << newCoeffs[i][j] << std::endl;
      }

   /*
if (thisBox.getIndex() == 365)
{
  std::vector<double> thisBoxgetC = thisBox.getC();
  std::cout << "from = " << thisBox.getCenter().coordToString() << std::endl;
  std::cout << "to = " << parentBox.getCenter().coordToString() << std::endl;
  for (unsigned int i=0; i<newCoeffs.size(); ++i)
    std::cout << "newCoeffs[" << i << "] = " << newCoeffs[i] << std::endl;
  std::cout << std::endl << std::endl;

  for (unsigned int i=0; i<thisBoxgetC.size(); ++i)
	std::cout << "thisBoxgetC[" << i << "] = " << thisBoxgetC[i] << std::endl;
}
    */
      numOpsIndirect+=pow(potential.getP(),2); // matrix-vector multiply
      numOpsIndirect+=potential.getP();        // vector addition
    }
  }
std::cout << "Completed Second Part of Upward Pass" << std::endl;
}



/* Explanation of downwardPass1
 *
 * This is the First Part of the Downward Pass.  The member function is a
 * nested loop.  The outer loop is over the levels of refinement and works
 * from level l = 2 (coarse refinement level) to the bottom level (most refined).
 * Once the level is set the middle loop is over the cells of the refinement level.
 * On each pass through the middle loop, the interactive list of the cell being
 * addressed is obtained.
 *
 * Recall that the E4 index is the interactive list of a cell.
 * For a given cell, the E4 list is the set of all peer cells that are the children
 * of the given cell's parent's nearest neighbors excluding the nearest neighbors of the cell.
 *
 * With the interactive list of the cell of interest obtained, the inner loop of the member
 * function is then started over the cells of the interactive list.  During a pass in the inner loop
 * the far-field coefficients of the interactive list cell set by the inner loop are S|R transformed
 * from the interactive list cell center to the center of the cell of interest set by the outer loop.
 * The new coefficients of the series due to sources in the interactive list are saved (added to) to
 * the Dtilde vector of the cell of interest.
 *
 *
 */
void FmmTree::downwardPass1()
{
  std::vector<int> thisBoxNeighborsE4Indexes;
  std::vector<double> from;
  std::vector<double> to;

  for (int el=2; el<numOfLevels; ++el)
  {

	for (unsigned int k=0; k<tree_structure[el].size(); ++k)
    {
      // getting the neighbor indices and interaction list index
      thisBoxNeighborsE4Indexes.resize(0);
      Box& thisBox = tree_structure[el][k];
      thisBox.getNeighborsE4Index(thisBoxNeighborsE4Indexes);

      // translating the far field series with old coefficients C to
      // a near field series with new coefficients Dtilde (see Main.cc notes)
      for (unsigned int j=0; j<thisBoxNeighborsE4Indexes.size(); ++j)
      {
        Box& thisBoxE4Neighbor = tree_structure[el][thisBoxNeighborsE4Indexes[j]];
        from = thisBoxE4Neighbor.getCenter().getCoord();
        ++numOpsIndirect;
        to = thisBox.getCenter().getCoord();
        ++numOpsIndirect;

        std::vector<std::vector<std::complex<double> > >
           newCoeffs = potential.getSR(from, to, thisBoxE4Neighbor.getC());
        thisBox.addToDtilde(newCoeffs);

        Point target(0.02641560818,0.02641560818,0.02641560818);
        Point source(0.90141560818,0.90141560818,0.90141560818);
        if (newCoeffs[0][0] != 0.0 && el == 2 && k == (unsigned int)target.getBoxIndex(2))
        {
            std::cout << "from (E4Neighbor center) = [" << from[0] << ", " << from[1] << ", " << from[2] << "]" << std::endl;
            std::cout << "to (parent center) = [" << to[0] << ", " << to[1] << ", " << to[2] << "]" << std::endl;

            std::vector<std::vector<std::complex<double> > > CofBox = thisBoxE4Neighbor.getC();
            for (unsigned int i=0; i<CofBox.size(); ++i)
            {
              for (unsigned int j=0; j<CofBox.size(); ++j)
                std::cout << " Cn[" << i << "][" << j << "] = " << CofBox[i][j];
              std::cout << std::endl;
            }

            for (unsigned int i=0; i<newCoeffs.size(); ++i)
            {
              for (unsigned int j=0; j<newCoeffs.size(); ++j)
                std::cout << " Cp[" << i << "][" << j << "] = " << newCoeffs[i][j];
              std::cout << std::endl;
            }

            std::vector<std::vector<std::complex<double> > > RVecP = potential.getRVector(target.getCoord(), thisBox.getCenter().getCoord());
            std::vector<std::vector<std::vector<std::complex<double> > > > RVecGradP = potential.getRVectorGrad(RVecP);
            for (unsigned int i=0; i<RVecP.size(); ++i)
            {
              for (unsigned int j=0; j<RVecP.size(); ++j)
                std::cout << " RVecP[" << i << "][" << j << "] = " << RVecP[i][j];
              std::cout << std::endl;
            }

            for (unsigned int i=0; i<RVecGradP[0].size(); ++i)
            {
              for (unsigned int j=0; j<RVecGradP[0].size(); ++j)
                std::cout << " RVGP[0][" << i << "][" << j << "] = " << RVecGradP[0][i][j];
              std::cout << std::endl;
            }

            for (unsigned int i=0; i<RVecGradP[1].size(); ++i)
            {
              for (unsigned int j=0; j<RVecGradP[1].size(); ++j)
                std::cout << " RVGP[1][" << i << "][" << j << "] = " << RVecGradP[1][i][j];
              std::cout << std::endl;
            }

            for (unsigned int i=0; i<RVecGradP[2].size(); ++i)
            {
              for (unsigned int j=0; j<RVecGradP[2].size(); ++j)
                std::cout << " RVGP[2][" << i << "][" << j << "] = " << RVecGradP[2][i][j];
              std::cout << std::endl;
            }

            std::complex<double> sum1 = 0.0;
            std::complex<double> sum2 = 0.0;
            std::complex<double> sum3 = 0.0;
            std::complex<double> sum4 = 0.0;
              for (unsigned int s = 0; s<RVecP.size(); ++s)
                for (unsigned int t = 0; t<RVecP.size(); ++t)
                {
                  sum1 += newCoeffs[s][t]*RVecP[s][t];
                  sum2 += newCoeffs[s][t]*RVecGradP[0][s][t];
                  sum3 += newCoeffs[s][t]*RVecGradP[1][s][t];
                  sum4 += newCoeffs[s][t]*RVecGradP[2][s][t];
                }
            std::cout << "parent near field sum for potential = " << sum1 << std::endl;
            std::cout << "parent near field sum for grad comp 1 = " << sum2 << std::endl;
            std::cout << "parent near field sum for grad comp 2 = " << sum3 << std::endl;
            std::cout << "parent near field sum for grad comp 3 = " << sum4 << std::endl;

          //std::vector<std::vector<std::complex<double> > > CofBox = thisBox.getC();
          //for (unsigned int i=0; i<CofBox.size(); ++i)
          //  for (unsigned int j=0; j<CofBox.size(); ++j)
          //    std::cout << "CofBox[" << i << "][" << j << "] = " << CofBox[i][j] << std::endl;

//          for (unsigned int i=0; i<newCoeffs.size(); ++i)
//            for (unsigned int j=0; j<newCoeffs.size(); ++j)
//              std::cout << "newCoeffs[" << i << "][" << j << "] = " << newCoeffs[i][j] << std::endl;
        }

/*
        if (newCoeffs[0][0] != 0.0)
        {
            std::cout << "from = [" << from[0] << ", " << from[1] << ", " << from[2] << "]" << std::endl;
            std::cout << "to = [" << to[0] << ", " << to[1] << ", " << to[2] << "]" << std::endl;

          //std::vector<std::vector<std::complex<double> > > CofBox = thisBox.getC();
          //for (unsigned int i=0; i<CofBox.size(); ++i)
          //  for (unsigned int j=0; j<CofBox.size(); ++j)
          //    std::cout << "CofBox[" << i << "][" << j << "] = " << CofBox[i][j] << std::endl;

          for (unsigned int i=0; i<newCoeffs.size(); ++i)
            for (unsigned int j=0; j<newCoeffs.size(); ++j)
              std::cout << "newCoeffs[" << i << "][" << j << "] = " << newCoeffs[i][j] << std::endl;
        }
*/

        numOpsIndirect+=pow(potential.getP(),2); // matrix-vector multiply
        numOpsIndirect+=potential.getP();        // vector addition
      }

    }

  }

  std::cout << "Completed downard pass part I" << std::endl;
}

/* Explanation of downwardPass2
 *
 * 											-------------------------------
 * 										   |   |   |   |   ||   |   |   |   |
 * 	c  - cell							   |--pi---|--pi---||--pi---|--pi---|
 * 	p  - parent							   |___|___|___|___||___|___|___|___|
 * 	n  - nearest neighbor				   |   |   |   |   ||   |   |   |   |
 * 	i  - interaction list cell 			   |--pi---|--pi---||--pi---|--pi---|
 *  pn - parent nearest neigh.  		   |   |   |   |   ||   |   |   |   |
 * 	pi - parent interaction list cell	   |===|===|===|===||===|===|===|===|
 * 										   | i | i | i | i || i | i |   |   |
 * 										   |--pn---|--pn---||--pn---|--pi---|
 * 										   |_i_|_n_|_n_|_n_||_i_|_i_|___|___|
 * 										   | i | n | c | n || i | i |   |   |
 * 										   |--pn---|---p---||--pn---|--pi---|
 * 										   |_i_|_n_|_n_|_n_||_i_|_i_|___|___|
 *
 *                                         l = 3
 *
 * This is the Second Part of the Downward Pass.  The main component of this
 * member function is a nested loop.  However the first step of the function is a loop
 * through refinement level l = 2 copying the Dtilde coefficients obtained in the First
 * Part of the Downward Pass to the D Coefficients.  Level l = 2 is special since the
 * interactive list makes up the rest of the domain when considered with a cell's
 * nearest neighbors (shown above).  The other levels do not have this property.
 *
 * However, as we move down another level and consider a particular cell at that level,
 * the rest of the domain that is outside the nearest neighbors of that cell can be expressed
 * as a union of the interactive lists of the parent and grandparent cells.  For example, consider
 * the cell c in the level l = 3 illustration above.  Outside of the nearest neighbors n of cell
 * c, the rest of the domain is the union of the interaction list of the cell with the interaction
 * list of the parent cell.
 *
 * More generally, since a child cell's interactive list are made up of children of the nearest neighbors
 * of the parent cell, and the parent cell's interactive list are peer cells beyond its nearest neighbors,
 * the interactive lists of a parent and child do not have any cells (or source points) in
 * common.  That is, the interactive lists of a parent and a child cell are mutually exlusive and
 * there is no overlap in accounting for source points in the union of a parent and child's
 * interactive lists.  Therefore, once the nearest neighbors of a cell are accounted for, the rest of
 * the domain is the union of the hierarchy of the parents' mutually exclusive interaction lists all the
 * way up to level l = 2.
 *
 * In the main part of the code, the outer loop is over the levels of refinement, working
 * from level l = 2 (coarse refinement level) to the bottom level (most refined level).
 * The middle loop is over the cells of the refinement level and collects the indexes of the
 * child cells of the cell fixed by the loop.  The inner loop loops over the collection
 * of child cells and performs an R|R translation of the series whose coefficients have been
 * collected in D.
 *
 * That is, in the 1st part of the DownwardPass, the far-field series of all
 * the sources in the interaction list of a cell were S|R translated to the center of the cell.
 * At el = 2 ( l = 2) in the outerloop of downwardPass2 these series are then translated
 * to the children cells of each cell at l = 2.  In other words, the sources that
 * are outside the nearest neighbors of a cell at l = 2 were translated to the child cells
 * and are contained in the vector D.  That means that a child cell at l = 3 had all sources taken
 * care of outside of its interaction list (outside of the parent's nearest neighbors).  We also know from the
 * First Part of the Downward Pass that the child cell contains in Dtilde the series for the sources
 * in its interaction list.  Therefore for each cell at l = 3, the only sources left to consider
 * are those in the cell's nearest neighbors.  This last step is taken care of in the last part of the
 * FmmTree member function FmmTree::solve().
 *
 * Once the R|R translation takes place of the series from the interaction list of a parent
 * to the child's center, the new coefficients are saved in the D vector of the child's cell.
 * The coefficients from the series of the sources in the interaction list of the child cell
 * are also added to the D vector.  Then all that is not accounted for are the sources in the
 * nearest neighbors of the child.
 *
 * We continue with the next pass through the outer loop el = 3 (refinement level = 3) to make
 * sure the code is clear.  At l = 3, the children are gathered and the series contained in D
 * of each cell at l = 3 is translated to the children.  That means the interaction lists
 * sources for the parent and the grandparent at l = 2 are translated to the child.  Therefore,
 * all sources are accounted for with the child except the interaction list and nearest neighbor
 * sources.  Before completing the loop, the interaction list is accounted for and all that is
 * left are the nearest neighbors sources, which are taken care of in the end of FmmTree::solve.
 *
 */

void FmmTree::downwardPass2()
{
  std::vector<double> from, to;
  std::vector<int> children_indexes;

  for (unsigned int i=0; i<tree_structure[2].size(); ++i)
  {
	std::vector<std::vector<std::complex<double> > > DtildeCoeffs = tree_structure[2][i].getDtilde();
	tree_structure[2][i].addToD(DtildeCoeffs);
    numOpsIndirect+=potential.getP();

	if (tree_structure[2][i].getIndex() == 0)
	{
 	  std::cout << "level = " << 2 << " and " << "tree_structure[2][" << i << "].getIndex = "
 			                  << tree_structure[2][i].getIndex() << std::endl;
      std::vector<std::vector<std::complex<double> > > D = tree_structure[2][i].getD();
      for (unsigned int i=0; i<D.size(); ++i)
        for (unsigned int j=0; j<D.size(); ++j)
          std::cout << "D[" << i << "][" << j << "] = " << D[i][j] << std::endl;
	}
  }

  for (int el=2; el<numOfLevels-1; ++el)
  {
    for (unsigned int m=0; m<tree_structure[el].size(); ++m)
    {
      Box& thisBox = tree_structure[el][m];
      from = thisBox.getCenter().getCoord();
      numOpsIndirect++;
      children_indexes.resize(0);
      thisBox.getChildrenIndex(children_indexes);
      for (unsigned int k=0; k<children_indexes.size(); ++k)
      {
    	Box& thisBoxChild = tree_structure[el+1][children_indexes[k]];
        to = thisBoxChild.getCenter().getCoord();
        numOpsIndirect++;
    	std::vector<std::vector<std::complex<double> > >
    	  newCoeffs = potential.getRR(from, to, thisBox.getD());

        thisBoxChild.addToD(newCoeffs);
        numOpsIndirect+=pow(potential.getP(),2); // matrix-vector multiply
        numOpsIndirect+=potential.getP();        // vector addition

    	std::vector<std::vector<std::complex<double> > >
    	  DtildeCoeffs = thisBoxChild.getDtilde();
        thisBoxChild.addToD(DtildeCoeffs);
        numOpsIndirect+=potential.getP();

        std::vector<std::vector<std::complex<double> > > DChild = thisBoxChild.getD();
        if (thisBoxChild.getIndex() == 0.0)
        {
         	std::cout << "level = " << el << " and " << "thisBox.getIndex = " << thisBox.getIndex()
         			  << " and thisBoxChild.getIndex = " << thisBoxChild.getIndex() << std::endl;
            std::cout << "from = [" << from[0] << ", " << from[1] << ", " << from[2] << "]" << std::endl;
            std::cout << "to = [" << to[0] << ", " << to[1] << ", " << to[2] << "]" << std::endl;

          for (unsigned int i=0; i<newCoeffs.size(); ++i)
            for (unsigned int j=0; j<newCoeffs.size(); ++j)
              std::cout << "newCoeffs[" << i << "][" << j << "] = " << newCoeffs[i][j] << std::endl;

          //std::vector<std::vector<std::complex<double> > > DtildeofChildBox = thisBoxChild.getDtilde();
          for (unsigned int i=0; i<DtildeCoeffs.size(); ++i)
            for (unsigned int j=0; j<DtildeCoeffs.size(); ++j)
              std::cout << "DtildeCoeffs[" << i << "][" << j << "] = " << DtildeCoeffs[i][j] << std::endl;

          for (unsigned int i=0; i<DChild.size(); ++i)
            for (unsigned int j=0; j<DChild.size(); ++j)
              std::cout << "DChild[" << i << "][" << j << "] = " << DChild[i][j] << std::endl;

        }

      }
    }
  }

}



std::vector<double> FmmTree::solveDirect(std::vector<double> &u)
{
  std::cout << "Starting Direct Solve..." << "\n";

  std::vector<double> v(y.size());
  std::complex<double> potential_direct_calculation;
  std::vector<double> yPoint, xPoint, diffXY(3);
  double lengthY, lengthX, maxXY, maxXYOne;


  for (unsigned int j=0; j<v.size(); ++j)
  {
    for (unsigned int i=0; i<x.size(); ++i)
    {
//      std::cout << "particle interaction = (" << j << "," << i << ")" << std::endl;
      // taking care of relative and absolute difference
      // issues for when x[i] and y[j] are both small
      // or both large (see explanation in FmmTree member function solve
      // above)
      yPoint = y[j].getCoord();  xPoint = x[i].getCoord();
      lengthY = pow(yPoint[0],2.0) + pow(yPoint[1],2.0) + pow(yPoint[2],2.0);
      lengthY = pow(lengthY,0.5);
      lengthX = pow(xPoint[0],2.0) + pow(xPoint[1],2.0) + pow(xPoint[2],2.0);
      lengthX = pow(lengthX,0.5);
      maxXY = std::max(lengthY, lengthX);
      maxXYOne = std::max(1.0,maxXY);

      diffXY[0] = y[j].getX() - x[i].getX();
      diffXY[1] = y[j].getY() - x[i].getY();
      diffXY[2] = y[j].getZ() - x[i].getZ();
      if ( pow(pow(diffXY[0],2.0) + pow(diffXY[1],2.0) + pow(diffXY[2],2.0), 0.5)
                <= std::numeric_limits<double>::epsilon()*maxXYOne)
      {
        // Do nothing - y[j] and x[i] are the same point
    	// (up to machine epsilon)
    	// This happens when the target and source points are the same.
    	// Specifically, this happens when y[j] and x[i] are a target
    	// point and a source point in the same box (and the target
    	// and source points are the same).
    	// Also, physically it does not make sense for a particle y[j] to act
    	// on itself
      }
      else // target and source points y[j] and x[i] are not the same
      {
        potential_direct_calculation
          = u[i] * potential.direct(y[j].getCoord(),
      		                        x[i].getCoord());
        v[j] += potential_direct_calculation.real();
        numOpsDirect++;
      }
    }
  }

  std::cout << "Completed Direct Solve..." << "\n";
//  std::cout << "numOpsDirect = " << numOpsDirect << std::endl;

  return v;
}




std::vector<std::vector<double> > FmmTree::solveDirectGrad(std::vector<double> &u)
{
  std::cout << "Starting Direct Solve..." << "\n";

  std::vector<std::vector<double> > v(3, std::vector<double>());
  v[0].resize(y.size());
  v[1].resize(y.size());
  v[2].resize(y.size());

  for (unsigned int i=0; i<y.size(); ++i)
  {
    v[0][i] = 0.0;
    v[1][i] = 0.0;
    v[2][i] = 0.0;
  }

  std::cout << "y.size() = " << y.size() << std::endl;
  std::cout << "x.size() = " << x.size() << std::endl;
  std::cout << "u.size() = " << u.size() << std::endl;

  std::vector<double> potential_direct_grad_calculation(3, 0.0);
  std::vector<double> yPoint, xPoint, diffXY(3);
  double lengthY, lengthX, maxXY, maxXYOne;


  for (unsigned int j=0; j<y.size(); ++j)
  {
    for (unsigned int i=0; i<x.size(); ++i)
    {
//      std::cout << "particle interaction = (" << j << "," << i << ")" << std::endl;
      // taking care of relative and absolute difference
      // issues for when x[i] and y[j] are both small
      // or both large (see explanation in FmmTree member function solve
      // above)
      yPoint = y[j].getCoord();  xPoint = x[i].getCoord();

      std::cout << "yPoint = [" << yPoint[0] << "," << yPoint[1] << "," << yPoint[2] << "]" << std::endl;
      std::cout << "xPoint = [" << xPoint[0] << "," << xPoint[1] << "," << xPoint[2] << "]" << std::endl;

      lengthY = pow(yPoint[0],2.0) + pow(yPoint[1],2.0) + pow(yPoint[2],2.0);
      lengthY = pow(lengthY,0.5);
      lengthX = pow(xPoint[0],2.0) + pow(xPoint[1],2.0) + pow(xPoint[2],2.0);
      lengthX = pow(lengthX,0.5);
      maxXY = std::max(lengthY, lengthX);
      maxXYOne = std::max(1.0,maxXY);

      diffXY[0] = y[j].getX() - x[i].getX();
      diffXY[1] = y[j].getY() - x[i].getY();
      diffXY[2] = y[j].getZ() - x[i].getZ();
      if ( pow(pow(diffXY[0],2.0) + pow(diffXY[1],2.0) + pow(diffXY[2],2.0), 0.5)
                <= std::numeric_limits<double>::epsilon()*maxXYOne)
      {
        // Do nothing - y[j] and x[i] are the same point
    	// (up to machine epsilon)
    	// This happens when the target and source points are the same.
    	// Specifically, this happens when y[j] and x[i] are a target
    	// point and a source point in the same box (and the target
    	// and source points are the same).
    	// Also, physically it does not make sense for a particle y[j] to act
    	// on itself
      }
      else // target and source points y[j] and x[i] are not the same
      {
    	std::cout << "u[" << i << "] = " << u[i] << std::endl;
        potential_direct_grad_calculation
          = potential.directGrad(y[j].getCoord(),
      		                        x[i].getCoord());
        potential_direct_grad_calculation[0] = u[i] * potential_direct_grad_calculation[0];
        potential_direct_grad_calculation[1] = u[i] * potential_direct_grad_calculation[1];
        potential_direct_grad_calculation[2] = u[i] * potential_direct_grad_calculation[2];
        //        v[0][j] += potential_direct_grad_calculation[0].real();
        //        v[1][j] += potential_direct_grad_calculation[0].real();
        //        v[2][j] += potential_direct_grad_calculation[0].real();
        v[0][j] += potential_direct_grad_calculation[0];
        v[1][j] += potential_direct_grad_calculation[1];
        v[2][j] += potential_direct_grad_calculation[2];
        numOpsDirect++;
      }
    }
  }

  std::cout << "Completed Direct Solve Grad..." << "\n";
//  std::cout << "numOpsDirect = " << numOpsDirect << std::endl;

  return v;
}
