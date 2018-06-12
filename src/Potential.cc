/*
 * Potential.cc
 *
 *  Created on: Jul 8, 2016
 *      Author: dbpc
 *
 *
 * The Potential class contains information on the potential function: derivatives, formula for
 * direct calculations, S|S, S|R, and R|R translations.  The approximations for the Taylor series
 * for the fast multipole method (FMM) are up order |alpha| = 3.  The S|R transformation uses double
 * this order to fill its p x p transformation matrix.  Therefore, Potential.cc contains calculations
 * of derivatives to order |alpha| = 6.
 *
 * The heart of the Potential function class are the S|S, R|R, and S|R translation matrices in the member functions
 * getSS, getSR, getRR.  The matrices are explained in the write-up in Main.cc.  The matrices coded below are
 * the transpose of the matrices in Main.cc, since C++ arrays run faster through rows than down columns.
 *
 * On the upward pass, the fast multipole method starts working through the cells at the lowest level
 * (highest refinement level) to form far-field approximations (just the coefficients of the approximation)
 * for each source in the cell.  The approximations are centered at the center of the cell.
 * The upward pass consists of moving up to the next coarser refinement level and translating these approximations
 * for a cell from the center of the cell to the center of the parent cell at the next coarser level.
 * This process is repeated until reaching level l = 2 (with the index count starting on l = 0 - the third refinement
 * level).  Since the powers of the translated series depend only on the new center, only the coefficients need
 * to be determined.  Each new coefficient is obtained from the S|S translation matrix and is a linear combination of
 * the coefficients at the old center.
 *
 * The first part of the downward pass starts at level l = 2.  For cells at l = 2 the domain is made up of only
 * nearest neighbors and an interaction list.  The sources in the interaction list of each cell at level l = 2 are
 * S|R translated from the center of the interaction list cell to the center of the cell of interest.  The S|R translation
 * transforms the far-field s-expansion to an r-expansion (going from S(y,x_*) and ||x - x_*|| < ||y - x_*|| to R(y,y*)
 * and ||x - y_*|| > ||y - y_*||).  Again, only the coefficients of the transformed series need be determined (the center
 * y_* of the powers are known) and are combinations of the previous coefficients given by the S|R transformation matrix.
 *
 * The process is repeated at the next refinement level (next lower level).
 * A cell at level l = 3 is a child cell of a parent cell at level l = 2.  The coefficients held by
 * the parent cell account for the sources of the parent's interaction list.  Those sources are R|R translated
 * to the child cell.  The only sources unaccounted for are in the parent and the nearest neighbors of the parent.
 * Some of the children of these parents are far enough away to be in an interaction list of the cell and
 * the sources can be approximated by an s-expansion.  The others are nearest neighbors of the cell.
 * The nearest neighbor sources are calculated directly.  The far-field expansions of the sources in the interaction
 * list are S|R translated as in level l = 2.
 *
 * For a cell at level l = 4, the sources in the interaction list of the parent cells at the previous levels are
 * incorporated into the coefficients of the parent cell.  An R|R translation to the child cell at l = 4 takes
 * care of those sources.  All that is left are the sources in the nearest neighbors of the parent.  Some of the
 * children of these parents are part of the interaction list of the cell and some are nearest neighbors.
 * The potential function is calculated directly for the near neighbor sources and the far field
 * approximations for the interaction list cells are S|R translated to the cell.  The process is continued down
 * to the lowest level (highest refinement level).
 *
 * A standard ordering of the derivatives is necessary for the different member functions to communicate.
 * We set the standard here and order first by level of the derivative (1st derivative, 2nd derivative,...).
 * We then order the multi-index right-to-left from highest to lowest respectively.
 * The order for 0 <= |alpha <= 6 are given below.
 * 1)     (0,0,0),  // zero order derivative     ( 2 \\ 2) = 1
 * 2)     (0,0,1),  // first order derivative    ( 3 \\ 2) = 3
 * 3)     (0,1,0),
 * 4)     (1,0,0),
 * 5)     (0,0,2),  // second order derivative   ( 4 \\ 2) = 6
 * 6)     (0,1,1),
 * 7)     (0,2,0),
 * 8)     (1,0,1),
 * 9)     (1,1,0),
 * 10)    (2,0,0),
 * 11)    (0,0,3),  // third order derivative    ( 5 \\ 2) = 10
 * 12)    (0,1,2),
 * 13)    (0,2,1),
 * 14)    (0,3,0),
 * 15)    (1,0,2),
 * 16)    (1,1,1),
 * 17)    (1,2,0),
 * 18)    (2,0,1),
 * 19)    (2,1,0),
 * 20)    (3,0,0),
 * 21)    (0,0,4),  // fourth order derivative  ( 6 \\ 2) = 15
 * 22)    (0,1,3),
 * 23)    (0,2,2),
 * 24)    (0,3,1),
 * 25)    (0,4,0),
 * 26)    (1,0,3),
 * 27)    (1,1,2),
 * 28)    (1,2,1),
 * 29)    (1,3,0),
 * 30)    (2,0,2),
 * 31)    (2,1,1),
 * 32)    (2,2,0),
 * 33)    (3,0,1),
 * 34)    (3,1,0),
 * 35)    (4,0,0),
 * 36)	  (0,0,5),  // fifth order derivative  ( 7 \\ 2) = 21
 * 37)    (0,1,4),
 * 38)    (0,2,3),
 * 39)    (0,3,2),
 * 40)    (0,4,1),
 * 41)    (0,5,0),
 * 42)    (1,0,4),
 * 43)    (1,1,3),
 * 44)    (1,2,2),
 * 45)    (1,3,1),
 * 46)    (1,4,0),
 * 47)    (2,0,3),
 * 48)    (2,1,2),
 * 49)    (2,2,1),
 * 50)    (2,3,0),
 * 51)    (3,0,2),
 * 52)    (3,1,1),
 * 53)    (3,2,0),
 * 54)    (4,0,1),
 * 55)    (4,1,0),
 * 56)    (5,0,0),
 * 57)    (0,0,6),  // sixth order derivatives ( 8 \\ 2) = 28
 * 58)    (0,1,5),
 * 59)    (0,2,4),
 * 60)    (0,3,3),
 * 61)    (0,4,2),
 * 62)    (0,5,1),
 * 63)    (0,6,0),
 * 64)    (1,0,5),
 * 65)    (1,1,4),
 * 66)    (1,2,3),
 * 67)    (1,3,2),
 * 68)    (1,4,1),
 * 69)    (1,5,0),
 * 70)    (2,0,4),
 * 71)    (2,1,3),
 * 72)    (2,2,2),
 * 73)    (2,3,1),
 * 74)    (2,4,0),
 * 75)    (3,0,3),
 * 76)    (3,1,2),
 * 77)    (3,2,1),
 * 78)    (3,3,0),
 * 79)    (4,0,2),
 * 80)    (4,1,1),
 * 81)    (4,2,0),
 * 82)    (5,0,1),
 * 83)    (5,1,0),
 * 84)    (6,0,0),
 * 85)     (0,0,7),  // seventh order derivatives ( 9 \\ 2) = 36
 * 86)     (0,1,6),
 * 87)     (0,2,5),
 * 88)     (0,3,4),
 * 89)     (0,4,3),
 * 90)     (0,5,2),
 * 91)     (0,6,1),
 * 92)     (0,7,0),
 * 93)     (1,0,6),
 * 94)     (1,1,5),
 * 95)     (1,2,4),
 * 96)     (1,3,3),
 * 97)     (1,4,2),
 * 98)     (1,5,1),
 * 99)     (1,6,0),
 * 100)    (2,0,5),
 * 101)    (2,1,4),
 * 102)    (2,2,3),
 * 103)    (2,3,2),
 * 104)    (2,4,1),
 * 105)    (2,5,0),
 * 106)    (3,0,4),
 * 107)    (3,1,3),
 * 108)    (3,2,2),
 * 109)    (3,3,1),
 * 110)    (3,4,0),
 * 111)    (4,0,3),
 * 112)    (4,1,2),
 * 113)    (4,2,1),
 * 114)    (4,3,0),
 * 115)    (5,0,2),
 * 116)    (5,1,1),
 * 117)    (5,2,0),
 * 118)    (6,0,1),
 * 119)    (6,1,0),
 * 120)    (7,0,0),
 * 121)    (0,0,8),  // eighth order derivatives ( 10 \\ 2) = 45
 * 122)    (0,1,7),
 * 123)    (0,2,6),
 * 124)    (0,3,5),
 * 125)    (0,4,4),
 * 126)    (0,5,3),
 * 127)    (0,6,2),
 * 128)    (0,7,1),
 * 129)    (0,8,0),
 * 130)    (1,0,7),
 * 131)    (1,1,6),
 * 132)    (1,2,5)
 * 133)    (1,3,4),
 * 134)    (1,4,3),
 * 135)    (1,5,2),
 * 136)    (1,6,1),
 * 137)    (1,7,0),
 * 138)    (2,0,6),
 * 139)    (2,1,5),
 * 140)    (2,2,4),
 * 141)    (2,3,3),
 * 142)    (2,4,2),
 * 143)    (2,5,1),
 * 144)    (2,6,0),
 * 145)    (3,0,5),
 * 146)    (3,1,4),
 * 147)    (3,2,3),
 * 148)    (3,3,2),
 * 149)    (3,4,1),
 * 150)    (3,5,0),
 * 151)    (4,0,4),
 * 152)    (4,1,3),
 * 153)    (4,2,2),
 * 154)    (4,3,1),
 * 155)    (4,4,0),
 * 156)    (5,0,3),
 * 157)    (5,1,2),
 * 158)    (5,2,1),
 * 159)    (5,3,0),
 * 160)    (6,0,2),
 * 161)    (6,1,1),
 * 162)    (6,2,0),
 * 163)    (7,0,1),
 * 164)    (7,1,0),
 * 165)    (8,0,0),
 */

#include <vector>
#include <complex>
#include <cmath>
#include <iostream>
#include <functional>  // std::modulus

#include "Potential.h"


/**
 * Header Interface for Class Point
 *
class Potential
{
  public:
	int p;
	unsigned int order_of_approximation;
	// p gives the number of terms used to approximate the Taylor series used in the Fast Multipole Method (FMM)
	// In Cartesian coordinates the number of terms corresponds to the summation of the derivatives at each order
	// which can be calculated using combinations
	// p = (0+2)! / 2!((0+2)-2)! + (1+2)! / 2!((1+2)-2)! + (2+2)! / 2!((2+2)-2)! + (3+2)! / 2!((3+2)-2)!
	//   =        1                       3              +        6              +       10               = 20
	//        |alpha| = 0        +    |alpha| = 1        +    |alpha| = 2        +    |alpha| = 3
	// In Dehnen's code the number of terms p used to approximate the Taylor series do not match up to the
	// derivatives in the same way since Dehnen's formulas incorporate the condition that defines harmonic
	// functions \partial_x^2 + \partial_y^2 + \partial_z^2 = 0 and use a different coordinate system to simplify
	// this incorporation.  Therefore, less terms are required in Dehnen's coordinate system to obtain the same
	// number of Taylor series terms encountered in Cartesian coordinates without using the property of harmonic functions.

    int DEFAULT_P = 16;
    int DEFAULT_ORDER_OF_APPROXIMATION = 3;


    Potential() { this->p = DEFAULT_P;  this->order_of_approximation = DEFAULT_ORDER_OF_APPROXIMATION; };
	Potential(int p, int order_of_approximation)
	 :
     p(p),
     order_of_approximation(order_of_approximation)
	{};
	int getP() { return p;};
	void setP(int p) { this->p = p; };
	int getOrderApprox() {return order_of_approximation;};
	void setOrderApprox(int order_of_approximation) { this->order_of_approximation = oder_of_approximation; };
	std::vector<std::vector<std::complex<double> > > getSS(std::vector<double> from,
			                                               std::vector<double> to,
			                                               const std::vector<std::vector<std::complex<double> > > &sCoeff_old);
	std::vector<std::vector<std::complex<double> > > getRR(std::vector<double> from,
			                                               std::vector<double> to,
			                                               const std::vector<std::vector<std::complex<double> > > &rCoeff_old);
	std::vector<std::vector<std::complex<double> > > getSR(std::vector<double> from, std::vector<double> to,
			                                               const std::vector<std::vector<std::complex<double> > > &sCoeff_old);

	std::vector<std::vector<std::complex<double> > > getRCoeff(std::complex<double> xi, std::complex<double> xstar);
	std::vector<std::vector<std::complex<double> > >   getSCoeff(std::vector<double> xi, std::vector<double> xstar);

    std::vector<std::vector<std::complex<double> > > getRVector(std::complex<double> y, std::complex<double> xstar);
    std::vector<std::vector<std::complex<double> > >   getSVector(std::complex<double> y, std::complex<double> xstar);

    double                             direct(std::vector<double> yj, std::vector<double> xi);

};
*/



/**
 * Explanation of getRCoeff:
 *
 * The member function Potential::getRCoeff returns the coefficients of the near-field expansion (r-expansion).
 * The member function Potential::getRCoeff returns the coefficients \theta_n^{m} (x_b - z_A) of the
 * For RR and SR translations we need theta_n^m out to twice the order_of_approximation.
 * This is due to the formula used for the RR and SR translation.  The SS translation formula is different
 * and getSS member function returns square matrices whose size equals the order_of_approximation.
 * Ex: If the order_of_approximation = 3 (higher order derivatives used to approximate Taylor series is 3)
 * Then we would need theta_n^m to go out to 0 <= n <= 6 and -6 <= m <= 6 with -n <= m <= n.
 * Therefore the matrices returned by getRCoeff and getSVec are twice the size of the matrices
 * returned by getSCoeff and getRVec
 *
 * The Taylor-series in Cartesian coordinates for the near-field expansion indicate the coefficients
 * are the derivatives D^{\alpha} f(z_b - x_b)^{\alpha} / (\alpha !) and the powers of the series are
 * (x_a - z_b)^{\alpha} in
 *
 *  f(x_b - x_a) = sum_{|alpha| >= 0} [ ( (1 / (alpha !)) (x_b - z_b)^{alpha} ) D^{alpha} f(z_b - x_a) ]
 *
 *  The series converges when || x_b - z_b || < || z_b - x_a || (or || y - x^t_* || < || x^t_* - x_i ||=|| y - x_* ||)
 *
 *  In Dehnen's coordinate system the same series has the form
 *
 *     \frac{1}{|x_b - x_a|} = \frac{1}{|(z_B - x_a) - (z_B - x_b)|}
 *                            = \sum_{n=0}^{\infty} \sum_{m=-n}^{n} \gamma_n^{m*} (z_B - x_b) \theta_n^{m} (z_B - x_a)
 *
 *
 * Recall the formulas for \theta_n^{m} (\mathbf{r})
 *
 *  (0)   \theta_0^{0*} (r) = \frac{ 1 }{ \|\mathbf{r}\| }
 *
 *  (1)   \theta_n^{n}  (r) = (2n-1) \frac{2\xi}{\|\mathbf{r}\|^2} \theta_{n-1}^{(n-1)}   for 1 \leq n
 *
 *  (2)   \theta_n^{-n} (r) = (2n-1) \frac{2\eta}{\|\mathbf{r}\|^2} \theta_{n-1}^{-(n-1)}  for 1 \leq n
 *
 *  (3)   \theta_n^{m}  (r) = ( (2n-1) z \theta_{n-1}^{m} -  ( (n - 1)^2 - m^2 ) \theta_{n-2}^{m} ) / \|\mathbf{r}\|^2  for 0 < m \leq n
 *
 *  (4)   \theta_n^{-m} (r) = ( (2n-1) z \theta_{n-1}^{-m} -  ( (n - 1)^2 - m^2 ) \theta_{n-2}^{-m} / \|\mathbf{r}\|^2  for -n \leq m < 0
 *
 * The coefficients \theta_n^{m*} (x_b - z_A) may be complex as can be seen from the formulas above
 * Therefore, the vector of coefficients that is returned - ans - is a complex vector
 * Further, the coefficients correspond to a full matrix shown below
 *
 * --                                                                                                                  --
 * | \theta_0^0     \theta_1^{-1}     \theta_2^{-2}     \theta_3^{-3}     \theta_4^{-4}     \hdots      \theta_n^{-n}    |
 * |                                                                                                                     |
 * | \theta_1^1     \theta_1^0        \theta_2^{-1}     \theta_3^{-2}     \theta_4^{-3}                 \theta_n^{-n+1}  |
 * |                                                                                                                     |
 * | \theta_2^2     \theta_2^1        \theta_2^0        \theta_3^{-1}     \theta_4^{-2}                 \theta_n^{-n+2}  |
 * |                                                                                                                     |
 * | \theta_3^3     \theta_3^2        \theta_3^1        \theta_3^0        \theta_4^{-1}                 \vdots           |
 * |                                                                                                                     |
 * | \theta_4^4     \theta_4^3        \theta_4^2        \theta_4^1        \theta_4^0                                     |
 * |                                                                                                                     |
 * | \vdots                                                                            \ddots                            |
 * |                                                                                                                     |
 * | \theta_n^n                                                                                         \theta_n^0       |
 * |                                                                                                                     |
 * --                                                                                                                  --
 */
std::vector<std::vector<std::complex<double> > > Potential::getRCoeff(std::vector<double> xi,
		                                                std::vector<double> xstar)
{
  // xi = x_a is the source point
  // xstar = z_B is the center of the target cell
////  p = 16;
  // order_of_approximation is the highest order derivatives we use to approximate a Taylor series
  //
////  order_of_approximation = 3;
  //      norm_r_squared - ||x_a - z_B||^2
  double norm_r_squared = std::pow(xi[0]-xstar[0],2.0) + std::pow(xi[1]-xstar[1],2.0) + std::pow(xi[2]-xstar[2],2.0);

  // Recall the coordinate \eta = n =  -0.5 * (x - iy)
  //        the coordinate \xi  = e =   0.5 * (x + iy)
  //        the target particle is x_b
  //        the target cell center is z_B
  std::complex<double> eta_x_a, xi_x_a;
  double x_a_z_component;
  eta_x_a.real(-0.5*xi[0]);  eta_x_a.imag(0.5*xi[1]);
  xi_x_a.real(0.5*xi[0]);    xi_x_a.imag(0.5*xi[1]);
  x_a_z_component = xi[2];

  std::complex<double> eta_z_B, xi_z_B;
  double z_B_z_component;
  eta_z_B.real(-0.5*xstar[0]);  eta_z_B.imag(0.5*xstar[1]);
  xi_z_B.real(0.5*xstar[0]);    xi_z_B.imag(0.5*xstar[1]);
  z_B_z_component = xstar[2];

  std::complex<double> difference_eta = eta_z_B - eta_x_a;
  std::complex<double> difference_xi  = xi_z_B - xi_x_a;
  double difference_z = z_B_z_component - x_a_z_component;

  // Calculating the Matrix
  //
  // The (0,0) position upper left corner \theta_0^{0*} = 1 / |r|
  // The first column and first row can be calculated using the recursive formulas (1) and (2)
  // Each diagonal is calculated using formulas (3) and (4)
  // Since we are using p = 16 terms (up to and including 3rd order derivatives), the matrix above will a 4x4 matrix
  // The count starts on 0, so we can set the matrix size using the highest order derivatives we are using to
  // approximate the series.
  // We will want to incorporate the order of approximation that we will be using as a class variable
  // May be a faster way using columns with C++

  // For RR, SR, and SS translations we need the square theta_n^m matrix the size of the order_of_approximation.
  std::vector<std::vector<std::complex<double> > >
     ans(order_of_approximation+1, std::vector<std::complex<double> >(order_of_approximation+1));

  // Set theta_0^{0*} to 1 / |r|
  //  (0)   \theta_0^{0*} (r) = \frac{ 1 }{ \|\mathbf{r}\| }
  ans[0][0] = 1.0 / std::pow(norm_r_squared, 0.5);
  //std::cout << "Set (0,0) position of getSVector matrix" << std::endl;

  // Calculate the first row and column
  //  (1)   \theta_n^{n}  (r) = (2n-1) \frac{2\xi}{\|\mathbf{r}\|^2} \theta_{n-1}^{(n-1)}   for 1 \leq n
  //  (2)   \theta_n^{-n} (r) = (2n-1) \frac{2\eta}{\|\mathbf{r}\|^2} \theta_{n-1}^{-(n-1)}  for 1 \leq n
  for (unsigned int i=1; i<=order_of_approximation; ++i)
  {
    ans[i][0] = ans[i-1][0] * (2.0*double(i) - 1.0) * 2.0 * difference_xi / norm_r_squared;
    ans[0][i] = ans[0][i-1] * (2.0*double(i) - 1.0) * 2.0 * difference_eta / norm_r_squared;
  }
  //std::cout << "Set 1st row and column of getSVector matrix" << std::endl;

  // Calculate each diagonal ********************************************* //
  //  (3)   \theta_n^{m}  (r) = ( (2n-1) z \theta_{n-1}^{m} -  ( (n - 1)^2 - m^2 ) \theta_{n-2}^{m} ) / \|\mathbf{r}\|^2  for 0 < m \leq n
  //  (4)   \theta_n^{-m} (r) = ( (2n-1) z \theta_{n-1}^{-m} -  ( (n - 1)^2 - m^2 ) \theta_{n-2}^{-m} / \|\mathbf{r}\|^2  for -n \leq m < 0

  // Main Diagonal - m = 0
  ans[1][1] = (2*1-1) * difference_z *ans[0][0] / norm_r_squared;
  //std::cout << "Set (1,1) position of getSVector matrix" << std::endl;
  for (unsigned int i=2; i<=order_of_approximation; ++i)
    ans[i][i] = ( (2.*double(i)-1.) * difference_z * ans[i-1][i-1] - (std::pow(double(i)-1,2.0)) * ans[i-2][i-2] ) / norm_r_squared;
  //std::cout << "Set main diagonal of getSVector matrix" << std::endl;

  // lower triangular part
  // i = 1: [i+j][j] = [2][1], [3][2], [4][3], [5][4] - 1st lower diagonal
  // i = 2: [i+j][j] = [3][1], [4][2], [5][3], [6][4] - 2nd lower diagonal
  // upper triangular part
  // i = 1: [j][i+j] = [1][2], [2][3], [3][4], [4][5]
  // i = 2: [j][i+j] = [1][3], [2][4], [3][5], [4][6]
  // Note: from matrix above - i coincides with m
  //                         - i+j coincides with n
  //       used in recursive formulas below
  for (unsigned int i=1; i<=order_of_approximation; ++i)  // main diagonal - i = 0
    for (unsigned int j=1; j<=order_of_approximation-i; ++j)
    {
      if (j == 1) // 2nd column or row - second term in recursive formula is zero
      {
        //std::cout << "Nested Loop if statement lower triangular part iteration i = " << i << ", j = " << j << std::endl;
        // lower triangular part
        ans[i+j][j] = ( (2.*double(i+j)-1.) * difference_z * ans[i+j-1][j-1] )
        		               / norm_r_squared ;
        //std::cout << "Nested Loop if statement upper triangular part iteration i = " << i << ", j = " << j << std::endl;
        // upper triangular part
        ans[j][i+j] = ( (2.*double(i+j)-1.) * difference_z * ans[j-1][i+j-1] )
        		               / norm_r_squared ;
      }
      else
      {
        //std::cout << "Nested Loop else statement lower triangular part iteration i = " << i << ", j = " << j << std::endl;
        // lower triangular part
        ans[i+j][j] = ( (2.*double(i+j)-1.) * difference_z * ans[i+j-1][j-1]
      		              - ( std::pow(double(i+j-1),2.0) - std::pow(double(i),2.0) ) * ans[i+j-2][j-2]
        		      )
        		         / norm_r_squared;
        //std::cout << "Nested Loop else statement upper triangular part iteration i = " << i << ", j = " << j << std::endl;
        // upper triangular part
        ans[j][i+j] = ( (2.*double(i+j)-1.) * difference_z * ans[j-1][i+j-1]
        		          - ( std::pow(double(i+j-1),2.0) - std::pow(double(i),2.0) ) * ans[j-2][i+j-2]
        	          )
                         / norm_r_squared;

      }
    }

  return ans;

}




/**
 * Explanation of getSCoeff

 * The far-field coefficients of the s-expansion are the powers
 *
 *    (x_a - z_a)^{alpha} / (\alpha !)  =  (xi - xstar)^{alpha} / (alpha !)
 *
 * where z_a is the center of the source cell and x_a is the source particle (point).
 * Recall the series has the form
 *
 *    psi(x - y)  =  sum_{alpha >= 0}  frac{ (x - x_*)^{alpha} }{ alpha! }  D^{alpha} psi(x_* - y)
 *
 * and the powers of the series are the derivatives
 *
 *    D^{alpha} psi(x_* - y)  =  D^{\alpha} f(z_a - x_b)  =  D^{alpha} psi(xstar - yj).
 *
 * In Dehnen's paper the far-field expansion is
 *
 *     \frac{1}{|r-x|} = \sum_{n=0}^{\infty} \gamma_n^{m*} (x) \theta_{n}^{m} (r)
 *
 * The coefficients are \gamma_n^{m*} (x) seen in
 *
 *     \frac{1}{|x_a - x_b|} = \frac{1}{|(x_a - z_A) - (x_b - z_A)|}
 *                           = \sum_{n=0}^{\infty} \sum_{m=-n}^{n} \gamma_n^{m*} (x_a - z_A) \theta_n^{m*} (x_b - z_A)
 *
 *
 * In Potential::getSCoeff
 *   xi    - is the source point x_a,
 *   xstar - is the center of the source cell z_a
 *
 *
 */
std::vector<std::vector<std::complex<double> > >
                                   Potential::getSCoeff(std::vector<double> xi,
		                                                std::vector<double> xstar)
{
  // Here xi             - x_a is the source particle
  //      xstar          - z_A is the source cell center

  // terms used to approximation the Taylor series
////  p = 16;
  // highest order of derivative used to approximate Taylor series, to which p corresponds
////  order_of_approximation = 3;
  // matrix of coefficients that will be returned
  // note - Taylor series includes zeroth derivative - therefore adding 1 to highest order
  std::vector<std::vector<std::complex<double> > > ans(order_of_approximation+1, std::vector<std::complex<double> >(order_of_approximation+1));
  // norm_r_squared - ||x_a - z_A||^2
  double norm_r_squared = std::pow(xi[0]-xstar[0],2.0) + std::pow(xi[1]-xstar[1],2.0) + std::pow(xi[2]-xstar[2],2.0);

  // **********************************************************************************************************
  // Setting the coordinates xi and eta in coordinate system xi, eta, z **************************************
  // **********************************************************************************************************
  // **********************************************************************************************************
  // Recall the coordinate \eta = n =  -0.5 * (x - iy)
  //        the coordinate \xi  = e =   0.5 * (x + iy)
  //        the target particle is x_b
  //        the target cell center is z_B
  std::complex<double> eta_x_a, xi_x_a;
  double x_a_z_component;
  eta_x_a.real(-0.5*xi[0]);  eta_x_a.imag(0.5*xi[1]);
  xi_x_a.real(0.5*xi[0]);    xi_x_a.imag(0.5*xi[1]);
  x_a_z_component = xi[2];

  std::complex<double> eta_z_A, xi_z_A;
  double z_A_z_component;
  eta_z_A.real(-0.5*xstar[0]);  eta_z_A.imag(0.5*xstar[1]);
  xi_z_A.real(0.5*xstar[0]);    xi_z_A.imag(0.5*xstar[1]);
  z_A_z_component = xstar[2];

  std::complex<double> difference_eta = eta_x_a - eta_z_A;
  std::complex<double> difference_xi  = xi_x_a - xi_z_A;
  double difference_z = x_a_z_component - z_A_z_component;

  // **********************************************************************************************************
  // Recalling the information on the formulas for coefficients gamma_n^{m*} **********************************
  // **********************************************************************************************************
  // **********************************************************************************************************
  // We use p = 16 to obtain all derivatives up to and including order 3.
  // The member function Potential::getSCoeff returns the coefficients \gamma_n^{m*} (x_a - z_A) of the
  // far-field expansion \frac{1}{|r-x|} = \sum_{n=0}^{\infty} \gamma_n^{m*} (x) \theta_{n}^{m} (r)
  // Then
  //                     \frac{1}{|x_a - x_b|} = \frac{1}{|(x_a - z_A) - (x_b - z_A)|}
  //                                           = \sum_{n=0}^{\infty} \sum_{m=-n}^{n} \gamma_n^{m*} (x_a - z_A) \theta_n^{m*} (x_b - z_A)
  //
  // Recall the formulas for \gamma_n^{m*} (\mathbf{r})
  //
  //  (0)   \gamma_0^{0*}  (r) = 1
  //
  //  (1)   \gamma_n^{n*}  (r) = (-1) \frac{\eta}{n} \gamma_{n-1}^{(n-1)*}   for 1 \leq n
  //
  //  (2)   \gamma_n^{-n*} (r) = (-1) \frac{\xi}{n}  \gamma_{n-1}^{-(n-1)*}  for 1 \leq n
  //
  //  (3)   \gamma_n^{m*}  (r) = \frac{(2n-1)z}{n^2 - m^2} \gamma_{n-1}^{m*} -  \frac{|r|^2}{n^2 - m^2} \gamma_{n-2}^{m*}   for m \geq 0
  //
  //  (4)   \gamma_n^{-m*} (r) = \frac{(2n-1)z}{n^2 - m^2} \gamma_{n-1}^{-m*} - \frac{|r|^2}{n^2 - m^2} \gamma_{n-2}^{-m*}  for m \geq 0

  // The coefficients \gamma_n^{m*} (x_a - z_A) may be complex as can be seen from the formulas above
  // Therefore, the vector of coefficients that is returned - ans - is a complex vector
  // Further, the coefficients correspond to a full matrix shown below
  //
  // --                                                                                                                  --
  // | \gamma_0^0*    \gamma_1^{-1*}    \gamma_2^{-2*}    \gamma_3^{-3*}    \gamma_4^{-4*}    \hdots      \gamma_n^{-n*}   |
  // |                                                                                                                     |
  // | \gamma_1^1*    \gamma_1^0*       \gamma_2^{-1*}    \gamma_3^{-2*}    \gamma_4^{-3*}                \gamma_n^{-n+1*} |
  // |                                                                                                                     |
  // | \gamma_2^2*    \gamma_2^1*       \gamma_2^0*       \gamma_3^{-1*}    \gamma_4^{-2*}                \gamma_n^{-n+2*} |
  // |                                                                                                                     |
  // | \gamma_3^3*    \gamma_3^2*       \gamma_3^1*       \gamma_3^0*       \gamma_4^{-1*}                \vdots           |
  // |                                                                                                                     |
  // | \gamma_4^4*    \gamma_4^3*       \gamma_4^2*       \gamma_4^1*       \gamma_4^0*                                    |
  // |                                                                                                                     |
  // | \vdots                                                                                 \ddots                       |
  // |                                                                                                                     |
  // | \gamma_n^n*                                                                                        \gamma_n^0*      |
  // |                                                                                                                     |
  // --                                                                                                                  --
  //
  // **********************************************************************************************************
  // Calculating the Matrix ***********************************************************************************
  // **********************************************************************************************************
  // **********************************************************************************************************
  //
  // The (0,0) position upper left corner \gamma_0^{0*} = 1
  // The first column and first row can be calculated using the recursive formulas (1) and (2)
  // Each diagonal is calculated using formulas (3) and (4)
  // Since we are using p = 16 terms (up to and including 3rd order derivatives), the matrix above will a 4x4 matrix
  // The count starts on 0, so we can set the matrix size using the highest order derivatives we are using to
  // approximate the series.
  // We will want to incorporate the order of approximation that we will be using as a class variable
  // May be a faster way using columns with C++


  // Set gamma_0^{0*} to 1
  //  (0)   \gamma_0^{0*}  (r) = 1
  ans[0][0] = 1.0;
  //std::cout << "Set (0,0) position of getSCoeff matrix" << std::endl;

  // Calculate the first row and column
  //  (1)   \gamma_n^{n*}  (r) = (-1) \frac{\eta}{n} \gamma_{n-1}^{(n-1)*}   for 1 \leq n
  //  (2)   \gamma_n^{-n*} (r) = (-1) \frac{\xi}{n}  \gamma_{n-1}^{-(n-1)*}  for 1 \leq n
  for (unsigned int i=1; i<=order_of_approximation; ++i)
  {
	ans[i][0] = ans[i-1][0] * (-1.0) * difference_eta / double(i);
    ans[0][i] = ans[0][i-1] * (-1.0) * difference_xi / double(i);
 //   std::cout << "ans[" << 0 << "][" << i << "] = " << ans[0][i] << std::endl;
  }
  //std::cout << "Set 1st row and column of getSCoeff matrix" << std::endl;

  // Calculate each diagonal ********************************************* //
  //  (3)   \gamma_n^{m*}  (r) = \frac{(2n-1)z}{n^2 - m^2} \gamma_{n-1}^{m*} -  \frac{|r|^2}{n^2 - m^2} \gamma_{n-2}^{m*}   for m \geq 0
  //  (4)   \gamma_n^{-m*} (r) = \frac{(2n-1)z}{n^2 - m^2} \gamma_{n-1}^{-m*} - \frac{|r|^2}{n^2 - m^2} \gamma_{n-2}^{-m*}  for m \geq 0

  // Main Diagonal - here m = 0 so power (n^2 - m^2) = n^2 in formula
  ans[1][1] = (2*1-1) * difference_z *ans[0][0] / pow(double(1),2.0);
  //std::cout << "Set (1,1) position of getSCoeff matrix" << std::endl;
  for (unsigned int i=2; i<=order_of_approximation; ++i)
    ans[i][i] = ( (2.*double(i)-1.) * difference_z * ans[i-1][i-1] - norm_r_squared * ans[i-2][i-2] )
                                 / std::pow(double(i),2.0);
  //std::cout << "Set main diagonal of getSCoeff matrix" << std::endl;

  // Stepping Through Nested Loop Below for Matrix Coefficients Below *************
  //
  // Note: i and j start their count on 1 since the first row and column and main diagonal are complete
  //       For an (n,n) matrix of coefficients there are n upper and lower diagonals
  //         Above, we have set the last upper and lower diagonal that is the lower left hand and upper
  //         right hand corner of the matrix and falls into the first row and column.
  //         Therefore, the bound order_of_approxiation on i ensures that i+j runs from 1st upper and lower
  //         off diagonal to the second to last upper and lower off diagonal
  //       The bound on j ensures that the inner nested loop does not go beyond the matrix size
  // Example: matrix is 4x4 according to order of approximation n = 3
  //          The coefficients on the main diagonal and 1st row and column have been taken care of so far
  //                           --                             --    --                                                          --
  //                           |   *       *       *       *    |   |   *       *               *                    *            |
  //                           |                                |   |                                                             |
  //                           |   *       *     [1][2]  [1][3] |   |   *       *             \gamma_2^{-1*}       \gamma_3^{-2*} |
  //  Matrix of Coefficients = |                                | = |                                                             |
  //                           |   *     [2][1]    *     [2][3] |   |   *     \gamma_2^{1*}     *                  \gamma_3^{-1*} |
  //                           |                                |   |                                                             |
  //                           |   *     [3][1]  [3][2]    *    |   |   *     \gamma_3^{2*}   \gamma_3^{1*}          *            |
  //                           |                                |   |                                                             |
  //                           --                             --    --                                                          --
  // Lower Triangular Part - \gamma_n^{m*} ***************************************************
  // i = 1:
  //   j = 1:  [i+j][j] = [2][1] = gamma_2^{1*} = gamma_{i+j}^{i*}         - 1st lower diagonal
  //   j = 2:  [i+j][j] =       [3][2] = gamma_3^{1*} = gamma_{i+j}^{i*}   - 1st lower diagonal
  // Stop: i+j = 1+2 = 3 = order_of_approximation
  // i = 2:
  //   j = 1:  [i+j][j] = [3][1] = gamma_3^{2*} = gamma_{i+j}^{i*}         - 2nd lower diagonal
  // Stop: i+j = 2+1 = 3 = order_of_approximation
  //
  // NOTE: from matrix and work above - i coincides with m
  //                                    and i+j coincides with n
  //                                    used in recursive formulas below
  //
  // Upper Triangular Part - \gamma_n^{m*} ***************************************************
  // i = 1:
  //   j = 1:  [j][i+j] = [1][2]= gamma_2^{-1*} = gamma_{i+j}^{-i*}          - 1st upper diagonal
  //   j = 2:  [j][i+j] =       [2][3]= gamma_{3}^{-1*} = gamma_{i+j}^{-i*}  - 1st upper diagonal
  // Stop: i+j = 1+2 = 3 = order_of_approximation
  // i = 2:
  //   j = 1:  [j][i+j] = [1][3]= gamma_3^{-2*} = gamma_{i+j}^{-i*}          - 2nd upper diagonal
  // Stop: i+j = 2+1 = 3 = order_of_approximation
  //
  // NOTE: from matrix and work above - i+j coincides with n
  //                                    i   coincides with m
  //                                    used in recursive formulas below
  //       Also, i indexes the diagonal begin worked on
  //
  for (unsigned int i=1; i<=order_of_approximation; ++i)
  {
	// * The diagonals of our matrix correspond to the recursive formulas for \gamma_n^{(m)*}
	// * Therefore, our nested i,j loop calculates the elements for each diagonal
	// * i locates the diagonal being worked on i = 1 is the first diagonal below the main diagonal
	//                                          i = 2 is the second diagonal below the main diagonal
	// * The main diagonal has already been calculated above
	//   Therefore we need i > 0 (i to start at 1) since i = 0 corresponds to the main diagonal (diagonal zero)
	//   Also, i <= order_of_approximation since we have order_of_approximation diagonals below the main diagonal
	// * Our first diagonal to calculate will be just below the main diagonal with i = 1
	// * j determines the column or row of the element we are calculating
	//   Therefore, we need j > 0 (j to start on 1) since column 0 and row 0 have already been calculated
	//   Also we need j <= order_of_approximation to not go beyond the number of rows of the matrix
	// * We increment by j down the diagonal we are working on: column i+j and row j
	//   Therefore, i+j <= order_of_approximation to not go beyond the number of columns of the matrix
	//   and cause a segmentation fault
	// * Further, once we set i, we know that we are working on the upper and lower diagonals i
	// * Therefore we are starting from position (i,0) or (0,i), incrementing j by 1 unit,
	//   and moving down that diagonal to (i+j,j) or (j,i+j)
	// * If we were working on the main diagonal, then we would know that there are n elements
	//   if the matrix size is n.  Each shift to a lower or upper diagonal has one unit less than
	//   the previous diagonal.  That is, if i = 1 then the number of elements in the diagonal is
	//   n_d_e = n - 1 = n - i, if i = 2 then the number of diagonal elements is n_d_e = n - 2 = n - i,
	//   and so on.  This number n_d_e tells us the upper bound for j = n - i
	//   j can therefore not have a fixed upper bound, and we must offset j (the number of diagonal elements)
	//   by i as shown below
	//  (3)   \gamma_n^{m*}  (r) = \frac{(2n-1)z}{n^2 - m^2} \gamma_{n-1}^{m*} -  \frac{|r|^2}{n^2 - m^2} \gamma_{n-2}^{m*}   for m \geq 0
	//  (4)   \gamma_n^{-m*} (r) = \frac{(2n-1)z}{n^2 - m^2} \gamma_{n-1}^{-m*} - \frac{|r|^2}{n^2 - m^2} \gamma_{n-2}^{-m*}  for m \geq 0
	for (unsigned int j=1; j<=order_of_approximation-i; ++j) //(order_of_approximation + 1)-i; ++j)
    {
      if (j == 1) // 2nd column or row - In 2nd column or row, second term in recursive formula is zero
    	          //                     and term would be located outside matrix index, throwing segmentation fault
    	          //                     Therefore we omit the second term
      {
        //std::cout << "Nested Loop if statement lower triangular part iteration i = " << i << ", j = " << j << std::endl;
        // Lower Triangular Part: i = m and i+j = n
        ans[i+j][j] = ( (2.*double(i+j)-1.) * difference_z * ans[i+j-1][j-1] )
        		               / ( std::pow(double(i+j),2.0) - std::pow(double(i),2.0) );
        //std::cout << "Nested Loop if statement upper triangular part iteration i = " << i << ", j = " << j << std::endl;
        // Upper Triangular Part: -i = m and i+j = n
        ans[j][i+j] = ( (2.*double(i+j)-1.) * difference_z * ans[j-1][i+j-1] )
        		               / ( std::pow(double(i+j),2.0) - std::pow(double(i),2.0) );
      }
      else
      {
        //std::cout << "Nested Loop else statement lower triangular part iteration i = " << i << ", j = " << j << std::endl;
        // Lower Triangular Part: i = m and i+j = n
        ans[i+j][j] = ( ((2.*double(i+j)-1.) * difference_z * ans[i+j-1][j-1]) - (norm_r_squared * ans[i+j-2][j-2]) )
                                              / ( std::pow(double(i+j),2.0) - std::pow(double(i),2.0) );
        //std::cout << "Nested Loop else statement upper triangular part iteration i = " << i << ", j = " << j << std::endl;
        // Upper Triangular Part: -i = m and i+j = n
        ans[j][i+j] = ( ((2.*double(i+j)-1.) * difference_z * ans[j-1][i+j-1]) - (norm_r_squared * ans[j-2][i+j-2]) )
                                              / ( std::pow(double(i+j),2.0) - std::pow(double(i),2.0) );

      }
    }
    //std::cout << "Outer Loop iteration i = " << i << std::endl;
  }
  //std::cout << "Finish Nested Loop" << std::endl;

  return ans;

}



/**
 * Explanation of getSR
 *
 * We are transforming the power series coefficient matrix sCoeff_old = \gamma_n^{(m)*} (x) = \gamma_n^{m*} (x_a - z_A^{(p)}}) matrix
 * to a theta_n^m coefficient matrix using the formula
 *
 *   \theta_n^{m} (x+y) = \sum_{k=0}^{order_of_approximation - n}
 *                                   \sum_{l=-k}^{k} \gamma_{k}^{(l)*} (-x)  theta_{n+k}^{m+l} (y)
 *                      = \sum_{k=0}^{order_of_approximation - n}
 *                                   \sum_{l=-k}^{k} (-1)^n \gamma_{k}^{(l)*} (x)  theta_{n+k}^{m+l} (y)
 *
 * Recall: \gamma_m^{(n)*} (-x) = (-1)^n \gamma_n^{(m)*} (x)
 *
 * The steps we take below are
 *
 * 1. We have sCoeff_old matrix of containing \gamma_n^{(m)*} (x)
 * 2. We evaluate theta_n^m(y) using getRCoeff member function
 * 3. We form theta_n^m(x+y) using
 *
 *   \theta_n^{m} (x+y) = \sum_{k=0}^{order_of_approximation - n}
 *                                   \sum_{l=-k}^{k} (-1)^n \gamma_{k}^{(l)*} (x)  theta_{n+k}^{m+l} (y)
 *
 *   where y = z_B^{(p)} - z_A^{(p)} = to - from
 *         x = x_a - z_A^{(p)}
 *   and
 *           x + y     = x_a - z_B^{(p)}
 *           x_a       = source point
 *           x_b       = target point
 *           z_A^{(p)} = parent source cell
 *           z_B^{(p)} = parent target cell
 *
 *   seen in
 *
 *   \mu_a \frac{1}{|x_b - x_a|} = \mu_a \frac{1}{|x_b - z_B + z_B - x_a|}
 *                               = \mu_a \frac{1}{|(z_B - x_a) - (z_B - x_b)|
 *                               = \mu_a \sum_{n=0}^{\infty} \sum_{m=-n}^{n} \gamma_n^{m*} (z_B - x_b) \theta_n^m(z_B - x_a)
 *                               = \mu_a \sum_{n=0}^{\infty} \sum_{m=-n}^{n} \gamma_n^{m*} (z_B - x_b)
 *                                                                                       \theta_n^m(z_B - z_A - (x_a - z_A))
 *                               = \mu_a \sum_{n=0}^{\infty} \sum_{m=-n}^{n} \gamma_n^{m*} (z_B - x_b)
 *                                                                  \sum_{k=0}^{\infty} \sum_{l=-k}^{k} \gamma_{k}^{l*}(x_a - z_A)
 *                                                                                                         theta_{n+k}^{m+l} (z_B - z_A)
 *
 */
std::vector<std::vector<std::complex<double> > > Potential::getSR(std::vector<double> from, std::vector<double> to,
		                                                          const std::vector<std::vector<std::complex<double> > > &sCoeff_old)
{
	  // from - refers to the center where the old coefficients were evaluated
	  //        z_A^{(p)} - the parent cell center z_A^{(p)} of the source point
	  // to   - refers to the center where the new coefficients will be evaluated
	  //        z_B^{(p)} - the parent cell center z_B^{(p)} of the target point
	  // t    - refers to y = from - to = z_A^{(p)} - z_B^{(p)}
	  // Here we have sCoeff_old = \gamma_n^{(m)*} (x) = \gamma_n^{m*} (x_a - z_A^{(p)}}) matrix

	  // need to assert from and to have 3 components

	  // terms used to approximate the Taylor series
////	  p = 16;
	  // highest order of derivative used to approximate Taylor series, to which p corresponds
////	  order_of_approximation = 3;

//	  for (unsigned int i=0; i<4; ++i)
//	    for (unsigned int j=0; j<4; ++j)
//	      std::cout << "sCoeff_old[" << i << "][" << j << "] = "
//	                << sCoeff_old[i][j] << std::endl;

	  // note - Taylor series includes zeroth derivative - therefore adding 1 to highest order derivative
	  std::vector<std::vector<std::complex<double> > > ans(order_of_approximation+1, std::vector<std::complex<double> >(order_of_approximation+1));

//	  // Here x = from and x+y = to and t = to - from = y
//	  std::vector<double>  t(3);
//	  t[0] = from[0] - to[0];  t[1] = from[1] - to[1]; t[2] = from[2] - to[2];


	  // Here we build the \theta_n^{(m)} (x) = \theta_n^{m} (z_B - z_A) matrix
	  // Recall the matrix dimensions
	  // std::vector<std::vector<std::complex<double> > >
	  //                  theta_y_matrix(order_of_approximation+1, std::vector<std::complex<double> >(order_of_approximation+1));
	  // Recall
	  // Here we are interested in \theta_n^m(z_B - x_a) = \theta_n^m(z_B - z_A - (x_a - z_A))
	  //                                                 = \sum_{k=0}^{\infty} \sum_{l=-k}^{k} \gamma_{k}^{l*}(x_a - z_A)
	  //                                                                                              \theta_{n+k}^{m+l} (z_B - z_A)
	  // Recall getRCoeff(xi, xstar) performs its operations using the difference (xstar - xi).
	  // Therefore, we let from = z_B = xstar and to = z_A = xi to
	  // We call getRCoeff(zA, zB) = getRCoeff(from, to)
	  std::vector<std::vector<std::complex<double> > > theta_y_matrix = this->getRCoeff(from, to);

//	  for (unsigned int i=0; i<4; ++i)
//	    for (unsigned int j=0; j<4; ++j)
//	      std::cout << "theta_y_matrix[" << i << "][" << j << "] = "
//	                << theta_y_matrix[i][j] << std::endl;


	  for (unsigned int j = 0; j<=order_of_approximation; ++j)
	    for (unsigned int i = 0; i<=order_of_approximation; ++i)
	      ans[i][j] = 0.0;

	  // Generating the matrix \gamma_n^{(m)*} (-y) = \theta_n^{m} (x+y) = theta_n^m (x_b - z_B^{(p)}})
	  // using the formula
	  //
	  //   \theta_n^{m} (x+y) = \sum_{k=0}^{order_of_approximation - n}
	  //                                  \sum_{l=-k}^{k} \gamma_{k}^{(l)*} (-y)  theta_{n+k}^{m+l} (x)
	  //
	  for (unsigned int n = 0; n <= order_of_approximation; ++n)
	    for (int m = -n; m<=int(n); ++m)
	    {
	      if (int(n) >= std::abs(m)) // assemble an element, otherwise n < m implies element is zero
	      {
	    	// theta_n^m matrix has entries only up to order_of_approximation
	    	// Therefore, n+k can only be up to order_of_approximation => n+k <= order_of_approximation
	    	//                    => k <= order_of_approximation - n
	    	// or else we will be looking for theta_{n+k}^{m+l} that are not there (indexing outside of the theta matrix)
	        for (unsigned int k=0; k<=(order_of_approximation - n); ++k)
	          for (int l=-k; l<=int(k); ++l)
	          {
	            if ( m >= 0 ) // setting lower triangular part of theta matrix ans -> theta_n^{(m)} (x+y)
	            {
	              if ( int(n+k) >= std::abs(m+l) )  // theta_{n+k}^{(m+l)*} not zero
	              {
	                if ( (m+l) >= 0)  // indexing lower triangular part of theta matrix theta_y_matrix -> theta_{n+k}^{(m+l)} (y)
	                {
	                  if ( l >= 0 )  // indexing lower triangular part of gamma matrix sCoeff_old -> gamma_k^{l*} (x)
	                	             // m >=0 and (m+l) >= 0 and l >= 0
	                	             // implies that m >= -l
	                  {
	            	    // ans[n][n-m] is the location of theta_n^{(m)*} (x+y)
	          	        // rCoeff_old[n+k][(n+k)-(m+l)] is the location of theta_{n+k}^{(m+l)} (x)
	          	        // theta_y_matrix[k][k-l] is the location of gamma_{k}^{(l)*} (-y)
	                    ans[n][n-m] += sCoeff_old[k][k-l] * theta_y_matrix[n+k][(n+k)-(m+l)];
	                  }
	                  else // l < 0 // indexing upper triangular part of gamma matrix sCoeff_old -> gamma_k^{l*} (x)
	                	            // (m + l) >= 0 and l < 0 => m >= -l > 0
	                  {
	                    // EX: ans[n][n-m] += rCoeff_old[n+k][(n+k)-(m+l)] * gamma_neg_y_matrix[k+l][k];
	                	//     n = 2, m = 1
	                	//       k = 1, l = -1
	                    //         theta_2^1(x+y) += theta_3^0(x) * gamma_1^{-1}(-y)
	                	//         theta_new[2][1] += theta_y_matrix[3][3] * gamma_[0][1]
	                	//         ans[n][2-1] += theta_y_matrix[n+k][n+k-(m+l)] * sCoeff_old[k+l][k]
	                  	//     n = 2, m = 2
	                  	//       k = 1, l = -1
	                      //         theta_2^2(x+y) += theta_3^1(y) * gamma_1^{-1}(-y)
	                  	//           theta_new[2][1] += theta_y_matrix[3][2] * sCoeff_old[0][1]
	                  	//           ans[n][2-1] += theta_y_matrix[n+k][n+k-(m+l)] * sCoeff_old[k+l][k]
	                    ans[n][n-m] += sCoeff_old[k+l][k] * theta_y_matrix[n+k][(n+k)-(m+l)];
	                  }
	                }
	                else // (m+l) < 0 and indexing upper triangular part of theta matrix theta_y_matrix -> theta_{n+k}^{(m+l)*} (y)
	                {
	                  if ( l >= 0 )  // indexing lower triangular part of gamma matrix sCoeff_old -> gamma_k^{l*} (x)
	                  {
	                	// do nothing: l >= 0 and (m+l) < 0  =>  m < -l < 0
	                	// m has to be negative and we are doing case m >= 0 (outer if condition)
	                  }
	                  else // l < 0 // indexing upper triangular part of gamma matrix sCoeff_old -> gamma_k^{l*} (x)
	                  {
	                    // EX: ans[n][n-m] += theta_y_matrix[n+k][(n+k)+(m+l)] * sCoeff_old[k+l][k];
	                  	//     n = 2, m = 1
	                  	//       k = 1, l = -1
	                    //         theta_2^1(x+y) += theta_3^{0}(y) * gamma_1^{(-1)*}(x)
	                  	//         theta_new[2][1] += theta_y_matrix[3][3] * sCoeff_old[0][1]
	                  	//         ans[n][n-m] += theta_y_matrix[n+k][n+k-(m+l)] * sCoeff_old[k+l][k]
	                    	//     n = 2, m = 0
	                    	//       k = 1, l = -1
	                        //         theta_2^0(x+y) += theta_3^{-1}(y) * gamma_1^{-1}(x)
	                    	//         theta_new[2][2] += theta_y_mamtrix[2][3] * sCoeff_old[0][1]
	                    	//         ans[n][2-1] += theta_y_mamtrix[n+k][n+k-(m+l)] * sCoeff_old[k+l][k]
	                    ans[n][n-m] += sCoeff_old[k+l][k] * theta_y_matrix[(n+k)+(m+l)][n+k];
	                  }
	                }
	              }
	              else // (n-k) < std::abs(m-l)
	              {
	                //  no need to do anything since gamma_{n-k}^{(m-l)*} = zero
	              }
	            }
	            else // m < 0 - setting lower triangular part of gamma matrix ans -> gamma_n^{(m)*} (x+y)
	            {
	              if ( int(n+k) >= std::abs(m+l) )  // gamma_{n-k}^{(m-l)*} not zero
	              {
	                if ( (m+l) >= 0)  // indexing lower triangular part of theta matrix theta_y_matrix -> thetea_{n+k}^{m+l} (y)
	                {
	                  if ( l >= 0 )  // indexing lower triangular part of gamma matrix rCoeff_old -> gamma_k^{l*} (x)
	                  {
	           	        // ans[n][n-m] is the location of theta_n^{m} (x+y)
	        	        // rCoeff_old[k][k-l] is the location of gamma_{k}^{l*} (x)
	        	        // theta_y_matrix[n+k][(n+k)-(m+l)] is the location of theta_{n+k}^{m+l} (y)
	                    ans[n+m][n] += sCoeff_old[k][k-l] * theta_y_matrix[n+k][(n+k)-(m+l)];
	                  }
	                  else // l < 0 - indexing upper triangular part of gamma matrix rCoefff_old -> gamma_k^{l*} (x)
	                  {
	                  	// do nothing: l < 0 and (m+l) >= 0  =>  m >= -l > 0
	                  	// m is positive but we are doing the case m < 0 (outer if condition)
	                  }
	                }
	                else // (m+l) < 0 and indexing upper triangular part of theta matrix theta_y_matrix -> theta_{n+k}^{(m+l)} (y)
	                {
	                  if ( l >= 0 )  // indexing lower triangular part of gamma matrix rCoeff_old -> gamma_k^{l*} (x)
	                  {
	            	    // ans[n][n-m] is the location of theta_n^{(m)*} (x+y)
	          	        // rCoeff_old[k][k-l] is the location of gamma_{k}^{l*} (y)
	          	        // theta_y_matrix[(n+k)+(m+l)][n+k] is the location of theta_{n+k}^{m+l} (y)
	                    ans[n+m][n] += sCoeff_old[k][k-l] * theta_y_matrix[(n+k)+(m+l)][n+k];
	                  }
	                  else // l < 0 // indexing upper triangular part of gamma matrix rCoeff_old -> gamma_k^{l*} (x)
	                  {
		            	// ans[n][n-m] is the location of theta_n^{(m)*} (x+y)
		          	    // rCoeff_old[k][k-l] is the location of gamma_{k}^{l*} (y)
		          	    // theta_y_matrix[(n+k)+(m+l)][n+k] is the location of theta_{n+k}^{m+l} (y)
	                    ans[n+m][n] += sCoeff_old[k+l][k] * theta_y_matrix[(n+k)+(m+l)][n+k];
	                  }
	                }
	              }
	              else // (n-k) < std::abs(m-l)
	              {
	                //  no need to do anything since gamma_{n-k}^{(m-l)*} = zero
	              }
	            }
	          }
	      }
	    }

	  return ans;

}

/*
std::vector<std::vector<std::vector<std::complex<double> > > > Potential::getSRgrad(std::vector<double> from, std::vector<double> to,
		                                                                            const std::vector<std::vector<std::complex<double> > > &sCoeff_old)
{
	  // from - refers to the center where the old coefficients were evaluated
	  //        z_A^{(p)} - the parent cell center z_A^{(p)} of the source point
	  // to   - refers to the center where the new coefficients will be evaluated
	  //        z_B^{(p)} - the parent cell center z_B^{(p)} of the target point
	  // t    - refers to y = from - to = z_A^{(p)} - z_B^{(p)}
	  // Here we have sCoeff_old = \gamma_n^{(m)*} (x) = \gamma_n^{m*} (x_a - z_A^{(p)}}) matrix

	  std::vector<std::vector<std::vector<std::complex<double> > > >
	      ans(order_of_approximation, std::vector<std::vector<std::complex<double> > >(order_of_approximation, std::vector<std::complex<double> >(3)));

//	  // Here x = from and x+y = to and t = to - from = y
//	  std::vector<double>  t(3);
//	  t[0] = from[0] - to[0];  t[1] = from[1] - to[1]; t[2] = from[2] - to[2];


	  // Here we build the \theta_n^{(m)} (x) = \theta_n^{m} (z_B - z_A) matrix
	  // Recall the matrix dimensions
	  // std::vector<std::vector<std::complex<double> > >
	  //                  theta_y_matrix(order_of_approximation+1, std::vector<std::complex<double> >(order_of_approximation+1));
	  // Recall
	  // Here we are interested in \theta_n^m(z_B - x_a) = \theta_n^m(z_B - z_A - (x_a - z_A))
	  //                                                 = \sum_{k=0}^{\infty} \sum_{l=-k}^{k} \gamma_{k}^{l*}(x_a - z_A)
	  //                                                                                              \theta_{n+k}^{m+l} (z_B - z_A)
	  // Recall getRCoeff(xi, xstar) performs its operations using the difference (xstar - xi).
	  // Therefore, we let from = z_B = xstar and to = z_A = xi to
	  // We call getRCoeff(zA, zB) = getRCoeff(from, to)
	  // --                                                                                                                  --
	  // | \theta_0^0     \theta_1^{-1}     \theta_2^{-2}     \theta_3^{-3}     \theta_4^{-4}     \hdots      \theta_n^{-n}    |
	  // |                                                                                                                     |
	  // | \theta_1^1     \theta_1^0        \theta_2^{-1}     \theta_3^{-2}     \theta_4^{-3}                 \theta_n^{-n+1}  |
	  // |                                                                                                                     |
	  // | \theta_2^2     \theta_2^1        \theta_2^0        \theta_3^{-1}     \theta_4^{-2}                 \theta_n^{-n+2}  |
	  // |                                                                                                                     |
	  // | \theta_3^3     \theta_3^2        \theta_3^1        \theta_3^0        \theta_4^{-1}                 \vdots           |
	  // |                                                                                                                     |
	  // | \theta_4^4     \theta_4^3        \theta_4^2        \theta_4^1        \theta_4^0                                     |
	  // |                                                                                                                     |
	  // | \vdots                                                                            \ddots                            |
	  // |                                                                                                                     |
	  // | \theta_n^n                                                                                         \theta_n^0       |
	  // |                                                                                                                     |
	  // --                                                                                                                  --

	  std::vector<std::vector<std::complex<double> > > theta_y_matrix = this->getRCoeff(from, to);

      // Here we will have three matrices based on theta_y_matrix
	  // theta_y_matrix_component_1 = [\hat{theta_1}_n^m] = [   0.5*(theta_{n+1}^{m+1} - theta_{n+1}^{m-1})]
	  // theta_y_matrix_component_2 = [\hat{theta_2}_n^m] = [-0.5*i*(theta_{n+1}^{m+1} + theta_{n+1}^{m-1})]
	  // theta_y_matrix_component_2 = [\hat{theta_3}_n^m] = [        theta_{n+1}^{m}                       ]

	  std::vector<std::vector<std::complex<double> > > theta_hat_y_matrix_1(theta_y_matrix.size()-1, std::vector<std::complex<double> >(theta_y_matrix.size()-1));
	  std::vector<std::vector<std::complex<double> > > theta_hat_y_matrix_2(theta_y_matrix.size()-1, std::vector<std::complex<double> >(theta_y_matrix.size()-1));
	  std::vector<std::vector<std::complex<double> > > theta_hat_y_matrix_3(theta_y_matrix.size()-1, std::vector<std::complex<double> >(theta_y_matrix.size()-1));

	  std::complex<double> imaginary_unit(0.0,1.0);

      // -----------------------------------------------------------------------------------------------------------
	  // Creating the theta_hat_y_matrix_1 -------------------------------------------------------------------------
	  // The elements are calculated using theta_y_matrix and 0.5*(theta_{n+1}^{m+1} - theta_{n+1}^{m-1}) ----------
      // -----------------------------------------------------------------------------------------------------------

	  // Creating main diagonal entries for theta_hat_y_matrix_1, theta_hat_y_matrix_2, and theta_hat_y_matrix_3
	  // In the calculation, 0.5*(theta_{n+1}^{m+1} - theta_{n+1}^{m-1}) or 0.5*i*(theta_{n+1}^{m+1} + theta_{n+1}^{m-1})
	  // theta_{n+1}^{m+1} is one row below and theta_{n+1}^{m-1} is one column to the right (in the theta_y_matrix)
	  // In the calculation, theta_{n+1}^{m} - theta_{n+1}^{m} is located on the same diagonal as theta_{n}^{m} just below
	  // (to the 1 column right and 1 row down
	  for (unsigned int i=0; i<theta_y_matrix.size()-1; ++i)    // columns
	  {
        theta_hat_y_matrix_1[i][i] = 0.5*(theta_y_matrix[i][i+1]-theta_y_matrix[i+1][i]);
        theta_hat_y_matrix_2[i][i] = -0.5 * imaginary_unit * (theta_y_matrix[i][i+1]+theta_y_matrix[i+1][i]);
        theta_hat_y_matrix_3[i][i] = theta_y_matrix[i+1][i+1];
	  }

	  // Creating lower triangular entries (below the main diagonal where row is greater than column)
	  // In the calculation, 0.5*(theta_{n+1}^{m+1} - theta_{n+1}^{m-1})
	  // theta_{n+1}^{m+1} is down 1 row and theta_{n+1}^{m-1} down 1 row and right 2 columns (in the theta_y_matrix)
	  for (unsigned int j=0; j<theta_hat_y_matrix_1.size(); ++j)    // columns
        for (unsigned int i=j+1; i<theta_hat_y_matrix_1.size(); ++i)  // rows
        {
          theta_hat_y_matrix_1[i][j] = 0.5 * (theta_y_matrix[i+1][j] - theta_y_matrix[i+1][j+2]);
          theta_hat_y_matrix_2[i][j] = -0.5 * imaginary_unit * (theta_y_matrix[i+1][j] - theta_y_matrix[i+1][j+2]);
          theta_hat_y_matrix_3[i][j] = theta_y_matrix[i+1][j+1];
        }

	  // Creating upper triangular entries (above the main diagonal where column is greater than row)
	  // In calculation, 0.5*(theta_{n+1}^{m+1} - theta_{n+1}^{m-1})
	  // theta_{n+1}^{m+1} is one column to the right and theta_{n+1}^{m-1} is to the right 2 columns and down one row (in the theta_y_matrix)
	  for (unsigned int i=0; i<theta_hat_y_matrix_1.size(); ++i)    // rows
        for (unsigned int j=i+1; j<theta_hat_y_matrix_1.size(); ++j)  // columns
        {
          theta_hat_y_matrix_1[i][j] = 0.5 * (theta_y_matrix[i+2][j+1] - theta_y_matrix[i][j+1]);
          theta_hat_y_matrix_2[i][j] = -0.5 * imaginary_unit * (theta_y_matrix[i+2][j+1] - theta_y_matrix[i][j+1]);
          theta_hat_y_matrix_3[i][j] = theta_y_matrix[i+1][j+1];
        }

	  for (unsigned int j = 0; j<=order_of_approximation-1; ++j)
	    for (unsigned int i = 0; i<=order_of_approximation-1; ++i)
	      for (unsigned int k = 0; k < 3; ++k)
  	        ans[i][j][k] = 0.0;

	  // Generating the matrix \gamma_n^{(m)*} (-y) = \theta_n^{m} (x+y) = theta_n^m (x_b - z_B^{(p)}})
	  // using the formula
	  //
	  //   \theta_n^{m} (x+y) = \sum_{k=0}^{order_of_approximation - n}
	  //                                  \sum_{l=-k}^{k} \gamma_{k}^{(l)*} (-y)  theta_{n+k}^{m+l} (x)
	  //
	  for (unsigned int n = 0; n <= order_of_approximation-1; ++n)
	    for (int m = -n; m<=int(n); ++m)
	    {
	      if (int(n) >= std::abs(m)) // assemble an element, otherwise n < m implies element is zero
	      {
	    	// theta_n^m matrix has entries only up to order_of_approximation
	    	// Therefore, n+k can only be up to order_of_approximation => n+k <= order_of_approximation
	    	//                    => k <= order_of_approximation - n
	    	// or else we will be looking for theta_{n+k}^{m+l} that are not there (indexing outside of the theta matrix)
	        for (unsigned int k=0; k<=(order_of_approximation-1 - n); ++k)
	          for (int l=-k; l<=int(k); ++l)
	          {
	            if ( m >= 0 ) // setting lower triangular part of theta matrix ans -> theta_n^{(m)} (x+y)
	            {
	              if ( int(n+k) >= std::abs(m+l) )  // theta_{n+k}^{(m+l)*} not zero
	              {
	                if ( (m+l) >= 0)  // indexing lower triangular part of theta matrix theta_y_matrix -> theta_{n+k}^{(m+l)} (y)
	                {
	                  if ( l >= 0 )  // indexing lower triangular part of gamma matrix sCoeff_old -> gamma_k^{l*} (x)
	                	             // m >=0 and (m+l) >= 0 and l >= 0
	                	             // implies that m >= -l
	                  {
	            	    // ans[n][n-m] is the location of theta_n^{(m)*} (x+y)
	          	        // rCoeff_old[n+k][(n+k)-(m+l)] is the location of theta_{n+k}^{(m+l)} (x)
	          	        // theta_y_matrix[k][k-l] is the location of gamma_{k}^{(l)*} (-y)
	                    ans[n][n-m][0] += sCoeff_old[k][k-l] * theta_hat_y_matrix_1[n+k][(n+k)-(m+l)];
	                    ans[n][n-m][1] += sCoeff_old[k][k-l] * theta_hat_y_matrix_2[n+k][(n+k)-(m+l)];
	                    ans[n][n-m][2] += sCoeff_old[k][k-l] * theta_hat_y_matrix_3[n+k][(n+k)-(m+l)];
	                  }
	                  else // l < 0 // indexing upper triangular part of gamma matrix sCoeff_old -> gamma_k^{l*} (x)
	                	            // (m + l) >= 0 and l < 0 => m >= -l > 0
	                  {
	                    // EX: ans[n][n-m] += rCoeff_old[n+k][(n+k)-(m+l)] * gamma_neg_y_matrix[k+l][k];
	                	//     n = 2, m = 1
	                	//       k = 1, l = -1
	                    //         theta_2^1(x+y) += theta_3^0(x) * gamma_1^{-1}(-y)
	                	//         theta_new[2][1] += theta_y_matrix[3][3] * gamma_[0][1]
	                	//         ans[n][2-1] += theta_y_matrix[n+k][n+k-(m+l)] * sCoeff_old[k+l][k]
	                  	//     n = 2, m = 2
	                  	//       k = 1, l = -1
	                      //         theta_2^2(x+y) += theta_3^1(y) * gamma_1^{-1}(-y)
	                  	//           theta_new[2][1] += theta_y_matrix[3][2] * sCoeff_old[0][1]
	                  	//           ans[n][2-1] += theta_y_matrix[n+k][n+k-(m+l)] * sCoeff_old[k+l][k]
	                    ans[n][n-m][0] += sCoeff_old[k+l][k] * theta_hat_y_matrix_1[n+k][(n+k)-(m+l)];
	                    ans[n][n-m][1] += sCoeff_old[k+l][k] * theta_hat_y_matrix_2[n+k][(n+k)-(m+l)];
	                    ans[n][n-m][2] += sCoeff_old[k+l][k] * theta_hat_y_matrix_3[n+k][(n+k)-(m+l)];
	                  }
	                }
	                else // (m+l) < 0 and indexing upper triangular part of theta matrix theta_y_matrix -> theta_{n+k}^{(m+l)*} (y)
	                {
	                  if ( l >= 0 )  // indexing lower triangular part of gamma matrix sCoeff_old -> gamma_k^{l*} (x)
	                  {
	                	// do nothing: l >= 0 and (m+l) < 0  =>  m < -l < 0
	                	// m has to be negative and we are doing case m >= 0 (outer if condition)
	                  }
	                  else // l < 0 // indexing upper triangular part of gamma matrix sCoeff_old -> gamma_k^{l*} (x)
	                  {
	                    // EX: ans[n][n-m] += theta_y_matrix[n+k][(n+k)+(m+l)] * sCoeff_old[k+l][k];
	                  	//     n = 2, m = 1
	                  	//       k = 1, l = -1
	                    //         theta_2^1(x+y) += theta_3^{0}(y) * gamma_1^{(-1)*}(x)
	                  	//         theta_new[2][1] += theta_y_matrix[3][3] * sCoeff_old[0][1]
	                  	//         ans[n][n-m] += theta_y_matrix[n+k][n+k-(m+l)] * sCoeff_old[k+l][k]
	                    	//     n = 2, m = 0
	                    	//       k = 1, l = -1
	                        //         theta_2^0(x+y) += theta_3^{-1}(y) * gamma_1^{-1}(x)
	                    	//         theta_new[2][2] += theta_y_mamtrix[2][3] * sCoeff_old[0][1]
	                    	//         ans[n][2-1] += theta_y_mamtrix[n+k][n+k-(m+l)] * sCoeff_old[k+l][k]
	                    ans[n][n-m][0] += sCoeff_old[k+l][k] * theta_hat_y_matrix_1[(n+k)+(m+l)][n+k];
	                    ans[n][n-m][1] += sCoeff_old[k+l][k] * theta_hat_y_matrix_2[(n+k)+(m+l)][n+k];
	                    ans[n][n-m][2] += sCoeff_old[k+l][k] * theta_hat_y_matrix_3[(n+k)+(m+l)][n+k];
	                  }
	                }
	              }
	              else // (n-k) < std::abs(m-l)
	              {
	                //  no need to do anything since gamma_{n-k}^{(m-l)*} = zero
	              }
	            }
	            else // m < 0 - setting lower triangular part of gamma matrix ans -> gamma_n^{(m)*} (x+y)
	            {
	              if ( int(n+k) >= std::abs(m+l) )  // gamma_{n-k}^{(m-l)*} not zero
	              {
	                if ( (m+l) >= 0)  // indexing lower triangular part of theta matrix theta_y_matrix -> thetea_{n+k}^{m+l} (y)
	                {
	                  if ( l >= 0 )  // indexing lower triangular part of gamma matrix rCoeff_old -> gamma_k^{l*} (x)
	                  {
	           	        // ans[n][n-m] is the location of theta_n^{m} (x+y)
	        	        // rCoeff_old[k][k-l] is the location of gamma_{k}^{l*} (x)
	        	        // theta_y_matrix[n+k][(n+k)-(m+l)] is the location of theta_{n+k}^{m+l} (y)
	                    ans[n+m][n][0] += sCoeff_old[k][k-l] * theta_hat_y_matrix_1[n+k][(n+k)-(m+l)];
	                    ans[n+m][n][1] += sCoeff_old[k][k-l] * theta_hat_y_matrix_2[n+k][(n+k)-(m+l)];
	                    ans[n+m][n][2] += sCoeff_old[k][k-l] * theta_hat_y_matrix_3[n+k][(n+k)-(m+l)];
	                  }
	                  else // l < 0 - indexing upper triangular part of gamma matrix rCoefff_old -> gamma_k^{l*} (x)
	                  {
	                  	// do nothing: l < 0 and (m+l) >= 0  =>  m >= -l > 0
	                  	// m is positive but we are doing the case m < 0 (outer if condition)
	                  }
	                }
	                else // (m+l) < 0 and indexing upper triangular part of theta matrix theta_y_matrix -> theta_{n+k}^{(m+l)} (y)
	                {
	                  if ( l >= 0 )  // indexing lower triangular part of gamma matrix rCoeff_old -> gamma_k^{l*} (x)
	                  {
	            	    // ans[n][n-m] is the location of theta_n^{(m)*} (x+y)
	          	        // rCoeff_old[k][k-l] is the location of gamma_{k}^{l*} (y)
	          	        // theta_y_matrix[(n+k)+(m+l)][n+k] is the location of theta_{n+k}^{m+l} (y)
	                    ans[n+m][n][0] += sCoeff_old[k][k-l] * theta_hat_y_matrix_1[(n+k)+(m+l)][n+k];
	                    ans[n+m][n][1] += sCoeff_old[k][k-l] * theta_hat_y_matrix_2[(n+k)+(m+l)][n+k];
	                    ans[n+m][n][2] += sCoeff_old[k][k-l] * theta_hat_y_matrix_3[(n+k)+(m+l)][n+k];
	                  }
	                  else // l < 0 // indexing upper triangular part of gamma matrix rCoeff_old -> gamma_k^{l*} (x)
	                  {
		            	// ans[n][n-m] is the location of theta_n^{(m)*} (x+y)
		          	    // rCoeff_old[k][k-l] is the location of gamma_{k}^{l*} (y)
		          	    // theta_y_matrix[(n+k)+(m+l)][n+k] is the location of theta_{n+k}^{m+l} (y)
	                    ans[n+m][n][0] += sCoeff_old[k+l][k] * theta_hat_y_matrix_1[(n+k)+(m+l)][n+k];
	                    ans[n+m][n][1] += sCoeff_old[k+l][k] * theta_hat_y_matrix_2[(n+k)+(m+l)][n+k];
	                    ans[n+m][n][2] += sCoeff_old[k+l][k] * theta_hat_y_matrix_3[(n+k)+(m+l)][n+k];
	                  }
	                }
	              }
	              else // (n-k) < std::abs(m-l)
	              {
	                //  no need to do anything since gamma_{n-k}^{(m-l)*} = zero
	              }
	            }
	          }
	      }
	    }

	  return ans;

}
*/



/**
 * Explanation of getSS
 *
 * The getSS function returns the new coefficients \gamma_n^{m*} (y) from the SS translation
 * applied to old coefficients \gamma_n^{m*} (x) using the translation formula
 *
 *   \gamma_n^{m*} (x+y) = \sum_{k=0}^n \sum_{l=-k}^{k} \gamma_{n-k}^{(m-l)*} (x) + gamma_{k}^{l} (y)
 *
 * The translation to obtain new coefficients is part of the Upward Pass Part II
 *
 * The Upward Pass Part II is the second step of the Fast Multipole Method.
 * In the first step of the upward pass, The Upward Pass Part I, the far-field expansions
 * for each source particle are determined.  The second step is the translation of the
 * far-field expansions from child source cell center x_c^{(c)} to parent source cell center x_c^{(p)}.
 * The translations from child to parent are done at each refinement level up to level l = 2.
 * Here we count the number of operations to complete this part of the Upward Pass.
 *
 * Here x is the location where the old coefficients were evaluated.  Looking at our far-field series
 * whose center is at the source point child cell with center z_A^{(c)}
 *
 *      \frac{1}{|x_a - x_b|} = \frac{1}{|(x_a - z_A^{(c)}) - (x_b - z_A^{(c)})|}
 *                            = \sum_{n=0}^{\infty} \sum_{m=-n}^{n} \gamma_n^{m*} (x_a - z_A^{(c)}) \theta_n^{m*} (x_b - z_A^{(c)})
 *
 * we have gamma_n^{(m)*} (x) = \gamma_n^{(m)*} (x_a - z_A^{(c)}  and  x = x_a - z_A^{(c)}
 *
 * We want to change the center to x + y = x_a - z_A^{(p)} which means y = x_a - z_A^{(p)} - (x_a - z_A^{(c)}
 * and then
 *
 *      y = z_A^{(c)} - z_A^{(p)} = from - to
 *
 * below we use t    = y
 *              from = z_A^{(c)} (from child center)
 *              to   = z_A^{(p)} (to parent center)
 *
 * From the translation formula
 *
 *   \gamma_n^{m*} (x+y) = \sum_{k=0}^n \sum_{l=-k}^{k} \gamma_{n-k}^{(m-l)*} (x) . gamma_{k}^{l*} (y)
 *
 * For each 0 <= m <= n <= 3 we will be doing a nested loop
 *
 * for k=0:n
 *   for l=-k:k
 *     if (n-k >= m-l)
 *       gamma_n^{(m)*} (x+y) += \sum_{k=0}^{n} \sum_{l=-k}^{k} gamma_{n-k}^{(m-l)*} (x) gamma_{k}^{l*} (y)
 *
 *
 * With respect to the matrices gamma_n^{(m)*} (x+y), gamma_{n-k}^{(m-l)*} (x), and \gamma_{k}^{(l)*}
 * we have the value of the subscript matching the row or column value exactly
 * However, the superscript value counts away from the main diagonal.  Therefore, if the subscript is n
 * Then we are working with either column or row n, depending on if the superscript is positive or negative.
 * If the superscript is m, then we count away from the main diagonal by subtracting n - |m|.  For example,
 * if m = 0 then we are on the main diagonal and n-m makes sense.
 *
 * // --                                                                                                                  --
 * // | \gamma_0^0*    \gamma_1^{-1*}    \gamma_2^{-2*}    \gamma_3^{-3*}    \gamma_4^{-4*}    \hdots      \gamma_n^{-n*}   |
 * // |                                                                                                                     |
 * // | \gamma_1^1*    \gamma_1^0*       \gamma_2^{-1*}    \gamma_3^{-2*}    \gamma_4^{-3*}                \gamma_n^{-n+1*} |
 * // |                                                                                                                     |
 * // | \gamma_2^2*    \gamma_2^1*       \gamma_2^0*       \gamma_3^{-1*}    \gamma_4^{-2*}                \gamma_n^{-n+2*} |
 * // |                                                                                                                     |
 * // | \gamma_3^3*    \gamma_3^2*       \gamma_3^1*       \gamma_3^0*       \gamma_4^{-1*}                \vdots           |
 * // |                                                                                                                     |
 * // | \gamma_4^4*    \gamma_4^3*       \gamma_4^2*       \gamma_4^1*       \gamma_4^0*                                    |
 * // |                                                                                                                     |
 * // | \vdots                                                                                 \ddots                       |
 * // |                                                                                                                     |
 * // | \gamma_n^n*                                                                                        \gamma_n^0*      |
 * // |                                                                                                                     |
 * // --                                                                                                                  --
 *
 * suppose n = 3 and m = 2
 *
 * let k = 0
 *       l = -0:0
 *       l = 0
 *         n - k = 3 - 0 = 3         |
 *         m - l = 2 - (0) = 2       |  n-k >= m-l True
 *           gamma_3^{(2)*} (x+y) += gamma_{3}^{(2)*} (x) gamma_{0}^{(0)*} (y)
 *              ans[3][1]   += sCoeff_old[3][1] * gamma_y_matrix[0][0]
 *              ans[n][n-m] += sCoeff_old[n-k][(n-k)-(m-l)] * gamma_y_matrix[k][k-l]
 * let k = 1
 *       l = -1:1
 *       l = -1
 *         n - k = 3 - 1 = 2         |
 *         m - l = 2 - (-1) = 3      |  n-k >= m-l False
 *         False
 *       l = 0
 *         n - k = 3 - 1 = 2         |
 *         m - l = 2 - (0) = 2       |  n-k >= m-l True
 *         True
 *           gamma_3^{(2)*} (x+y) += gamma_{2}^{(2)*} (x) \gamma_{1}^{(0)*} (y)
 *              ans[3][1]         += sCoeff_old[2][0] * gamma_y_matrix[1][1]
 *              ans[n][n-m]       += sCoeff_old[n-k][(n-k)-(m-l)] * gamma_y_matrix[k][k-l]
 *       l = 1
 *         n - k = 3 - 1 = 2         |
 *         m - l = 2 - (1) = 1       |  n-k >= m-l True
 *         True
 *           gamma_3^{(2)*} (x+y) += gamma_{2}^{(1)*} (x) \gamma_{1}^{(1)*} (y)
 *                            = sCoeff_old[2][1] * gamma_y_matrix[1][0]
 *              ans[3][1]         += sCoeff_old[2][1] * gamma_y_matrix[1][0]
 *              ans[n][n-m]       += sCoeff_old[n-k][(n-k)-(m-l)] * gamma_y_matrix[k][k-l]
 *
 * let k = 2
 *       l = -2:2
 *       l = -2
 *         n - k = 3 - 2 = 1         |
 *         m - l = 2 - (-2) = 4      |  n-k >= m-l False
 *         False
 *       l = -1
 *         n - k = 3 - 2 = 1         |
 *         m - l = 2 - (-1) = 3      |  n-k >= m-l False
 *         False
 *       l = 0
 *         n - k = 3 - 2 = 1         |
 *         m - l = 2 - (0) = 2       |  n-k >= m-l False
 *         False
 *       l = 1
 *         n - k = 3 - 2 = 1         |
 *         m - l = 2 - (1) = 1       |  n-k >= m-l True
 *         True
 *           gamma_3^{(2)*} (x+y) += gamma_{1}^{(1)*} (x) \gamma_{2}^{(1)*} (y)
 *                                 = sCoeff_old[1][0] * gamma_y_matrix[2][1]
 *       l = 2
 *         n - k = 3 - 2 = 1         |
 *         m - l = 2 - (2) = 0       |  n-k >= m-l True
 *         True
 *           gamma_3^{(2)*} (x+y) += gamma_{1}^{(0)*} (x) \gamma_{2}^{(2)*} (y)
 *                                 = sCoeff_old[1][1] * gamma_y_matrix[2][0]
 * let k = 3
 *       l = -3
 *         n - k = 3 - 3 = 0         |
 *         m - l = 2 - (-3) = 5      |  n-k >= m-l False
 *         False
 *       l = -2
 *         n - k = 3 - 3 = 0         |
 *         m - l = 2 - (-2) = 4      |  n-k >= m-l False
 *         False
 *       l = -1
 *         n - k = 3 - 3 = 0         |
 *         m - l = 2 - (-1) = 3      |  n-k >= m-l False
 *         False
 *       l = 0
 *         n - k = 3 - 3 = 0         |
 *         m - l = 2 - (0) = 2       |  n-k >= m-l False
 *         False
 *       l = 1
 *         n - k = 3 - 3 = 0         |
 *         m - l = 2 - (1) = 1       |  n-k >= m-l False
 *         False
 *       l = 2
 *         n - k = 3 - 2 = 0         |
 *         m - l = 2 - (2) = 0       |  n-k >= m-l True
 *         True
 *           gamma_3^{(2)*} (x+y) += gamma_{0}^{(0)*} (x) \gamma_{3}^{(2)*} (y)
 *                                 = sCoeff_old[0][0] * gamma_y_matrix[3][1]
 *
 *
 * We write the results in matrix notation
 *
 *   \gamma_n^{m*} (x+y) = \sum_{k=0}^n \sum_{l=-k}^{k} \gamma_{n-k}^{(m-l)*} (x)  gamma_{k}^{l*} (y)
 *
 *                                gamma_0^{( 0)*} (x)        gamma_1^{(-1)*} (x)           gamma_1^{( 0)*} (x)            gamma_1^{( 1)*} (y)           gamma_2^{(-2)*} (y)          gamma_2^{(-1)*} (y)    gamma_2^{( 0)*} (y)
 *                                     n-k   m-l                  n-k   m-l                     n-k   m-l                      n-k   m-l                        k    l                       k    l
 *
 *                                      k    l                     k    l                        k    l                         k    l                     n-k   m-l                    n-k   m-l
 *  --                --      --                                                                                                                                                                                                                                                                                                                   --      --                   --
 * | gamma_0^{( 0)*} (x+y) |    | gamma_0^{( 0)*} (y)      & 0                           & 0                            & 0                           & 0                           & 0                   & 0                   & 0           & 0            & 0            & 0            & 0            & 0            & 0           & 0           & 0           |     | gamma_0^{( 0)*) (y)  |
 * | gamma_1^{(-1)*} (x+y) |    | gamma_1^{(-1)*} (y)      & gamma_0^{( 0)*} (y)         & 0                            & 0                           & 0                           & 0                   & 0                   & 0           & 0            & 0            & 0            & 0            & 0            & 0           & 0           & 0           |     | gamma_1^{(-1)*} (y)  |
 * | gamma_1^{( 0)*} (x+y) |    | gamma_1^{( 0)*} (y)      & 0                           & gamma_0^{( 0)*} (y)          & 0                           & 0                           & 0                   & 0                   & 0           & 0            & 0            & 0            & 0            & 0            & 0           & 0           & 0           |     | gamma_1^{( 0)*} (y)  |
 * | gamma_1^{( 1)*} (x+y) |    | theta_1^{( 1)*} (y)      & 0                           & 0                            & gamma_0^{( 0)*} (y)         & 0                           & 0                   & 0                   & 0           & 0            & 0            & 0            & 0            & 0            & 0           & 0           & 0           |     | gamma_1^{( 1)*} (y)  |
 * | gamma_2^{(-2)*} (x+y) |    | gamma_2^{(-2)*} (y)      & gamma_1^{(-1)*} (y)         & 0                            & 0                           & gamma_0^{( 0)*} (y)         & 0                   & 0                   & 0           & 0            & 0            & 0            & 0            & 0            & 0           & 0           & 0           |     | gamma_2^{(-2)*} (y)  |
 * | gamma_2^{(-1)*} (x+y) |    | gamma_2^{(-1)*} (y)      & gamma_1^{( 0)*} (y)         & gamma_1^{(-1)*} (y)          & 0                           & 0                           & gamma_0^{( 0)*} (y) & 0                   & 0           & 0            & 0            & 0            & 0            & 0            & 0           & 0           & 0           |     | gamma_2^{(-1)*} (y)  |
 * | gamma_2^{( 0)*} (x+y) |    | gamma_2^{( 0)*} (y)      & gamma_1^{( 1)*} (y)         & gamma_1^{( 0)*} (y)          & gamma_1^{(-1)*} (y)         & 0                           & 0                   & gamma_0^{( 0)*} (y) & 0           & 0            & 0            & 0            & 0            & 0            & 0           & 0           & 0           |     | gamma_2^{( 0)*} (y)  |
 * | theta_2^{ 1} (x+y) |    | theta_2^{1}  & theta_3^{-0} & theta_3^{1}  & theta_3^{2}  & 0            & 0            & 0            & 0            & 0           & 0            & 0            & 0            & 0            & 0            & 0           & 0           |     | gamma_2^1    |
 * | theta_2^{ 2} (x+y) |    | theta_2^{2}  & theta_3^{1}  & theta_3^{2}  & theta_3^{3}  & 0            & 0            & 0            & 0            & 0           & 0            & 0            & 0            & 0            & 0            & 0           & 0           |     | gamma_2^2    |
 * | theta_3^{-3} (x+y) |  = | theta_3^{-3} & 0            & 0            & 0            & 0            & 0            & 0            & 0            & 0           & 0            & 0            & 0            & 0            & 0            & 0           & 0           |     | gamma_3^{-3} |
 * | theta_3^{-2} (x+y) |    | theta_3^{-2} & 0            & 0            & 0            & 0            & 0            & 0            & 0            & 0           & 0            & 0            & 0            & 0            & 0            & 0           & 0           |     | gamma_3^{-2} |
 * | theta_3^{-1} (x+y) |    | theta_3^{-1} & 0            & 0            & 0            & 0            & 0            & 0            & 0            & 0           & 0            & 0            & 0            & 0            & 0            & 0           & 0           |     | gamma_3^{-1} |
 * | theta_3^{0}  (x+y) |    | theta_3^{0}  & 0            & 0            & 0            & 0            & 0            & 0            & 0            & 0           & 0            & 0            & 0            & 0            & 0            & 0           & 0           |     | gamma_3^0    |
 * | theta_3^{1}  (x+y) |    | theta_3^{1}  & 0            & 0            & 0            & 0            & 0            & 0            & 0            & 0           & 0            & 0            & 0            & 0            & 0            & 0           & 0           |     | gamma_3^1    |
 * | theta_3^{2}  (x+y) |    | theta_3^{2}  & 0            & 0            & 0            & 0            & 0            & 0            & 0            & 0           & 0            & 0            & 0            & 0            & 0            & 0           & 0           |     | gamma_3^2    |
 * | theta_3^{3}  (x+y) |    | theta_3^{3}  & 0            & 0            & 0            & 0            & 0            & 0            & 0            & 0           & 0            & 0            & 0            & 0            & 0            & 0           & 0           |     | gamma_3^3    |
 *  --                --      --                                                                                                                                                                                                                                        --      --           _--
 *
 */


std::vector<std::vector<std::complex<double> > > Potential::getSS(std::vector<double> from, std::vector<double> to,
		                                                          const std::vector<std::vector<std::complex<double> > > &sCoeff_old)
{
  // from - refers to the center where the old coefficients were evaluated
  //        z_A^{(c)} - the child cell center z_A^{(c)} of the source point
  // to   - refers to the center where the new coefficients will be evaluated
  //        z_A^{(p)} - the parent cell center z_A^{(p)} of the source point
  // t    - refers to y = z_A^{(c)} - z_A^{(p)}

  // terms used to approximate the Taylor series
////  p = 16;
  // highest order of derivative used to approximate Taylor series, to which p corresponds
////  order_of_approximation = 3;
  // matrix of coefficients that will be returned gamma_n^{(m)*} (x+y)
  // The interchangibility of x and y in the formula
  //
  //   \gamma_n^{m*} (x+y) = \sum_{k=0}^n \sum_{l=-k}^{k} \gamma_{n-k}^{(m-l)*} (x)  gamma_{k}^{l} (y)
  //   \gamma_n^{m*} (x+y) = \sum_{k=0}^n \sum_{l=-k}^{k} \gamma_{n-k}^{(m-l)*} (y)  gamma_{k}^{l} (x)
  //
  // means that the we will get mirror terms and can account for both terms at once in the loop
  // Note: x = x - x_c^{(c)}
  //       y = x_c^{(c)} - x_c^{(p)}
  //
//  for (unsigned int i=0; i<4; ++i)
//    for (unsigned int j=0; j<4; ++j)
//      std::cout << "sCoeff_old[" << i << "][" << j << "] = "
//                << sCoeff_old[i][j] << std::endl;

  // note - Taylor series includes zeroth derivative - therefore adding 1 to highest order derivative
  std::vector<std::vector<std::complex<double> > > ans(order_of_approximation+1, std::vector<std::complex<double> >(order_of_approximation+1));

  // Translation formula for \gamma_n^{m*} is
  // Here x = from and x+y = to and t = to - from = y
  std::vector<double>  t(3);
  t[0] = from[0] - to[0];  t[1] = from[1] - to[1]; t[2] = from[2] - to[2];


  // Our first step is to build the \gamma_n^{(m)*} (y) = \gamma_n^{m*} (z_c^{(c) - z_c^{(p)}}) matrix
  // We do this by calling getSCoeff with from = xi = z_c^{(c)} and to = xstar = z_c^{(p)}
  // Looking at the getSCoeff code we see that a matrix for \gamma_n^{(m)*} (xi - xstart) = \gamma_n^{(m)*} (y)
  // will be generated
  // recall
  //std::vector<std::vector<std::complex<double> > >
  //                  gamma_y_matrix(order_of_approximation+1, std::vector<std::complex<double> >(order_of_approximation+1));
  std::vector<std::vector<std::complex<double> > > gamma_y_matrix = this->getSCoeff(from, to);

//  if (from[0] == 0.9375)
//    for (unsigned int i=0; i<(order_of_approximation+1); ++i)
//      for (unsigned int j=0; j<(order_of_approximation+1); ++j)
//        std::cout << "gamma_y_matrix[" << i << "][" << j << "] = "
//                  << gamma_y_matrix[i][j] << std::endl;


  for (unsigned int i = 0; i<=order_of_approximation; ++i)
    for (unsigned int j = 0; j<=order_of_approximation; ++j)
      ans[i][j] = 0.0;
  // Example loop
  // n = 1; m = -1
  // n >= std::abs(m) : True 1 >= |-1|
  //   k = 0
  //     l = 0 => m - l < 0
  //       ans[n+m][n] = ans[0][1]
  //                 =
  //   k = 1
  //   \gamma_n^{m*} (x+y) = \sum_{k=0}^n \sum_{l=-k}^{k} \gamma_{n-k}^{(m-l)*} (x)  gamma_{k}^{l} (y)
  // next we do the nested for loop through
  for (int n = 0; n<=int(order_of_approximation); ++n)
    for (int m = -n; m<=n; ++m)
    {
      if (n >= std::abs(m))
        for (int k=0; k<=n; ++k)
          for (int l=-k; l<=k; ++l)
          {
            if ( m >= 0 ) // setting lower triangular part of gamma matrix ans -> gamma_n^{(m)*} (x+y)
            {
              if ( (n-k) >= std::abs(m-l) )  // gamma_{n-k}^{(m-l)*} not zero
              {
                if ( (m-l) >= 0)  // indexing lower triangular part of gamma matrix sCoeff_old -> gamma_{n-k}^{(m-l)*} (x)
                {
                  if ( l >= 0 )  // indexing lower triangular part of gamma matrix gamma_y_matrix -> gamma_k^{l*} (y)
                  {
            	    // ans[n][n-m] is the location of gamma_n^{(m)*} (x+y)
          	        // sCoeff_old[n-k][(n-k)-(m-l)] is the location of gamma_{n-k}^{(m-l)*} (x)
          	        // gamma_y_matrix[k][k-l] is the location of gamma_{k}^{(l)*} (y)
                    ans[n][n-m] += sCoeff_old[n-k][(n-k)-(m-l)] * gamma_y_matrix[k][k-l];
                  }
                  else // l < 0 // indexing upper triangular part of gamma matrix gamma_y_matrix -> gamma_k^{l*} (y)
                  {
                    ans[n][n-m] += sCoeff_old[n-k][(n-k)-(m-l)] * gamma_y_matrix[k+l][k];
                  }
                }
                else // (m-l) < 0 and indexing upper triangular part of gamma matrix sCoeff_old -> gamma_{n-k}^{(m-l)*} (x)
                {
                  if ( l >= 0 )  // indexing lower triangular part of gamma matrix gamma_y_matrix -> gamma_k^{l*} (y)
                  {
            	      // ans[n][n-m] is the location of gamma_n^{(m)*} (x+y)
          	      // sCoeff_old[n-k][(n-k)-(m-l)] is the location of gamma_{n-k}^{(m-l)*} (x)
          	      // gamma_y_matrix[k][k-l] is the location of gamma_{k}^{(l)*} (y)
                    ans[n][n-m] += sCoeff_old[(n-k)+(m-l)][n-k] * gamma_y_matrix[k][k-l];
                  }
                  else // l < 0 // indexing upper triangular part of gamma matrix gamma_y_matrix -> gamma_k^{l*} (y)
                  {
                    ans[n][n-m] += sCoeff_old[(n-k)+(m-l)][n-k] * gamma_y_matrix[k+l][k];
                  }
                }
              }
              else // (n-k) < std::abs(m-l)
              {
                //  no need to do anything since gamma_{n-k}^{(m-l)*} = zero
              }
            }
            else // m < 0 - setting lower triangular part of gamma matrix ans -> gamma_n^{(m)*} (x+y)
            {
              if ( (n-k) >= std::abs(m-l) )  // gamma_{n-k}^{(m-l)*} not zero
              {
                if ( (m-l) >= 0)  // indexing lower triangular part of gamma matrix sCoeff_old -> gamma_{n-k}^{(m-l)*} (x)
                {
                  if ( l >= 0 )  // indexing lower triangular part of gamma matrix gamma_y_matrix -> gamma_k^{l*} (y)
                  {
           	        // ans[n][n-m] is the location of gamma_n^{(m)*} (x+y)
        	        // sCoeff_old[n-k][(n-k)-(m-l)] is the location of gamma_{n-k}^{(m-l)*} (x)
        	        // gamma_y_matrix[k][k-l] is the location of gamma_{k}^{(l)*} (y)
                    ans[n+m][n] += sCoeff_old[n-k][(n-k)-(m-l)] * gamma_y_matrix[k][k-l];
                  }
                  else // l < 0 // indexing upper triangular part of gamma matrix gamma_y_matrix -> gamma_k^{l*} (y)
                  {
                    ans[n+m][n] += sCoeff_old[n-k][(n-k)-(m-l)] * gamma_y_matrix[k+l][k];
                  }
                }
                else // (m-l) < 0 and indexing upper triangular part of gamma matrix sCoeff_old -> gamma_{n-k}^{(m-l)*} (x)
                {
                  if ( l >= 0 )  // indexing lower triangular part of gamma matrix gamma_y_matrix -> gamma_k^{l*} (y)
                  {
            	    // ans[n][n-m] is the location of gamma_n^{(m)*} (x+y)
          	        // sCoeff_old[n-k][(n-k)-(m-l)] is the location of gamma_{n-k}^{(m-l)*} (x)
          	        // gamma_y_matrix[k][k-l] is the location of gamma_{k}^{(l)*} (y)
                    ans[n+m][n] += sCoeff_old[(n-k)+(m-l)][n-k] * gamma_y_matrix[k][k-l];
//                    std::cout << "n = " << n << std::endl;
//                    std::cout << "m = " << m << std::endl;
//                    std::cout << "l = " << l << std::endl;
//                    std::cout << "ans[" << n+m << "][" << n << "] = " << ans[n+m][n] << std::endl;
                  }
                  else // l < 0 // indexing upper triangular part of gamma matrix gamma_y_matrix -> gamma_k^{l*} (y)
                  {
                    ans[n+m][n] += sCoeff_old[(n-k)+(m-l)][n-k] * gamma_y_matrix[k+l][k];
                  }
                }
              }
              else // (n-k) < std::abs(m-l)
              {
                //  no need to do anything since gamma_{n-k}^{(m-l)*} = zero
              }
            }
          }
//      std::cout << "ans[" << n << "][" << n-m << "] = " << ans[n][n-m] << std::endl;
    }

  return ans;
}


/* Explanation of member function Potential::getRR
 *
 * Recall the near-field series for target y and source x and target parent cell center y_c^{(p)}
 * and target child cell center y_c^{(c)}
 *
 * \frac{1}{|y-x|} = \frac{1}{|(y-y_c^{(p)}) -(x-y_c^{(p)})|}
 *                 = \sum_{n=0}^{\infty} \sum_{m=-n}^{n} \gamma_n^{m*} (y-y_c^{(p)}) theta_m^n (x-y_c^{(p)})
 *
 * \frac{1}{|y-x|} = \frac{1}{|(y-y_c^{(c)}) -(x-y_c^{(c)})|}
 *                 = \sum_{n=0}^{\infty} \sum_{m=-n}^{n} \gamma_n^{m*} (y-y_c^{(c)}) theta_m^n (x-y_c^{(c)})
 *
 * Potential::getRR member function performs an RR transformation to determine new coefficients
 * theta_m^n (x-y_c^{(c)}) from old coefficients theta_m^n (x-y_c^{(p)})
 *
 * The function returns a coefficientt matrix - the set of near-field coefficients theta_n^{(m)} (x+y)
 * The coefficients correspond to a change of center from a parent target cell center to a child target
 * cell center.  The coefficients are determined by the formula
 *
 *   \theta_n^{m} (x+y) = \sum_{k=0}^{\infty} \sum_{l=-k}^{k} \gamma_{k}^{(l)*} (-y)  theta_{n+k}^{m+l} (x)
 *
 * We set x and y to satisfy the formulation below
 *
 *   \theta_n^{m} (x+y) = \sum_{k=0}^{order_of_approximation - n}
 *                                   \sum_{l=-k}^{k} (-1)^n \gamma_{k}^{(l)*} (-y)  theta_{n+k}^{m+l} (x)
 *
 *   where y = z_B^{(p)} - z_A^{(p)} = to - from
 *         x = x_a - z_A^{(p)}
 *   and
 *           x + y     = x_a - z_B^{(p)}
 *           x_a       = source point
 *           x_b       = target point
 *           z_A^{(p)} = parent source cell
 *           z_B^{(p)} = parent target cell
 *
 *   seen in
 *
 *   \mu_a \frac{1}{|x_b - x_a|} = \mu_a \frac{1}{|x_b - z_B + z_B - x_a|}
 *                               = \mu_a \frac{1}{|(z_B - x_a) - (z_B - x_b)|
 *                               = \mu_a \sum_{n=0}^{\infty} \sum_{m=-n}^{n} \gamma_n^{m*} (z_B - x_b) \theta_n^m(z_B - x_a)
 *                               = \mu_a \sum_{n=0}^{\infty} \sum_{m=-n}^{n} \gamma_n^{m*} (z_B - x_b)
 *                                                                                       \theta_n^m(z_B - z_A - (x_a - z_A))
 *                               = \mu_a \sum_{n=0}^{\infty} \sum_{m=-n}^{n} \gamma_n^{m*} (z_B - x_b)
 *                                                                  \sum_{k=0}^{\infty} \sum_{l=-k}^{k} \gamma_{k}^{l*}(x_a - z_A)
 *                                                                                                         theta_{n+k}^{m+l} (z_B - z_A)
 * That means the parent (old) coefficients for far-field series that we will be using
 * have the form
 *
 *    \theta_n^m(z_B^{(p)} - x_a)
 *
 * We return this form to be consistent.  Therefore our transformation formula
 *
 *   \theta_n^{m} (x+y) = \sum_{k=0}^{order_of_approximation - n}
 *                                   \sum_{l=-k}^{k} (-1)^n \gamma_{k}^{(l)*} (x)  theta_{n+k}^{m+l} (y)
 *
 * has the form
 *
 *   \theta_n^{m} (z_B^{(c)} - x_a) = \theta_n^{m} (z_B^{(c)} - z_B^{(p)} + z_B^{(p)} - x_a)
 *   = \sum_{k=0}^{order_of_approximation - n}
 *                                          \sum_{l=-k}^{k} (-1)^n \gamma_{k}^{(l)*} (-(z_B^{(c)} - z_B^{(p)}))  theta_{n+k}^{m+l} (z_B^{(p)} - x_a)
 *   = \sum_{k=0}^{order_of_approximation - n}
 *                                          \sum_{l=-k}^{k} (-1)^n \gamma_{k}^{(l)*} (z_B^{(p)} - z_B^{(c)})  theta_{n+k}^{m+l} (z_B^{(p)} - x_a)
 *
 *   where z_B^{(c)} - z_B^{(p)} = to - from
 *         z_B^{(p)} - x_a       = from - source (old coefficients center)
 *
 * We note the difference between this formula and the formula for the SS transformation
 *
 *   \gamma_n^{m*} (x+y) = \sum_{k=0}^n \sum_{l=-k}^{k} \gamma_{n-k}^{(m-l)*} (x)  gamma_{k}^{l} (y)
 *
 * Specifically (i) the subtraction in the subscript and superscript of \gamma_{n-k}^{(m-l)*} (x) indicating
 * no coefficients from the previous center beyond n and (ii) the fact that the best SS transformation we can
 * do for a given order_of_approximation is the same order of approximation which is due k only incrementable
 * up to n (else n-k <= 0 making \gamma_{n-k}^{(m-l)*} (x) = 0).  Therefore, we can use the formula for
 * \gamma_n^{m*} (x+y) with (x+y) being the new center coefficients, x the old center coefficients and
 * y the shift.
 *
 * In contrast, the formula for \theta_n^{m}
 *
 *   \theta_n^{m} (x+y) = \sum_{k=0}^{\infty} \sum_{l=-k}^{k} \gamma_{k}^{(l)*} (-y)  theta_{n+k}^{m+l} (x)
 *                      ~ \sum_{k=0}^{p-n} \sum_{l=-k}^{k} \gamma_{k}^{(l)*} (-y)  theta_{n+k}^{m+l} (x)
 *
 * where p is the order_of_approximation
 *
 * The equality indicates no limit on how well we can do the RR transformation for a given order_of_approximation
 * of the series we are transforming which can be seen by the upper limit \infty of the sum on k.  Due to the
 * summation in the superscript and subscript of theta_{n+k}^{m+l} we see that if we try to use the formula directly with
 * x+y the new center, x the old center, and y the shift, then we would need more theta evaluations at the old
 * center and have less theta evaluations at the new center.  In the downward pass we would need larger matrices
 * at the higher levels in order for the RR shifts down the levels to return the desired order_of_approximation
 * at the lowest level.  Rather than do this we keep all series at the same order and truncate at the order_of_approximation
 * to p-n as shown above
 *
 *
 * Recall theta_m^n (r) = Delta_m^n \frac{1}{|r|} = \frac{1}{|r|} if m = n = 0
 *        gamma_m^n (r) = |r| \frac{1}{|r|} = 1                   if m = n = 0
 *
 * Let the order_of_approximation = 3
 *
 * Examples
 *
 * theta_0^0 (x+y)    = \sum_{k=0}^{3-0} \sum_{l=-k}^{k} \gamma_k^{l*} (-y) theta_{n+k}^{m+l} (x)
 *                    = gamma_0^{( 0)*} (-y)  theta_{0}^{0} (x) +    (1)
 *                      gamma_1^{(-1)*} (-y) theta_{1}^{-1} (x) +    (2)
 *                      gamma_1^{( 0)*} (-y)  theta_{1}^{0} (x) +    (3)
 *                      gamma_1^{( 1)*} (-y)  theta_{1}^{1} (x) +    (4)
 *                      gamma_2^{(-2)*} (-y) theta_{2}^{-2} (x) +    (5)
 *                      gamma_2^{(-1)*} (-y) theta_{2}^{-1} (x) +    (6)
 *                      gamma_2^{( 0)*} (-y)  theta_{2}^{0} (x) +    (7)
 *                      gamma_2^{( 1)*} (-y)  theta_{2}^{1} (x) +    (8)
 *                      gamma_2^{( 2)*} (-y)  theta_{2}^{2} (x) +    (9)
 *                      gamma_3^{(-3)*} (-y) theta_{3}^{-3} (x) +   (10)
 *                      gamma_3^{(-2)*} (-y) theta_{3}^{-2} (x) +   (11)
 *                      gamma_3^{(-1)*} (-y) theta_{3}^{-1} (x) +   (12)
 *                      gamma_3^{( 0)*} (-y)  theta_{3}^{0} (x) +   (13)
 *                      gamma_3^{( 1)*} (-y)  theta_{3}^{1} (x) +   (14)
 *                      gamma_3^{( 2)*} (-y)  theta_{3}^{2} (x) +   (15)
 *                      gamma_3^{( 3)*} (-y)  theta_{3}^{3} (x)     (16)
 *
 * n = 0, m = 0
 *   n >= |m| - continue (theta_n^m != 0)
 *     m >= 0
 *       n + k >= |m + l| - continue (theta_{n+k}^{m+l} != 0)
 *         m + l >= 0
 *           l >= 0
 *             1) theta_0^0 (x)  *  gamma_0^{0*} (-y)
 *             2) theta_1^0 (x)  *  gamma_1^{0*} (-y)
 *             3) theta_1^1 (x)  *  gamma_1^{1*} (-y)
 *             4) theta_2^0 (x)  *  gamma_2^{0*} (-y)
 *             5) theta_2^1 (x)  *  gamma_2^{1*} (-y)
 *             6) theta_2^2 (x)  *  gamma_2^{2*} (-y)
 *             7) theta_3^0 (x)  *  gamma_3^{0*} (-y)
 *             8) theta_3^1 (x)  *  gamma_3^{1*} (-y)
 *             9) theta_3^2 (x)  *  gamma_3^{2*} (-y)
 *            10) theta_3^3 (x)  *  gamma_3^{3*} (-y)
 *           l < 0
 *             none
 *         m + l < 0
 *           l >= 0
 *             none
 *           l < 0
 *             1) theta_1^{-1} (x)  *  gamma_1^{(-1)*} (-y)
 *             2) theta_2^{-2} (x)  *  gamma_2^{(-2)*} (-y)
 *             3) theta_2^{-1} (x)  *  gamma_2^{(-1)*} (-y)
 *             4) theta_3^{-3} (x)  *  gamma_3^{(-3)*} (-y)
 *             5) theta_3^{-2} (x)  *  gamma_3^{(-2)*} (-y)
 *             6) theta_3^{-1} (x)  *  gamma_3^{(-1)*} (-y)
 *
 * theta_1^{-1} (x+y) = \sum_{k=0}^{3-1} \sum_{l=-k}^{k} \gamma_k^{l*} (-y) theta_{n+k}^{m+l} (x)
 *                    = gamma_0^{( 0)*} (-y)  theta_{1}^{-1} (x) +  (1)
 *                      gamma_1^{(-1)*} (-y)  theta_{2}^{-2} (x) +  (2)
 *                      gamma_1^{( 0)*} (-y)  theta_{2}^{-1} (x) +  (3)
 *                      gamma_1^{( 1)*} (-y)  theta_{2}^{ 0} (x) +  (4)
 *                      gamma_2^{(-2)*} (-y)  theta_{3}^{-3} (x) +  (5)
 *                      gamma_2^{(-1)*} (-y)  theta_{3}^{-2} (x) +  (6)
 *                      gamma_2^{( 0)*} (-y)  theta_{3}^{-1} (x) +  (7)
 *                      gamma_2^{( 1)*} (-y)  theta_{3}^{ 0} (x) +  (8)
 *                      gamma_2^{( 2)*} (-y)  theta_{3}^{ 1} (x) +  (9)
 *
 * n = 1, m = -1
 *   n >= |m| - continue (theta_n^m != 0)
 *     m >= 0
 *       n + k >= |m + l| - continue (theta_{n+k}^{m+l} != 0)
 *         m + l >= 0
 *           l >= 0
 *             none
 *           l < 0
 *             none
 *         m + l < 0
 *           l >= 0
 *             do nothing
 *             m + l < 0 and l >= 0 implies m < -l <= 0
 *             implies that m < 0 but we are doing the case m >= 0
 *           l < 0
 *       n + K < |m + l|
 *         do nothing (theta_{n+k}^{m+l} = 0)
 *     m < 0
 *       n + k >= |m + l| - continue (theta_{n+k}^{m+l} != 0)
 *         m + l >= 0
 *           l >= 0
 *             1) k = 1, l =  1 => m + l = (-1) + ( 1) =  0 >= 0 :  theta_2^{ 0} (x)  *  gamma_1^{( 1)*} (-y)
 *                theta_{n+k}^{m+l} = theta_2^{0} (x) = theta_old[2][2] = theta_old[n+k][(n+k)-(m+l)]
 *                gamma_{k}^{l*}    = gamma_1^{( 1)*} (-y) = gamma_neg_y[1][0] = gamma_neg_y[k][k-l]
 *             2) k = 2, l =  1 => m + l = (-1) + ( 1) =  0 >= 0 :  theta_3^{ 0) (x)  *  gamma_2^{( 1)*} (-y)
 *                theta_{n+k}^{m+l} = theta_3^{0} (x) = theta_old[3][3] = theta_old[n+k][(n+k)-(m+l)]
 *                gamma_{k}^{l*}    = gamma_2^{( 1)*} (-y) = gamma_neg_y[2][1] = gamma_neg_y[k][k-l]
 *             3) k = 2, l =  2 => m + l = (-1) + ( 2) =  1 >= 0 :  theta_3^{ 1) (x)  *  gamma_2^{( 2)*} (-y)
 *                theta_{n+k}^{m+l} = theta_3^{1} (x) = theta_old[3][2] = theta_old[n+k][(n+k)-(m+l)]
 *                gamma_{k}^{l*}    = gamma_2^{( 2)*} (-y) = gamma_neg_y[2][0] = gamma_neg_y[k][k-l]
 *           l < 0
 *             do nothing
 *             m + l >= 0 and l< 0 implies m >= -l > 0
 *             implies that m > 0 but we are doing the case m < 0
 *         m + l < 0
 *           l >= 0
 *             1) k = 0, l =  0 => m + l = (-1) + ( 0) = -1 < 0
 *                theta_{n+k}^{m+l} = theta_1^{-1} (x) = theta_old[0]1] = theta_old[(n+k)+(m+l)][n+k]
 *                gamma_{k}^{l*}    = gamma_0^{( 0)*} (-y) = gamma_neg_y[0][0] = gamma_neg_y[k][k-l]
 *             2) k = 1, l =  0 => m + l = (-1) + ( 0) = -1 < 0 :  theta_2^{-1} (x)  *  gamma_1^{0*} (-y)
 *                theta_{n+k}^{m+l} = theta_2^{-1} (x) = theta_old[1][2] = theta_old[(n+k)+(m+l)][n+k]
 *                gamma_{k}^{l*}    = gamma_1^{( 0)*} (-y) = gamma_neg_y[1][1] = gamma_neg_y[k][k-l]
 *             3) k = 2, l =  0 => m + l = (-1) + ( 0) = -1 < 0 :  theta_3^{-1) (x)  *  gamma_2^{( 0)*} (-y)
 *                theta_{n+k}^{m+l} = theta_3^{-1} (x) = theta_old[2][3] = theta_old[(n+k)+(m+l)][n+k]
 *                gamma_{k}^{l*}    = gamma_2^{( 0)*} (-y) = gamma_neg_y[2][2] = gamma_neg_y[k][k-l]
 *           l < 0
 *             1) k = 1, l = -1 => m + l = (-1) + (-1) = -2 < 0 :  theta_2^{-2) (x)  *  gamma_1^{(-1)*} (-y)
 *                theta_{n+k}^{m+l} = theta_2^{-2} (x) = theta_old[0][2] = theta_old[(n+k)+(m+l)][n+k]
 *                gamma_{k}^{l*}    = gamma_1^{(-1)*} (-y) = gamma_neg_y[0][1] = gamma_neg_y[k+l][k]
 *             2) k = 2, l = -2 => m + l = (-1) + (-2) = -3 < 0 :  theta_3^{-3) (x)  *  gamma_2^{(-2)*} (-y)
 *                theta_{n+k}^{m+l} = theta_3^{-3} (x) = theta_old[0][3] = theta_old[(n+k)+(m+l)][n+k]
 *                gamma_{k}^{l*}    = gamma_2^{(-2)*} (-y) = gamma_neg_y[0][2] = gamma_neg_y[k+l][k]
 *             3) k = 2, l = -1 => m + l = (-1) + (-1) = -2 < 0 :  theta_3^{-2) (x)  *  gamma_2^{(-1)*} (-y)
 *                theta_{n+k}^{m+l} = theta_3^{-2} (x) = theta_old[1][3] = theta_old[(n+k)+(m+l)][n+k]
 *                gamma_{k}^{l*}    = gamma_2^{(-1)*} (-y) = gamma_neg_y[1][2] = gamma_neg_y[k+l][k]
 *
 *
 * theta_2^1 (x+y)    = \sum_{k=0}^{3-2} \sum_{l=-k}^{k} \gamma_k^{l*} (-y) theta_{n+k}^{m+l} (x)
 *                    = gamma_0^{(0)*} (-y)  theta_{2}^{1} (x) +
 *                      gamma_1^{(-1)*} (-y) theta_{3}^{0} (x) +
 *                      gamma_1^{(0)*} (-y)  theta_{3}^{1} (x) +
 *                      gamma_1^{(1)*} (-y)  theta_{3}^{2} (x)
 *
 * n = 2, m =  1
 *   n >= |m| - continue (theta_n^m != 0)
 *     m >= 0
 *       n + k >= |m + l| - continue (theta_{n+k}^{m+l} != 0)
 *         m + l >= 0
 *           l >= 0
 *             1) k = 0, l =  0 => m + l = ( 1) + ( 0) =  1 >= 0 :  theta_2^{ 1} (x)  *  gamma_0^{( 0)*} (-y)
 *                theta_{n+k}^{m+l} = theta_2^{1} (x) = theta_old[2][1] = theta_old[n+k][(n+k)-(m+l)]
 *                gamma_{k}^{l*}    = gamma_0^{( 0)*} (-y) = gamma_neg_y[0][0] = gamma_neg_y[k][k-l]
 *             2) k = 1, l =  0 => m + l = ( 1) + ( 0) =  1 >= 0 :  theta_3^{ 1} (x)  *  gamma_1^{( 0)*} (-y)
 *                theta_{n+k}^{m+l} = theta_3^{1} (x) = theta_old[3][2] = theta_old[n+k][(n+k)-(m+l)]
 *                gamma_{k}^{l*}    = gamma_1^{( 0)*} (-y) = gamma_neg_y[1][1] = gamma_neg_y[k][k-l]
 *             3) k = 1, l =  1 => m + l = ( 1) + ( 1) =  2 >= 0 :  theta_3^{ 2} (x)  *  gamma_1^{( 1)*} (-y)
 *                theta_{n+k}^{m+l} = theta_3^{2} (x) = theta_old[3][1] = theta_old[n+k][(n+k)-(m+l)]
 *                gamma_{k}^{l*}    = gamma_1^{( 1)*} (-y) = gamma_neg_y[1][0] = gamma_neg_y[k][k-l]
 *           l < 0
 *             1) k = 1, l = -1 => m + l = ( 1) + (-1) =  0 >= 0 :  theta_2^{ 1} (x)  *  gamma_1^{(-1)*} (-y)
 *                theta_{n+k}^{m+l} = theta_3^{0} (x) = theta_old[3][3] = theta_old[n+k][(n+k)-(m+l)]
 *                gamma_{k}^{l*}    = gamma_1^{(-1)*} (-y) = gamma_neg_y[0][1] = gamma_neg_y[k+l][k]
 *         m + l < 0
 *           l >= 0
 *             do nothing
 *             m + l < 0 and l >= 0 implies m < -l <= 0
 *             implies that m < 0 but we are doing the case m >= 0
 *           l < 0
 *             none
 *       n + K < |m + l|
 *         do nothing (theta_{n+k}^{m+l} = 0)
 *     m < 0
 *       n + k >= |m + l| - continue (theta_{n+k}^{m+l} != 0)
 *         m + l >= 0
 *           l >= 0
 *             none
 *           l < 0
 *             do nothing
 *             m + l >= 0 and l< 0 implies m >= -l > 0
 *             implies that m > 0 but we are doing the case m < 0
 *         m + l < 0
 *           l >= 0
 *             none
 *           l < 0
 *             none
 *
 * theta_2^1 (x+y)    = \sum_{k=0}^{3-2} \sum_{l=-k}^{k} \gamma_k^{l*} (-y) theta_{n+k}^{m+l} (x)
 *                    = gamma_0^{(0)*} (-y)  theta_{2}^{1} (x) +
 *                      gamma_1^{(-1)*} (-y) theta_{3}^{0} (x) +
 *                      gamma_1^{(0)*} (-y)  theta_{3}^{1} (x) +
 *                      gamma_1^{(1)*} (-y)  theta_{3}^{2} (x)
 *
 * n = 3, m = -2
 *   n >= |m| - continue (theta_n^m != 0)
 *     m >= 0
 *       n + k >= |m + l| - continue (theta_{n+k}^{m+l} != 0)
 *         m + l >= 0
 *           l >= 0
 *             none
 *           l < 0
 *             none
 *         m + l < 0
 *           l >= 0
 *             do nothing
 *             m + l < 0 and l >= 0 implies m < -l <= 0
 *             implies that m < 0 but we are doing the case m >= 0
 *           l < 0
 *             none
 *       n + K < |m + l|
 *         do nothing (theta_{n+k}^{m+l} = 0)
 *     m < 0
 *       n + k >= |m + l| - continue (theta_{n+k}^{m+l} != 0)
 *         m + l >= 0
 *           l >= 0
 *             1) k = 0, l =  0 => m + l = (-2) + ( 0) = -2 < 0 :  theta_3^{-2} (x)  *  gamma_0^{( 0)*} (-y)
 *           l < 0
 *             do nothing
 *             m + l >= 0 and l< 0 implies m >= -l > 0
 *             implies that m > 0 but we are doing the case m < 0
 *         m + l < 0
 *           l >= 0
 *             none
 *           l < 0
 *             none
 *
 * We write the results in matrix notation
 *  --                --      --                                                                                                                                                                                                                                        --      --           --
 * | theta_0^{0}  (x+y) |    | theta_0^{0}  & theta_1^{-1} & theta_1^{0}  & theta_1^{1}  & theta_2^{-2} & theta_2^{-1} & theta_2^{0}  & theta_2^{1}  & theta_2^{2} & theta_3^{-3} & theta_3^{-2} & theta_3^{-1} & theta_3^{0}  & theta_3^{1}  & theta_3^{2} & theta_3^{3} |     | gamma_0^0    |
 * | theta_1^{-1} (x+y) |    | theta_1^{-1} & theta_2^{-2} & theta_2^{-1} & theta_2^{0}  & theta_3^{-3} & theta_3^{-2} & theta_3^{-1} & theta_3^{0}  & theta_3^{1} & 0            & 0            & 0            & 0            & 0            & 0           & 0 		  |     | gamma_1^{-1} |
 * | theta_1^{0}  (x+y) |    | theta_1^{0}  & theta_2^{-1} & theta_2^{0}  & theta_2^{1}  & theta_3^{-2} & theta_3^{-1} & theta_3^{0}  & theta_3^{1}  & theta_3^{2} & 0            & 0            & 0            & 0            & 0            & 0           & 0           |     | gamma_1^0    |
 * | theta_1^{1}  (x+y) |    | theta_1^{1}  & theta_2^{0}  & theta_2^{1}  & theta_2^{2}  & theta_3^{-1} & theta_3^{0}  & theta_3^{1}  & theta_3^{2}  & theta_3^{3} & 0            & 0            & 0            & 0            & 0            & 0           & 0           |     | gamma_1^1    |
 * | theta_2^{-2} (x+y) |    | theta_2^{-2} & theta_3^{-3} & theta_3^{-2} & theta_3^{-1} & 0            & 0            & 0            & 0            & 0           & 0            & 0            & 0            & 0            & 0            & 0           & 0           |     | gamma_2^{-2} |
 * | theta_2^{-1} (x+y) |    | theta_2^{-1} & theta_3^{-2} & theta_3^{-1} & theta_3^{0}  & 0            & 0            & 0            & 0            & 0           & 0            & 0            & 0            & 0            & 0            & 0           & 0           |     | gamma_2^{-1} |
 * | theta_2^{ 0} (x+y) |    | theta_2^{0}  & theta_3^{-1} & theta_3^{0}  & theta_3^{1}  & 0            & 0            & 0            & 0            & 0           & 0            & 0            & 0            & 0            & 0            & 0           & 0           |     | gamma_2^0    |
 * | theta_2^{ 1} (x+y) |    | theta_2^{1}  & theta_3^{-0} & theta_3^{1}  & theta_3^{2}  & 0            & 0            & 0            & 0            & 0           & 0            & 0            & 0            & 0            & 0            & 0           & 0           |     | gamma_2^1    |
 * | theta_2^{ 2} (x+y) |    | theta_2^{2}  & theta_3^{1}  & theta_3^{2}  & theta_3^{3}  & 0            & 0            & 0            & 0            & 0           & 0            & 0            & 0            & 0            & 0            & 0           & 0           |     | gamma_2^2    |
 * | theta_3^{-3} (x+y) |  = | theta_3^{-3} & 0            & 0            & 0            & 0            & 0            & 0            & 0            & 0           & 0            & 0            & 0            & 0            & 0            & 0           & 0           |     | gamma_3^{-3} |
 * | theta_3^{-2} (x+y) |    | theta_3^{-2} & 0            & 0            & 0            & 0            & 0            & 0            & 0            & 0           & 0            & 0            & 0            & 0            & 0            & 0           & 0           |     | gamma_3^{-2} |
 * | theta_3^{-1} (x+y) |    | theta_3^{-1} & 0            & 0            & 0            & 0            & 0            & 0            & 0            & 0           & 0            & 0            & 0            & 0            & 0            & 0           & 0           |     | gamma_3^{-1} |
 * | theta_3^{0}  (x+y) |    | theta_3^{0}  & 0            & 0            & 0            & 0            & 0            & 0            & 0            & 0           & 0            & 0            & 0            & 0            & 0            & 0           & 0           |     | gamma_3^0    |
 * | theta_3^{1}  (x+y) |    | theta_3^{1}  & 0            & 0            & 0            & 0            & 0            & 0            & 0            & 0           & 0            & 0            & 0            & 0            & 0            & 0           & 0           |     | gamma_3^1    |
 * | theta_3^{2}  (x+y) |    | theta_3^{2}  & 0            & 0            & 0            & 0            & 0            & 0            & 0            & 0           & 0            & 0            & 0            & 0            & 0            & 0           & 0           |     | gamma_3^2    |
 * | theta_3^{3}  (x+y) |    | theta_3^{3}  & 0            & 0            & 0            & 0            & 0            & 0            & 0            & 0           & 0            & 0            & 0            & 0            & 0            & 0           & 0           |     | gamma_3^3    |
 *  --                --      --                                                                                                                                                                                                                                        --      --           _--
 *
 * Let's go through the cases 0 <= n <= 3 and -3 <= m <= 3
 *
 * The previous location of the series that is being RR shifted is x = y - y_c^{(c)} = x_b - z_B^{(c)}
 * Therefore, the old coefficients theta_n^m (x) are passed into the function as the argument
 * const std::vector<std::vector<std::double> > > &rCoeff_old
 *
 *   \theta_n^{m} (x+y) = \sum_{k=0}^{\infty} \sum_{l=-k}^{k} \gamma_{k}^{(l)*} (-y)  theta_{n+k}^{m+l} (x)
 *                      ~ \sum_{k=0}^{p-n} \sum_{l=-k}^{k} \gamma_{k}^{(l)*} (-y)  theta_{n+k}^{m+l} (x)
 *
 * We will need to calculate the coefficients \gamma_{k}^{(l)*} (-y) for the value -y = z_B^{(c)} - z_B^{(p)}
 * As in the member function getSS where we use getSCoeff, we will use getRCoeff to calculate the matrix
 * for \gamma_{k}^{(l)*} (-y) where -y = z_B^{(c)} - z_B^{(p)}
 */
std::vector<std::vector<std::complex<double> > > Potential::getRR(std::vector<double> from, std::vector<double> to,
		                                                          const std::vector<std::vector<std::complex<double> > > &rCoeff_old)
{
  // from - refers to the center where the old coefficients were evaluated
  //        z_B^{(c)} - the child cell center z_B^{(c)} of the target point
  // to   - refers to the center where the new coefficients will be evaluated
  //        z_B^{(p)} - the parent cell center z_B^{(p)} of the target point
  // t    - refers to y = z_B^{(c)} - z_B^{(p)}

  // need to assert from and to have 3 components

  // terms used to approximate the Taylor series
////  p = 16;
  // highest order of derivative used to approximate Taylor series, to which p corresponds
////  order_of_approximation = 3;

//  for (unsigned int i=0; i<4; ++i)
//    for (unsigned int j=0; j<4; ++j)
//      std::cout << "rCoeff_old[" << i << "][" << j << "] = "
//                << rCoeff_old[i][j] << std::endl;

  // note - Taylor series includes zeroth derivative - therefore adding 1 to highest order derivative
  std::vector<std::vector<std::complex<double> > > ans(order_of_approximation+1, std::vector<std::complex<double> >(order_of_approximation+1));

// Translation formula for \gamma_n^{m*} is
// Here x = from and x+y = to and t = to - from = y
//  std::vector<double>  t(3);
//  t[0] = from[0] - to[0];  t[1] = from[1] - to[1]; t[2] = from[2] - to[2];


  // Here we build the \gamma_n^{(m)*} (-y) = \gamma_n^{m*} (z_B^{(p) - z_B^{(c)}}) matrix
  // Recall getSCoeff(xi, xstar) performs its operations using the difference (xi - xstar).
  // Therefore, we let to = z_B^{(c)} = xstar and from = z_B^{(p)} = xi
  // Recall the matrix dimensions
  // std::vector<std::vector<std::complex<double> > >
  //                  gamma_y_matrix(order_of_approximation+1, std::vector<std::complex<double> >(order_of_approximation+1));
  std::vector<std::vector<std::complex<double> > > gamma_y_matrix = this->getSCoeff(from, to);

//  for (unsigned int i=0; i<4; ++i)
//    for (unsigned int j=0; j<4; ++j)
//      std::cout << "gamma_y_matrix[" << i << "][" << j << "] = "
//                << gamma_y_matrix[i][j] << std::endl;


  for (unsigned int j = 0; j<=order_of_approximation; ++j)
    for (unsigned int i = 0; i<=order_of_approximation; ++i)
      ans[i][j] = 0.0;

  // Generating the matrix \gamma_n^{(m)*} (-y) = \theta_n^{m} (x+y) = theta_n^m (x_b - z_B^{(p)}})
  // using the formula
  //
  //   \theta_n^{m} (x+y) = \sum_{k=0}^{order_of_approximation - n}
  //                                  \sum_{l=-k}^{k} \gamma_{k}^{(l)*} (-y)  theta_{n+k}^{m+l} (x)
  //
  for (unsigned int n = 0; n <= order_of_approximation; ++n)
    for (int m = -n; m<=int(n); ++m)
    {
      if (int(n) >= std::abs(m)) // assemble an element, otherwise n < m implies element is zero
      {
    	// theta_n^m matrix has entries only up to order_of_approximation
    	// Therefore, n+k can only be up to order_of_approximation => n+k <= order_of_approximation
    	//                    => k <= order_of_approximation - n
    	// or else we will be looking for theta_{n+k}^{m+l} that are not there (indexing outside of the theta matrix)
        for (unsigned int k=0; k<=(order_of_approximation - n); ++k)
          for (int l=-k; l<=int(k); ++l)
          {
            if ( m >= 0 ) // setting lower triangular part of gamma matrix ans -> gamma_n^{(m)*} (x+y)
            {
              if ( int(n+k) >= std::abs(m+l) )  // gamma_{n+k}^{(m+l)*} not zero
              {
                if ( (m+l) >= 0)  // indexing lower triangular part of theta matrix rCoeff_old -> theta_{n+k}^{(m+l)} (x)
                {
                  if ( l >= 0 )  // indexing lower triangular part of gamma matrix gamma_y_matrix -> gamma_k^{l*} (-y)
                	             // m >=0 and (m+l) >= 0 and l >= 0
                	             // implies that m >= -l
                  {
            	    // ans[n][n-m] is the location of theta_n^{(m)*} (x+y)
          	        // rCoeff_old[n+k][(n+k)-(m+l)] is the location of theta_{n+k}^{(m+l)} (x)
          	        // gamma_y_matrix[k][k-l] is the location of gamma_{k}^{(l)*} (-y)
                    ans[n][n-m] += rCoeff_old[n+k][(n+k)-(m+l)] * gamma_y_matrix[k][k-l];
//if(m==0 && n==0)
//{
//    std::cout << "m = 0, n = 0, k = " << k << " l = " << l << std::endl;
//    std::cout << "gamma_y_matrix = " << gamma_y_matrix[k][k-l] << std::endl;
//    std::cout << "rCoeff_old = " << rCoeff_old[n+k][(n+k)-(m+l)] << std::endl;
//}
                  }
                  else // l < 0 // indexing upper triangular part of gamma matrix gamma_y_matrix -> gamma_k^{l*} (-y)
                	            // (m + l) >= 0 and l < 0 => m >= -l > 0
                  {
                    // EX: ans[n][n-m] += rCoeff_old[n+k][(n+k)-(m+l)] * gamma_y_matrix[k+l][k];
                	//     n = 2, m = 1
                	//       k = 1, l = -1
                    //         theta_2^1(x+y) += theta_3^0(x) * gamma_1^{-1}(-y)
                	//         theta_new[2][1] += theta_old[3][3] * gamma_[0][1]
                	//         ans[n][2-1] += rCoeff_old[n+k][n+k-(m+l)] * gamma_y_matrix[k+l][k]
                  	//     n = 2, m = 2
                  	//       k = 1, l = -1
                      //         theta_2^2(x+y) += theta_3^1(x) * gamma_1^{-1}(-y)
                  	//         theta_new[2][1] += theta_old[3][2] * gamma_[0][1]
                  	//         ans[n][2-1] += rCoeff_old[n+k][n+k-(m+l)] * gamma_y_matrix[k+l][k]
                    ans[n][n-m] += rCoeff_old[n+k][(n+k)-(m+l)] * gamma_y_matrix[k+l][k];
//if(m==0 && n==0)
//{
//                      std::cout << "m = 0, n = 0, k = " << k << " l = " << l << std::endl;
//                      std::cout << "gamma_y_matrix = " << gamma_y_matrix[k+l][k] << std::endl;
//                      std::cout << "rCoeff_old = " << rCoeff_old[n+k][(n+k)-(m+l)] << std::endl;
//}
                  }
                }
                else // (m+l) < 0 and indexing upper triangular part of theta matrix rCoeff_old -> theta_{n+k}^{(m+l)*} (x)
                {
                  if ( l >= 0 )  // indexing lower triangular part of gamma matrix gamma_y_matrix -> gamma_k^{l*} (-y)
                  {
                	// do nothing: l >= 0 and (m+l) < 0  =>  m < -l < 0
                	// m has to be negative and we are doing case m >= 0 (outer if condition)
              	        // ans[n-m][n] is the location of gamma_n^{(m)*} (x+y)
            	        // sCoeff_old[n-k][(n-k)-(m-l)] is the location of gamma_{n-k}^{(m-l)*} (x)
            	        // gamma_y_matrix[k][k-l] is the location of gamma_{k}^{(l)*} (y)
                        // ans[n-m][n] += rCoeff_old[(n+k)+(m+l)][n+k] * gamma_y_matrix[k][k-l];
                  }
                  else // l < 0 // indexing upper triangular part of gamma matrix gamma_y_matrix -> gamma_k^{l*} (-y)
                  {
                    // EX: ans[n][n-m] += rCoeff_old[n+k][(n+k)+(m+l)] * gamma_y_matrix[k+l][k];
                  	//     n = 2, m = 1
                  	//       k = 1, l = -2
                    //         theta_2^1(x+y) += theta_3^0(x) * gamma_1^{-1}(-y)
                  	//         theta_new[2][1] += theta_old[3][3] * gamma_[0][1]
                  	//         ans[n][2-1] += rCoeff_old[n+k][n+k-(m+l)] * gamma_y_matrix[k+l][k]
                    	//     n = 2, m = 0
                    	//       k = 1, l = -1
                        //         theta_2^0(x+y) += theta_3^{-1}(x) * gamma_1^{-1}(-y)
                    	//         theta_new[2][2] += theta_old[2][3] * gamma_[0][1]
                    	//         ans[n][2-1] += rCoeff_old[n+k][n+k-(m+l)] * gamma_y_matrix[k+l][k]
                    ans[n][n-m] += rCoeff_old[(n+k)+(m+l)][n+k] * gamma_y_matrix[k+l][k];
//if(m==0 && n==0)
//{
//    std::cout << "m = 0, n = 0, k = " << k << " l = " << l << std::endl;
//    std::cout << "gamma_y_matrix = " << gamma_y_matrix[k+l][k] << std::endl;
//    std::cout << "rCoeff_old = " << rCoeff_old[(n+k)+(m+l)][n+k] << std::endl;
//}
                  }
                }
              }
              else // (n-k) < std::abs(m-l)
              {
                //  no need to do anything since gamma_{n-k}^{(m-l)*} = zero
              }
            }
            else // m < 0 - setting lower triangular part of gamma matrix ans -> gamma_n^{(m)*} (x+y)
            {
              if ( int(n+k) >= std::abs(m+l) )  // gamma_{n-k}^{(m-l)*} not zero
              {
                if ( (m+l) >= 0)  // indexing lower triangular part of gamma matrix sCoeff_old -> gamma_{n-k}^{(m-l)*} (x)
                {
                  if ( l >= 0 )  // indexing lower triangular part of gamma matrix gamma_y_matrix -> gamma_k^{l*} (-y)
                  {
           	        // ans[n][n-m] is the location of theta_n^{m} (x+y)
        	        // rCoeff_old[n+k][(n+k)-(m+l)] is the location of theta_{n+k}^{(m+l)*} (x)
        	        // gamma_y_matrix[k][k-l] is the location of gamma_{k}^{(l)*} (-y)
                    ans[n+m][n] += rCoeff_old[n+k][(n+k)-(m+l)] * gamma_y_matrix[k][k-l];
//if(m==0 && n==0)
//{
//    std::cout << "m = 0, n = 0, k = " << k << " l = " << l << std::endl;
//    std::cout << "gamma_y_matrix = " << gamma_y_matrix[k][k-l] << std::endl;
//    std::cout << "rCoeff_old = " << rCoeff_old[n+k][(n+k)-(m+l)] << std::endl;
//}
                  }
                  else // l < 0 - indexing upper triangular part of gamma matrix gamma_y_matrix -> gamma_k^{l*} (-y)
                  {
                  	// do nothing: l < 0 and (m+l) >= 0  =>  m >= -l > 0
                  	// m is positive but we are doing the case m < 0 (outer if condition)
                        // ans[n+m][n] += rCoeff_old[n+k][(n+k)-(m+l)] * gamma_y_matrix[k+l][k];
                  }
                }
                else // (m+l) < 0 and indexing upper triangular part of theta matrix rCoeff_old -> theta_{n+k}^{(m+l)} (x)
                {
                  if ( l >= 0 )  // indexing lower triangular part of gamma matrix gamma_y_matrix -> gamma_k^{l*} (-y)
                  {
            	    // ans[n][n-m] is the location of gamma_n^{(m)*} (x+y)
          	        // sCoeff_old[n-k][(n-k)-(m-l)] is the location of gamma_{n-k}^{(m-l)*} (x)
          	        // gamma_y_matrix[k][k-l] is the location of gamma_{k}^{(l)*} (y)
                    ans[n+m][n] += rCoeff_old[(n+k)+(m+l)][n+k] * gamma_y_matrix[k][k-l];
//if(m==0 && n==0)
//{
//    std::cout << "m = 0, n = 0, k = " << k << " l = " << l << std::endl;
//    std::cout << "gamma_y_matrix = " << gamma_y_matrix[k][k-l] << std::endl;
//    std::cout << "rCoeff_old = " << rCoeff_old[(n+k)+(m+l)][n+k] << std::endl;
//}
                  }
                  else // l < 0 // indexing upper triangular part of gamma matrix gamma_y_matrix -> gamma_k^{l*} (-y)
                  {
                    ans[n+m][n] += rCoeff_old[(n+k)+(m+l)][n+k] * gamma_y_matrix[k+l][k];
//if(m==0 && n==0)
//{
//    std::cout << "m = 0, n = 0, k = " << k << " l = " << l << std::endl;
//    std::cout << "gamma_y_matrix = " << gamma_y_matrix[k+l][k] << std::endl;
//    std::cout << "rCoeff_old = " << rCoeff_old[(n+k)+(m+l)][n+k] << std::endl;
//}
                  }
                }
              }
              else // (n-k) < std::abs(m-l)
              {
                //  no need to do anything since gamma_{n-k}^{(m-l)*} = zero
              }
            }
          }
      }
    }

  return ans;

}



/*
std::vector<std::vector<std::vector<std::complex<double> > > > Potential::getRRgrad(std::vector<double> from, std::vector<double> to,
 		                                                               const std::vector<std::vector<std::vector<std::complex<double> > > > &rCoeff_old)
 {
   // from - refers to the center where the old coefficients were evaluated
   //        z_B^{(c)} - the child cell center z_B^{(c)} of the target point
   // to   - refers to the center where the new coefficients will be evaluated
   //        z_B^{(p)} - the parent cell center z_B^{(p)} of the target point
   // t    - refers to y = z_B^{(c)} - z_B^{(p)}

   // need to assert from and to have 3 components

   // terms used to approximate the Taylor series
 ////  p = 16;
   // highest order of derivative used to approximate Taylor series, to which p corresponds
 ////  order_of_approximation = 3;

 //  for (unsigned int i=0; i<4; ++i)
 //    for (unsigned int j=0; j<4; ++j)
 //      std::cout << "rCoeff_old[" << i << "][" << j << "] = "
 //                << rCoeff_old[i][j] << std::endl;

   // note - Taylor series includes zeroth derivative - therefore adding 1 to highest order derivative
   std::vector<std::vector<std::vector<std::complex<double> > > >
       ans(order_of_approximation+1, std::vector<std::vector<std::complex<double> > >(order_of_approximation+1, std::vector<std::complex<double> >(3)));

 // Translation formula for \gamma_n^{m*} is
 // Here x = from and x+y = to and t = to - from = y
 //  std::vector<double>  t(3);
 //  t[0] = from[0] - to[0];  t[1] = from[1] - to[1]; t[2] = from[2] - to[2];


   // Here we build the \gamma_n^{(m)*} (-y) = \gamma_n^{m*} (z_B^{(p) - z_B^{(c)}}) matrix
   // Recall getSCoeff(xi, xstar) performs its operations using the difference (xi - xstar).
   // Therefore, we let to = z_B^{(c)} = xstar and from = z_B^{(p)} = xi
   // Recall the matrix dimensions
   // std::vector<std::vector<std::complex<double> > >
   //                  gamma_y_matrix(order_of_approximation+1, std::vector<std::complex<double> >(order_of_approximation+1));
   std::vector<std::vector<std::complex<double> > > gamma_y_matrix = this->getSCoeff(from, to);

 //  for (unsigned int i=0; i<4; ++i)
 //    for (unsigned int j=0; j<4; ++j)
 //      std::cout << "gamma_y_matrix[" << i << "][" << j << "] = "
 //                << gamma_y_matrix[i][j] << std::endl;


   for (unsigned int j = 0; j<=order_of_approximation; ++j)
     for (unsigned int i = 0; i<=order_of_approximation; ++i)
       for (unsigned int k = 0; k<3; ++k)
         ans[i][j][k] = 0.0;

   // Generating the matrix \gamma_n^{(m)*} (-y) = \theta_n^{m} (x+y) = theta_n^m (x_b - z_B^{(p)}})
   // using the formula
   //
   //   \theta_n^{m} (x+y) = \sum_{k=0}^{order_of_approximation - n}
   //                                  \sum_{l=-k}^{k} \gamma_{k}^{(l)*} (-y)  theta_{n+k}^{m+l} (x)
   //
   for (unsigned int n = 0; n <= order_of_approximation; ++n)
     for (int m = -n; m<=int(n); ++m)
     {
       if (int(n) >= std::abs(m)) // assemble an element, otherwise n < m implies element is zero
       {
     	// theta_n^m matrix has entries only up to order_of_approximation
     	// Therefore, n+k can only be up to order_of_approximation => n+k <= order_of_approximation
     	//                    => k <= order_of_approximation - n
     	// or else we will be looking for theta_{n+k}^{m+l} that are not there (indexing outside of the theta matrix)
         for (unsigned int k=0; k<=(order_of_approximation - n); ++k)
           for (int l=-k; l<=int(k); ++l)
           {
             if ( m >= 0 ) // setting lower triangular part of gamma matrix ans -> gamma_n^{(m)*} (x+y)
             {
               if ( int(n+k) >= std::abs(m+l) )  // gamma_{n+k}^{(m+l)*} not zero
               {
                 if ( (m+l) >= 0)  // indexing lower triangular part of theta matrix rCoeff_old -> theta_{n+k}^{(m+l)} (x)
                 {
                   if ( l >= 0 )  // indexing lower triangular part of gamma matrix gamma_y_matrix -> gamma_k^{l*} (-y)
                 	             // m >=0 and (m+l) >= 0 and l >= 0
                 	             // implies that m >= -l
                   {
             	    // ans[n][n-m] is the location of theta_n^{(m)*} (x+y)
           	        // rCoeff_old[n+k][(n+k)-(m+l)] is the location of theta_{n+k}^{(m+l)} (x)
           	        // gamma_y_matrix[k][k-l] is the location of gamma_{k}^{(l)*} (-y)
                     ans[n][n-m][0] += rCoeff_old[n+k][(n+k)-(m+l)][0] * gamma_y_matrix[k][k-l];
                     ans[n][n-m][1] += rCoeff_old[n+k][(n+k)-(m+l)][1] * gamma_y_matrix[k][k-l];
                     ans[n][n-m][2] += rCoeff_old[n+k][(n+k)-(m+l)][2] * gamma_y_matrix[k][k-l];
 //if(m==0 && n==0)
 //{
 //    std::cout << "m = 0, n = 0, k = " << k << " l = " << l << std::endl;
 //    std::cout << "gamma_y_matrix = " << gamma_y_matrix[k][k-l] << std::endl;
 //    std::cout << "rCoeff_old = " << rCoeff_old[n+k][(n+k)-(m+l)] << std::endl;
 //}
                   }
                   else // l < 0 // indexing upper triangular part of gamma matrix gamma_y_matrix -> gamma_k^{l*} (-y)
                 	            // (m + l) >= 0 and l < 0 => m >= -l > 0
                   {
                     // EX: ans[n][n-m] += rCoeff_old[n+k][(n+k)-(m+l)] * gamma_y_matrix[k+l][k];
                 	//     n = 2, m = 1
                 	//       k = 1, l = -1
                     //         theta_2^1(x+y) += theta_3^0(x) * gamma_1^{-1}(-y)
                 	//         theta_new[2][1] += theta_old[3][3] * gamma_[0][1]
                 	//         ans[n][2-1] += rCoeff_old[n+k][n+k-(m+l)] * gamma_y_matrix[k+l][k]
                   	//     n = 2, m = 2
                   	//       k = 1, l = -1
                       //         theta_2^2(x+y) += theta_3^1(x) * gamma_1^{-1}(-y)
                   	//         theta_new[2][1] += theta_old[3][2] * gamma_[0][1]
                   	//         ans[n][2-1] += rCoeff_old[n+k][n+k-(m+l)] * gamma_y_matrix[k+l][k]
                     ans[n][n-m][0] += rCoeff_old[n+k][(n+k)-(m+l)][0] * gamma_y_matrix[k+l][k];
                     ans[n][n-m][1] += rCoeff_old[n+k][(n+k)-(m+l)][1] * gamma_y_matrix[k+l][k];
                     ans[n][n-m][2] += rCoeff_old[n+k][(n+k)-(m+l)][2] * gamma_y_matrix[k+l][k];
 //if(m==0 && n==0)
 //{
 //                      std::cout << "m = 0, n = 0, k = " << k << " l = " << l << std::endl;
 //                      std::cout << "gamma_y_matrix = " << gamma_y_matrix[k+l][k] << std::endl;
 //                      std::cout << "rCoeff_old = " << rCoeff_old[n+k][(n+k)-(m+l)] << std::endl;
 //}
                   }
                 }
                 else // (m+l) < 0 and indexing upper triangular part of theta matrix rCoeff_old -> theta_{n+k}^{(m+l)*} (x)
                 {
                   if ( l >= 0 )  // indexing lower triangular part of gamma matrix gamma_y_matrix -> gamma_k^{l*} (-y)
                   {
                 	// do nothing: l >= 0 and (m+l) < 0  =>  m < -l < 0
                 	// m has to be negative and we are doing case m >= 0 (outer if condition)
               	        // ans[n-m][n] is the location of gamma_n^{(m)*} (x+y)
             	        // sCoeff_old[n-k][(n-k)-(m-l)] is the location of gamma_{n-k}^{(m-l)*} (x)
             	        // gamma_y_matrix[k][k-l] is the location of gamma_{k}^{(l)*} (y)
                         // ans[n-m][n] += rCoeff_old[(n+k)+(m+l)][n+k] * gamma_y_matrix[k][k-l];
                   }
                   else // l < 0 // indexing upper triangular part of gamma matrix gamma_y_matrix -> gamma_k^{l*} (-y)
                   {
                     // EX: ans[n][n-m] += rCoeff_old[n+k][(n+k)+(m+l)] * gamma_y_matrix[k+l][k];
                   	//     n = 2, m = 1
                   	//       k = 1, l = -2
                     //         theta_2^1(x+y) += theta_3^0(x) * gamma_1^{-1}(-y)
                   	//         theta_new[2][1] += theta_old[3][3] * gamma_[0][1]
                   	//         ans[n][2-1] += rCoeff_old[n+k][n+k-(m+l)] * gamma_y_matrix[k+l][k]
                     	//     n = 2, m = 0
                     	//       k = 1, l = -1
                         //         theta_2^0(x+y) += theta_3^{-1}(x) * gamma_1^{-1}(-y)
                     	//         theta_new[2][2] += theta_old[2][3] * gamma_[0][1]
                     	//         ans[n][2-1] += rCoeff_old[n+k][n+k-(m+l)] * gamma_y_matrix[k+l][k]
                     ans[n][n-m][0] += rCoeff_old[(n+k)+(m+l)][n+k][0] * gamma_y_matrix[k+l][k];
                     ans[n][n-m][1] += rCoeff_old[(n+k)+(m+l)][n+k][1] * gamma_y_matrix[k+l][k];
                     ans[n][n-m][2] += rCoeff_old[(n+k)+(m+l)][n+k][2] * gamma_y_matrix[k+l][k];
 //if(m==0 && n==0)
 //{
 //    std::cout << "m = 0, n = 0, k = " << k << " l = " << l << std::endl;
 //    std::cout << "gamma_y_matrix = " << gamma_y_matrix[k+l][k] << std::endl;
 //    std::cout << "rCoeff_old = " << rCoeff_old[(n+k)+(m+l)][n+k] << std::endl;
 //}
                   }
                 }
               }
               else // (n-k) < std::abs(m-l)
               {
                 //  no need to do anything since gamma_{n-k}^{(m-l)*} = zero
               }
             }
             else // m < 0 - setting lower triangular part of gamma matrix ans -> gamma_n^{(m)*} (x+y)
             {
               if ( int(n+k) >= std::abs(m+l) )  // gamma_{n-k}^{(m-l)*} not zero
               {
                 if ( (m+l) >= 0)  // indexing lower triangular part of gamma matrix sCoeff_old -> gamma_{n-k}^{(m-l)*} (x)
                 {
                   if ( l >= 0 )  // indexing lower triangular part of gamma matrix gamma_y_matrix -> gamma_k^{l*} (-y)
                   {
            	        // ans[n][n-m] is the location of theta_n^{m} (x+y)
         	        // rCoeff_old[n+k][(n+k)-(m+l)] is the location of theta_{n+k}^{(m+l)*} (x)
         	        // gamma_y_matrix[k][k-l] is the location of gamma_{k}^{(l)*} (-y)
                     ans[n+m][n][0] += rCoeff_old[n+k][(n+k)-(m+l)][0] * gamma_y_matrix[k][k-l];
                     ans[n+m][n][1] += rCoeff_old[n+k][(n+k)-(m+l)][1] * gamma_y_matrix[k][k-l];
                     ans[n+m][n][2] += rCoeff_old[n+k][(n+k)-(m+l)][2] * gamma_y_matrix[k][k-l];
 //if(m==0 && n==0)
 //{
 //    std::cout << "m = 0, n = 0, k = " << k << " l = " << l << std::endl;
 //    std::cout << "gamma_y_matrix = " << gamma_y_matrix[k][k-l] << std::endl;
 //    std::cout << "rCoeff_old = " << rCoeff_old[n+k][(n+k)-(m+l)] << std::endl;
 //}
                   }
                   else // l < 0 - indexing upper triangular part of gamma matrix gamma_y_matrix -> gamma_k^{l*} (-y)
                   {
                   	// do nothing: l < 0 and (m+l) >= 0  =>  m >= -l > 0
                   	// m is positive but we are doing the case m < 0 (outer if condition)
                         // ans[n+m][n] += rCoeff_old[n+k][(n+k)-(m+l)] * gamma_y_matrix[k+l][k];
                   }
                 }
                 else // (m+l) < 0 and indexing upper triangular part of theta matrix rCoeff_old -> theta_{n+k}^{(m+l)} (x)
                 {
                   if ( l >= 0 )  // indexing lower triangular part of gamma matrix gamma_y_matrix -> gamma_k^{l*} (-y)
                   {
             	    // ans[n][n-m] is the location of gamma_n^{(m)*} (x+y)
           	        // sCoeff_old[n-k][(n-k)-(m-l)] is the location of gamma_{n-k}^{(m-l)*} (x)
           	        // gamma_y_matrix[k][k-l] is the location of gamma_{k}^{(l)*} (y)
                     ans[n+m][n][0] += rCoeff_old[(n+k)+(m+l)][n+k][0] * gamma_y_matrix[k][k-l];
                     ans[n+m][n][1] += rCoeff_old[(n+k)+(m+l)][n+k][1] * gamma_y_matrix[k][k-l];
                     ans[n+m][n][2] += rCoeff_old[(n+k)+(m+l)][n+k][2] * gamma_y_matrix[k][k-l];
 //if(m==0 && n==0)
 //{
 //    std::cout << "m = 0, n = 0, k = " << k << " l = " << l << std::endl;
 //    std::cout << "gamma_y_matrix = " << gamma_y_matrix[k][k-l] << std::endl;
 //    std::cout << "rCoeff_old = " << rCoeff_old[(n+k)+(m+l)][n+k] << std::endl;
 //}
                   }
                   else // l < 0 // indexing upper triangular part of gamma matrix gamma_y_matrix -> gamma_k^{l*} (-y)
                   {
                     ans[n+m][n][0] += rCoeff_old[(n+k)+(m+l)][n+k][0] * gamma_y_matrix[k+l][k];
                     ans[n+m][n][1] += rCoeff_old[(n+k)+(m+l)][n+k][1] * gamma_y_matrix[k+l][k];
                     ans[n+m][n][2] += rCoeff_old[(n+k)+(m+l)][n+k][2] * gamma_y_matrix[k+l][k];
 //if(m==0 && n==0)
 //{
 //    std::cout << "m = 0, n = 0, k = " << k << " l = " << l << std::endl;
 //    std::cout << "gamma_y_matrix = " << gamma_y_matrix[k+l][k] << std::endl;
 //    std::cout << "rCoeff_old = " << rCoeff_old[(n+k)+(m+l)][n+k] << std::endl;
 //}
                   }
                 }
               }
               else // (n-k) < std::abs(m-l)
               {
                 //  no need to do anything since gamma_{n-k}^{(m-l)*} = zero
               }
             }
           }
       }
     }

   return ans;

 }
*/

/**
 *  Explanation of getRVector
 */
// getRVector member function returns the powers of the R-expansion power series.
// The member function returns the powers \gamma_n^{m*} (x_b - z_B) of the near-field R-expansion
// The coefficients for the near-field expansion in Cartesian coordinates are the derivatives
// D^{\alpha} f(x_a - z_a) / {alpha !} and the powers are (z_a - x_b)^{\alpha} in
//
//  f(x_b - x_a) = sum_{|alpha| >= 0} [ ( (1 / (alpha !)) (x_b - z_b)^{alpha} ) D^{alpha} f(z_b - x_a) ]
//
//  The series converges when || x_b - z_b || < || z_b - x_a || (or || y - x^t_* || < || x^t_* - x_i ||=|| y - x_* ||)
//
//  In Dehnen's coordinate system the same series has the form
//
//     \frac{1}{|x_a - x_b|} = \frac{1}{|(x_b - z_B) - (x_a - z_B)|}
//                            = \sum_{n=0}^{\infty} \sum_{m=-n}^{n} \gamma_n^{m*} (x_b - z_B) \theta_n^{m} (x_a - z_B)
//
//  and the coefficients of the series are \theta_n^{m} (x_a - z_B) and the powers are \gamma_n^{m*} (x_b - z_B)

// Here x_b = y is the target point and z_B = xstar is the target cell center.  We use p = 16 to obtain all
// derivatives up to and including order 3.
//

 std::vector<std::vector<std::complex<double> > >  Potential::getRVector(std::vector<double> y, std::vector<double> xstar)
{
  // Here y              - x_b is the target particle
  //      xstar          - z_B is the target cell center

  // terms used to approximation the Taylor series
////  p = 16;
  // highest order of derivative used to approximate Taylor series, to which p corresponds
////  order_of_approximation = 3;
  // matrix of coefficients that will be returned
  // note - Taylor series includes zeroth derivative - therefore adding 1 to highest order
  std::vector<std::vector<std::complex<double> > > ans(order_of_approximation+1, std::vector<std::complex<double> >(order_of_approximation+1));
  // norm_r_squared - ||x_b - z_B||^2
  double norm_r_squared = std::pow(y[0]-xstar[0],2.0) + std::pow(y[1]-xstar[1],2.0) + std::pow(y[2]-xstar[2],2.0);

  // **********************************************************************************************************
  // Setting the coordinates xi and eta in coordinate system xi, eta, z **************************************
  // **********************************************************************************************************
  // **********************************************************************************************************
  // Recall the coordinate \eta = n =  -0.5 * (x - iy)
  //        the coordinate \xi  = e =   0.5 * (x + iy)
  //        the target particle is x_b
  //        the target cell center is z_B
  std::complex<double> eta_x_b, xi_x_b;
  double x_b_z_component;
  eta_x_b.real(-0.5*y[0]);  eta_x_b.imag(0.5*y[1]);
  xi_x_b.real(0.5*y[0]);    xi_x_b.imag(0.5*y[1]);
  x_b_z_component = y[2];

  std::complex<double> eta_z_B, xi_z_B;
  double z_B_z_component;
  eta_z_B.real(-0.5*xstar[0]);  eta_z_B.imag(0.5*xstar[1]);
  xi_z_B.real(0.5*xstar[0]);    xi_z_B.imag(0.5*xstar[1]);
  z_B_z_component = xstar[2];

  std::complex<double> difference_eta = eta_z_B - eta_x_b;
  std::complex<double> difference_xi  = xi_z_B - xi_x_b;
  double difference_z = z_B_z_component - x_b_z_component;

  // **********************************************************************************************************
  // Recalling the information on the formulas for coefficients gamma_n^{m*} **********************************
  // **********************************************************************************************************
  // **********************************************************************************************************
  // We use p = 16 to obtain all derivatives up to and including order 3.
  // The member function Potential::getRVector returns the powers \gamma_n^{m*} (x_b - z_B) of the
  // near-field expansion \frac{1}{|r-x|} = \sum_{n=0}^{\infty} \gamma_n^{m*} (r) \theta_{n}^{m} (x)
  // Then
  //                     \frac{1}{|x_a - x_b|} = \frac{1}{|(x_b - z_B) - (x_a - z_B)|}
  //                                           = \sum_{n=0}^{\infty} \sum_{m=-n}^{n} \gamma_n^{m*} (x_b - z_B) \theta_n^{m*} (x_a - z_B)
  //
  // Recall the formulas for \gamma_n^{m*} (\mathbf{r})
  //
  //  (0)   \gamma_0^{0*}  (r) = 1
  //
  //  (1)   \gamma_n^{n*}  (r) = (-1) \frac{\eta}{n} \gamma_{n-1}^{(n-1)*}   for 1 \leq n
  //
  //  (2)   \gamma_n^{-n*} (r) = (-1) \frac{\xi}{n}  \gamma_{n-1}^{-(n-1)*}  for 1 \leq n
  //
  //  (3)   \gamma_n^{m*}  (r) = \frac{(2n-1)z}{n^2 - m^2} \gamma_{n-1}^{m*} -  \frac{|r|^2}{n^2 - m^2} \gamma_{n-2}^{m*}   for m \geq 0
  //
  //  (4)   \gamma_n^{-m*} (r) = \frac{(2n-1)z}{n^2 - m^2} \gamma_{n-1}^{-m*} - \frac{|r|^2}{n^2 - m^2} \gamma_{n-2}^{-m*}  for m \geq 0

  // The coefficients \gamma_n^{m*} (x_a - z_A) may be complex as can be seen from the formulas above
  // Therefore, the vector of coefficients that is returned - ans - is a complex vector
  // Further, the coefficients correspond to a full matrix shown below
  //
  // --                                                                                                                  --
  // | \gamma_0^0*    \gamma_1^{-1*}    \gamma_2^{-2*}    \gamma_3^{-3*}    \gamma_4^{-4*}    \hdots      \gamma_n^{-n*}   |
  // |                                                                                                                     |
  // | \gamma_1^1*    \gamma_1^0*       \gamma_2^{-1*}    \gamma_3^{-2*}    \gamma_4^{-3*}                \gamma_n^{-n+1*} |
  // |                                                                                                                     |
  // | \gamma_2^2*    \gamma_2^1*       \gamma_2^0*       \gamma_3^{-1*}    \gamma_4^{-2*}                \gamma_n^{-n+2*} |
  // |                                                                                                                     |
  // | \gamma_3^3*    \gamma_3^2*       \gamma_3^1*       \gamma_3^0*       \gamma_4^{-1*}                \vdots           |
  // |                                                                                                                     |
  // | \gamma_4^4*    \gamma_4^3*       \gamma_4^2*       \gamma_4^1*       \gamma_4^0*                                    |
  // |                                                                                                                     |
  // | \vdots                                                                                 \ddots                       |
  // |                                                                                                                     |
  // | \gamma_n^n*                                                                                        \gamma_n^0*      |
  // |                                                                                                                     |
  // --                                                                                                                  --
  //
  // **********************************************************************************************************
  // Calculating the Matrix ***********************************************************************************
  // **********************************************************************************************************
  // **********************************************************************************************************
  //
  // The (0,0) position upper left corner \gamma_0^{0*} = 1
  // The first column and first row can be calculated using the recursive formulas (1) and (2)
  // Each diagonal is calculated using formulas (3) and (4)
  // Since we are using p = 16 terms (up to and including 3rd order derivatives), the matrix above will a 4x4 matrix
  // The count starts on 0, so we can set the matrix size using the highest order derivatives we are using to
  // approximate the series.
  // We will want to incorporate the order of approximation that we will be using as a class variable
  // May be a faster way using columns with C++


  // Set gamma_0^{0*} to 1
  //  (0)   \gamma_0^{0*}  (r) = 1
  ans[0][0] = 1.0;
  //std::cout << "Set (0,0) position of getSCoeff matrix" << std::endl;

  // Calculate the first row and column
  //  (1)   \gamma_n^{n*}  (r) = (-1) \frac{\eta}{n} \gamma_{n-1}^{(n-1)*}   for 1 \leq n
  //  (2)   \gamma_n^{-n*} (r) = (-1) \frac{\xi}{n}  \gamma_{n-1}^{-(n-1)*}  for 1 \leq n
  for (unsigned int i=1; i<=order_of_approximation; ++i)
  {
	ans[i][0] = ans[i-1][0] * (-1.0) * difference_eta / double(i);
    ans[0][i] = ans[0][i-1] * (-1.0) * difference_xi / double(i);
  }
  //std::cout << "Set 1st row and column of getSCoeff matrix" << std::endl;

  // Calculate each diagonal ********************************************* //
  //  (3)   \gamma_n^{m*}  (r) = \frac{(2n-1)z}{n^2 - m^2} \gamma_{n-1}^{m*} -  \frac{|r|^2}{n^2 - m^2} \gamma_{n-2}^{m*}   for m \geq 0
  //  (4)   \gamma_n^{-m*} (r) = \frac{(2n-1)z}{n^2 - m^2} \gamma_{n-1}^{-m*} - \frac{|r|^2}{n^2 - m^2} \gamma_{n-2}^{-m*}  for m \geq 0

  // Main Diagonal - here m = 0 so power (n^2 - m^2) = n^2 in formula
  ans[1][1] = (2*1-1) * difference_z *ans[0][0] / pow(double(1),2.0);
  //std::cout << "Set (1,1) position of getSCoeff matrix" << std::endl;
  for (unsigned int i=2; i<=order_of_approximation; ++i)
    ans[i][i] = ( (2.*double(i)-1.) * difference_z * ans[i-1][i-1] - norm_r_squared * ans[i-2][i-2] ) / std::pow(double(i),2.0);
  //std::cout << "Set main diagonal of getSCoeff matrix" << std::endl;

  // Stepping Through Nested Loop Below for Matrix Coefficients Below *************
  //
  // Note: i and j start their count on 1 since the first row and column and main diagonal are complete
  //       For an (n,n) matrix of coefficients there are n upper and lower diagonals
  //         Above, we have set the last upper and lower diagonal that is the lower left hand and upper
  //         right hand corner of the matrix and falls into the first row and column.
  //         Therefore, the bound order_of_approxiation on i ensures that i+j runs from 1st upper and lower
  //         off diagonal to the second to last upper and lower off diagonal
  //       The bound on j ensures that the inner nested loop does not go beyond the matrix size
  // Example: matrix is 4x4 according to order of approximation n = 3
  //          The coefficients on the main diagonal and 1st row and column have been taken care of so far
  //                           --                             --    --                                                          --
  //                           |   *       *       *       *    |   |   *       *               *                    *            |
  //                           |                                |   |                                                             |
  //                           |   *       *     [1][2]  [1][3] |   |   *       *             \gamma_2^{-1*}       \gamma_3^{-2*} |
  //  Matrix of Coefficients = |                                | = |                                                             |
  //                           |   *     [2][1]    *     [2][3] |   |   *     \gamma_2^{1*}     *                  \gamma_3^{-1*} |
  //                           |                                |   |                                                             |
  //                           |   *     [3][1]  [3][2]    *    |   |   *     \gamma_3^{2*}   \gamma_3^{1*}          *            |
  //                           |                                |   |                                                             |
  //                           --                             --    --                                                          --
  // Lower Triangular Part - \gamma_n^{m*} ***************************************************
  // i = 1:
  //   j = 1:  [i+j][j] = [2][1] = gamma_2^{1*} = gamma_{i+j}^{i*}         - 1st lower diagonal
  //   j = 2:  [i+j][j] =       [3][2] = gamma_3^{1*} = gamma_{i+j}^{i*}   - 1st lower diagonal
  // Stop: i+j = 1+2 = 3 = order_of_approximation
  // i = 2:
  //   j = 1:  [i+j][j] = [3][1] = gamma_3^{2*} = gamma_{i+j}^{i*}         - 2nd lower diagonal
  // Stop: i+j = 2+1 = 3 = order_of_approximation
  //
  // NOTE: from matrix and work above - i coincides with m
  //                                    and i+j coincides with n
  //                                    used in recursive formulas below
  //
  // Upper Triangular Part - \gamma_n^{m*} ***************************************************
  // i = 1:
  //   j = 1:  [j][i+j] = [1][2]= gamma_2^{-1*} = gamma_{i+j}^{-i*}          - 1st upper diagonal
  //   j = 2:  [j][i+j] =       [2][3]= gamma_{3}^{-1*} = gamma_{i+j}^{-i*}  - 1st upper diagonal
  // Stop: i+j = 1+2 = 3 = order_of_approximation
  // i = 2:
  //   j = 1:  [j][i+j] = [1][3]= gamma_3^{-2*} = gamma_{i+j}^{-i*}          - 2nd upper diagonal
  // Stop: i+j = 2+1 = 3 = order_of_approximation
  //
  // NOTE: from matrix and work above - i+j coincides with n
  //                                    and -i coincides with m
  //                                    used in recursive formulas below
  for (unsigned int i=1; i<=order_of_approximation; ++i)
  {
	// * Working down column i+j and across row i+j
	// * Since there are only (order_of_approximation + 1) rows and columns
	//   and the count for the index into the matrices starts on zero,
	//   we need i+j < order_of_approximation + 1 to not go beyond the matrix index
	// * j must also be less than the size of the matrix since j indexes a row or column
	//   j < order_of_approximation + 1
	// * Further, once we set i, we know that we are working on the upper and lower diagonals i
	// * Therefore we are starting from position (i,0) or (0,i), incrementing j by 1 unit,
	//   and moving down that diagonal to (i+j,j) or (j,i+j)
	// * If we were working on the main diagonal, then we would know that there are n elements
	//   if the matrix size is n.  Each shift to a lower or upper diagonal has one unit less than
	//   the previous diagonal.  That is, if i = 1 then the number of elements in the diagonal is
	//   n_d_e = n - 1 = n - i, if i = 2 then the number of diagonal elements is n_d_e = n - 2 = n - i,
	//   and so on.
	// * j can therefore not have a fixed upper bound, and we must offset j (the number of diagonal elements)
	//   by i as shown below
	//  (3)   \gamma_n^{m*}  (r) = \frac{(2n-1)z}{n^2 - m^2} \gamma_{n-1}^{m*} -  \frac{|r|^2}{n^2 - m^2} \gamma_{n-2}^{m*}   for m \geq 0
	//  (4)   \gamma_n^{-m*} (r) = \frac{(2n-1)z}{n^2 - m^2} \gamma_{n-1}^{-m*} - \frac{|r|^2}{n^2 - m^2} \gamma_{n-2}^{-m*}  for m \geq 0
	for (unsigned int j=1; j<=order_of_approximation-i; ++j) //(order_of_approximation + 1)-i; ++j)
    {
      if (j == 1) // 2nd column or row - In 2nd column or row, second term in recursive formula is zero
    	          //                     and would be outside matrix index
      {
        //std::cout << "Nested Loop if statement lower triangular part iteration i = " << i << ", j = " << j << std::endl;
        // Lower Triangular Part: i = m and i+j = n
        ans[i+j][j] = ( (2.*double(i+j)-1.) * difference_z * ans[i+j-1][j-1] )
        		               / ( std::pow(double(i+j),2.0) - std::pow(double(i),2.0) );
        //std::cout << "Nested Loop if statement upper triangular part iteration i = " << i << ", j = " << j << std::endl;
        // Upper Triangular Part: -i = m and i+j = n
        ans[j][i+j] = ( (2.*double(i+j)-1.) * difference_z * ans[j-1][i+j-1] )
        		               / ( std::pow(double(i+j),2.0) - std::pow(double(i),2.0) );
      }
      else
      {
        //std::cout << "Nested Loop else statement lower triangular part iteration i = " << i << ", j = " << j << std::endl;
        // Lower Triangular Part: i = m and i+j = n
        ans[i+j][j] = ( ((2.*double(i+j)-1.) * difference_z * ans[i+j-1][j-1]) - (norm_r_squared * ans[i+j-2][j-2]) )
                                              / ( std::pow(double(i+j),2.0) - std::pow(double(i),2.0) );
        //std::cout << "Nested Loop else statement upper triangular part iteration i = " << i << ", j = " << j << std::endl;
        // Upper Triangular Part: -i = m and i+j = n
        ans[j][i+j] = ( ((2.*double(i+j)-1.) * difference_z * ans[j-1][i+j-1]) - (norm_r_squared * ans[j-2][i+j-2]) )
                                              / ( std::pow(double(i+j),2.0) - std::pow(double(i),2.0) );

      }
    }
    //std::cout << "Outer Loop iteration i = " << i << std::endl;
  }
  //std::cout << "Finish Nested Loop" << std::endl;

  return ans;

}



// | \gamma_0^0*    \gamma_1^{-1*}    \gamma_2^{-2*}    \gamma_3^{-3*}    \gamma_4^{-4*}   \gamma_5^{-5*}    \hdots      \gamma_n^{-n*}   |
// |                                                                                                                                      |
// | \gamma_1^1*    \gamma_1^0*       \gamma_2^{-1*}    \gamma_3^{-2*}    \gamma_4^{-3*}   \gamma_5^{-4*}                \gamma_n^{-n+1*} |
// |                                                                                                                                      |
// | \gamma_2^2*    \gamma_2^1*       \gamma_2^0*       \gamma_3^{-1*}    \gamma_4^{-2*}   \gamma_5^{-3*}                \gamma_n^{-n+2*} |
// |                                                                                                                                      |
// | \gamma_3^3*    \gamma_3^2*       \gamma_3^1*       \gamma_3^0*       \gamma_4^{-1*}   \gamma_5^{-2*}                \vdots           |
// |                                                                                                                                      |
// | \gamma_4^4*    \gamma_4^3*       \gamma_4^2*       \gamma_4^1*       \gamma_4^0*      \gamma_5^{-1*}                                 |
// |                                                                                                                                             |
// | \gamma_5^5*    \gamma_5^4*       \gamma_5^3*       \gamma_5^2*       \gamma_5^1*      \gamma_5^0*                                           |
// |                                                                                                                                             |
// | \vdots                             \vdots                                                           \ddots                                  |
// |                                                                                                                                             |
// | \gamma_n^n*    \gamma_n^{(n-1)*}  \hdots                                                            \gamma_n^{(n-(n-1))*}  \gamma_n^0*      |
// |                                                                                                                     |

// gamma_hat_n^{m*}(r) = -0.5 * (Delta_1^{-1} gamma_n^{m*}(r) - Delta_1^1 gamma_n^{m*}(r)
/*     = -0.5 * ( Delta_1^{-1} gamma_n^{m*}(r) - Delta_1^1 gamma_n^{m*}(r) )                             */
/*     = -0.5 * ( (-1)^{1+(-1)} (-1)^m gamma_{n-1}^{-m+(-1)} - (-1)^{1+1} (-1)^m gamma_{n-1}^{-m+1} )    */
/*     = -0.5 * ( (-1)*(-1)^{m+1} gamma_{n-1}^{-(m+1)} - (-1)(-1)^{m-1} gamma_{n-1}^{-(m-1)}             */
/*     = -0.5 * ( -1 gamma_{n-1}^{(m+1)*} + gamma_{n-1}^{(m-1)*} )                                       */
std::vector<std::vector<std::vector<std::complex<double> > > > Potential::getRVectorGrad(std::vector<std::vector<std::complex<double> > > & RVec)
{
  //const std::complex<double> imaginary_unit(0.0,1.0);
	  //std::vector<std::vector<std::complex<double> > > ans(order_of_approximation+1, std::vector<std::complex<double> >(order_of_approximation+1));
	  std::vector<std::vector<std::vector<std::complex<double> > > >
	      ans(3, std::vector<std::vector<std::complex<double> > > (order_of_approximation+1, std::vector<std::complex<double> >(order_of_approximation+1)));

	  ans[0][0][0] = 0.0;
	  ans[1][0][0] = 0.0;
	  ans[2][0][0] = 0.0;
	  //std::cout << "Set (0,0) position of getSCoeff matrix" << std::endl;

	  // Calculate the first row and column
	  //  (1)   \gamma_n^{n*}  (r) = (-1) \frac{\eta}{n} \gamma_{n-1}^{(n-1)*}   for 1 \leq n
	  //  (2)   \gamma_n^{-n*} (r) = (-1) \frac{\xi}{n}  \gamma_{n-1}^{-(n-1)*}  for 1 \leq n
	  for (unsigned int i=1; i<=order_of_approximation; ++i)
	  {
		ans[0][i][0] =  RVec[i-1][0] * (0.5);
	    ans[0][0][i] = -RVec[0][i-1] * (0.5);

	    ans[1][i][0] =  RVec[i-1][0] * (0.5) * std::complex<double>{0.0,-1.0};
	    ans[1][0][i] =  RVec[0][i-1] * (0.5) * std::complex<double>{0.0,-1.0};

	    ans[2][i][0] =  0.0;
	    ans[2][0][i] =  0.0;
	  }

	  // main diagonal
	  ans[0][1][1] = 0.0;
	  ans[1][1][1] = 0.0;
	  ans[2][1][1] = -RVec[0][0] * std::complex<double>{-1.0, 0.0};
	  for (unsigned int i=2; i<=order_of_approximation; ++i)
	  {
	    ans[0][i][i] = (0.5) * (-RVec[i-1][i-2] + RVec[i-2][i-1]);
	    ans[1][i][i] = (0.5) * ( RVec[i-1][i-2] + RVec[i-2][i-1]) * std::complex<double>{0.0,-1.0};
	    ans[2][i][i] = -RVec[i-1][i-1] * std::complex<double>{-1.0, 0.0};
	  }

	  // working along rows (for speed - vector of vectors)
	  for (unsigned int i=1; i<=order_of_approximation; ++i)
       for (unsigned int j=1; j<=order_of_approximation; ++j)
	    {
         if (i > j)    // working on lower triangular part  e.g., (i,j) = (4,1)
         {
           if (j == 1)  // 2nd column where gamma_{n-1}^{(m+1)*} = 0 (e.g. (i,j) = (3,1) has gamma_hat_3^{2*}
           {
             ans[0][i][j] = (0.5) * (RVec[i-1][j]);
             ans[1][i][j] = (0.5) * (RVec[i-1][j]) * std::complex<double>{0.0,-1.0};
             ans[2][i][j] = -RVec[i-1][j-1] * std::complex<double>{-1.0, 0.0};
           }
           else // (j > 1)
           {
             ans[0][i][j] = (0.5) * (-RVec[i-1][j-2] + RVec[i-1][j]);
             ans[1][i][j] = (0.5) * ( RVec[i-1][j-2] + RVec[i-1][j]) * std::complex<double>{0.0,-1.0};
             ans[2][i][j] = -RVec[i-1][j-1] * std::complex<double>{-1.0, 0.0};
           }
         }
         else if (i == j)
         {} // do nothing - main diagonal already done
         else // (i < j) - working on upper triangular part e.g., (i,j) = (1,3)
		 {
           if (i == 1)  // 2nd row where gamma_{n-1}^{(m-1)*} = 0 (e.g. (i,j) = (1,3) has gamma_hat_3^{-2*}
           {
             ans[0][i][j] = (0.5) * (-RVec[i][j-1]);
             ans[1][i][j] = (0.5) * ( RVec[i][j-1]) * std::complex<double>{0.0,-1.0};
             ans[2][i][j] = -RVec[i-1][j-1] * std::complex<double>{-1.0, 0.0};
           }
           else // (i > 1)
           {
             ans[0][i][j] = (0.5) * (-RVec[i][j-1] + RVec[i-2][j-1]);
             ans[1][i][j] = (0.5) * ( RVec[i][j-1] + RVec[i-2][j-1]) * std::complex<double>{0.0,-1.0};
             ans[2][i][j] = -RVec[i-1][j-1] * std::complex<double>{-1.0, 0.0};
           }
		 }
	    }

	  return ans;

}




/**
 *  Explanation of getSVector
 */
// getSVector member function returns the powers of the S-expansion power series.
// The function is programmed similarly to the getSCoeff function.  See getSCoeff function for more details.
// We use p = 16 to obtain all derivatives up to and including order 3.
// The member function Potential::getSCoeff returns the coefficients \theta_n^{m} (x_b - z_A) of the
// far-field expansion \frac{1}{|r-x|} = \sum_{n=0}^{\infty} \gamma_n^{m*} (x) \theta_{n}^{m} (r)
// where
//                     \frac{1}{|x_a - x_b|} = \frac{1}{|(x_a - z_A) - (x_b - z_A)|}
//                                           = \sum_{n=0}^{\infty} \sum_{m=-n}^{n} \gamma_n^{m*} (x_a - z_A) \theta_n^{m*} (x_b - z_A)
//
// Recall the formulas for \theta_n^{m} (\mathbf{r})
//
//  (0)   \theta_0^{0*} (r) = \frac{ 1 }{ \|\mathbf{r}\| }
//
//  (1)   \theta_n^{n}  (r) = (2n-1) \frac{2\xi}{\|\mathbf{r}\|^2} \theta_{n-1}^{(n-1)}   for 1 \leq n
//
//  (2)   \theta_n^{-n} (r) = (2n-1) \frac{2\eta}{\|\mathbf{r}\|^2} \theta_{n-1}^{-(n-1)}  for 1 \leq n
//
//  (3)   \theta_n^{m}  (r) = ( (2n-1) z \theta_{n-1}^{m} -  ( (n - 1)^2 - m^2 ) \theta_{n-2}^{m} ) / \|\mathbf{r}\|^2  for 0 < m \leq n
//
//  (4)   \theta_n^{-m} (r) = ( (2n-1) z \theta_{n-1}^{-m} -  ( (n - 1)^2 - m^2 ) \theta_{n-2}^{-m} / \|\mathbf{r}\|^2  for -n \leq m < 0
//
// The coefficients \theta_n^{m*} (x_b - z_A) may be complex as can be seen from the formulas above
// Therefore, the vector of coefficients that is returned - ans - is a complex vector
// Further, the coefficients correspond to a full matrix shown below
//
// --                                                                                                                  --
// | \theta_0^0     \theta_1^{-1}     \theta_2^{-2}     \theta_3^{-3}     \theta_4^{-4}     \hdots      \theta_n^{-n}    |
// |                                                                                                                     |
// | \theta_1^1     \theta_1^0        \theta_2^{-1}     \theta_3^{-2}     \theta_4^{-3}                 \theta_n^{-n+1}  |
// |                                                                                                                     |
// | \theta_2^2     \theta_2^1        \theta_2^0        \theta_3^{-1}     \theta_4^{-2}                 \theta_n^{-n+2}  |
// |                                                                                                                     |
// | \theta_3^3     \theta_3^2        \theta_3^1        \theta_3^0        \theta_4^{-1}                 \vdots           |
// |                                                                                                                     |
// | \theta_4^4     \theta_4^3        \theta_4^2        \theta_4^1        \theta_4^0                                     |
// |                                                                                                                     |
// | \vdots                                                                            \ddots                            |
// |                                                                                                                     |
// | \theta_n^n                                                                                         \theta_n^0       |
// |                                                                                                                     |
// --                                                                                                                  --
//
std::vector<std::vector<std::complex<double> > > Potential::getSVector(std::vector<double> y, std::vector<double> xstar)
{
  // y = x_b is the target point
  // xstar = z_A is the center of the source cell
////  p = 16;
////  order_of_approximation = 3;
  //      norm_r_squared - ||x_b - z_A||^2
  double norm_r_squared = std::pow(y[0]-xstar[0],2.0) + std::pow(y[1]-xstar[1],2.0) + std::pow(y[2]-xstar[2],2.0);
  //
  // Recall the coordinate \eta = n =  -0.5 * (x - iy)
  //        the coordinate \xi  = e =   0.5 * (x + iy)
  //        the target particle is x_b
  //        the target cell center is z_B
  std::complex<double> eta_x_b, xi_x_b;
  double x_b_z_component;
  eta_x_b.real(-0.5*y[0]);  eta_x_b.imag(0.5*y[1]);
  xi_x_b.real(0.5*y[0]);    xi_x_b.imag(0.5*y[1]);
  x_b_z_component = y[2];

  std::complex<double> eta_z_A, xi_z_A;
  double z_A_z_component;
  eta_z_A.real(-0.5*xstar[0]);  eta_z_A.imag(0.5*xstar[1]);
  xi_z_A.real(0.5*xstar[0]);    xi_z_A.imag(0.5*xstar[1]);
  z_A_z_component = xstar[2];

  std::complex<double> difference_eta = eta_x_b - eta_z_A;
  std::complex<double> difference_xi  = xi_x_b - xi_z_A;
  double difference_z = x_b_z_component - z_A_z_component;

  // Calculating the Matrix
  //
  // The (0,0) position upper left corner \theta_0^{0*} = 1 / |r|
  // The first column and first row can be calculated using the recursive formulas (1) and (2)
  // Each diagonal is calculated using formulas (3) and (4)
  // Since we are using p = 16 terms (up to and including 3rd order derivatives), the matrix above will a 4x4 matrix
  // The count starts on 0, so we can set the matrix size using the highest order derivatives we are using to
  // approximate the series.
  // We will want to incorporate the order of approximation that we will be using as a class variable
  // May be a faster way using columns with C++

  std::vector<std::vector<std::complex<double> > > ans(order_of_approximation+1, std::vector<std::complex<double> >(order_of_approximation+1));

  // Set theta_0^{0*} to 1 / |r|
  //  (0)   \theta_0^{0*} (r) = \frac{ 1 }{ \|\mathbf{r}\| }
  ans[0][0] = 1.0 / std::pow(norm_r_squared, 0.5);
  //std::cout << "Set (0,0) position of getSVector matrix" << std::endl;

  // Calculate the first row and column
  //  (1)   \theta_n^{n}  (r) = (2n-1) \frac{2\xi}{\|\mathbf{r}\|^2} \theta_{n-1}^{(n-1)}   for 1 \leq n
  //  (2)   \theta_n^{-n} (r) = (2n-1) \frac{2\eta}{\|\mathbf{r}\|^2} \theta_{n-1}^{-(n-1)}  for 1 \leq n
  for (unsigned int i=1; i<=order_of_approximation; ++i)
  {
    ans[i][0] = ans[i-1][0] * (2.0*double(i) - 1.0) * 2.0 * difference_xi / norm_r_squared;
    ans[0][i] = ans[0][i-1] * (2.0*double(i) - 1.0) * 2.0 * difference_eta / norm_r_squared;
  }
  //std::cout << "Set 1st row and column of getSVector matrix" << std::endl;

  // Calculate each diagonal ********************************************* //
  //  (3)   \theta_n^{m}  (r) = ( (2n-1) z \theta_{n-1}^{m} -  ( (n - 1)^2 - m^2 ) \theta_{n-2}^{m} ) / \|\mathbf{r}\|^2  for 0 < m \leq n
  //  (4)   \theta_n^{-m} (r) = ( (2n-1) z \theta_{n-1}^{-m} -  ( (n - 1)^2 - m^2 ) \theta_{n-2}^{-m} / \|\mathbf{r}\|^2  for -n \leq m < 0

  // Main Diagonal - m = 0
  ans[1][1] = (2*1-1) * difference_z *ans[0][0] / norm_r_squared;
  //std::cout << "Set (1,1) position of getSVector matrix" << std::endl;
  for (unsigned int i=2; i<=order_of_approximation; ++i)
    ans[i][i] = ( (2.*double(i)-1.) * difference_z * ans[i-1][i-1] - (std::pow(double(i)-1,2.0)) * ans[i-2][i-2] ) / norm_r_squared;
  //std::cout << "Set main diagonal of getSVector matrix" << std::endl;

  // lower triangular part
  // i = 1: [i+j][j] = [2][1], [3][2], [4][3], [5][4] - 1st lower diagonal
  // i = 2: [i+j][j] = [3][1], [4][2], [5][3], [6][4] - 2nd lower diagonal
  // upper triangular part
  // i = 1: [j][i+j] = [1][2], [2][3], [3][4], [4][5]
  // i = 2: [j][i+j] = [1][3], [2][4], [3][5], [4][6]
  // Note: from matrix above - i coincides with m
  //                         - i+j coincides with n
  //       used in recursive formulas below
  //  (3)   \theta_n^{m}  (r) = ( (2n-1) z \theta_{n-1}^{m} -  ( (n - 1)^2 - m^2 ) \theta_{n-2}^{m} ) / \|\mathbf{r}\|^2  for 0 < m \leq n
  //  (4)   \theta_n^{-m} (r) = ( (2n-1) z \theta_{n-1}^{-m} -  ( (n - 1)^2 - m^2 ) \theta_{n-2}^{-m} / \|\mathbf{r}\|^2  for -n \leq m < 0
  for (unsigned int i=1; i<=order_of_approximation; ++i)  // main diagonal - i = 0
    for (unsigned int j=1; j<=order_of_approximation-i; ++j)
    {
      if (j == 1) // 2nd column or row - second term in recursive formula is zero and outside matrix so not including
      {
        //std::cout << "Nested Loop if statement lower triangular part iteration i = " << i << ", j = " << j << std::endl;
        // lower triangular part
        ans[i+j][j] = ( (2.*double(i+j)-1.) * difference_z * ans[i+j-1][j-1] )
        		               / norm_r_squared ;
        //std::cout << "Nested Loop if statement upper triangular part iteration i = " << i << ", j = " << j << std::endl;
        // upper triangular part
        ans[j][i+j] = ( (2.*double(i+j)-1.) * difference_z * ans[j-1][i+j-1] )
        		               / norm_r_squared ;
      }
      else // j > 1 - include both terms of recursive formula
      {
        //std::cout << "Nested Loop else statement lower triangular part iteration i = " << i << ", j = " << j << std::endl;
        // lower triangular part
        ans[i+j][j] = ( (2.*double(i+j)-1.) * difference_z * ans[i+j-1][j-1]
      		              - ( std::pow(double(i+j-1),2.0) - std::pow(double(i),2.0) ) * ans[i+j-2][j-2]
        		      )
        		         / norm_r_squared;
        //std::cout << "Nested Loop else statement upper triangular part iteration i = " << i << ", j = " << j << std::endl;
        // upper triangular part
        ans[j][i+j] = ( (2.*double(i+j)-1.) * difference_z * ans[j-1][i+j-1]
        		          - ( std::pow(double(i+j-1),2.0) - std::pow(double(i),2.0) ) * ans[j-2][i+j-2]
        	          )
                         / norm_r_squared;

      }
    }

//  std::cout << "norm_r_squared from SVec is |r|^2 = " << norm_r_squared << std::endl;
//  std::cout << "ans[1][2] = " << ans[1][2] << std::endl;
//  std::cout << "ans[0][1] = " << ans[0][1] << std::endl;
  return ans;

}


// direct (exact) calculation of the potential acting on yj by xi
double Potential::direct(std::vector<double> yj, std::vector<double> xi)
{
  double norm_yj_minus_xi = pow(yj[0]-xi[0],2.0) + pow(yj[1]-xi[1],2.0) + pow(yj[2]-xi[2],2.0);
  norm_yj_minus_xi = pow(norm_yj_minus_xi,0.5);

  double ans = 1.0 / pow(norm_yj_minus_xi,1.0);
  return ans;
}


// direct (exact) calculation of the potential acting on yj by xi
//std::vector<double>   Potential::directGrad(std::vector<double> yj, std::vector<double> xi)
std::vector<double>   Potential::directGrad(std::vector<double> yj, std::vector<double> xi)
{
  double norm_yj_minus_xi = pow(yj[0]-xi[0],2.0) + pow(yj[1]-xi[1],2.0) + pow(yj[2]-xi[2],2.0);
  norm_yj_minus_xi = pow(norm_yj_minus_xi,0.5);

  std::vector<double> ans(3,0.0);
//  double ans;
//  ans[0] = xi[0]-yj[0] / pow(norm_yj_minus_xi,3.0);
//  ans[1] = xi[1]-yj[1] / pow(norm_yj_minus_xi,3.0);
//  ans[2] = xi[2]-yj[2] / pow(norm_yj_minus_xi,3.0);
  ans[0] = (yj[0]-xi[0]) / pow(norm_yj_minus_xi,3.0);
  ans[1] = (yj[1]-xi[1]) / pow(norm_yj_minus_xi,3.0);
  ans[2] = (yj[2]-xi[2]) / pow(norm_yj_minus_xi,3.0);
  return ans;
}






