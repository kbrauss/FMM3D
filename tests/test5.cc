/*
 * test5.cc
 *
 *  Created on: Apr 21, 2017
 *      Author: dbpc
 */


int main()
{
	  // note that these refinement levels counts start on l = 1
	  int DEFAULT_NUM_LEVEL = 3;                   // default refinement level
	  // p represents the number of derivatives to approximate Taylor series to order |alpha|
	  // The derivatives at order |alpha| = n are the combinations (n+2)_C_2 = ((n +2) \\ 2) = ((n+2)*(n+1))/(2*1)
	  // We set p = 20 to approximate series to the third derivatives (1 + 3 + 6 + 10 = 20)
	  //                                                               0th 1st 2nd 3rd  derivatives
	  // The count is different with Dehnen's formulas.  Each level of approximation in terms of order of derivative
	  // relates to the formula
	  //            \frac{1}{n!} (x \cdot \nabla)^n = \sum_{m=-n}^{n} \gamma_n^{m*} (-x) \Delta_n^m
	  // Therefore, for the zero order derivatives (1) we have one term of the summation m = 0
	  //            for the 1st order derivatives  (3) we have three terms for the summation m = -1, 0, 1
	  //            for the 2nd order derivatives  (5) we have 5 terms for the summation     m = -2, -1, 0, 1, 2
	  //            for the 3rd order deirvatives  (7) we have 7 terms for the summation     m = -3, -2, -1, 0, 1, 2, 3
	  // In Dehnen's coordinate system, if we we approximate the far-field or near-field series
	  // to all the third derivatives where all 0, 1st, 2nd, and 3rd order derivatives add up to  1 + 3 + 5 + 7 = 16)
	  //                                                                                          0th 1st 2nd 3rd derivatives
	  // The harmonic property appears to reduce the number of terms in comparison to the Euclidean FMM algorithm
	  // We determine p by setting the level of approximation (highers order of derivatives use in Taylor series)
	  // abs_alpha - below set to 3 to use third order derivatives
	  unsigned int abs_alpha = 3;
	  // p is the number of terms used in Taylor Series approximation
	  // Setting p to 1 and adding number of derivatives at each order below
	  // using (n+2)_C_2 combinations formula.
	  int p = 0;
	  int combination = 1; // number of derivatives for order |alpha| = 0
	  p = p + combination;
	  for (unsigned int i = 1; i <= abs_alpha; ++i)
	  {
	    combination = 2*i+1;  // number derivatives for order |alpha| = i
	    p = p + combination;
	  }

	  if (p == 16)
	    std::cout << "p = " << p << " correct" << std::endl;
	  else // p \neq 20
	    std::cout << "p != " << 20 << "incorrect" << std::endl;




	  /*******************************************************************************
	   * Target and Source Points
	   *******************************************************************************/

	  // for this the source particle is x[0] and the target y[0]
	  // the interaction list cell that will be used contains point x[1] used below
	  // to located the cell.
	  std::vector<Point>  x;
	  std::vector<Point>  y;
	  std::vector<double> u;
	  x.resize(4);
	  y.resize(1);
	  u.resize(4);
	  for (unsigned int i=0; i<u.size(); ++i)
	    u[i] = 1.0;

	  // target point in cell of interest - located in cell n = 5 at l = 3
	  y[0].setX(0.15);
	  y[0].setY(0.03);
	  y[0].setZ(0.20);

	  // point in interaction list cell - located in cell n = 45 at l = 3
	  x[1].setX(0.38);
	  x[1].setY(0.06);
	  x[1].setZ(0.40);

	  // source - located in cell n = 511 at l = 3
	  x[0].setX(0.90);
	  x[0].setY(0.925);
	  x[0].setZ(0.95);


	  /*******************************************************************************
	   * Direct Calculation of the Potential Function f(x) = 1/||x-y||
	   *******************************************************************************/

	  std::cout << "**************************************************" << std::endl;
	  std::cout << "Direct Calculation" << std::endl;
	  std::cout << "**************************************************" << std::endl;
	  std::cout << std::endl;

	  Potential potential(p,abs_alpha);

	  // performing direct calculation for potential of source x acting on target y below and printing value.
	  double direct_calculation = potential.direct(y[0].getCoord(), x[0].getCoord());
	  std::cout << "direct calculation = " << direct_calculation << std::endl;

	  std::cout << std::endl;
	  std::cout << std::endl;
	  std::cout << std::endl;
	  std::cout << "**************************************************" << std::endl;
	  std::cout << "**************************************************" << std::endl;

	  /*******************************************************************************
	   * far field expansion - s-expansion
	   *                                                                                                           _______________________
	   *                                                                                                          |           |           |
	   *                                                                                                          |           | * s = x_b |
	   *                                                                                                          |           |           |
	   *                                                                                                          |           |     *     |
	   *                                                                                                          |           |     z_b,c |
	   *                                                                                                          |-----------*-----------|
	   *                                                                                                          |           | z_b,p     |
	   *                                                                                                          |           |           |
	   *                                                                                                          |           |           |
	   *                                                                                                          |           |           |
	   *                                                                                                          |___________|___________|
	   *
	   *
	   *
	   *
	   *
	   *
	   *
	   *
	   *
	   *
	   *                                                   _______________________
	   *                                                  |           |           |
	   *                                                  |           | * o = x_b |
	   *                                                  |           |     *     |
	   *                                                  |           |    z_b    |
	   *                                                  |           |           |
	   *                                                  |-----------i-----------|
	   *                                                  |           |           |
	   *                                                  |           |           |
	   *                                                  |           |           |
	   *                                                  |           |           |
	   *                           ___________ ___________|___________|___________|
	   *                          |           |           |
	   *                          |           |           |
	   *                          |           |           |
	   *                          |           |           |
	   *                          |           |           |
	   *                          |-----------n-----------|
	   *                          |           |           |
	   *                          |           |           |
	   *                          |           |           |
	   *                          |           |           |
	   *   _______________________|___________|___________|
	   *  |   y = x_a |           |
	   *  |  *        |           |
	   *  |     *     |           |
	   *  |     z_a,c |           |
	   *  |           |           |
	   *  |-----------*-----------|
	   *  |           | z_a,p     |
	   *  |           |           |
	   *  |           |           |
	   *  |           |           |
	   *  |___________|___________|
	   *
	   *
	   *
	   * Recall that for the far-field expansion sum_{i=1}^{p} ( b_m (x_i, x_*) S_m(y,x_*) )
	   * the powers (and factorial) (1/alpha) (x_b - z_b)^{alpha} are the coefficients b_m
	   * and the derivatives D^{alpha} f(z_b - x_a) are the powers S_m.  The expansion is
	   * the Taylor series
	   *
	   *      f(x_b - x_a) = sum_{|alpha| >= 0} [ ( (1 / (alpha !)) (x_b - z_b)^{alpha} ) D^{alpha} f(z_b - x_a) ]
	   *
	   * We recall that the series converges when || x_b - z_b || < || z_b - x_a ||
	   * (or || x - x_* || < || x_* - y ||=|| y - x_* ||)
	   * All the particles of a cell then have the same powers and we can combine the series of all
	   * the particles into a single series.  The coefficients of this new series are a summation of
	   * the coefficients of all the old series.
	   *
	   * The code below performs a far-field approximation for the potential function of the source x_b = x[0]
	   * acting on target x_a = y[0] and prints the value.  The center z_a is center_yBox
	   *
	   *******************************************************************************/

	  std::cout << "**************************************************" << std::endl;
	  std::cout << "Far-Field Expansion" << std::endl;
	  std::cout << "**************************************************" << std::endl;
	  std::cout << std::endl;

	  std::complex<double> r_expansion = 0.0;

	  // refinement level count starts at l = 0 for call to getBoxIndex
	  // checking that points located in expected boxes
	  int xBoxIndex = x[0].getBoxIndex(DEFAULT_NUM_LEVEL);
	  int xBoxParentIndex = x[0].getBoxIndex(DEFAULT_NUM_LEVEL-1);
	  int yBoxIndex = y[0].getBoxIndex(DEFAULT_NUM_LEVEL);
	  int yBoxParentIndex = y[0].getBoxIndex(DEFAULT_NUM_LEVEL-1);

	  Box xBox(DEFAULT_NUM_LEVEL, xBoxIndex, p);
	  Box xBoxParent(DEFAULT_NUM_LEVEL-1, xBoxParentIndex, p);
	  Box yBox(DEFAULT_NUM_LEVEL, yBoxIndex, p);
	  Box yBoxParent(DEFAULT_NUM_LEVEL-1, yBoxParentIndex, p);

	  std::cout << "x[0] = [" << x[0].getX() << "," << x[0].getY() << "," << x[0].getZ() << "]"
			    << std::endl;
	  std::cout << "x[0] is located in Box " << xBox.getIndex()
			    << " at level l = 3."<< std::endl;
	  std::cout << "x[0] is located in Box " << xBox.getParentIndex()
			    << " at level l = 2."<< std::endl;
	  std::cout << "y[0] = [" << y[0].getX() << "," << y[0].getY() << "," << y[0].getZ() << "]"
			    << std::endl;
	  std::cout << "y[0] is located in Box " << yBox.getIndex()
			    << " at level l = 3."<< std::endl;
	  std::cout << "y[0] is located in Box " << yBox.getParentIndex()
			    << " at level l = 2."<< std::endl;


	  FmmTree fmmtree(DEFAULT_NUM_LEVEL, x, y, potential);

	  Point center_xBox = xBox.getCenter();
	  Point center_xBoxParent = xBoxParent.getCenter();
	  std::cout << "center_xBox(x,y,z) = " << "center_xBox(" << center_xBox.getX() << ","
			    << center_xBox.getY() << "," << center_xBox.getZ() << ")" << std::endl;
	  std::cout << "center_xBoxParent(x,y,z) = " << "center_xBoxParent(" << center_xBoxParent.getX() << ","
			    << center_xBoxParent.getY() << "," << center_xBoxParent.getZ() << ")" << std::endl;
	  Point center_yBox = yBox.getCenter();
	  Point center_yBoxParent = yBoxParent.getCenter();
	  std::cout << "center_yBox(x,y,z) = " << "center_yBox(" << center_yBox.getX() << ","
			    << center_yBox.getY() << "," << center_yBox.getZ() << ")" << std::endl;
	  std::cout << "center_yBoxParent(x,y,z) = " << "center_yBoxParent(" << center_yBoxParent.getX() << ","
			    << center_yBoxParent.getY() << "," << center_yBoxParent.getZ() << ")" << std::endl;

	  // Obtaining far-field coefficients (s-expansion coefficients)
	  // The far-field coefficients of the s-expansion are the powers
	  //
	  //    (x_a - z_a)^{alpha} / (\alpha !)  =  (xi - xstar)^{alpha} / (alpha !)
	  //
	  // where z_a is the center of the source cell and x_a is the source particle (point).
	  // Recall the series has the form
	  //
	  //    psi(x - y)  =  sum_{alpha >= 0}  frac{ (x - x_*)^{alpha} }{ alpha! }  D^{alpha} psi(x_* - y)
	  //
	  // and the powers of the series are the derivatives
	  //
	  //    D^{alpha} psi(x_* - y)  =  D^{\alpha} f(z_a - x_b)  =  D^{alpha} psi(xstar - yj).
	  //
	  // In Dehnen's paper the far-field expansion is
	  //
	  //     \frac{1}{|r-x|} = \sum_{n=0}^{\infty} \gamma_n^{m*} (x) \theta_{n}^{m} (r)
	  //
	  // The coefficients are \gamma_n^{m*} (x) seen in
	  //
	  //     \frac{1}{|x_a - x_b|} = \frac{1}{|(x_a - z_A) - (x_b - z_A)|}
	  //                           = \sum_{n=0}^{\infty} \sum_{m=-n}^{n} \gamma_n^{m*} (x_a - z_A) \theta_n^{m*} (x_b - z_A)

	  // In Dehnen's paper
	  std::cout << "We are here" << std::endl;
	  std::vector<std::vector<std::complex<double> > > y_SCoeff = potential.getSCoeff(x[0].getCoord(), center_xBox.getCoord());

	  std::vector<std::vector<std::complex<double> > > RCoeffTranslated = potential.getSR(center_xBox.getCoord(), center_yBox.getCoord(),
	 		                                                                              y_SCoeff);
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
	  for (unsigned int i=0; i<4; ++i)
	    for (unsigned int j=0; j<4; ++j)
	      std::cout << "RCoeffTranslated[" << i << "][" << j << "] = "
	                << RCoeffTranslated[i][j] << std::endl;
	  // obtaining far-field powers (s-expansion powers)
	  // far field powers are D^{\alpha} f (x_b - z_a) = D^{\alpha} f (x[0] - center_yBox)
	  // Note: due to the formula
	  //
	  //   \mu_a \frac{1}{|x_b - x_a|} = \mu_a \frac{1}{|x_b - z_B + z_B - x_a|}
	  //                              = \mu_a \frac{1}{|(z_B - x_a) - (z_B - x_b)|
	  //                              = \mu_a \sum_{n=0}^{\infty} \sum_{m=-n}^{n} \gamma_n^{m*} (z_B - x_b) \theta_n^m(z_B - x_a)
	  //                              = \mu_a \sum_{n=0}^{\infty} \sum_{m=-n}^{n} \gamma_n^{m*} (z_B - x_b)
	  //                                                                                      \theta_n^m(z_B - z_A - (x_a - z_A))
	  //                              = \mu_a \sum_{n=0}^{\infty} \sum_{m=-n}^{n} \gamma_n^{m*} (z_B - x_b)
	  //                                                                 \sum_{k=0}^{\infty} \sum_{l=-k}^{k} \gamma_{k}^{l*}(x_a - z_A)
	  //                                                                                                        theta_{n+k}^{m+l} (z_B - z_A)
	  //
	  // We have the order z_B - x_b to calculate the powers \gamma_n^{m*} (z_B - x_b) using the
	  // Potential class member function Potential::getRVector(y, xstar) with Potential::getRVector(x_b, z_B)
	  std::vector<std::vector<std::complex<double> > > y_RVec = potential.getRVector(y[0].getCoord(), center_yBox.getCoord());

	  for (unsigned int i=0; i<=abs_alpha; ++i)
	    for (unsigned int j=0; j<=abs_alpha; ++j)
	      std::cout << "y_RVec[" << i << "][" << j << "] = " << y_RVec[i][j] << std::endl;

	  for (unsigned int i=0; i<=abs_alpha; ++i)
	    for (unsigned int j=0; j<=abs_alpha; ++j)
	    {
	      r_expansion = r_expansion + RCoeffTranslated[j][i] * y_RVec[j][i];
	    }

	  std::cout << "r_expansion = " << r_expansion.real() << " + i" << r_expansion.imag() << std::endl;

	  std::cout << std::endl;
	  std::cout << std::endl;
	  std::cout << std::endl;
	  std::cout << "**************************************************" << std::endl;
	  std::cout << "**************************************************" << std::endl;

}


