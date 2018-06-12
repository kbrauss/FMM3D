/*
 * test4.cc
 *
 *  Created on: Apr 14, 2017
 *      Author: dbpc
 */
/**
 *  test4.cc
 *  Created on: April 14, 2017
 *      Author: dbpc
 *
 *  The program tests the RR transformation from a near-field series approximation to
 *  a near-field series approximation.
 *
 *  Recall the formula
 *
 *  \gamma_n^{m*} (x+y) = \sum_{k=0}^n \sum_{l=-k}^k \gamma_k^{l*} (y) \gamma_{n-k}^{(m-l)*}(x)
 *
 *  Let x denote the child center
 *      y denote the parent center
 *
 *  The gamma_k^{l*} (y) can be calculated once and used for the entire upward pass for a given level
 *  since the distance y will be the same for uniform partitions.
 *
 *  There is a good chance the gamma_{n-k}^{*m-l)*} (x) have already been calculated in the setup
 *  of the first part of the upward pass where the coefficients \gamma_n^{m*} are calculated for all
 *  source particles at the lowest refinement level l = L.
 *
 *  For example, take p = 16.  From our notes on the getSCoeff member function of Potential.cc we know
 *  that we would have generated a matrix of coefficients for each source particle that is 4x4 and contains
 *  the coefficients \gamma_0^0(x), \gamma_1^{-1}(x), \gamma_1^0(x), \gamma_1^1(x), ... , \gamma_3^3(x).
 *  That is \gamma_i^j(x) with 0 \leq i \leq 3  and -3 \leq j \leq 3 with |j| \leq i.
 *
 *  Using
 *
 *  \gamma_n^{m*} (x+y) = \sum_{k=0}^n \sum_{l=-k}^k \gamma_k^{l*} (y) \gamma_{n-k}^{(m-l)*}(x)
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
	  // Therefore, for the zero order derivatives (1) we have one term of the summation     m = 0
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
	  else // p \neq 16
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
	   * near field expansion - r-expansion
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
	   * At each level the particles of a cell have the same powers and can be combined to form a
	   * single series.  The coefficients of this new series are a summation of the coefficients of all the old series.
	   *
	   * The code below performs a near-field approximation for the potential function of the source x_b = x[0]
	   * acting on target x_a = y[0] and prints the value.  The center z_b is center_yBox
	   *
	   *******************************************************************************/

	  std::cout << "**************************************************" << std::endl;
	  std::cout << "Near-Field Expansion" << std::endl;
	  std::cout << "**************************************************" << std::endl;
	  std::cout << std::endl;

	  std::complex<double> r_expansion = 0.0;
	  std::complex<double> r_expansion_translated = 0.0;

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

	  Point center_yBox = yBox.getCenter();
	  Point center_yBoxParent = yBoxParent.getCenter();
	  std::cout << "center_yBox(x,y,z) = " << "center_yBox(" << center_yBox.getX() << ","
			    << center_yBox.getY() << "," << center_yBox.getZ() << ")" << std::endl;
	  std::cout << "center_yBoxParent(x,y,z) = " << "center_yBoxParent(" << center_yBoxParent.getX() << ","
			    << center_yBoxParent.getY() << "," << center_yBoxParent.getZ() << ")" << std::endl;

	  // obtaining far-field coefficients (s-expansion coefficients)
	  // far field coefficients are (z_a - x_a)^{\alpha} = (center_yBox - y[0])^{\alpha}
	  // Note: y[0] is the target ( in Main.cc it is x )
	  std::cout << "We are here" << std::endl;
	  std::vector<std::vector<std::complex<double> > > x_RCoeff = potential.getRCoeff(x[0].getCoord(), center_yBoxParent.getCoord());

	  for (unsigned int i=0; i<=abs_alpha; ++i)
	    for (unsigned int j=0; j<=abs_alpha; ++j)
	        std::cout << "x_RCoeff[" << i << "][" << j << "] = "
	                  << x_RCoeff[i][j] << std::endl;

	  // obtaining far-field powers (s-expansion powers)
	  // far field powers are D^{\alpha} f (x_b - z_a) = D^{\alpha} f (x[0] - center_yBox)
	  std::vector<std::vector<std::complex<double> > > x_RVecParent = potential.getRVector(y[0].getCoord(), center_yBoxParent.getCoord());
	  std::vector<std::vector<std::complex<double> > > x_RVecChild = potential.getRVector(y[0].getCoord(), center_yBox.getCoord());

	  for (unsigned int i=0; i<4; ++i)
	    for (unsigned int j=0; j<4; ++j)
	      std::cout << "x_RVecChild[" << i << "][" << j << "] = "
	                << x_RVecChild[i][j] << std::endl;

	  std::vector<std::vector<std::complex<double> > > RCoeffTranslated = potential.getRR(center_yBox.getCoord(), center_yBoxParent.getCoord(),
	 		                                                                              x_RCoeff);

	  for (unsigned int i=0; i<=abs_alpha; ++i)
	    for (unsigned int j=0; j<=abs_alpha; ++j)
	    {
	      r_expansion = r_expansion + x_RCoeff[j][i] * x_RVecChild[j][i];
	    }

	  std::cout << "r_expansion = " << r_expansion.real() << " + i" << r_expansion.imag() << std::endl;

	  for (unsigned int i=0; i<=abs_alpha; ++i)
	    for (unsigned int j=0; j<=abs_alpha; ++j)
	        std::cout << "RCoeffTranslated[" << i << "][" << j << "] = "
	                  << RCoeffTranslated[i][j] << std::endl;

	  for (unsigned int i=0; i<=abs_alpha; ++i)
	    for (unsigned int j=0; j<=abs_alpha; ++j)
	        std::cout << "x_RVecParent[" << i << "][" << j << "] = "
	                  << x_RVecParent[i][j] << std::endl;

	  for (unsigned int i=0; i<=abs_alpha; ++i)
	    for (unsigned int j=0; j<=abs_alpha; ++j)
	    {
	      r_expansion_translated = r_expansion_translated + RCoeffTranslated[j][i] * x_RVecChild[j][i];
	    }

	  std::cout << "r_expansion_translated = " << r_expansion_translated.real() << " + i" << r_expansion_translated.imag() << std::endl;

	  std::cout << std::endl;
	  std::cout << std::endl;
	  std::cout << std::endl;
	  std::cout << "**************************************************" << std::endl;
	  std::cout << "**************************************************" << std::endl;

}

