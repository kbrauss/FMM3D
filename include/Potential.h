/*
 * Potential.h
 *
 *  Created on: Jul 13, 2016
 *      Author: dbpc
 */

#ifndef POTENTIAL_H_
#define POTENTIAL_H_

#include <vector>
#include <complex>

/* We note that the only class that is passed in p - the number of terms in the far-field
 * and near field Taylor series expansions is the class Potential
 * However, in the constructor for FmmTree class, the instance of the Potential class
 * is also passed.
 * This is an important point as far as p - the number of terms of the series is concerned -
 * since in this pass boxes that make up the fmmtree array of class FmmTree are passed
 * p, and p is used to set the size of the vectors D, C, and DTilde.
 * Therefore, there are only two classes that actually use the number of terms of the
 * Taylor series of the FMM code - classes Potential and Box.
 * The importance of this is in the modularity of the code.  The only parts of code that
 * one has to be concerned about when increasing accuracy of approximation by increasing p
 * are the classes Potential and Box.  However, looking at the two classes, the only class
 * where code changes are made to implement more terms in the series approximations is the
 * Potential class.
 *
 */

class Potential
{
  public:
	int p;
	unsigned int order_of_approximation;
	// p gives the number of derivatives - here p = 20 means up to and including third order derivatives
    //int DEFAULT_P = 16;
    //unsigned int DEFAULT_ORDER_OF_APPROXIMATION = 3;
    unsigned int DEFAULT_ORDER_OF_APPROXIMATION = 4;
    int DEFAULT_P = (DEFAULT_ORDER_OF_APPROXIMATION+1)*(DEFAULT_ORDER_OF_APPROXIMATION+1);

    Potential() { this->p = DEFAULT_P;  this->order_of_approximation = DEFAULT_ORDER_OF_APPROXIMATION; };
	Potential(int p, int order_of_approximation)
	 :
     p(p),
     order_of_approximation(order_of_approximation)
//     multiIndexOrders(164,std::vector<unsigned int>(3))
	{};
	int getP() { return p;};
	void setP(int p) { this->p = p; };
	int getOrderApprox() {return order_of_approximation;};
	void setOrderApprox(int order_of_approximation) { this->order_of_approximation = order_of_approximation; };
	std::vector<std::vector<std::complex<double> > > getSS(std::vector<double> from,
			                                               std::vector<double> to,
			                                               const std::vector<std::vector<std::complex<double> > > &sCoeff_old);
	std::vector<std::vector<std::complex<double> > > getRR(std::vector<double> from,
			                                               std::vector<double> to,
			                                               const std::vector<std::vector<std::complex<double> > > &rCoeff_old);
//	std::vector<std::vector<std::vector<std::complex<double> > > > getRRgrad(std::vector<double> from,
//			                                                                 std::vector<double> to,
//			                                                                 const std::vector<std::vector<std::vector<std::complex<double> > > > &rCoeff_old);
	std::vector<std::vector<std::complex<double> > > getSR(std::vector<double> from, std::vector<double> to,
			                                               const std::vector<std::vector<std::complex<double> > > &sCoeff_old);
//	std::vector<std::vector<std::vector<std::complex<double> > > > getSRgrad(std::vector<double> from, std::vector<double> to,
//			                                                                 const std::vector<std::vector<std::complex<double> > > &sCoeff_old);

	std::vector<std::vector<std::complex<double> > > getRCoeff(std::vector<double> xi, std::vector<double> xstar);
	std::vector<std::vector<std::complex<double> > > getSCoeff(std::vector<double> xi, std::vector<double> xstar);

	std::vector<std::vector<std::complex<double> > > getRVector(std::vector<double> y, std::vector<double> xstar);
	std::vector<std::vector<std::vector<std::complex<double> > > > getRVectorGrad(std::vector<std::vector<std::complex<double> > > &RVec);
    std::vector<std::vector<std::complex<double> > > getSVector(std::vector<double> y, std::vector<double> xstar);

    double                             direct(std::vector<double> yj, std::vector<double> xi);
    std::vector<double>              directGrad(std::vector<double> yj, std::vector<double> xi);
    //double                             directGrad(std::vector<double> yj, std::vector<double> xi);

};




#endif /* POTENTIAL_H_ */
