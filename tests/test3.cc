/*
 * test3.cc
 *
 *  Created on: Apr 2, 2017
 *      Author: dbpc
 */

/**
 *  test3.cc
 *  Created on: April 2, 2017
 *      Author: dbpc
 *
 *  The program tests the SS transformation from a far-field series approximation to
 *  a far-field series approximation.
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
 *
 *  Let n = 0, m = 0
 *  Then for \gamma_0^0(x+y) = \sum_{k=0}^n sum_{l=-k}^k \gamma_{n-k}^{(m-l)*} (x) \gamma_k^{l*} (y)
 *                           = \gamma_{0}^{0*} (x) \gamma_{0}^{0*} (y)
 *
 *  Let n = 1, m = 0
 *  Then for \gamma_1^0(x+y) = \sum_{k=0}^n sum_{l=-k}^k [ \gamma_{n-k}^{(m-l)*} (x) \gamma_k^{l*} (y) ]
 *                           = (     k = 0,     l = 0   =>   n - k = 1 - 0 = 1 and m -l = 0 - 0 = 0 =>
 *                               n - k = 1, m - l = 0
 *                             )
 *                              \gamma_{1}^{0*} (x) \gamma_{0}^{0*} (y)
 *                             (     k = 1,     l = -1  =>   n - k = 1 - 1 = 0 and m - l = 0 - (-1) = 1 =>
 *                               n - k = 0, m - l = 1
 *                             )
 *                             + \gamma_{0}^{(1)*} (x) \gamma_{1}^{(-1)*} (y)
 *                             (     k = 1,     l = 0   =>   n - k = 1 - 1 = 0 and m - l = 0 - (0) = 0 =>
 *                               n - k = 0, m - l = 0
 *                             )
 *                             + \gamma_{0}^{(0)*} (x) \gamma_{1}^{(0)*} (y)
 *                             (     k = 1,     l = 1   =>   n - k = 1 - 1 = 0 and m - l = 0 - (1) = -1 =>
 *                               n - k = 0, m - l = -1
 *                             )
 *                             + \gamma_{0}^{(-1)*} (x) \gamma_{1}^{(1)*} (y)
 *                           = \gamma_{1}^{0*} (x) \gamma_{0}^{0*} (y)
 *                             + \gamma_{0}^{(0)*} (x) \gamma_{1}^{(0)*} (y)
 *
 *  Note that \gamma_n^{m*} = 0 when |m| > n as in \gamma_{0}^{(1)*} (x) and \gamma_{0}^{(-1)*} (x)
 *
 *  Let n = 1, m = 1
 *  Then for \gamma_1^1(x+y) = \sum_{k=0}^n sum_{l=-k}^k [ \gamma_{n-k}^{(m-l)*} (x) \gamma_k^{l*} (y) ]
 *                           = (     k = 0,     l = 0   =>   n - k = 1 - 0 = 1 and m -l = 1 - 0 = 1 =>
 *                               n - k = 1, m - l = 1
 *                             )
 *                              \gamma_{1}^{1*} (x) \gamma_{0}^{0*} (y)
 *                             (     k = 1,     l = -1  =>   n - k = 1 - 1 = 0 and m - l = 1 - (-1) = 2 =>
 *                               n - k = 0, m - l = 2
 *                             )
 *                             + \gamma_{0}^{(2)*} (x) \gamma_{1}^{(-1)*} (y)
 *                             (     k = 1,     l = 0   =>   n - k = 1 - 1 = 0 and m - l = 1 - (0) = 1 =>
 *                               n - k = 0, m - l = 1
 *                             )
 *                             + \gamma_{0}^{(1)*} (x) \gamma_{1}^{(0)*} (y)
 *                             (     k = 1,     l = 1   =>   n - k = 1 - 1 = 0 and m - l = 1 - (1) = 0 =>
 *                               n - k = 0, m - l = 0
 *                             )
 *                             + \gamma_{0}^{(0)*} (x) \gamma_{1}^{(1)*} (y)
 *                           = \gamma_{1}^{1*} (x) \gamma_{0}^{0*} (y)
 *                             + \gamma_{0}^{(0)*} (x) \gamma_{1}^{(1)*} (y)
 *
 *  Let n = 2, m = 0
 *  Then for \gamma_2^0(x+y) = \sum_{k=0}^2 sum_{l=-k}^k [ \gamma_{n-k}^{(m-l)*} (x) \gamma_k^{l*} (y) ]
 *                           = \gamma_{2}^{(0)*} (x) \gamma_0^{0*} (y)              // k = 0, l = 0
 *                             + \gamma_{1}^{(1)*}  (x)  \gamma_1^{(-1)*} (y)         // k = 1, l = -1
 *                             + \gamma_{1}^{(0)*}  (x)  \gamma_1^{(0)*}  (y)         // k = 1, l =  0
 *                             + \gamma_{1}^{(-1)*} (x)  \gamma_1^{(1)*}  (y)         // k = 1, l =  1
 *                             + \gamma_{0}^{(2)*}  (x)  \gamma_2^{(-2)*} (y)           // k = 2, l = -2
 *                             + \gamma_{0}^{(1)*}  (x)  \gamma_2^{(-1)*} (y)           // k = 2, l = -1
 *                             + \gamma_{0}^{(0)*}  (x)  \gamma_2^{(0)*} (y)            // k = 2, l =  0
 *                             + \gamma_{0}^{(-1)*} (x)  \gamma_2^{(1)*} (y)            // k = 2, l =  1
 *                             + \gamma_{0}^{(-2)*} (x)  \gamma_2^{(2)*} (y)            // k = 2, l =  2
 *                           = \gamma_{2}^{(0)*} (x) \gamma_0^{0*} (y)              // k = 0, l = 0
 *                             + \gamma_{1}^{(1)*}  (x)  \gamma_1^{(-1)*} (y)         // k = 1, l = -1
 *                             + \gamma_{1}^{(0)*}  (x)  \gamma_1^{(0)*}  (y)         // k = 1, l =  0
 *                             + \gamma_{1}^{(-1)*} (x)  \gamma_1^{(1)*}  (y)         // k = 1, l =  1
 *                             + \gamma_{0}^{(0)*}  (x)  \gamma_2^{(0)*} (y)            // k = 2, l =  0
 *
 *  Let n = 2, m = 1
 *  Then for \gamma_2^1(x+y) = \sum_{k=0}^2 sum_{l=-k}^k [ \gamma_{n-k}^{(m-l)*} (x) \gamma_k^{l*} (y) ]
 *                           = \gamma_{2}^{(1)*} (x) \gamma_0^{0*} (y)              // k = 0, l = 0
 *                             + \gamma_{1}^{(2)*}  (x)  \gamma_1^{(-1)*} (y)         // k = 1, l = -1
 *                             + \gamma_{1}^{(1)*}  (x)  \gamma_1^{(0)*}  (y)         // k = 1, l =  0
 *                             + \gamma_{1}^{(0)*} (x)  \gamma_1^{(1)*}  (y)         // k = 1, l =  1
 *                             + \gamma_{0}^{(3)*}  (x)  \gamma_2^{(-2)*} (y)           // k = 2, l = -2
 *                             + \gamma_{0}^{(2)*}  (x)  \gamma_2^{(-1)*} (y)           // k = 2, l = -1
 *                             + \gamma_{0}^{(1)*}  (x)  \gamma_2^{(0)*} (y)            // k = 2, l =  0
 *                             + \gamma_{0}^{(0)*} (x)  \gamma_2^{(1)*} (y)            // k = 2, l =  1
 *                             + \gamma_{0}^{(-1)*} (x)  \gamma_2^{(2)*} (y)            // k = 2, l =  2
 *                           = \gamma_{2}^{(1)*} (x) \gamma_0^{0*} (y)              // k = 0, l = 0
 *                             + \gamma_{1}^{(1)*}  (x)  \gamma_1^{(0)*}  (y)         // k = 1, l =  0
 *                             + \gamma_{0}^{(0)*} (x)  \gamma_2^{(1)*} (y)            // k = 2, l =  1
 *
 *  Let n = 2, m = 2
 *  Then for \gamma_2^2(x+y) = \sum_{k=0}^2 sum_{l=-k}^k [ \gamma_{n-k}^{(m-l)*} (x) \gamma_k^{l*} (y) ]
 *                           = \gamma_{2}^{(2)*} (x) \gamma_0^{0*} (y)              // k = 0, l = 0
 *                             + \gamma_{1}^{(3)*}  (x)  \gamma_1^{(-1)*} (y)         // k = 1, l = -1
 *                             + \gamma_{1}^{(2)*}  (x)  \gamma_1^{(0)*}  (y)         // k = 1, l =  0
 *                             + \gamma_{1}^{(1)*} (x)  \gamma_1^{(1)*}  (y)         // k = 1, l =  1
 *                             + \gamma_{0}^{(4)*}  (x)  \gamma_2^{(-2)*} (y)           // k = 2, l = -2
 *                             + \gamma_{0}^{(3)*}  (x)  \gamma_2^{(-1)*} (y)           // k = 2, l = -1
 *                             + \gamma_{0}^{(2)*}  (x)  \gamma_2^{(0)*} (y)            // k = 2, l =  0
 *                             + \gamma_{0}^{(1)*} (x)  \gamma_2^{(1)*} (y)             // k = 2, l =  1
 *                             + \gamma_{0}^{(0)*} (x)  \gamma_2^{(2)*} (y)             // k = 2, l =  2
 *                           = \gamma_{2}^{(2)*} (x) \gamma_0^{0*} (y)              // k = 0, l = 0
 *                             + \gamma_{1}^{(1)*} (x)  \gamma_1^{(1)*}  (y)          // k = 1, l =  1
 *                             + \gamma_{0}^{(0)*} (x)  \gamma_2^{(2)*} (y)             // k = 2, l =  2
 *
 *  Let n = 3, m = 0
 *  Then for \gamma_3^0(x+y) = \sum_{k=0}^3 sum_{l=-k}^k [ \gamma_{n-k}^{(m-l)*} (x) \gamma_k^{l*} (y) ]
 *                           = \gamma_{3}^{(0)*} (x) \gamma_0^{0*} (y)              // k = 0, l = 0
 *                             + \gamma_{2}^{(1)*}  (x)  \gamma_1^{(-1)*} (y)         // k = 1, l = -1
 *                             + \gamma_{2}^{(0)*}  (x)  \gamma_1^{(0)*}  (y)         // k = 1, l =  0
 *                             + \gamma_{2}^{(-1)*} (x)  \gamma_1^{(1)*}  (y)         // k = 1, l =  1
 *                               + \gamma_{1}^{(2)*}  (x)  \gamma_2^{(-2)*} (y)         // k = 2, l = -2
 *                               + \gamma_{1}^{(1)*}  (x)  \gamma_2^{(-1)*} (y)         // k = 2, l = -1
 *                               + \gamma_{1}^{(0)*}  (x)  \gamma_2^{(0)*} (y)          // k = 2, l =  0
 *                               + \gamma_{1}^{(-1)*} (x)  \gamma_2^{(1)*} (y)          // k = 2, l =  1
 *                               + \gamma_{1}^{(-2)*} (x)  \gamma_2^{(2)*} (y)          // k = 2, l =  2
 *                                 + \gamma_{0}^{(3)*}  (x)  \gamma_3^{(-3)*} (y)         // k = 3, l = -3
 *                                 + \gamma_{0}^{(2)*}  (x)  \gamma_3^{(-2)*} (y)         // k = 3, l = -2
 *                                 + \gamma_{0}^{(1)*}  (x)  \gamma_3^{(-1)*} (y)         // k = 3, l = -1
 *                                 + \gamma_{0}^{(0)*} (x)   \gamma_3^{(0)*} (y)          // k = 3, l =  0
 *                                 + \gamma_{0}^{(-1)*} (x)  \gamma_3^{(1)*} (y)          // k = 3, l =  1
 *                                 + \gamma_{0}^{(-2)*} (x)  \gamma_3^{(2)*} (y)          // k = 3, l =  2
 *                                 + \gamma_{0}^{(-3)*} (x)  \gamma_3^{(3)*} (y)          // k = 3, l =  3
 *                           = \gamma_{3}^{(0)*} (x) \gamma_0^{0*} (y)              // k = 0, l = 0
 *                             + \gamma_{2}^{(1)*}  (x)  \gamma_1^{(-1)*} (y)         // k = 1, l = -1
 *                             + \gamma_{2}^{(0)*}  (x)  \gamma_1^{(0)*}  (y)         // k = 1, l =  0
 *                             + \gamma_{2}^{(-1)*} (x)  \gamma_1^{(1)*}  (y)         // k = 1, l =  1
 *                               + \gamma_{1}^{(1)*}  (x)  \gamma_2^{(-1)*} (y)         // k = 2, l = -1
 *                               + \gamma_{1}^{(0)*}  (x)  \gamma_2^{(0)*} (y)          // k = 2, l =  0
 *                               + \gamma_{1}^{(-1)*} (x)  \gamma_2^{(1)*} (y)          // k = 2, l =  1
 *                                 + \gamma_{0}^{(0)*} (x)   \gamma_3^{(0)*} (y)          // k = 3, l =  0
 *
 *  Let n = 3, m = 1
 *  Then for \gamma_3^1(x+y) = \sum_{k=0}^3 sum_{l=-k}^k [ \gamma_{n-k}^{(m-l)*} (x) \gamma_k^{l*} (y) ]
 *                           = \gamma_{3}^{(1)*} (x) \gamma_0^{0*} (y)              // k = 0, l = 0
 *                             + \gamma_{2}^{(2)*}  (x)  \gamma_1^{(-1)*} (y)         // k = 1, l = -1
 *                             + \gamma_{2}^{(1)*}  (x)  \gamma_1^{(0)*}  (y)         // k = 1, l =  0
 *                             + \gamma_{2}^{(0)*}  (x)  \gamma_1^{(1)*}  (y)         // k = 1, l =  1
 *                               + \gamma_{1}^{(3)*}  (x)  \gamma_2^{(-2)*} (y)         // k = 2, l = -2
 *                               + \gamma_{1}^{(2)*}  (x)  \gamma_2^{(-1)*} (y)         // k = 2, l = -1
 *                               + \gamma_{1}^{(1)*}  (x)  \gamma_2^{(0)*} (y)          // k = 2, l =  0
 *                               + \gamma_{1}^{(0)*}  (x)  \gamma_2^{(1)*} (y)          // k = 2, l =  1
 *                               + \gamma_{1}^{(-1)*} (x)  \gamma_2^{(2)*} (y)          // k = 2, l =  2
 *                                 + \gamma_{0}^{(4)*}  (x)  \gamma_3^{(-3)*} (y)         // k = 3, l = -3
 *                                 + \gamma_{0}^{(3)*}  (x)  \gamma_3^{(-2)*} (y)         // k = 3, l = -2
 *                                 + \gamma_{0}^{(2)*}  (x)  \gamma_3^{(-1)*} (y)         // k = 3, l = -1
 *                                 + \gamma_{0}^{(1)*}  (x)  \gamma_3^{(0)*} (y)          // k = 3, l =  0
 *                                 + \gamma_{0}^{(0)*}  (x)  \gamma_3^{(1)*} (y)          // k = 3, l =  1
 *                                 + \gamma_{0}^{(-1)*} (x)  \gamma_3^{(2)*} (y)          // k = 3, l =  2
 *                                 + \gamma_{0}^{(-2)*} (x)  \gamma_3^{(3)*} (y)          // k = 3, l =  3
 *                           = \gamma_{3}^{(1)*} (x) \gamma_0^{0*} (y)              // k = 0, l = 0
 *                             + \gamma_{2}^{(2)*}  (x)  \gamma_1^{(-1)*} (y)         // k = 1, l = -1
 *                             + \gamma_{2}^{(1)*}  (x)  \gamma_1^{(0)*}  (y)         // k = 1, l =  0
 *                             + \gamma_{2}^{(0)*}  (x)  \gamma_1^{(1)*}  (y)         // k = 1, l =  1
 *                               + \gamma_{1}^{(1)*}   (x)  \gamma_2^{(0)*}  (y)         // k = 2, l = 0
 *                               + \gamma_{1}^{(0)*}   (x)  \gamma_2^{(1)*}  (y)         // k = 2, l = 1
 *                               + \gamma_{1}^{(-1)*}  (x)  \gamma_2^{(2)*}  (y)         // k = 2, l =  0
 *                                 + \gamma_{0}^{(0)*} (x)  \gamma_3^{(1)*}  (y)         // k = 3, l =  0
 *
 *  Let n = 3, m = 2
 *  Then for \gamma_3^2(x+y) = \sum_{k=0}^3 sum_{l=-k}^k [ \gamma_{n-k}^{(m-l)*} (x) \gamma_k^{l*} (y) ]
 *                           = \gamma_{3}^{(2)*} (x) \gamma_0^{0*} (y)              // k = 0, l = 0
 *                             + \gamma_{2}^{(3)*}  (x)  \gamma_1^{(-1)*} (y)         // k = 1, l = -1
 *                             + \gamma_{2}^{(2)*}  (x)  \gamma_1^{(0)*}  (y)         // k = 1, l =  0
 *                             + \gamma_{2}^{(1)*}  (x)  \gamma_1^{(1)*}  (y)         // k = 1, l =  1
 *                               + \gamma_{1}^{(4)*}  (x)  \gamma_2^{(-2)*} (y)         // k = 2, l = -2
 *                               + \gamma_{1}^{(3)*}  (x)  \gamma_2^{(-1)*} (y)         // k = 2, l = -1
 *                               + \gamma_{1}^{(2)*}  (x)  \gamma_2^{(0)*} (y)          // k = 2, l =  0
 *                               + \gamma_{1}^{(1)*}  (x)  \gamma_2^{(1)*} (y)          // k = 2, l =  1
 *                               + \gamma_{1}^{(0)*}  (x)  \gamma_2^{(2)*} (y)          // k = 2, l =  2
 *                                 + \gamma_{0}^{(5)*}  (x)  \gamma_3^{(-3)*} (y)         // k = 3, l = -3
 *                                 + \gamma_{0}^{(4)*}  (x)  \gamma_3^{(-2)*} (y)         // k = 3, l = -2
 *                                 + \gamma_{0}^{(3)*}  (x)  \gamma_3^{(-1)*} (y)         // k = 3, l = -1
 *                                 + \gamma_{0}^{(2)*}  (x)  \gamma_3^{(0)*} (y)          // k = 3, l =  0
 *                                 + \gamma_{0}^{(1)*}  (x)  \gamma_3^{(1)*} (y)          // k = 3, l =  1
 *                                 + \gamma_{0}^{(0)*}  (x)  \gamma_3^{(2)*} (y)          // k = 3, l =  2
 *                                 + \gamma_{0}^{(-1)*} (x)  \gamma_3^{(3)*} (y)          // k = 3, l =  3
 *                           = \gamma_{3}^{(2)*} (x) \gamma_0^{0*} (y)              // k = 0, l = 0
 *                             + \gamma_{2}^{(2)*}  (x)  \gamma_1^{(0)*} (y)         // k = 1, l =  0
 *                             + \gamma_{2}^{(1)*}  (x)  \gamma_1^{(1)*}  (y)         // k = 1, l =  1
 *                               + \gamma_{1}^{(1)*}   (x)  \gamma_2^{(1)*}  (y)        // k = 2, l = 1
 *                               + \gamma_{1}^{(0)*}   (x)  \gamma_2^{(2)*}  (y)        // k = 2, l = 2
 *                                 + \gamma_{0}^{(0)*} (x)  \gamma_3^{(2)*}  (y)          // k = 3, l =  2
 *
 *  Let n = 3, m = 3
 *  Then for \gamma_3^3(x+y) = \sum_{k=0}^3 sum_{l=-k}^k [ \gamma_{n-k}^{(m-l)*} (x) \gamma_k^{l*} (y) ]
 *                           = \gamma_{3}^{(3)*} (x) \gamma_0^{0*} (y)              // k = 0, l = 0
 *                             + \gamma_{2}^{(4)*}  (x)  \gamma_1^{(-1)*} (y)         // k = 1, l = -1
 *                             + \gamma_{2}^{(3)*}  (x)  \gamma_1^{(0)*}  (y)         // k = 1, l =  0
 *                             + \gamma_{2}^{(2)*}  (x)  \gamma_1^{(1)*}  (y)         // k = 1, l =  1
 *                               + \gamma_{1}^{(5)*}  (x)  \gamma_2^{(-2)*} (y)         // k = 2, l = -2
 *                               + \gamma_{1}^{(4)*}  (x)  \gamma_2^{(-1)*} (y)         // k = 2, l = -1
 *                               + \gamma_{1}^{(3)*}  (x)  \gamma_2^{(0)*} (y)          // k = 2, l =  0
 *                               + \gamma_{1}^{(2)*}  (x)  \gamma_2^{(1)*} (y)          // k = 2, l =  1
 *                               + \gamma_{1}^{(1)*}  (x)  \gamma_2^{(2)*} (y)          // k = 2, l =  2
 *                                 + \gamma_{0}^{(6)*}  (x)  \gamma_3^{(-3)*} (y)         // k = 3, l = -3
 *                                 + \gamma_{0}^{(5)*}  (x)  \gamma_3^{(-2)*} (y)         // k = 3, l = -2
 *                                 + \gamma_{0}^{(4)*}  (x)  \gamma_3^{(-1)*} (y)         // k = 3, l = -1
 *                                 + \gamma_{0}^{(3)*}  (x)  \gamma_3^{(0)*} (y)          // k = 3, l =  0
 *                                 + \gamma_{0}^{(3)*}  (x)  \gamma_3^{(1)*} (y)          // k = 3, l =  1
 *                                 + \gamma_{0}^{(1)*}  (x)  \gamma_3^{(2)*} (y)          // k = 3, l =  2
 *                                 + \gamma_{0}^{(0)*}  (x)  \gamma_3^{(3)*} (y)          // k = 3, l =  3
 *                           = \gamma_{3}^{(3)*} (x) \gamma_0^{0*} (y)              // k = 0, l = 0
 *                             + \gamma_{2}^{(2)*}  (x)  \gamma_1^{(1)*} (y)          // k = 1, l =  0
 *                               + \gamma_{1}^{(1)*}   (x)  \gamma_2^{(2)*}  (y)        // k = 2, l = 1
 *                                 + \gamma_{0}^{(0)*} (x)  \gamma_3^{(3)*}  (y)          // k = 3, l =  2
 *
 *  For all n,m we always get the first summation (k=0,l=0)
 *  Looking at n=3, m=0 with difference n - m = 3
 *    We always get the first summation (k=0,l=0)
 *    For the second summation (k=1) we get l = -1, l =  0, l =  1                      since n-k = 2 (-1 <= l <= 1) and  1 >= m-l >= -1
 *    For the third  summation (k=2) we get l = -1, l =  0, l =  1                      since n-k = 1 (-2 <= l <= 2) and  2 >= m-l >= -2
 *    For the fourth summation (k=3) we get l =  0                                      since n-k = 0 (-3 <= l <= 3) and  3 >= m-l >= -3
 *  Looking at n=3, m=1 with difference n - m = 2
 *    m = 1 implies the terms we get back same number of the terms to the left and right of the center of each summation (0 <= k <= n, l = 0)
 *    For the second summation (k=1) we get l = -1, l =  0, l =  1                      since n-k = 2 (-1 <= l <= 1) and  2 >= m-l >=  0
 *    For the third  summation (k=2) we get l =  0, l =  1, l =  2                      since n-k = 1 (-2 <= l <= 2) and  3 >= m-l >= -1
 *    For the fourth summation (k=3) we get l =  1                                      since n-k = 0 (-3 <= l <= 3) and  4 >= m-l >= -2
 *  Looking at n=3, m=2 with difference n - m = 1
 *    m = 1 implies the terms we get back same number of the terms to the left and right of the center of each summation (0 <= k <= n, l = 0)
 *    For the second summation (k=1) we get l =  0, l =  1                              since n-k = 2 (-1 <= l <= 1) and  3 >= m-l >=  1
 *    For the third  summation (k=2) we get l =  1, l =  2                              since n-k = 1 (-2 <= l <= 2) and  4 >= m-l >=  0
 *    For the fourth summation (k=3) we get l =  2                                      since n-k = 0 (-3 <= l <= 3) and  5 >= m-l >= -1
 *  Looking at n=3, m=3 with difference n - m = 1
 *    m = 1 implies the terms we get back same number of the terms to the left and right of the center of each summation (0 <= k <= n, l = 0)
 *    For the second summation (k=1) we get l =  1                                      since n-k = 2 (-1 <= l <= 1) and  4 >= m-l >=  2
 *    For the third  summation (k=2) we get l =  2                                      since n-k = 1 (-2 <= l <= 2) and  5 >= m-l >=  1
 *    For the fourth summation (k=3) we get l =  3                                      since n-k = 0 (-3 <= l <= 3) and  6 >= m-l >=  0
 *
 *  Looking at n=4, m=0 with difference n - m = 4
 *    We always get the first summation (k=0,l=0)
 *      \gamma_4^{(0)*} (x) \gamma_0^{(0)*} (y)
 *    For the second summation (k=1) we get l = -1, l =  0, l =  1,                      since n-k = 3 (-1 <= l <= 1) and  1 >= m-l >= -1
 *      \gamma_3^{(1)*} (x) \gamma_1^{(-1)*} (y)
 *      \gamma_3^{(0)*} (x) \gamma_1^{(0)*} (y)
 *      \gamma_3^{(-1)*} (x) \gamma_1^{(1)*} (y)
 *    For the third  summation (k=2) we get l = -2, l = -1, l =  0, l =  1, l =  2       since n-k = 2 (-2 <= l <= 2) and  2 >= m-l >= -2
 *      \gamma_2^{(2)*} (x) \gamma_2^{(-2)*} (y)
 *      \gamma_2^{(1)*} (x) \gamma_2^{(-1)*} (y)
 *      \gamma_2^{(0)*} (x) \gamma_2^{(0)*} (y)
 *      \gamma_2^{(-1)*} (x) \gamma_2^{(1)*} (y)
 *      \gamma_2^{(-2)*} (x) \gamma_2^{(2)*} (y)
 *    For the fourth summation (k=3) we get l = -1, l =  0, l =  1                       since n-k = 1 (-3 <= l <= 3) and  3 >= m-l >= -3
 *      \gamma_1^{(1)*} (x) \gamma_3^{(-1)*} (y)
 *      \gamma_1^{(0)*} (x) \gamma_3^{(0)*} (y)
 *      \gamma_1^{(-1)*} (x) \gamma_3^{(1)*} (y)
 *    For the fifth  summation (k=4) we get l =  0                                       since n-k = 0 (-4 <= l <= 4) and  4 >= m-l >= -4
 *      \gamma_0^{(0)*} (x) \gamma_4^{(0)*} (y)
 *
 *  Looking at n=4, m=1 with difference n - m = 3
 *    We always get the first summation (k=0,l=0)
 *      \gamma_4^{(1)*} (x) \gamma_0^{(0)*} (y)
 *    For the second summation (k=1) we get l = -1, l =  0, l =  1                      since n-k = 3 (-1 <= l <= 1) and  2 >= m-l >=  0
 *      \gamma_3^{(2)*} (x) \gamma_1^{(-1)*} (y)
 *      \gamma_3^{(1)*} (x) \gamma_1^{(0)*} (y)
 *      \gamma_3^{(0)*} (x) \gamma_1^{(1)*} (y)
 *    For the third  summation (k=2) we get l = -1, l =  0, l =  1, l =  2              since n-k = 2 (-2 <= l <= 2) and  3 >= m-l >= -1
 *      \gamma_2^{(2)*} (x) \gamma_2^{(-1)*} (y)
 *      \gamma_2^{(1)*} (x) \gamma_2^{(0)*} (y)
 *      \gamma_2^{(0)*} (x) \gamma_2^{(1)*} (y)
 *      \gamma_2^{(-1)*} (x) \gamma_2^{(2)*} (y)
 *    For the fourth summation (k=3) we get l = 0, l =  1, l =  2                      since n-k = 1 (-3 <= l <= 3) and  4 >= m-l >= -2
 *      \gamma_1^{(1)*} (x) \gamma_3^{(0)*} (y)
 *      \gamma_1^{(0)*} (x) \gamma_3^{(1)*} (y)
 *      \gamma_1^{(-1)*} (x) \gamma_3^{(2)*} (y)
 *    For the fifth summation (k=4) we get l =  1                                      since n-k = 0 (-3 <= l <= 3) and  4 >= m-l >= -2
 *      \gamma_0^{(0)*} (x) \gamma_4^{(1)*} (y)
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
 * From the 1st part of the upward pass, the matrix above has been created for \gamma_n^{m*} (x)
 * We first need to create the matrix above for \gamma_n^{m*} (y)
 * We then form the matrix for \gamma_n^{m*} (x+y)
 *
 * To obtain each \gamma_n^{m*} (x+y)
 * We nested loop over
 *   0 <= k <= n
 *     -k <= l <= k
 *        check to make sure that n-k >= m - l
 *        Use the two matrices above to form the product
 *          \gamma_{(n-k)}^{(m-l)*} (x) \gamma_{k}^{l*} (y)
 *
 * The resulting \gamma_n^{m*} (x+y) is our coefficient matrix for the power series centered at (x+y)
 * having powers \theta_n^{m*} (x_b - z_A^{(p)})
 * where z_A^{(p)} is the center of the parent cell containing the source x_a
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

  std::complex<double> s_expansion = 0.0;

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

  // obtaining far-field coefficients (s-expansion coefficients)
  // far field coefficients are (z_a - x_a)^{\alpha} = (center_yBox - y[0])^{\alpha}
  // Note: y[0] is the target ( in Main.cc it is x )
  std::cout << "We are here" << std::endl;
  std::vector<std::vector<std::complex<double> > > y_SCoeff = potential.getSCoeff(x[0].getCoord(), center_xBox.getCoord());

  std::vector<std::vector<std::complex<double> > > SCoeffTranslated = potential.getSS(center_xBox.getCoord(), center_xBoxParent.getCoord(),
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
      std::cout << "SCoeffTranslated[" << i << "][" << j << "] = "
                << SCoeffTranslated[i][j] << std::endl;
  // obtaining far-field powers (s-expansion powers)
  // far field powers are D^{\alpha} f (x_b - z_a) = D^{\alpha} f (x[0] - center_yBox)
  std::vector<std::vector<std::complex<double> > > y_SVec = potential.getSVector(y[0].getCoord(), center_xBoxParent.getCoord());

  for (unsigned int i=0; i<=abs_alpha; ++i)
    for (unsigned int j=0; j<=abs_alpha; ++j)
      std::cout << "y_SVec[" << i << "][" << j << "] = " << y_SVec[i][j] << std::endl;

  for (unsigned int i=0; i<=abs_alpha; ++i)
    for (unsigned int j=0; j<=abs_alpha; ++j)
    {
      s_expansion = s_expansion + SCoeffTranslated[j][i] * y_SVec[j][i];
    }

  std::cout << "s_expansion = " << s_expansion.real() << " + i" << s_expansion.imag() << std::endl;

  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "**************************************************" << std::endl;
  std::cout << "**************************************************" << std::endl;

}
