/**
 * Main.cc
 *
 *  @date Created on: Dec 19, 2016
 *  @author Keith D Brauss
 */

/**
 * @mainpage
 *
 * @section section_toc Table of Contents
 * <ul>
 *   <li> @ref error_in_one_dimension </li>
 *     <ul>
 *       <li> @ref motivation </li>
 *       <li> @ref single_variable_potential_function </li>
 *       <li> @ref fmm_expansions </li>
 *         <ul>
 *           <li> @ref far_field_expansion </li>
 *           <li> @ref near_field_expansion </li>
 *         </ul>
 *       <li> @ref approximation_error </li>
 *         <ul>
 *           <li> @ref rate_of_convergence </li>
 *           <li> @ref upward_pass_part_one </li>
 *           <li> @ref upward_pass_part_two </li>
 *           <li> @ref downward_pass_part_one </li>
 *           <li> @ref downward_pass_part_two </li>
 *         </ul>
 *       <li>
 *     </ul>
 * </ul>
 *
 * @section error_in_one_dimension Error Analysis of the Fast Multipole Method for Single-Variable Potential Function
 *
 * @subsection motivation Motivation: MHD Problem Goal
 *
 * We wish to approximate the Biot-Savart magnetic field calculation used in the numerical solution of
 * the MHD velocity-current formulation.  A direction calculation of the Biot-Savart integral for all \f$N\f$ quadrature
 * points in the domain requires \f$O(N^2)\f$ operations.  The Fast Multipole Method can bring the operation
 * count down to \f$O(N \log(N))\f$.  We wish to maintain the same rates of convergence for the unknowns that are obtained
 * using the direct calculation.
 *
 * The only unknown encountered in the Biot-Savart integral is the current density.  Therefore, the rate of convergence
 * of the magnetic field will affect the convergence rate of the current density which will in turn affect the other unknowns.
 * If we can maintain the same convergence rate for the magnetic field seen in the direction calculations of the Biot-Savart
 * integral, then we expect the convergence rate of the current density to hold as well.
 *
 * The Biot-Savart integral maps current densities in \f$L_2(\Omega)\f$ to \f$L_2(\Omega)\f$.
 * Consider the convergence of the magnetic
 * field when calculating the Biot-Savart integral directly.  We use the fact that the operator \f$\mathbf{B}(\mathbf{J})(\mathbf{x})\f$
 * is a bounded linear operator with \f$\|\mathbf{B}(\mathbf{J})\|_{L_2(\Omega)} \leq C \|\mathbf{J}\|_{L_2(\Omega)}\f$.  Then
 *
 * \f{eqnarray*}{
 *   \| \mathbf{B}(\mathbf{J}) - \mathbf{B}(\mathbf{J}^h) \|_{L_2(\Omega)}
 *      =
 *         \| \mathbf{B}(\mathbf{J} - \mathbf{J}^h) \|_{L_2(\Omega)}
 *      \leq
 *         C \| \mathbf{J} - \mathbf{J}^h \|_{L_2(\Omega)}
 * \f}
 *
 * The finite elements we have chosen for \f$\mathbf{J}\f$ result in a quadratic rate of convergence.
 * Therefore, the magnetic field also experiences quadratic convergence according to the result above.
 *
 * We wish to show that approximating the Biot-Savart integral with the Fast Multipole Method will not
 * affect this rate of convergence.  Let \f$B_{FMM}\f$ denote the Fast Multipole Method approximation to the
 * Biot-Savart integral and \f$T_n\f$ the Taylor series approximation to
 * \f$\frac{\mathbf{x}-\mathbf{y}}{\|\mathbf{x}-\mathbf{y}\|^3}\f$ that accompanies the method.  Then
 *
 * \f{eqnarray*}{
 *   \|B(J)-B_{FMM}(J^h)\|
 * \f}
 *
 *
 *
 *
 * @subsection single_variable_potential_function Single Variable Potential Function
 *
 * The single variable potential function may be the easiest case to examine to understand the error
 * encountered in the Fast Multipole Method.  The single-variable potential function requires a smaller
 * amount of bookkeeping.  For our potential function, we have
 *
 * \f[ \psi(z) = \frac{1}{z} \f]
 *
 * We are interested in the interaction of a source particle \f$x\f$ with a target particle \f$y\f$
 * defined by the potential function.  The expression of interest is the potential function in the
 * form
 *
 * \f[ \psi(y-x) = \frac{1}{y-x} \f]
 *
 *
 * @subsection fmm_expansions FMM Expansions of the Potential Function
 *
 * Recall that FMM expresses the potential function as a series (far-field and near-field series) to
 * minimize the operation count required to calculate the potential function for multipole source and
 * target particles by expressing the terms of the series as a product of two expressions having only
 * x and y respectively (a separation of variables).
 *
 * To write the single-variable potential function
 * above as a series, recall the geometric series \f$\frac{1}{1-x} = \sum_{n=0}^{\infty} x^n\f$.
 *
 *
 * @subsubsection far_field_expansion Far-Field Expansion of the Potential Function
 *
 * Factoring the potential function, we can apply the formula for the series
 *
 * \f[ \psi(y-x) = \frac{1}{y-x} = \frac{1}{y} \frac{1}{1-\frac{x}{y}} = \frac{1}{y} \sum_{n=0}^{\infty} \left( \frac{x}{y} \right)^{n}
 *               = \sum_{n=0}^{\infty} \frac{x^n}{y^{n+1}} \f]
 *
 * Following the Fast Multipole Method, we introduce a center \f$c\f$
 *
 * \f[ \psi(y-x) = \psi((y-c)-(x-c)) = \frac{1}{(y-c)-(x-c)} = \sum_{n=0}^{\infty} \frac{(x-c)^n}{(y-c)^{n+1}}  \f]
 *
 * Since the geometric series \f$ \frac{1}{1-x} = \sum_{n=0}^{\infty} x^n \f$ converges only if \f$ |x| < 1 \f$, we have
 * convergence of the series above if
 *
 * \f[ \left| \frac{x-c}{y-c} \right| < 1  \quad \mbox{ or } \quad   |x-c| < |y-c| \f].
 *
 * If \f$x\f$ is closer to \f$c\f$ than \f$y\f$, the series converges.  We refer to this expansion (series) as the far-field expansion.
 * In the Fast Multipole Method \f$c\f$ will be the center of the cell containing the source particle \f$x\f$.  \f$y\f$ will therefore be far from
 * this cell, center and source particle.
 *
 *
 * @subsubsection near_field_expansion Near-Field Expansion of the Potential Function
 *
 * With a little algebra we can write another series representation for the potential function
 *
 * \f[ \psi(y-x) = \frac{1}{y-x} = \frac{1}{y-c-(x-c)} = \frac{-1}{x-c} \frac{1}{1-\frac{y-c}{x-c}} = \frac{-1}{x} \sum_{n=0}^{\infty} \left( \frac{y-c}{x-c} \right)^{n}
 *               = \sum_{n=0}^{\infty} \frac{-(y-c)^n}{(x-c)^{n+1}} \f]
 *
 * We can see the series converges if
 *
 * \f[ \left| \frac{y-c}{x-c} \right| < 1  \quad \mbox{ or } \quad   |y-c| < |x-c| \f].
 *
 * If \f$y\f$ is closer to \f$c\f$ than \f$x\f$, the series converges.  We refer to this expansion (series) as the near-field expansion.
 * In our fast multipole method \f$c\f$ will be the center of the cell containing the target particle \f$y\f$.  \f$y\f$ will be near \f$c\f$.
 *
 * We have several levels of refinement for our domain using the Fast Multipole Method.  We take the domain to be the unit line [0,1].
 * The first level of refinement (level l = 0) will consist of one cell - covering the domain.  In one dimension the cells will be
 * intervals.  The single cell at l = 0 is the interval (0,1).  As we work down the levels the domain becomes more refined with
 * number of cells covering the domain increasing.
 *
 * At level 2 (level l = 1), there will be \f$2^1 = 2\f$ cells (0,1/2) and (1/2,1).
 * At level \f$l=3\f$ there will be \f$2^2 = 4\f$ cells: \f$(0,1/4)\f$,\f$(1/4,1/2)\f$,\f$(1/2,3/4)\f$, and \f$(3/4,1)\f$.
 * At level \f$l = 4\f$ there will be \f$2^3 = 8\f$ cells:
 *
 * \f[ (0,1/8),(1/8,2/8),(2/8,3/8),(3/8,4/8),(4/8,5/8),(5/8,6/8),(6/8,7/8), \quad \mbox{ and } \quad (7/8,8/8) \f].
 *
 * The length of a cell at the lowest level we denote as \f$h\f$.
 *
 *
 * @subsection approximation_error Error Due to Approximation
 *
 * The far-field and near-field expansions are truncations of an infinite series and they will differ from the direct calculation
 * of the potential function.  The FMM is an approximation method and errors are encountered.
 *
 * Further approximations occur in the Fast Multipole Method during
 * translations of the center for near and far-field expansions to other locations.  The two translations are called R|R and S|S translations,
 * respectively near-field to near-field and far-field to far-field.  The FMM also has an S|R transformation from a far-field expansion to
 * a near-field expansion.  At the same time the S|R transformation translates the center from being near a source point to being near a target
 * point.
 *
 * All the translations above consist are constructed using infinite series (an infinte matrix when put together).  For calculation purposes all series
 * are truncated and the translations have approximation error.
 *
 * As we increase the number of refinement levels we use in the Fast Multipole Method, we would like the approximation error to decrease with respect
 * to the lowest level \f$l = L\f$ (grid most refined).
 *
 * Again, let the size of a cell at the lowest FMM level be denoted by \f$h\f$. Level \f$l = 0\f$ has \f$h_0 = 1\f$, level \f$l = 1\f$ has
 * \f$h_1 = \frac{1}{2} h_0 = \frac{1}{2}\f$, level \f$l = 2\f$ has \f$h_2 = \frac{1}{2} h_1 = \left( \frac{1}{2} \right)^2\f$, level
 * \f$l = 3\f$ has \f$h_3 = \frac{1}{2} h_2 = \left( \frac{1}{2} \right)^3\f$, and in general level
 *
 * \f[ l = L \quad \mbox{ has } \quad h_L = \frac{1}{2} h_{L-1} = \left( \frac{1}{2} \right)^L. \f]
 *
 *
 * @subsubsection rate_of_convergence A Convergence Rate
 *
 * We would like the rate at which the approximation error decreases with respect to the cell size \f$h\f$ to be some power of the cell size \f$O(h^n)\f$.
 * Let our desired rate of convergence be \f$O(h^2)\f$.  Since the near-field and far-field expansions are Taylor series expansions
 *
 * \f[ f(x) = \sum_{n=0}^{\infty} \frac{f^{(n)}(a)}{n!} (x-a)^n  = T_n(x) + R_n(x) \f]
 *
 * where the Taylor polynomial \f$T_n\f$ and remainder \f$R_n\f$ are defined as
 *
 * \f[ T_n(x) = \sum_{i=0}^{n} \frac{f^{(i)}(a)}{i!} (x-a)^i \quad \mbox{ and } \quad
 *     R_n(x) = \sum_{i=n+1}^{\infty} \frac{f^{(i)}(a)}{i!} (x-a)^i
 * \f]
 *
 * and \f$ |f^{(n+1)} (x) | \leq M \f$ for \f$|x-a| \leq d\f$ (\f$f\f$ bounded in a neighborhood around its center \f$a\f$) implies
 *
 * \f[ | R_n(x) | \leq  \frac{ \max \left( | f^{(n+1)}(x) | \right) }{(n+1)!} |x-a|^{n+1} \quad \mbox{ for } a - d \leq x \leq a + d \f]
 *
 * The remainder \f$|R_n(x)|\f$ is the error in approximating a Taylor series by its \f$n\f$th Taylor polynomial.  We wish the error
 * \f$|R_n(x)|\f$ to be a function of the lowest level mesh size \f$h\f$.  Then as we increase the number of FMM levels and the
 * lowest level cell size decreases, the error will decrease proportionally.
 *
 *
 * @subsubsection upward_pass_part_one Upward Pass Part One
 *
 * To construct the Fast Multipole Method so that it converges at a rate of \f$O(h^2)\f$ we start with the first step of the FMM - the Upward Pass.
 * The lowest FMM level for our example is \f$l=5\f$ (highest refinement level).
 *
 * The first part of the Upward Pass determines far-field expansions (just the coefficients \f$a_n = (x-x_c)^n\f$
 * of the power series \f$\sum_{n=0}^{\infty} a_n S_n \f$ )
 *
 * \f[ \sum_{n=0}^{\infty} a_n S_n =  \sum_{n=0}^{\infty} (x-x_c)^n \frac{1}{(y-x_c)^{n+1}} \f]
 *
 * associated with each source particle \f$x\f$, where \f$x_c\f$ is the center of the cell containing the source particle \f$x\f$.
 * Previously factoring the terms of the series into two expressions that separate the source \f$y\f$ and target \f$x\f$ allows us to
 * form the coefficients of the series without knowing the source particle and only knowing the target particle.
 *
 * Let \f$z = y-x_c\f$ and \f$t = x_c - x\f$.  Then \f$z+t = y-x\f$ and the Taylor series for \f$\phi(z+t) = \frac{1}{z+t}\f$ is
 *
 * \f{eqnarray*}{
 *     \phi(z+t) & = & \sum_{n=0}^{\infty} \frac{\phi^{(n)}(z)}{n!} t^n
 *                 = \frac{1}{0! z} t^0 + \frac{(-1)^1}{1! z^2} t^1 + \frac{(-1)^2 2!}{2! z^3} t^2 + \frac{(-1)^3 3!}{3! z^4} t^3
 *                   + \hdots + \frac{(-1)^n n!}{n! z^{n+1}} t^n + R_n(z+t) \\
 *               & = &
 *                   \frac{1}{z} + \frac{(-1)^1}{z^2} t + \frac{(-1)^2}{z^3} t^2 + \frac{(-1)^3}{z^4} t^3
 *                   + \hdots + \frac{(-1)^n}{z^{n+1}} t^n + R_n(z+t) \\
 *               & = &
 *                   \frac{1}{y-x_c} + \frac{(-1)^1}{(y-x_c)^2} (x_c-x) + \frac{(-1)^2}{(y-x_c)^3} (x_c-x)^2 + \frac{(-1)^3}{(y-x_c)^4} (x_c-x)^3
 *                   + \hdots + \frac{(-1)^n}{(y-x_c)^{n+1}} (x_c-x)^n + R_n(y-x) \\
 *               & = &
 *                   \frac{1}{y-x_c} + \frac{(-1)^1}{(y-x_c)^2} (-1)^1(x-x_c) + \frac{(-1)^2}{(y-x_c)^3} (-1)^2(x-x_c)^2 + \frac{(-1)^3}{(y-x_c)^4} (-1)^3(x-x_c)^3
 *                   + \hdots + \frac{(-1)^n}{(y-x_c)^{n+1}} (-1)^n(x-x_c)^n + R_n(y-x) \\
 *               & = &
 *                   \frac{1}{y-x_c} + \frac{(x-x_c)}{(y-x_c)^2}  + \frac{(x-x_c)^2}{(y-x_c)^3}  + \frac{(x-x_c)^3}{(y-x_c)^4}
 *                   + \hdots + \frac{(x-x_c)^n}{(y-x_c)^{n+1}}  + R_n(y-x) \\
 *               & = &
 *                   \sum_{k=0}^{n} \frac{(x-x_c)^k}{(y-x_c)^{k+1}} + R_n(y-x)
 *                 =
 *                   \sum_{k=0}^{n} a_k S_k + R_n(y-x)
 * \f}
 *
 * We see \f$\phi^{(n)}(z) = \frac{(-1)^n n!}{z^{n+1}}\f$.  Replacing \f$f = \phi\f$, \f$x = z + t \f$, \f$a = z \f$ and \f$x-a = t\f$ in
 *
 * \f[ | R_n(x) | \leq \frac{ max \left( | f^{(n+1)}(x) | \right) }{(n+1)!} |x-a|^{n+1} \f]
 *
 * we have for \f$|t| \leq d\f$
 *
 * \f{eqnarray*}{
 * | R_n(z+t) |
 *              & \leq &
 *                        \frac{ max \left( | \phi^{(n+1)}(z+t) | \right) }{(n+1)!} |t|^{n+1}.
 * \f}
 *
 * We switch to the source cell center \f$x_c\f$, target \f$y\f$, and source \f$x\f$ variables using
 * \f$z = y-x_c\f$, \f$t = x_c - x\f$, and \f$z+t = y-x\f$.  For \f$|t| = |x-x_c| < d\f$
 *
 * \f{eqnarray*}{
 * | R_n(y-x) |
 *              & \leq &
 *                        \frac{ max \left( | \phi^{(n+1)}(y-x) | \right) }{(n+1)!} |x_c - x|^{n+1}.
 * \f}
 *
 * Since \f$|t| = |x-x_c| < \frac{1}{2} h\f$, define \f$d := \frac{1}{2} h\f$.  We note that \f$|x-x_c| < \frac{\sqrt{2}}{2}\f$
 * in two dimensions and \f$|x-x_c| < \frac{\sqrt{3}}{2}\f$ in three dimensions.
 *
 * The far-field expansions are only used when the source cell is outside the near neighbors of the target cell
 * and in the interactive list or beyond.  This means
 *
 * \f[
 *      |y-x| > 1.0h, \quad |y-x_c| > 1.5h, \quad \mbox{ and } \quad  |x-x_c| < \frac{1}{2} h \quad \implies \quad
 *      \frac{|x-x_c|}{|y-x_c|} < \frac{\frac{1}{2} h}{1.5h} < \frac{0.5h}{1.5h} = \frac{1}{3}.
 * \f]
 *
 * Therefore,
 *
 * \f{eqnarray*}{
 * | R_n(y-x) |
 *              & \leq &
 *                        \frac{ max \left( | \phi^{(n+1)}(y-x) | \right) }{(n+1)!} |x_c - x|^{n+1} \\
 *              & = &
 *                        max \left( \left| \frac{(-1)^{n+1} (n+1)!}{(y-x)^{n+2}} \right| \right) \frac{1}{(n+1)!} |x_c - x|^{n+1} \\
 *              & \leq &
 *                        \frac{1}{h^{n+2}} \left( \frac{1}{2} h \right)^{n+1} \\
 *              & = &
 *                        \frac{1}{h} \left( \frac{1}{2} \right)^{n+1}
 * \f}
 *
 * Recall that we want the error in the far-field approximation to decrease as we increase the number of levels we
 * use in the FMM.  Again, \f$h\f$ represents the length of a cell at the lowest FMM level (highest refinement level \f$l = L\f$)
 * and we want the error to decrease by \f$O(h^2)\f$ as the number of FMM levels increases.
 *
 * For the unit line \f$[0,1]\f$, we have \f$h=1\f$ at the highest FMM level (coarsest grid \f$l = 0\f$),
 *
 * \f[ h=\frac{1}{2} \mbox{ at } l = 1, \quad h = \left( \frac{1}{2} \right)^2 \mbox{ at } l = 2,
 *     \quad \hdots \quad h = \left( \frac{1}{2} \right)^L \mbox{ at } l = L
 * \f]
 *
 * We would like
 *
 * \f[  | R_n(y-x) | \leq O(h^2) \f]
 *
 * Therefore, for \f$l = k\f$, we have \f$h = \left( \frac{1}{2} \right)^{k}\f$
 *
 * \f{eqnarray*}{
 *    |R_n(y-x)|
 *    & \leq &
 *                 \frac{1}{h} \left( \frac{1}{2} \right)^{n+1} \\
 *    &  =   &
 *                 \frac{1}{h} \left( \frac{1}{2} \right)^{k} \left( \left( \frac{1}{2} \right)^{k} \right)^2 \left( \frac{1}{2} \right)^{n+1-k-2k} \\
 *    &  =   &
 *                 h^2 \left( \frac{1}{2} \right)^{n - (3k-1)}
 * \f}
 *
 * If we choose \f$n = 3k-1\f$ then \f$|R_n(y-x)| \leq h^2\f$, and we have the desired rate of convergence for the far-field expansion with respect to
 * increased refinement levels.  The fast multipole method can be used when the lowest level is l = 2.  The number of terms needed for the far-field
 * expansion is \f$n = 3k-1 = 3(2) - 1 = 5\f$.
 *
 * \f{eqnarray*}{
 *   l & = & 1 \quad \implies \quad n = 3k-1 = 3(1) - 1 = 2  \\
 *   l & = & 2 \quad \implies \quad n = 3k-1 = 3(2) - 1 = 5  \\
 *   l & = & 3 \quad \implies \quad n = 3k-1 = 3(3) - 1 = 8  \\
 *   l & = & 4 \quad \implies \quad n = 3k-1 = 3(4) - 1 = 11  \\
 *   l & = & 5 \quad \implies \quad n = 3k-1 = 3(5) - 1 = 14  \\
 *     & \hdots &  \\
 *   l & = & L \quad \implies \quad n = 3k-1 = 3L - 1
 * \f}
 *
 * The terms needed for \f$O(h^2)\f$ convergence of the far-field expansion
 * increases by 3 with each increment in number of refinement levels used by the Fast Multipole Method.
 *
 * Since the lowest level \f$l = 5\f$, then 15 terms (\f$n = 14\f$ with index starting on zero) are needed for a far-field expansion
 * to have error of order \f$O(h^2)\f$.
 *
 *
 * @subsubsection upward_pass_part_two Upward Pass Part Two
 *
 * We continue to the next step in the Fast Multipole Method - the second part of the Upward Pass.  In the
 * second part of the Upward Pass the far-field expansions are translated from the lowest level cell centers
 * upward to the cell centers at the next level.  The translation process is continued until reaching level \f$l = 2\f$.
 *
 * We translate the far-field expansion \f$\sum_{n=0}^{\infty} a_n(x-x_c^{(c)}) S_n(y-x_c^{(c)})\f$ where
 * \f$S_k(y-x_c^{(c)}) = \frac{1}{(y-x_c^{(c)})^{k+1}}\f$ from the child center \f$x_c^{(c)}\f$
 * to the far-field expansion \f$\sum_{n=0}^{\infty} a_n(x-x_c^{(p)}) S_n(y-x_c^{(p)})\f$ at the parent center \f$x_c^{(p)}\f$
 * by performing a Taylor series for each power \f$S_n(y-x_c^{(c)})\f$ of the far-field expansion to be translated.  We write
 * the resulting series in terms of the powers \f$S_n(y-x_c^{(p)})\f$.
 *
 * Let \f$z+t = y-x_c^{(c)}\f$, \f$z = y-x_c^{(p)}\f$, and \f$t = z+t - z = x_c^{(p)} - x_c^{(c)}\f$.
 * Note that \f$|y-x_c^{(c)}| > \frac{3}{2} h\f$ and \f$|x_c^{(p)} - x_c^{(c)}| \leq \frac{1}{2} h\f$.
 *
 * \f$\mathbf{S_0} - \mathbf{\mbox{1st Power}}\f$
 *
 * \f{eqnarray*}{
 *     S_0(z+t) & = & \frac{1}{z+t} = \sum_{n=0}^{\infty} \frac{\phi^{(n)}(z)}{n!} t^n \\
 *              & = &
 *                   \frac{1}{0! z} t^0 + \frac{(-1)^1}{1! z^2} t^1 + \frac{(-1)^2 2!}{2! z^3} t^2 + \frac{(-1)^3 3!}{3! z^4} t^3
 *                   + \hdots + \frac{(-1)^n n!}{n! z^{n+1}} t^n + R_n(z+t) \\
 *               & = &
 *                   S_0(z) + (-1)^1 t S_1(z) + (-1)^2 t^2 S_2(z) + (-1)^3 t^3 S_3(z)
 *                   + \hdots + (-1)^n t^n S_n(z) + R_n(z+t) \\
 *               & = & \sum_{k=0}^{n} (-1)^{k} t^k S_k(z) + R_n(z+t)
 * \f}
 *
 * We have \f$\phi^{(n)}(z) = \frac{(-1)^n n!}{z^{n+1}}\f$.  For \f$l = k\f$ and \f$h = \left( \frac{1}{2} \right)^k\f$ we have
 *
 * \f{eqnarray*}{
 *    |R_n(z+t)| & = & |R_n(y-x_c^{(c)})|
 *                       \leq
 *                                 \frac{ max \left( | \phi^{(n+1)}(y-x_c^{(c)}) | \right) }{(n+1)!} |x_c^{(p)} - x_c^{(c)}|^{n+1} \\
 *                       & = &  \max \left( \left| \frac{(-1)^{n+1} (n+1)!}{(y-x_c^{(c)})^{n+2}} \right| \right)  \frac{1}{(n+1)!}
 *                                                        |x_c^{(p)} - x_c^{(c)}|^{n+1} \\
 *                       & \leq & \frac{1}{(1.5h)^{n+2}} \left( \frac{1}{2} h \right)^{n+1} \\
 *                       & = & \frac{2}{3h} \left( \frac{1}{3} \right)^{n+1} \\
 *                       & = &
 *                             \frac{2}{3h} \left( \frac{1}{3} \right)^{k} \left( \left( \frac{1}{3} \right)^{k} \right)^2 \left( \frac{1}{3} \right)^{n+1-k-2k} \\
 *                       & \leq &
 *                             \frac{2}{3} h^2 \left( \frac{1}{3} \right)^{n-(3k-1)}
 * \f}
 *
 * We take the lowest FMM level to be \f$l = 5\f$ with \f$h = h_5 = \left( \frac{1}{2} \right)^5\f$ and choose \f$n = 14\f$ such that
 *
 * \f{eqnarray*}{
 *    |R_n(y-x_c^{(c)})|
 *               & \leq &
 *                        \frac{2}{3} h_3^2 \left( \frac{1}{3} \right)^{n-(3k-1)} \\
 *               & = &
 *                        \frac{2}{3} h_3^2 \left( \frac{1}{3} \right)^{14-(3(5)-1)} \\
 *               & = &    \frac{2}{3} h_3^2   \quad \quad ( = C h_3^2 )
 * \f}
 *
 * The bounding constant \f$C\f$ will grow with each power (see below).  We want to keep \f$C\f$ the
 * same for each power, else the error depends on the number of terms used to approximate a power which
 * depends on the lowest FMM level and therefore the size of \f$h\f$.
 *
 *
 * \f$\mathbf{S_1} - \mathbf{\mbox{Second Power}}\f$
 *
 * For the second power \f$S_1\f$ of the far-field expansion
 *
 * \f{eqnarray*}{
 *     S_1(z+t) & = & \frac{1}{(z+t)^2} = \sum_{n=0}^{\infty} \frac{S_1^{(n)}(z)}{n!} t^n \\
 *              & = &
 *                   \frac{1}{0! z^2} t^0 + \frac{(-1)^1 2!}{1! z^3} t^1 + \frac{(-1)^2 3!}{2! z^4} t^2 + \frac{(-1)^3 4!}{3! z^5} t^3
 *                   + \hdots + \frac{(-1)^n (n+1)!}{n! z^{n+2}} t^n + R_n(x) \\
 *               & = &
 *                   S_1(z) + \frac{(-1)^1 2!}{1!} t S_1(z) + \frac{(-1)^2 3!}{2!} t^2 S_2(z) + \frac{(-1)^3 4!}{3!} t^3 S_3(z)
 *                   + \hdots + \frac{(-1)^n (n+1)!}{n!} t^n S_{n+1} (z) + R_n(x) \\
 *               & = & \sum_{k=0}^{n} (-1)^{k} \begin{pmatrix} k+1 \\ 1 \end{pmatrix} t^k S_{k+1} (z) + R_n(z+t)
 * \f}
 *
 * where \f$l = k\f$ and \f$h = \left(\frac{1}{2}\right)^k\f$ implies
 *
 * \f{eqnarray*}{
 * | R_n(y-x_c^{(c)}) |
 *              & \leq &
 *                        \frac{ max \left( | S_1^{(n+1)}(y-x_c^{(c)}) | \right) }{(n+1)!} |x_c^{(p)} - x_c^{(c)}|^{n+1} \\
 *              & = &
 *                        max \left( \left| \frac{(-1)^{n+1} (n+2)!}{1!(y-x_c^{(c)})^{n+3}} \right| \right) \frac{1}{(n+1)!} |x_c^{(p)} - x_c^{(c)}|^{n+1} \\
 *              & \leq &
 *                        \frac{(n+2)}{1! \left( \frac{3}{2} h \right)^{n+3}} \left( \frac{1}{2} h \right)^{n+1} \\
 *              & = &
 *                        \frac{(n+2)}{1!} \left(\frac{2}{3}\right)^2 \left(\frac{1}{h}\right)^2 \left( \frac{1}{3} \right)^{n+1} \\
 *              & \leq &
 *                        \frac{(n+2)}{1!} \left(\frac{1}{h}\right)^2 \left( \frac{1}{3} \right)^{n+1} \\
 *              & = &
 *                        \frac{(n+2)}{1!} \left(\frac{1}{h}\right)^2
 *                            \left( \left( \frac{1}{3} \right)^{k} \right)^2 \left( \left( \frac{1}{3} \right)^{k} \right)^2 \left( \frac{1}{3} \right)^{n+1-2k-2k} \\
 *              & = &
 *                        \begin{pmatrix} n+2 \\ 1 \end{pmatrix}  h^2 \left( \frac{1}{3} \right)^{n-(4k-1)} \\
 * \f}
 *
 * The lowest FMM level for our example is \f$l = k = 5\f$ and this implies \f$h = h_5 = \left( \frac{1}{2} \right)^5\f$
 *
 * \f{eqnarray*}{
 *      | R_n(y-x_c^{(c)}) |
 *               & \leq &
 *                        (n+2) h_5^2  \left( \frac{1}{3} \right)^{n-(4(5)-1)} \\
 *               & = &     (n+2) h_5^2  \left( \frac{1}{3} \right)^{n-19}
 * \f}
 *
 * If \f$n = 19\f$
 *
 * \f[    | R_n(y-x_c^{(c)}) | \leq  21 h_3^2 = C h_3^2 \f]
 *
 * If we wish to keep the bounding constant \f$C\f$ fixed  with \f$ C \approx 1 \f$
 * then we will have to determine how many terms \f$m\f$ to add to \f$n\f$.  Setting \f$n = 19\f$, we solve
 *
 * \f[ \frac{ n+2 }{ 3^m } = 1 \quad \implies \quad  ln \left( n+2 \right) = ln(1) + m ln(3)
 *           \quad \implies \quad  m = \frac{ln(19)}{ln(3)} \approx 2.68
 * \f]
 *
 * Therefore, we take \f$n + m = 19 + 3 = 22\f$.
 *
 * \f{eqnarray*}{
 *      | R_n(y-x_c^{(c)}) |
 *               & \leq &
 *                         (22+2) h_3^2  \left( \frac{1}{3} \right)^{22-19} \\
 *               & = &
 *                         \frac{24}{3^3} h_3^2
 *             \approx     0.89 h_3^2
 * \f}
 *
 *
 * \f$\mathbf{S_2} - \mathbf{\mbox{Third Power}}\f$
 *
 * For the third power \f$S_2\f$
 *
 * \f{eqnarray*}{
 *     S_2(z+t) & = & \frac{1}{(z+t)^3} = \sum_{n=0}^{\infty} \frac{S_2^{(n)}(z)}{n!} t^n \\
 *              & = &
 *                   \frac{1}{0! z^3} t^0 + \frac{(-1)^1 3}{1! z^4} t^1 + \frac{(-1)^2 (4)(3)}{2! z^5} t^2 + \frac{(-1)^3 (5)(4)(3)}{3! z^6} t^3
 *                   + \hdots + \frac{(-1)^n (n+2)!}{2! n! z^{n+3}} t^n + R_n(z+t) \\
 *               & = &
 *                   S_2(z) + \frac{(-1)^1 3!}{2! 1!} t S_3(z) + \frac{(-1)^2 4!}{2! 2!} t^2 S_4(z) + \frac{(-1)^3 5!}{2! 3!} t^3 S_5(z)
 *                   + \hdots + \frac{(-1)^n (n+2)!}{2! n!} t^n S_{n+2} (z) + R_n(z+t) \\
 *               & = & \sum_{k=0}^{n} (-1)^{k} \begin{pmatrix} k+2 \\ 2 \end{pmatrix} t^k S_{k+2} (z) + R_n(z+t)
 * \f}
 *
 * where \f$l = k\f$ and \f$h_k = \left(\frac{1}{2}\right)^k\f$ implies
 *
 * \f{eqnarray*}{
 * | R_n(y-x_c^{(c)}) |
 *              & \leq &
 *                        \frac{ max \left( | S_2^{(n+1)}(y-x_c^{(c)}) | \right) }{(n+1)!} |x_c - x|^{n+1} \\
 *              & = &
 *                        max \left( \left| \frac{(-1)^{n+1} (n+3)!}{2! (y-x_c^{(c)})^{n+4}} \right| \right) \frac{1}{(n+1)!} |x_c^{(p)} - x_c^{(c)}|^{n+1} \\
 *              & \leq &
 *                        \frac{(n+3)(n+2)}{2! \left( \frac{3}{2} h_k \right)^{n+4}} \left( \frac{1}{2} h_k \right)^{n+1} \\
 *              & = &
 *                        \frac{(n+3)(n+2)}{2! h_k^3} \left( \frac{2}{3} \right)^{3} \left( \frac{1}{3} \right)^{n+1} \\
 *              & = &
 *                        \frac{(n+3)(n+2)}{2! h_k^3} \left( \left( \frac{1}{3} \right)^{k} \right)^3
 *                                            \left( \left( \frac{1}{3} \right)^{k} \right)^2
 *                                            \left( \frac{1}{3} \right)^{n+1-3k-2k} \\
 *              & \leq &
 *                        \left( \frac{2}{3} \right)^3 \frac{(n+3)(n+2)}{2!} h_k^2  \left( \frac{1}{3} \right)^{n-(5k-1)} \\
 *              & = &
 *                        \left( \frac{2}{3} \right)^3 \begin{pmatrix} n+3 \\ 2 \end{pmatrix} h_k^2  \left( \frac{1}{3} \right)^{n-(5k-1)}
 * \f}
 *
 * For \f$k = 5\f$
 *
 * \f{eqnarray*}{
 *      | R_n(y-x_c^{(c)}) |
 *               & \leq &
 *                        \left( \frac{2}{3} \right)^3 \frac{(n+3)(n+2)}{2!} h_3^2  \left( \frac{1}{3} \right)^{n-(5(5)-1)} \\
 *               & = &    \left( \frac{2}{3} \right)^3 \frac{(n+3)(n+2)}{2!} h_3^2  \left( \frac{1}{3} \right)^{n-24} \\
 * \f}
 *
 * If \f$n = 24\f$, then \f$(n+3)(n+2) = (27)(26) = 702\f$ and \f$2! = 2\f$.
 *
 * \f{eqnarray*}{
 *     | R_n(y-x_c^{(c)}) | & \leq &  \left( \frac{2}{3} \right)^3 \frac{702}{2} h_3^2 \left( \frac{1}{3} \right)^{24-24} \\
 *                  & = &  104 h_3^2
 * \f}
 *
 *
 * We note that \f$C = 104 > 1\f$.  Let's not be concerned with the number of terms \f$m\f$ we will need to add to \f$n\f$ just yet.
 * Our interest at this point will be the pattern as we step through the levels.
 *
 *
 *
 * \f$\mathbf{S_3} - \mathbf{\mbox{Fourth Power}}\f$
 *
 * For the fourth power
 *
 * \f{eqnarray*}{
 *     S_3(z+t) & = & \frac{1}{(z+t)^4} = \sum_{n=0}^{\infty} \frac{S_3^{(n)}(z)}{n!} t^n
 *                =   \sum_{k=0}^{n} (-1)^{k} \begin{pmatrix} k+3 \\ 3 \end{pmatrix} t^k S_{k+3} (z) + R_n(x) \\
 * \f}
 *
 *
 * For \f$l = k\f$ and \f$h_k = \left(\frac{1}{2}\right)^k\f$
 *
 * \f{eqnarray*}{
 * | R_n(y-x_c^{(c)}) |
 *              & \leq &
 *                        \frac{ max \left( | S_3^{(n+1)}(y-x_c^{(c)}) | \right) }{(n+1)!} |x_c^{(p)} - x_c^{(c)}|^{n+1}
 *                =
 *                        \left( \frac{2}{3} \right)^{4} \begin{pmatrix} n+4 \\ 3 \end{pmatrix} h_k^2  \left( \frac{1}{2} \right)^{n-(6k-1)}
 * \f}
 *
 * For \f$k = 5\f$
 *
 * \f{eqnarray*}{
 *      | R_n(y-x^{(c)}) |
 *               & \leq &
 *                        \left( \frac{2}{3} \right)^{4} \frac{(n+4)(n+3)(n+2)}{3!} h_3^2  \left( \frac{1}{2} \right)^{n-(6(5)-1)} \\
 *               & = &    \left( \frac{2}{3} \right)^{4} \frac{(n+4)(n+3)(n+2)}{3!} h_3^2  \left( \frac{1}{2} \right)^{n-29}
 * \f}
 *
 *
 * \f$\mathbf{S_4} - \mathbf{\mbox{Fifth Power}}\f$
 *
 * For the fifth power \f$S_4\f$
 *
 * \f{eqnarray*}{
 *     S_4(z+t) & = & \frac{1}{(z+t)^4} = \sum_{n=0}^{\infty} \frac{S_4^{(n)}(z)}{n!} t^n \\
 *              & = & \sum_{k=0}^{n} (-1)^{k} \begin{pmatrix} k+4 \\ 4 \end{pmatrix} t^k S_{k+4} (z) + R_n(x)
 * \f}
 *
 * \f{eqnarray*}{
 * | R_n(y-x^{(c)}) |
 *              & \leq &
 *                       \left( \frac{2}{3} \right)^5 \begin{pmatrix} n+5 \\ 4 \end{pmatrix} h_k^2  \left( \frac{1}{2} \right)^{n-(7k-1)}
 * \f}
 *
 * For \f$k = 5\f$
 *
 * \f{eqnarray*}{
 *      | R_n(y-x^{(c)}) |
 *               & \leq &
 *                        \left( \frac{2}{3} \right)^5 \frac{(n+5)(n+4)(n+3)(n+2)}{4!} h_3^2  \left( \frac{1}{3} \right)^{n-(7(5)-1)} \\
 *               & = &    \left( \frac{2}{3} \right)^5 \frac{(n+5)(n+4)(n+3)(n+2)}{4!} h_3^2  \left( \frac{1}{3} \right)^{n-34} \\
 * \f}
 *
 *
 * \f$\mathbf{S_m} - \mathbf{\mbox{mth Power}}\f$
 *
 * For the mth power \f$S_m\f$
 *
 * \f{eqnarray*}{
 *     S_m(z+t) & = & \frac{1}{(z+t)^m} = \sum_{n=0}^{\infty} \frac{S_m^{(n)}(z)}{n!} t^n
 *                = \sum_{k=0}^{n} (-1)^{k} \begin{pmatrix} k+m \\ m \end{pmatrix} t^k S_{k+m} (z) + R_n(x)
 * \f}
 *
 * \f{eqnarray*}{
 * | R_n(y-x^{(c)}) |
 *              & \leq &
 *                      \left( \frac{2}{3} \right)^{m+1} \begin{pmatrix} n+m+1 \\ m \end{pmatrix} h_k^2  \left( \frac{1}{2} \right)^{n-((m+3)k-1)}
 * \f}
 *
 * From the first part of the upward pass we know that at the lowest level we will need to approximate the far-field expansion out to the fourteenth term.
 * Setting \f$m = 14\f$ we have
 *
 * \f{eqnarray*}{
 * | R_n(y-x^{(c)}) |
 *              & \leq &
 *                      \left( \frac{2}{3} \right)^{m+1} \begin{pmatrix} n+m+1 \\ m \end{pmatrix} h_k^2  \left( \frac{1}{2} \right)^{n-(17k-1)}
 *
 * \f}
 *
 *
 * For \f$k = 5\f$
 *
 * \f{eqnarray*}{
 *      | R_n(y-x^{(c)}) |
 *               & \leq &
 *                       \left( \frac{2}{3} \right)^{m+1} \begin{pmatrix} n+15 \\ 14 \end{pmatrix} h_3^2  \left( \frac{1}{2} \right)^{n-(17(5)-1)} \\
 *               & = &   \left( \frac{2}{3} \right)^{m+1} \begin{pmatrix} n+15 \\ 14 \end{pmatrix} h_3^2  \left( \frac{1}{2} \right)^{n-84}
 * \f}
 *
 * we can see that we will need at least 84 terms.
 *
 *
 *
 * \f$\mathbf{Level 4} - \mathbf{S_m - \mbox{ mth Power}}\f$
 *
 * We continue up at level \f$l = 4\f$ to translate the series at the child cells at to level \f$l = 4\f$ to the
 * parent cells at level \f$l = 3\f$.
 *
 *  |                                                               .                                                               |            level l = 0
 *  0                                                                                                                               1
 *
 *  |                               .                               |                               .                               |            level l = 1
 * 0/2                                                             1/2                                                             2/2
 *
 *  |      8h       .       8h      |               .               |               .               |               .               |            level l = 2
 * 0/4                             1/4                             2/4                             3/4                             4/4
 *
 *  |   2h  .   2h  |   2h  .  2h   |       .       |       .       |       .       |       .       |       .       |       .       |            level l = 3
 * 0/8             1/8             2/8             3/8             4/8             5/8             6/8             7/8             8/8
 *
 *  | h . h | h . h | h . h |   .   |   .   |   .   |   .   |   .   |   .   |   .   |   .   |   .   |   .   |   .   |   .   |   .   |            level l = 4
 * 0/16    1/16    2/16    3/16    4/16    5/16    6/16    7/16    8/16    9/16   10/16   11/16   12/16   13/16   14/16   15/16   16/16
 *
 *  | . | . | . | . | . | . | . | . | . | . | . | . | . | . | . | . | . | . | . | . | . | . | . | . | . | . | . | . | . | . | . | . |            level l = 5
 * 0/32    2/32    4/32    6/32    8/32   10/32   12/32   14/32   16/32   18/32   20/32   22/32   24/32   26/32   28/32   30/32   32/32
 *
 *  |-h-|
 *
 * Let \f$z + t = y - x_c^{(c)}\f$, \f$z = y - x_c^{(p)}\f$, and \f$t = z+t-z = x_c^{(p)} - x_c^{(c)}\f$.
 * For the child cells at level \f$l = 4\f$, \f$|y-x_c^{(c)}| > 3h\f$ and \f$|x_c^{(p)} - x_c^{(c)}| \leq h\f$.
 *
 *
 *
 * For the second power \f$S_1\f$ of the far-field expansion
 *
 *
 * \f{eqnarray*}{
 *     S_0(z+t) & = & \frac{1}{z+t} = \sum_{n=0}^{\infty} \frac{\phi^{(n)}(z)}{n!} t^n \\
 *               & = & \sum_{k=0}^{n} (-1)^{k} t^k S_k(z) + R_n(z+t)
 * \f}
 *
 * We have \f$\phi^{(n)}(z) = \frac{(-1)^n n!}{z^{n+1}}\f$.  For \f$l = k\f$ and \f$h = \left( \frac{1}{2} \right)^k\f$ we have
 *
 * \f{eqnarray*}{
 *    |R_n(z+t)| & = & |R_n(y-x_c^{(c)})|
 *                       \leq
 *                                 \frac{ max \left( | \phi^{(n+1)}(y-x_c^{(c)}) | \right) }{(n+1)!} |x_c^{(p)} - x_c^{(c)}|^{n+1} \\
 *                       & = &  \max \left( \left| \frac{(-1)^{n+1} (n+1)!}{(y-x_c^{(c)})^{n+2}} \right| \right)  \frac{1}{(n+1)!}
 *                                                        |x_c^{(p)} - x_c^{(c)}|^{n+1} \\
 *                       & \leq & \frac{1}{(3h)^{n+2}} \left( h \right)^{n+1} \\
 *                       & = & \frac{1}{3h} \left( \frac{1}{3} \right)^{n+1} \\
 *                       & = &
 *                             \frac{1}{h} \left( \frac{1}{3} \right)^{k} \left( \left( \frac{1}{3} \right)^{k} \right)^2 \left( \frac{1}{3} \right)^{n+2-k-2k} \\
 *                       & \leq &
 *                             h^2 \left( \frac{1}{3} \right)^{n-(3k-2)}
 * \f}
 *
 *
 * For the second power \f$S_1\f$ of the far-field expansion
 *
 * \f{eqnarray*}{
 *     S_1(z+t) & = & \frac{1}{(z+t)^2} = \sum_{n=0}^{\infty} \frac{S_1^{(n)}(z)}{n!} t^n \\
 *               & = & \sum_{k=0}^{n} (-1)^{k} \begin{pmatrix} k+1 \\ 1 \end{pmatrix} t^k S_{k+1} (z) + R_n(z+t)
 * \f}
 *
 * where \f$l = k\f$ and \f$h = \left(\frac{1}{2}\right)^k\f$ implies
 *
 * \f{eqnarray*}{
 * | R_n(y-x_c^{(c)}) |
 *              & \leq &
 *                        \frac{ max \left( | S_1^{(n+1)}(y-x_c^{(c)}) | \right) }{(n+1)!} |x_c^{(p)} - x_c^{(c)}|^{n+1} \\
 *              & = &
 *                        max \left( \left| \frac{(-1)^{n+1} (n+2)!}{1!(y-x_c^{(c)})^{n+3}} \right| \right) \frac{1}{(n+1)!} |x_c^{(p)} - x_c^{(c)}|^{n+1} \\
 *              & \leq &
 *                        \frac{(n+2)}{1! \left( 3 h \right)^{n+3}} \left(  h \right)^{n+1} \\
 *              & = &
 *                        \frac{(n+2)}{1!} \left(\frac{1}{3}\right)^2 \left(\frac{1}{h}\right)^2 \left( \frac{1}{3} \right)^{n+1} \\
 *              & \leq &
 *                        \frac{(n+2)}{1!} \left(\frac{1}{h}\right)^2 \left( \frac{1}{3} \right)^{n+3} \\
 *              & = &
 *                        \frac{(n+2)}{1!} \left(\frac{1}{h}\right)^2
 *                            \left( \left( \frac{1}{3} \right)^{k} \right)^2 \left( \left( \frac{1}{3} \right)^{k} \right)^2 \left( \frac{1}{3} \right)^{n+3-2k-2k} \\
 *              & = &
 *                        \begin{pmatrix} n+2 \\ 1 \end{pmatrix}  h^2 \left( \frac{1}{3} \right)^{n-(4k-3)} \\
 * \f}
 *
 *
 * For the third power \f$S_2\f$
 *
 * \f{eqnarray*}{
 *     S_2(z+t) & = & \frac{1}{(z+t)^3} = \sum_{n=0}^{\infty} \frac{S_2^{(n)}(z)}{n!} t^n \\
 *               & = & \sum_{k=0}^{n} (-1)^{k} \begin{pmatrix} k+2 \\ 2 \end{pmatrix} t^k S_{k+2} (z) + R_n(z+t)
 * \f}
 *
 * where \f$l = k\f$ and \f$h = h_k = \left(\frac{1}{2}\right)^k\f$ implies
 *
 * \f{eqnarray*}{
 * | R_n(y-x_c^{(c)}) |
 *              & \leq &
 *                        \frac{ max \left( | S_2^{(n+1)}(y-x_c^{(c)}) | \right) }{(n+1)!} |x_c - x|^{n+1} \\
 *              & = &
 *                        max \left( \left| \frac{(-1)^{n+1} (n+3)!}{2! (y-x_c^{(c)})^{n+4}} \right| \right) \frac{1}{(n+1)!} |x_c^{(p)} - x_c^{(c)}|^{n+1} \\
 *              & \leq &
 *                        \frac{(n+3)(n+2)}{2! \left( 3 h \right)^{n+4}} \left( h \right)^{n+1} \\
 *              & = &
 *                        \frac{(n+3)(n+2)}{2! h^3} \left( \frac{1}{3} \right)^{3} \left( \frac{1}{3} \right)^{n+1} \\
 *              & = &
 *                        \frac{(n+3)(n+2)}{2! h^3} \left( \left( \frac{1}{3} \right)^{k} \right)^3
 *                                            \left( \left( \frac{1}{3} \right)^{k} \right)^2
 *                                            \left( \frac{1}{3} \right)^{n+4-3k-2k} \\
 *              & \leq &
 *                         \begin{pmatrix} n+3 \\ 2 \end{pmatrix} h^2  \left( \frac{1}{3} \right)^{n-(5k-4)}
 * \f}
 *
 *
 * For the mth power \f$S_m\f$
 *
 * \f{eqnarray*}{
 *     S_m(z+t) & = & \frac{1}{(z+t)^m} = \sum_{n=0}^{\infty} \frac{S_m^{(n)}(z)}{n!} t^n
 *                = \sum_{k=0}^{n} (-1)^{k} \begin{pmatrix} k+m \\ m \end{pmatrix} t^k S_{k+m} (z) + R_n(x)
 * \f}
 *
 * \f{eqnarray*}{
 * | R_n(y-x^{(c)}) |
 *              & \leq &
 *                       \begin{pmatrix} n+m+1 \\ m \end{pmatrix} h^2  \left( \frac{1}{3} \right)^{n-((m+3)k-1)}
 * \f}
 *
 * From the first part of the upward pass we know that at the lowest level we will need to approximate the far-field expansion out to the fourteenth term.
 * Setting \f$m = 14\f$ we have
 *
 * \f{eqnarray*}{
 * | R_n(y-x^{(c)}) |
 *              & \leq &
 *                      \left( \frac{2}{3} \right)^{m+1} \begin{pmatrix} n+m+1 \\ m \end{pmatrix} h_k^2  \left( \frac{1}{2} \right)^{n-(17k-1)}
 *
 * \f}
 *
 *
 * For \f$k = 5\f$
 *
 * \f{eqnarray*}{
 *      | R_n(y-x^{(c)}) |
 *               & \leq &
 *                       \left( \frac{2}{3} \right)^{m+1} \begin{pmatrix} n+15 \\ 14 \end{pmatrix} h_3^2  \left( \frac{1}{2} \right)^{n-(17(5)-1)} \\
 *               & = &   \left( \frac{2}{3} \right)^{m+1} \begin{pmatrix} n+15 \\ 14 \end{pmatrix} h_3^2  \left( \frac{1}{2} \right)^{n-84}
 * \f}
 *
 * we can see that we will need at least 84 terms.
 *
 *
 *
 *
 *
 *
 * \f$\mathbf{S|S} - \mathbf{\mbox{Transformation Matrix}}\f$
 *
 * Looking back at the approximations, \f$n=8\f$ for \f$S_0\f$, \f$n=13\f$ for \f$S_1\f$, \f$n=18\f$ for \f$S_2\f$, \f$n=23\f$ for \f$S_3\f$,
 * \f$n=28\f$ for \f$S_4\f$, \f$n=33\f$ for \f$S_5\f$, \f$n=38\f$ for \f$S_6\f$, \f$n=43\f$ for \f$S_7\f$, and \f$n=47\f$ for \f$S_8\f$.
 * The increments in \f$n\f$ are by \f$5\f$ as we increment through the powers except for the last power (if we incremented by 5 for \f$S_8\f$
 * we would have \f$n=48\f$ and a constant closer to 1).
 *
 * Let's go ahead and use \f$n=48\f$ for \f$S_8\f$.  We will need the resulting far-field approximation centered at the parent cell at level
 * \f$l=2\f$ to have 49 terms (need to count the zero term for \f$n=48\f$).  Our next step in the Fast Multipole Method is the translation of
 * the 49-term far-field series at level \f$l=2\f$ to near-field series in the first part of the Downward Pass.
 *
 *
 *
 *
 *
 *
 *
 * @subsubsection downward_pass_part_one Downward Pass Part One
 *
 * Once all the source terms at the lowest FMM level have been translated up to level \f$l = 2\f$,
 * we are ready to start the first part of the downward pass.
 * We note that each cell at level \f$l = 2\f$ has a far-field expansion \f$\sum_{n=0}^{\infty} a_n(x-x_c^{(p)}) S_n(y-x_c^{(p)})\f$
 * for all the source particles that the cell contains.
 *
 * Level \f$l=2\f$ is unique since there are no cells in the domain that are beyond a cell's interaction list.
 * The first step of the downward pass loops through all the cells at \f$l = 2\f$ and performs an S|R translation of the far-field expansions
 *  that are associated with the cell's interaction list to near-field expansions
 * .
 *
 * Fix one of the cell's of the loop.  The far-field expansion centered at a cell in the interaction list is translated to from
 * the interaction list cell center \f$x_c\f$ (cell contains sources) to the loop cell's center \f$y_c\f$ (cell contains targets).
 * During the S|R translation the series changes from a far-field series \f$\sum_{n=0}^{\infty} a_n (x-x_c) S_n(y-x_c)\f$ to a
 * near-field series \f$\sum_{n=0}^{\infty} b_n (x-y_c) R_n(y-y_c)\f$ where the initial center is \f$z+t = y-x_c^{(p)}\f$, and the
 * final center is \f$z = y-y_c\f$, and \f$t = (z+t)-z = y_c-x_c\f$.
 *
 * To construct the translation we use Taylor series of the form
 *
 * \f{eqnarray*}{
 *     \phi(z+t) & = & \frac{1}{z+t} = \sum_{n=0}^{\infty} \frac{\phi^{(n)} (t)}{n!} z^n \\
 *               & = & \frac{1}{0! t} z^0 + \frac{(-1)^1 1!}{1! t^2} z^1 + \frac{(-1)^2 2!}{2! t^3} z^2
 *                         + \frac{(-1)^3 3!}{3! t^4} z^3 + \hdots + \frac{(-1)^n n!}{n! t^{n+1}} z^n + R_n(z+t) \\
 *               & = & \sum_{k=0}^{n} \frac{(-1)^k k!}{k! t^{k+1}} z^k + R_n(z+t) \\
 *     \phi(y-x_c^{(p)}) & = & \sum_{k=0}^{n} \frac{(-1)^k}{(y_c^{(p)}-x_c^{(p)})^{k+1}} (y-y_c^{(p)})^n + R_n(y-x_c^{(p)}) \\
 *                       & = & \sum_{k=0}^{n} b_k(y_c^{(p)}-x_c^{(p)}) V_k(y-y_c^{(p)}) + R_n(y-x_c^{(p)})
 * \f}
 *
 * where \f$z + t = y - x_c^{(p)}\f$, \f$z = y - y_c^{(p)}\f$, \f$t = y_c^{(p)} - x_c^{(p)}\f$, and
 * \f$V_n(z) = z^n = (y-y_c^{(p)})^n\f$.
 *
 * Recall
 *
 * \f[ f(x) = \sum_{n=0}^{\infty} \frac{f^{(n)}(a)}{n!} (x-a)^n
 *          = \sum_{k=0}^{n} \frac{f^{(k)}(a)}{k!} (x-a)^k + \sum_{k=n+1}^{\infty} \frac{f^{(k)}(a)}{k!} (x-a)^k
 *          = T_n(x) + R_n(x)
 * \f]
 *
 * and when \f$|f^{(n+1)}(x)| \leq M\f$ for \f$|x-a| \leq d\f$
 *
 * \f[ |R_n(x)| \leq \frac{ \max(|f^{(n+1)}(x)|)}{(n+1)!} |x-a|^{n+1} \quad \mbox{for } a-d \leq x \leq a+d. \f]
 *
 * Here \f$f = \phi\f$, \f$a = t = y_c^{(p)}-x_c^{(p)}\f$, \f$x = z+t = y - x_c^{(p)}\f$, and \f$x-a =y-y_c^{(p)}\f$.  Therefore
 *
 * \f{eqnarray*}{
 *   |R_n(y - x_c^{(p)})| & \leq & \max \left( \left| \frac{(-1)^{n+1}(n+1)!}{(y_c^{(p)}-x_c^{(p)})^{n+2}} \right| \right) \frac{1}{(n+1)!} |y-y_c^{(p)}|^{n+1} \\
 *                        & \leq & \frac{1}{(2h)^{n+2}} \left( \frac{1}{2}h \right)^{n+1} \\
 *                        & =    &  \frac{1}{2} \frac{1}{h} \left( \frac{1}{2} \right)^{2n+2}
 * \f}
 *
 * where we use \f$|y_c^{(p)}-x_c^{(p)}| > 2h\f$ and \f$|y-y_c^{(p)}| < \frac{1}{2}h\f$.
 *
 * If we are at level \f$l = k\f$ and \f$h = h_k = \left( \frac{1}{2} \right)^{k}\f$ (\f$h = h_0 = \left( \frac{1}{2} \right)^0 = 1\f$ at level \f$l=0\f$)
 *
 * \f{eqnarray*}{
 *    |R_n(y - x_c^{(p)})| & \leq & \frac{1}{2} \frac{1}{h_k} \left( \frac{1}{2} \right)^{2n+2} \\
 *                         & =    & \frac{1}{2} \frac{1}{h_k} \left( \frac{1}{2} \right)^k \left( \left( \frac{1}{2} \right)^k \right)^2 \left( \frac{1}{2} \right)^{2n+2-k-2k} \\
 *                         & \leq & \frac{1}{2} h_k^2 \left( \frac{1}{4} \right)^{2n-(3k-2)}
 * \f}
 *
 * We are ready to start constructing the \f$S|R\f$ transformation matrix.  The \f$S|R\f$ transformation translates and transforms the powers \f$S_n(z+t)\f$
 * of the far-field expansion using Taylor series like the one above.  Since the \f$S|S\f$ transformation
 * matrix returns far-field expansions having 49 terms, we know that we will be approximating all 49 powers: \f$S_0, S_1, S_2, \ldots, S_{48}\f$.
 *
 *
 * \f$\mathbf{S_0} - \mathbf{\mbox{First Power}}\f$
 *
 * Since \f$S_0 = \phi\f$, the Taylor series for the first power \f$S_0\f$ is
 *
 * \f{eqnarray*}{
 *    S_0(z+t) & = & \frac{1}{z+t}
 *               =   \frac{1}{0! t} z^0 + \frac{(-1)^1 1!}{1! t^2} z^1 + \frac{(-1)^2 2!}{2! t^3} z^2
 *                         + \frac{(-1)^3 3!}{3! t^4} z^3 + \hdots + \frac{(-1)^n n!}{n! t^{n+1}} z^n + R_n(z+t) \\
 *             & = & \frac{1}{t} + \frac{(-1)}{t^2} z + \frac{(-1)^2}{t^3} z^2 + \frac{(-1)^3}{t^4} z^3
 *                        + \hdots + \frac{(-1)^n}{t^{n+1}} z^n + R_n(z+t) \\
 *             & = & \frac{1}{t} V_0(z) + \frac{(-1)}{t^2} V_1(z) + \frac{(-1)^2}{t^3} V_2(z) + \frac{(-1)^3}{t^4} V_3(z)
 *                        \hdots + \frac{(-1)^n}{t^{n+1}} V_n(z) + R_n(z+t) \\
 *             & = & \sum_{k=0}^{n} \frac{(-1)^k}{(y_c^{(p)}-x_c^{(p)})^{k+1}} (y-y_c^{(p)})^k + R_n(y-x_c^{(p)}) \\
 * \f}
 *
 * and
 *
 * \f{eqnarray*}{
 * | R_n(y-x_c^{(p)}) |
 *                         & \leq & \frac{1}{2} h_k^2 \left( \frac{1}{2} \right)^{2n-(3k-2)}
 * \f}
 *
 * For \f$l = 2\f$ and \f$h = \left( \frac{1}{2} \right)^2\f$ we have
 *
 * \f{eqnarray*}{
 * | R_n(y-x_c^{(p)}) |
 *                         & \leq & \frac{1}{2} h_2^2 \left( \frac{1}{2} \right)^{2n-(3(2)-2)} \\
 *                         & =    & \frac{1}{2} h_2^2 \left( \frac{1}{2} \right)^{2n-4}
 * \f}
 *
 * If we choose \f$n = 2\f$ we have \f$C = \frac{1}{2}\f$ and
 *
 * \f{eqnarray*}{
 * | R_n(y-x_c^{(p)}) |
 *                         & \leq & \frac{1}{2} h_2^2  = C h_2^2.
 * \f}
 *
 *
 * \f$\mathbf{S_1} - \mathbf{\mbox{Second Power}}\f$
 *
 * For \f$S_1\f$ the Taylor series is
 *
 * \f{eqnarray*}{
 *    S_1(z+t) & = & \frac{1}{(z+t)^2}
 *               =   \frac{1}{0! t^2} z^0 + \frac{(-1)^1 2!}{1!1! t^3} z^1 + \frac{(-1)^2 3!}{1!2! t^4} z^2
 *                         + \frac{(-1)^3 4!}{1! 3! t^5} z^3 + \hdots + \frac{(-1)^n (n+1)!}{1! n! t^{n+2}} z^n + R_n(z+t) \\
 *             & = & \frac{1}{t^2} + (-1)^1 \begin{pmatrix} 2 \\ 1 \end{pmatrix} \frac{1}{t^3} z
 *                                 + (-1)^2 \begin{pmatrix} 3 \\ 2 \end{pmatrix} \frac{1}{t^4} z^2
 *                                 + (-1)^3 \begin{pmatrix} 4 \\ 1 \end{pmatrix} \frac{1}{t^5} z^3
 *                        + \hdots + (-1)^n \begin{pmatrix} (n+1) \\ 1 \end{pmatrix} \frac{1}{t^{n+2}} z^n
 *                        + R_n(z+t) \\
 *             & = & (-1)^0 \begin{pmatrix} 1 \\ 1 \end{pmatrix} \frac{1}{t^2} V_0(z)
 *                   + (-1)^1 \begin{pmatrix} 2 \\ 1 \end{pmatrix} \frac{1}{t^3} V_1(z)
 *                   + (-1)^2 \begin{pmatrix} 3 \\ 1 \end{pmatrix} \frac{1}{t^4} V_2(z)
 *                   + (-1)^3 \begin{pmatrix} 4 \\ 1 \end{pmatrix} \frac{1}{t^5} V_3(z)
 *                   + \hdots
 *                   + (-1)^n \begin{pmatrix} n+1 \\ 1 \end{pmatrix} \frac{1}{t^{n+2}} V_n(z)
 *                   + R_n(z+t) \\
 *             & = & \sum_{k=0}^{n} (-1)^k \begin{pmatrix} k+1 \\ 1 \end{pmatrix} \frac{1}{(y_c^{(p)}-x_c^{(p)})^{k+2}} (y-y_c^{(p)})^k
 *                   + R_n(y-x_c^{(p)}) \\
 * \f}
 *
 * and if \f$l = k\f$ and \f$h = \left( \frac{1}{2} \right)^{k}\f$
 *
 * \f{eqnarray*}{
 * | R_n(y-x_c^{(p)}) |
 *                         & \leq &  \max_{y_c^{(p)}} \left( \left| (-1)^{n+1} \begin{pmatrix} n+2 \\ 1 \end{pmatrix} \frac{1}{(y_c^{(p)}-x_c^{(p)})^{n+3}} \right| \right)
 *                                           (y-y_c^{(p)})^{n+1} \\
 *                         & \leq &  \begin{pmatrix} n+2 \\ 1 \end{pmatrix} \frac{1}{(2h)^{n+3}}
 *                                           \left( \frac{1}{2} h \right)^{n+1}                                         \\
 *                         & =    &  \left( \frac{1}{2} \right)^2 \begin{pmatrix} n+2 \\ 1 \end{pmatrix}  \frac{1}{h^2}
 *                                           \left( \frac{1}{2} \right)^{2n+2}                                          \\
 *                         & =    &  \left( \frac{1}{2} \right)^2 \begin{pmatrix} n+2 \\ 1 \end{pmatrix}  \frac{1}{h^2}
 *                                           \left( \left( \frac{1}{2} \right)^{k} \right)^2 \left( \left( \frac{1}{2} \right)^{k} \right)^2
 *                                           \left( \frac{1}{2} h \right)^{2n+2-2k-2k}                                  \\
 *                         & \leq &  \left( \frac{1}{2} \right)^2 \begin{pmatrix} n+2 \\ 1 \end{pmatrix}  h^2 \left( \frac{1}{2} \right)^{2n-(4k-2)}
 * \f}
 *
 * For \f$l = 2\f$ and \f$h = \left( \frac{1}{2} \right)^2\f$ we have
 *
 * \f{eqnarray*}{
 * | R_n(y-x_c^{(p)}) |
 *                         & \leq & \left( \frac{1}{2} \right)^2 \begin{pmatrix} n+2 \\ 1 \end{pmatrix} h_2^2 \left( \frac{1}{2} \right)^{2n-(4(2)-2)} \\
 *                         & =    & \left( \frac{1}{2} \right)^2 \begin{pmatrix} n+2 \\ 1 \end{pmatrix} h_2^2 \left( \frac{1}{2} \right)^{2n-6}
 * \f}
 *
 * If we choose \f$n = 3\f$ we have \f$C = \left( \frac{1}{2} \right)^2 \begin{pmatrix} 5 \\ 1 \end{pmatrix} = \frac{5}{4} \f$
 *
 * \f{eqnarray*}{
 * | R_n(y-x_c^{(p)}) |
 *                         & \leq & \left( \frac{1}{2} \right)^2 \begin{pmatrix} 3+2 \\ 1 \end{pmatrix} h_2^2  = C h_2^2.
 * \f}
 *
 *
 *
 * \f$\mathbf{S_2} - \mathbf{\mbox{Third Power}}\f$
 *
 * For \f$S_2\f$ the Taylor series is
 *
 * \f{eqnarray*}{
 *    S_2(z+t) & = & \frac{1}{(z+t)^3}
 *               =   \frac{1}{0! t^3} z^0 + \frac{(-1)^1 3!}{2!1! t^4} z^1 + \frac{(-1)^2 4!}{2!2! t^5} z^2
 *                         + \frac{(-1)^3 5!}{2! 3! t^6} z^3 + \hdots + \frac{(-1)^n (n+2)!}{2! n! t^{n+3}} z^n + R_n(z+t) \\
 *             & = & (-1)^n \begin{pmatrix} (n+2) \\ 2 \end{pmatrix} \frac{1}{t^{n+3}} z^n
 *                        + R_n(z+t) \\
 *             & = & \sum_{k=0}^{n} (-1)^k \begin{pmatrix} k+2 \\ 2 \end{pmatrix} \frac{1}{(y_c^{(p)}-x_c^{(p)})^{k+3}} (y-y_c^{(p)})^k
 *                   + R_n(y-x_c^{(p)}) \\
 * \f}
 *
 * and if \f$l = k\f$ and \f$h = \left( \frac{1}{2} \right)^{k}\f$
 *
 * \f{eqnarray*}{
 * | R_n(y-x_c^{(p)}) |
 *                         & \leq &  \max_{y_c^{(p)}} \left( \left| (-1)^{n+1} \begin{pmatrix} n+3 \\ 2 \end{pmatrix} \frac{1}{(y_c^{(p)}-x_c^{(p)})^{n+4}} \right| \right)
 *                                           (y-y_c^{(p)})^{n+1} \quad \quad |y-y_c^{(p)}| \leq h \\
 *                         & \leq &  \begin{pmatrix} n+3 \\ 2 \end{pmatrix} \frac{1}{(2h)^{n+4}}
 *                                           \left( \frac{1}{2} h \right)^{n+1}                                         \\
 *                         & =    &  \left( \frac{1}{2} \right)^3 \begin{pmatrix} n+3 \\ 2 \end{pmatrix}  \frac{1}{h^3}
 *                                           \left( \frac{1}{2} \right)^{2n+2}                                          \\
 *                         & =    &  \left( \frac{1}{2} \right)^3 \begin{pmatrix} n+3 \\ 2 \end{pmatrix}  \frac{1}{h^3}
 *                                           \left( \left( \frac{1}{2} \right)^{k} \right)^3 \left( \left( \frac{1}{2} \right)^{k} \right)^2
 *                                           \left( \frac{1}{2} h \right)^{2n+2-3k-2k}                                  \\
 *                         & \leq &  \left( \frac{1}{2} \right)^3 \begin{pmatrix} n+3 \\ 2 \end{pmatrix}  h^2 \left( \frac{1}{2} \right)^{2n-(5k-2)}
 * \f}
 *
 * For \f$l = 2\f$ and \f$h = \left( \frac{1}{2} \right)^2\f$ we have
 *
 * \f{eqnarray*}{
 * | R_n(y-x_c^{(p)}) |
 *                         & \leq & \left( \frac{1}{2} \right)^3 \begin{pmatrix} n+3 \\ 2 \end{pmatrix} h_2^2 \left( \frac{1}{2} \right)^{2n-(5(2)-2)} \\
 *                         & =    & \left( \frac{1}{2} \right)^3 \begin{pmatrix} n+3 \\ 2 \end{pmatrix} h_2^2 \left( \frac{1}{2} \right)^{2n-8}
 * \f}
 *
 * If we choose \f$n = 4\f$ we have \f$C = \left( \frac{1}{2} \right)^3 \begin{pmatrix} 4+3 \\ 2 \end{pmatrix} = \frac{21}{8} = 2.625 \f$
 *
 * \f{eqnarray*}{
 * | R_n(y-x_c^{(p)}) |
 *                         & \leq & \left( \frac{1}{2} \right)^3 \begin{pmatrix} 4+3 \\ 2 \end{pmatrix} h_2^2  = C h_2^2.
 * \f}
 *
 *
 * \f$\mathbf{S_3} - \mathbf{\mbox{Fourth Power}}\f$
 *
 * For \f$S_3\f$ the Taylor series is
 *
 * \f{eqnarray*}{
 *    S_3(z+t) & = & \frac{1}{(z+t)^4}
 *               =   \sum_{k=0}^{n} (-1)^k \begin{pmatrix} k+3 \\ 3 \end{pmatrix} \frac{1}{(y_c^{(p)}-x_c^{(p)})^{k+4}} (y-y_c^{(p)})^k
 *                   + R_n(y-x_c^{(p)}) \\
 * \f}
 *
 * and if \f$l = k\f$ and \f$h = \left( \frac{1}{2} \right)^{k}\f$
 *
 * \f{eqnarray*}{
 * | R_n(y-x_c^{(p)}) |
 *                         & \leq &  \max_{y_c^{(p)}} \left( \left| (-1)^{n+1} \begin{pmatrix} n+4 \\ 3 \end{pmatrix} \frac{1}{(y_c^{(p)}-x_c^{(p)})^{n+5}} \right| \right)
 *                                           (y-y_c^{(p)})^{n+1} \quad \quad |y-y_c^{(p)}| \leq h \\
 *                         & \leq &  \left( \frac{1}{2} \right)^4 \begin{pmatrix} n+4 \\ 3 \end{pmatrix}  h^2 \left( \frac{1}{2} \right)^{2n-(6k-2)}
 * \f}
 *
 * For \f$l = 2\f$ and \f$h = \left( \frac{1}{2} \right)^2\f$ we have
 *
 * \f{eqnarray*}{
 * | R_n(y-x_c^{(p)}) |
 *                         & =    & \left( \frac{1}{2} \right)^4 \begin{pmatrix} n+4 \\ 3 \end{pmatrix} h_2^2 \left( \frac{1}{2} \right)^{2n-10}
 * \f}
 *
 * If we choose \f$n = 5\f$ we have \f$C = \left( \frac{1}{2} \right)^4 \begin{pmatrix} 5+4 \\ 3 \end{pmatrix} = \frac{84}{16} = 5.25 \f$
 *
 * \f{eqnarray*}{
 * | R_n(y-x_c^{(p)}) |
 *                         & \leq & \left( \frac{1}{2} \right)^4 \begin{pmatrix} 5+4 \\ 3 \end{pmatrix} h_2^2  = C h_2^2.
 * \f}
 *
 *
 * \f$\mathbf{S_4} - \mathbf{\mbox{Fifth Power}}\f$
 *
 * For \f$S_4\f$ the Taylor series is
 *
 * \f{eqnarray*}{
 *    S_4(z+t) & = & \frac{1}{(z+t)^5}
 *               =   \sum_{k=0}^{n} (-1)^k \begin{pmatrix} k+4 \\ 4 \end{pmatrix} \frac{1}{(y_c^{(p)}-x_c^{(p)})^{k+5}} (y-y_c^{(p)})^k
 *                   + R_n(y-x_c^{(p)}) \\
 * \f}
 *
 * and if \f$l = k\f$ and \f$h = \left( \frac{1}{2} \right)^{k}\f$
 *
 * \f{eqnarray*}{
 * | R_n(y-x_c^{(p)}) |
 *                         & \leq &  \max_{y_c^{(p)}} \left( \left| (-1)^{n+1} \begin{pmatrix} n+5 \\ 4 \end{pmatrix} \frac{1}{(y_c^{(p)}-x_c^{(p)})^{n+6}} \right| \right)
 *                                           (y-y_c^{(p)})^{n+1} \quad \quad |y-y_c^{(p)}| \leq h \\
 *                         & \leq &  \left( \frac{1}{2} \right)^5 \begin{pmatrix} n+5 \\ 4 \end{pmatrix}  h^2 \left( \frac{1}{2} \right)^{2n-(7k-2)}
 * \f}
 *
 * For \f$l = 2\f$ and \f$h = \left( \frac{1}{2} \right)^2\f$ we have
 *
 * \f{eqnarray*}{
 * | R_n(y-x_c^{(p)}) |
 *                         & =    & \left( \frac{1}{2} \right)^5 \begin{pmatrix} n+5 \\ 4 \end{pmatrix} h_2^2 \left( \frac{1}{2} \right)^{2n-12}
 * \f}
 *
 * If we choose \f$n = 6\f$ we have \f$C = \left( \frac{1}{2} \right)^5 \begin{pmatrix} 6+5 \\ 4 \end{pmatrix} = \frac{330}{32} = 10.3125 \f$
 *
 * \f{eqnarray*}{
 * | R_n(y-x_c^{(p)}) |
 *                         & \leq & \left( \frac{1}{2} \right)^5 \begin{pmatrix} 6+5 \\ 4 \end{pmatrix} h_2^2  = C h_2^2.
 * \f}
 *
 * We wish to keep \f$C\f$ constant at \f$C \approx 1\f$.  We increase \f$n\f$ by 2 so that
 *
 * \f{eqnarray*}{
 * | R_n(y-x_c^{(p)}) |
 *                         & \leq & \left( \frac{1}{2} \right)^5 \begin{pmatrix} n+5 \\ 4 \end{pmatrix} h_2^2 \left( \frac{1}{2} \right)^{2n-12} \\
 *                         & =    & \left( \frac{1}{2} \right)^5 \begin{pmatrix} 8+4 \\ 4 \end{pmatrix} h_2^2 \left( \frac{1}{2} \right)^{2(8)-12} \\
 *                         & \approx & 0.97
 * \f}
 *
 * So we set \f$n = 6 + 2 = 8\f$
 *
 *
 * \f$\mathbf{S_5} - \mathbf{\mbox{Sixth Power}}\f$
 *
 * For \f$S_5\f$ the Taylor series is
 *
 * \f{eqnarray*}{
 *    S_5(z+t) & = & \frac{1}{(z+t)^6}
 *               =   \sum_{k=0}^{n} (-1)^k \begin{pmatrix} k+5 \\ 5 \end{pmatrix} \frac{1}{(y_c^{(p)}-x_c^{(p)})^{k+6}} (y-y_c^{(p)})^k
 *                   + R_n(y-x_c^{(p)}) \\
 * \f}
 *
 * and if \f$l = k\f$ and \f$h = \left( \frac{1}{2} \right)^{k}\f$
 *
 * \f{eqnarray*}{
 * | R_n(y-x_c^{(p)}) |
 *                         & \leq &  \max_{y_c^{(p)}} \left( \left| (-1)^{n+1} \begin{pmatrix} n+6 \\ 5 \end{pmatrix} \frac{1}{(y_c^{(p)}-x_c^{(p)})^{n+7}} \right| \right)
 *                                           (y-y_c^{(p)})^{n+1} \quad \quad |y-y_c^{(p)}| \leq h \\
 *                         & \leq &  \left( \frac{1}{2} \right)^6 \begin{pmatrix} n+6 \\ 5 \end{pmatrix}  h^2 \left( \frac{1}{2} \right)^{2n-(8k-2)}
 * \f}
 *
 * For \f$l = 2\f$ and \f$h = \left( \frac{1}{2} \right)^2\f$ we have
 *
 * \f{eqnarray*}{
 * | R_n(y-x_c^{(p)}) |
 *                         & =    & \left( \frac{1}{2} \right)^6 \begin{pmatrix} n+6 \\ 5 \end{pmatrix} h_2^2 \left( \frac{1}{2} \right)^{2n-14}
 * \f}
 *
 * If we choose \f$n = 7\f$ we have \f$C = \left( \frac{1}{2} \right)^6 \begin{pmatrix} 7+6 \\ 5 \end{pmatrix} = \frac{1287}{64} = 20.109375 \f$
 *
 * \f{eqnarray*}{
 * | R_n(y-x_c^{(p)}) |
 *                         & \leq & \left( \frac{1}{2} \right)^6 \begin{pmatrix} 7+6 \\ 5 \end{pmatrix} h_2^2  = C h_2^2.
 * \f}
 *
 * We increase \f$n\f$ by 3 so that
 *
 * \f{eqnarray*}{
 * | R_n(y-x_c^{(p)}) |
 *                         & \leq & \left( \frac{1}{2} \right)^6 \begin{pmatrix} n+6 \\ 5 \end{pmatrix} h_2^2 \left( \frac{1}{2} \right)^{2n-14} \\
 *                         & =    & \left( \frac{1}{2} \right)^6 \begin{pmatrix} 10+6 \\ 5 \end{pmatrix} h_2^2 \left( \frac{1}{2} \right)^{2(10)-14} \\
 *                         & \approx & 1.07
 * \f}
 *
 * So we set \f$n = 7 + 3 = 10\f$
 *
 * \f$\mathbf{S_6} - \mathbf{\mbox{Seventh Power}}\f$
 *
 * \f$\vdots\f$
 *
 * \f$\mathbf{S_{48}} - \mathbf{\mbox{Forty-Ninth Power}}\f$
 *
 *
 * For \f$S_{48}\f$ the Taylor series is
 *
 * \f{eqnarray*}{
 *    S_{48}(z+t) & = & \frac{1}{(z+t)^6}
 *                =   \sum_{k=0}^{n} (-1)^k \begin{pmatrix} k+48 \\ 48 \end{pmatrix} \frac{1}{(y_c^{(p)}-x_c^{(p)})^{k+49}} (y-y_c^{(p)})^k
 *                   + R_n(y-x_c^{(p)}) \\
 * \f}
 *
 * and if \f$l = k\f$ and \f$h = \left( \frac{1}{2} \right)^{k}\f$
 *
 * \f{eqnarray*}{
 * | R_n(y-x_c^{(p)}) |
 *                         & \leq &  \max_{y_c^{(p)}} \left( \left| (-1)^{n+1} \begin{pmatrix} n+49 \\ 48 \end{pmatrix} \frac{1}{(y_c^{(p)}-x_c^{(p)})^{n+50}} \right| \right)
 *                                           (y-y_c^{(p)})^{n+1} \quad \quad |y-y_c^{(p)}| \leq h \\
 *                         & \leq &  \left( \frac{1}{2} \right)^{49} \begin{pmatrix} n+49 \\ 48 \end{pmatrix}  h^2 \left( \frac{1}{2} \right)^{2n-(51k-2)}
 * \f}
 *
 * For \f$l = 2\f$ and \f$h = \left( \frac{1}{2} \right)^2\f$ we have
 *
 * \f{eqnarray*}{
 * | R_n(y-x_c^{(p)}) |
 *                         & =    & \left( \frac{1}{2} \right)^{49} \begin{pmatrix} n+49 \\ 48 \end{pmatrix} h_2^2 \left( \frac{1}{2} \right)^{2n-100}
 * \f}
 *
 * If we choose \f$n = 50\f$ we have \f$C = \left( \frac{1}{2} \right)^{49} \begin{pmatrix} 50+49 \\ 48 \end{pmatrix}
 *                                                   = \frac{2.472 \times 10^{30}}{5.629 \times 10^{14}} = 8.609542 \times 10^{13} \f$
 *
 * We can predict \f$C\f$ based on the previous series.  Since \f$C\f$ nearly doubles with each power and for \f$S_2\f$ we have \f$C = 2.625\f$,
 * we expect \f$C = 2.625 2^{48-2} = 1.85 \times 10^{14}\f$ which is approximately off by a factor of 5.
 *
 * To keep the bounding constant \f$C \approx 1\f$ we solve
 *
 * \f[ \left( \frac{1}{2} \right)^{49} \begin{pmatrix} 50+49 \\ 48 \end{pmatrix} 2^{-2m} = 1
 *     \quad \implies \quad
 *     2m = \frac{\ln(8.609542 \times 10^{13})}{\ln(2)} \approx 46.29
 *     \quad \implies \quad
 *     m \approx 23.145
 * \f]
 *
 * We check to see how far off \f$C\f$ is from 1 using \f$n=50+23\f$.
 *
 * \f[ C = \left( \frac{1}{2} \right)^{49} \begin{pmatrix} (50+23)+49 \\ 48 \end{pmatrix} \left( \frac{1}{2} \right)^{2(50+23)-100}
 *      \approx 607,058 \f]
 *
 * Solving again
 *
 * \f[ \left( \frac{1}{2} \right)^{49} \begin{pmatrix} (50+23)+49 \\ 48 \end{pmatrix} \left( \frac{1}{2} \right)^{2(50+23)-100} 2^{-2m} = 1
 *     \quad \implies \quad
 *     2m = \frac{\ln(607,058)}{\ln(2)} \approx 19.21
 *     \quad \implies \quad
 *     m \approx 9.61
 * \f]
 *
 * We check C again
 *
 * \f[ C = \left( \frac{1}{2} \right)^{49} \begin{pmatrix} (50+23+10)+49 \\ 48 \end{pmatrix} \left( \frac{1}{2} \right)^{2(50+23+10)-100}
 *       \approx 65.43 \f]
 *
 *
 * \f[ \left( \frac{1}{2} \right)^{49} \begin{pmatrix} (50+23+10)+49 \\ 48 \end{pmatrix} \left( \frac{1}{2} \right)^{2(50+23+10)-100} 2^{-2m} = 1
 *     \quad \implies \quad
 *     2m = \frac{\ln(65.43)}{\ln(2)} \approx 6.03
 *     \quad \implies \quad
 *     m \approx 3.02
 * \f]
 *
 * We check to see how far off \f$C\f$ is from 1 using \f$n=50+23+10+3\f$.
 *
 * \f{eqnarray*}{
 *     C & = & \left( \frac{1}{2} \right)^{49} \begin{pmatrix} (50+23+10+3)+49 \\ 48 \end{pmatrix} \left( \frac{1}{2} \right)^{2(50+23+10+3)-100} \\
 *       & = & \left( \frac{1}{2} \right)^{49} \begin{pmatrix} (86)+49 \\ 48 \end{pmatrix} \left( \frac{1}{2} \right)^{2(86)-100} \\
 *       & = & \left( \frac{1}{2} \right)^{49} \begin{pmatrix} 135 \\ 48 \end{pmatrix} \left( \frac{1}{2} \right)^{172-100} \\
 *       \approx 6.87
 * \f}
 *
 * Solving again
 *
 * \f[ \left( \frac{1}{2} \right)^49 \begin{pmatrix} (50+23+10+3)+49 \\ 48 \end{pmatrix} \left( \frac{1}{2} \right)^{2(50+23+10+3)-100} 2^{-2m} = 1
 *     \quad \implies \quad
 *     2m = \frac{\ln(6.87)}{\ln(2)} \approx 2.78
 *     \quad \implies \quad
 *     m \approx 1.39
 * \f]
 *
 * We check to see how far off \f$C\f$ is from 1 using \f$n=50+23+10+3+1\f$.
 *
 * \f{eqnarray*}
 *     C & = & \left( \frac{1}{2} \right)^{49} \begin{pmatrix} (50+23+10+3+1)+49 \\ 48 \end{pmatrix} \left( \frac{1}{2} \right)^{2(50+23+10+3+1)-100} \\
 *       & = & \left( \frac{1}{2} \right)^{49} \begin{pmatrix} 87+49 \\ 48 \end{pmatrix} \left( \frac{1}{2} \right)^{2(87)-100} \\
 *       & = & \left( \frac{1}{2} \right)^{49} \begin{pmatrix} 136 \\ 48 \end{pmatrix} \left( \frac{1}{2} \right)^{2(87)-100} \\
 *       & = & \left( \frac{1}{2} \right)^{49} \frac{136!}{48!(136-48)!} \left( \frac{1}{2} \right)^{2(85)-100} \left( \frac{1}{2} \right)^{2} \\
 *       & = & \left( \frac{1}{2} \right)^{49} \frac{(136)135!}{48!(136-48)(135-48)!} \left( \frac{1}{2} \right)^{2(86)-100} \left( \frac{1}{2} \right)^{2} \\
 *       & = & \left( \frac{1}{2} \right)^{49} \frac{135!}{48!(135-48)!} \left( \frac{1}{2} \right)^{2(86)-100} \frac{136}{(136-48)} \left( \frac{1}{2} \right)^{2} \\
 *       & = & \left( \frac{1}{2} \right)^{49} \begin{pmatrix} 135 \\ 48 \end{pmatrix} \left( \frac{1}{2} \right)^{172-100} \frac{136}{(136-48)} \left( \frac{1}{2} \right)^{2} \\
 *       & \approx & (6.87)(0.39) \approx 2.68
 * \f}
 *
 * Adding 2 more powers of \f$\frac{1}{2}\f$ will bring the result below \f$1\f$.  Checking
 *
 * \f{eqnarray*}{
 *     C & = & \left( \frac{1}{2} \right)^49 \begin{pmatrix} (50+23+10+3+1+2)+49 \\ 48 \end{pmatrix} \left( \frac{1}{2} \right)^{2(50+23+10+3+1+2)-100} \\
 *       & = & \left( \frac{1}{2} \right)^49 \begin{pmatrix} 138 \\ 48 \end{pmatrix} \left( \frac{1}{2} \right)^{2(50+23+10+3+1+2)-100} \\
 *       & = & \left( \frac{1}{2} \right)^49 \frac{138!}{48!(138-48)!} \left( \frac{1}{2} \right)^{178-100} \\
 *       & = & \left( \frac{1}{2} \right)^49 \frac{(138)(137)136!}{48!(138-48)(137-48)(136-48)!} \left( \frac{1}{2} \right)^{174+4-100} \\
 *       & = & \left( \frac{1}{2} \right)^49 \frac{136!}{48!(136-48)!} \left( \frac{1}{2} \right)^{174-100} \frac{(138)(137)}{48!(138-48)(137-48)} \left( \frac{1}{2} \right)^{4} \\
 *       & = & \left( \frac{1}{2} \right)^49 \begin{pmatrix} 136 \\ 48 \end{pmatrix} \left( \frac{1}{2} \right)^{174-100} \frac{(138)(137)}{(138-48)(137-48)} \left( \frac{1}{2} \right)^{4} \\
 *       & \approx & (2.68)(0.15) \approx 0.40
 * \f}
 *
 * Therefore we choose \f$2n = (2)(50+23+10+3+1+2) = 2(50+39)\f$ \f$\implies\f$ \f$n = 50+39 = 89\f$
 *
 * \f{eqnarray*}{
 * | R_n(y-x_c^{(p)}) |
 *                         & \leq & \left( \frac{1}{2} \right)^{49} \begin{pmatrix} 50+39+49 \\ 48 \end{pmatrix} \left( \frac{1}{2} \right)^{(2)(89)-100} h_2^2 = C h_2^2.
 * \f}
 *
 *
 *
 *
 *
 * @subsubsection downward_pass_part_two Downward Pass Part Two
 *
 * The last part of the downward pass is the translation of near-field series from parent centers \f$y_c^{(p)}\f$
 * to child centers \f$y_c^{(c)}\f$.  Here we will be translating from a parent center at level \f$l=2\f$ to
 * a child center at level \f$l=3\f$.
 *
 * Recall the near-field expansion has the form \f$\sum_{n=0}^{\infty} b_n V_n\f$.
 * Let \f$z+t = y-y_c^{(p)}\f$, \f$z = y-y_c^{(c)}\f$ and \f$t = y_c^{(c)} - y_c^{(p)}\f$.
 * To construct the R|R translation we form the Taylor series expansions for each \f$V_n\f$.  From the
 * previous section we know that we will use Taylor series for \f$V_0\f$ to \f$V_{89}\f$.
 *
 *
 * \f$\mathbf{V_0} - \mathbf{\mbox{First Power}}\f$
 *
 * For the first power \f$V_0\f$
 *
 * \f{eqnarray*}{
 *   V_0 (z+t) = (z+t)^{0} = 1 = \sum_{n=0}^{\infty} \frac{V_0^{(n)}(z)}{n!} t^n  = 1 = z^0 = V_0(z)
 * \f}
 *
 * We have no error for the Taylor series of \f$V_0(z+t)\f$ since it consists of just one term.
 *
 *
 * \f$\mathbf{V_1} - \mathbf{\mbox{Second Power}}\f$
 *
 * For the second power \f$V_1\f$
 *
 * \f{eqnarray*}{
 *   V_1 (z+t) & = & (z+t)^{1} = \sum_{n=0}^{\infty} \frac{V_1^{(n)}(z)}{n!} t^n  \\
 *             & = & \frac{z}{0!} t^0 + \frac{1}{1!} t^1  \\
 *             & = &  V_1(z) + \frac{t}{1!} V_0(z) \\
 * \f}
 *
 * We have no error for the the second power if we include both terms that make up the Taylor series.
 *
 * \f$\mathbf{V_2} - \mathbf{\mbox{Third Power}}\f$
 *
 * For the third power \f$V_2\f$
 *
 * \f{eqnarray*}{
 *   V_2 (z+t) & = & (z+t)^{2} = \sum_{n=0}^{\infty} \frac{V_2^{(n)}(z)}{n!} t^n  \\
 *             & = & \frac{z^2}{0!} t^0 + \frac{2z}{1!} t^1 + \frac{2}{2!} t^2  \\
 *             & = &  V_2(z) + \frac{2t}{1!} V_1(z) + \frac{2t^2}{2!} V_0(z) \\
 * \f}
 *
 * We have no error for the the third power if we include all three terms that make up the Taylor series.
 *
 * \f$\mathbf{V_3} - \mathbf{\mbox{Fourth Power}}\f$
 *
 * For the fourth power \f$V_3\f$
 *
 * \f{eqnarray*}{
 *   V_3 (z+t) & = & (z+t)^{3} = \sum_{n=0}^{\infty} \frac{V_3^{(n)}(z)}{n!} t^n  \\
 *             & = & \frac{z^3}{0!} t^0 + \frac{3z^2}{1!} t^1 + \frac{(3)(2)z}{2!} t^2 + \frac{(3)(2)}{3!} t^3  \\
 *             & = &  V_3(z) + \frac{3t}{1!} V_2(z) + \frac{(3)(2)t^2}{2!} V_1(z) + \frac{(3)(2)(1)t^3}{3!} V_0(z) \\
 *             & = & \begin{pmatrix} 3 \\ 0 \end{pmatrix} V_3(z)
 *                   + \begin{pmatrix} 3 \\ 1 \end{pmatrix} V_2(z)
 *                   + \begin{pmatrix} 3 \\ 2 \end{pmatrix} V_1(z)
 *                   + \begin{pmatrix} 3 \\ 3 \end{pmatrix} V_0(z) \\
 *             & = & \sum_{n=0}^{3} \begin{pmatrix} 3 \\ n \end{pmatrix} V_{3-n}(z)
 * \f}
 *
 * We have no error for the the third power
 *
 *
 *
 *
 *
 *
 *
 *
 *
 */


/*
 * test1.cc
 *
 *  Created on: Mar 26, 2017
 *      Author: dbpc
 */



#include<iostream>
#include<vector>
#include<cmath>

#include "Point.h"
#include "Util.h"
//#include "Example1.h"
#include "Potential.h"
#include "FmmTree.h"
//#include "PotentialGradient.h"
//#include "FmmTreeGradient.h"


/**
 *  test1.cc
 *  Created on: March 27, 2017
 *      Author: dbpc
 *
 *  The program below tests the S|S, S|R, and R|R transformations of Potential.cc.
 *  The program performs the translations for a target as done in the Fast
 *  Multipole Method (FMM) at the refinement level = 3.  The series are truncated to
 *  order |alpha| = 3 resulting in p = 20 derivatives.  We note the S|S, S|R, and R|R
 *  transformation matrices are p x p matrices, and therefore the S|R matrix utilizes derivatives
 *  up to order |alpha| = 6.  For reference, the refinement levels l = 2 and l = 3 are shown below.
 *
  * level l = 2 refinement
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
  *
  *
  *
  * level l = 3 refinement
  *                                                                                                                                       y = 1.00
  *                                                                                                                             1.00 _______ _______ _______ _______
  *                                                                                                                                 |219|223|251|255|475|479|507|511|
  *                                                                                                 y = 0.875                       |---27--|---31--|---59--|---63--|
  *                                                                                                                             0.75|218|2221250|254|474|4785506|510|
  *                                                                                       1.00 _______ _______ _______ _______      |211|215|243|247|467|471|499|503|
  *               |                                                                           |217|221|249|253|473|477|505|509|     |---26--|---30--|---58--|---62--|
  *                                                           y = 0.75                        |---27--|---31--|---59--|---63--| 0.50|210|214|242|246|466|470|498|502|
  *                                                                                       0.75|216|2201248|252|472|4765504|508|     |155|159|187|191|411|415|443|447|
  *                                                 1.00 _______ _______ _______ _______      |209|213|241|245|465|469|497|501|     |---19--|---23--|---51--|---55--|
  *               z                                     |203|207|235|239|459|463|491|495|     |---26--|---30--|---58--|---62--| 0.25|154|1580186|190|410|4144442|446|
  *               |      y = 0.625                      |---25--|---29--|---57--|---61--| 0.50|208|212|240|244|464|468|496|500|     |147|151|179|183|403|407|435|439|
  *               |                                 0.75|202|2061234|238|458|4625490|494|     |153|157|185|189|409|413|441|445|     |---18--|---22--|---50--|---54--|
  *          1.00 |_______________________________      |195|199|227|231|451|455|483|487|     |---19--|---23--|---51--|---55--|     |146|150|178|182|402|406|434|438|
  *               |201|205|233|237|457|461|489|493|     |---24--|---28--|---56--|---60--| 0.25|152|1560184|188|408|4124440|444|    0.0     0.25    0.5    0 .75    1.0       x
  *               |---25--|---29--|---57--|---61--| 0.50|194|198|226|230|450|454|482|486|     |145|149|177|181|401|405|433|437|
  *          0.75 |200|2041232|236|456|4605488|492|     |139|143|171|175|395|399|427|431|     |---18--|---22--|---50--|---54--|
  *               |193|197|225|229|449|453|481|485|     |---17--|---21--|---49--|---53--|     |144|148|176|180|400|404|432|436|
  *               |---24--|---28--|---56--|---60--| 0.25|138|1420170|174|394|3984426|430|    0.0     0.25    0.5     0.75    1.0       x
  *          0.50 |192|196|224|228|448|452|480|484|     |131|135|163|167|387|391|419|423|
  *               |137|141|169|173|393|397|425|429|     |---16--|---20--|---48--|---52--|
  *               |---17--|---21--|---49--|---53--|     |130|134|162|166|386|390|418|422|
  *          0.25 |136|1400168|172|392|3964424|428|    0.0     0.25    0.5     0.75    1.0       x
  *               |129|133|161|165|385|389|417|421|
  *               |---16--|---20--|---48--|---52--|
  *               |128|132|160|164|384|388|416|420|______
  *              0.0     0.25    0.5     0.75    1.0       x
  *
  *
  *                                                                                                                               y = 0.50
  *                                                                                                                     1.00 _______ _______ _______ _______
  *                                                                                                                         |91 |95 |123|127|347|351|379|383|
  *                                                                                         y = 0.375                       |---11--|---15--|---43--|---47--|
  *                                                                                                                     0.75|90_|94_1122|126|346|3505378|382|
  *                                                                               1.00 _______ _______ _______ _______      |83 |87 |115|119|339|343|371|375|
  *       |                                                                           |89 |93 |121|125|345|349|377|381|     |---10--|---14--|---42--|---46--|
  *                                                   y = 0.25                        |---11--|---15--|---43--|---47--| 0.50|82_|86_|114|118|338|342|370|374|
  *                                                                               0.75|88_|92_1120|124|344|3485376|380|     |27 |31 | 59|63 |283|287|315|319|
  *                                         1.00 _______ _______ _______ _______      |81 |85 |113|117|337|341|369|373|     |---3---|---7---|---35--|---39--|
  *      z                                      |75 |79 |107|111|331|335|363|367|     |---10--|---14--|---42--|---46--| 0.25|26_|30_0_58|62_|282|2864314|318|
  *       |      y = 0.125                      |---9---|---13--|---41--|---45--| 0.50|80_|84_|112|116|336|340|368|372|     |19 |23 |51 |55 |275|279|307|311|
  *       |                                 0.75|74_|78_1106|110|330|3345362|366|     |25 |29 | 57| 61|281|285|313|317|     |---2---|---6---|---34--|---38--|
  *  1.00 |_______________________________      |67 |71 |99 |103|323|327|355|359|     |---3---|---7---|---35--|---39--|     |18_|22_|50_|54_|274|278|306|310|
  *       |73 |77 |105|109|329|333|361|365|     |---8---|---12--|---40--|---44--| 0.25|24_|28_0_56|_60|280|2844312|316|    0.0     0.25    0.5    0 .75    1.0       x
  *       |---9---|---13--|---41--|---45--| 0.50|66_|70_|98_|102|322|326|354|358|     |17 |21 |49 |53 |273|277|305|309|
  *  0.75 |72_|76_1104|108|328|3325360|364|     |11 |15 | 43|47 |267|271|299|303|     |---2---|---6---|---34--|---38--|
  *       |65 |69 |97 |101|321|325|353|357|     |---1---|---5---|---33--|---37--|     |16_|20_|48_|52_|272|276|304|308|
  *       |---8---|---12--|---40--|---44--| 0.25|10_|14_0_42|46_|266|2704298|302|    0.0     0.25    0.5     0.75    1.0       x
  *  0.50 |64_|68_|96_|100|320|324|352|356|     | 3 | 7 |35 |39 |259|263|291|295|
  *       | 9 |13 | 41|45 |265|269|297|301|     |---0---|---4---|---32--|---36--|
  *       |---1---|---5---|---33--|---37--|     |_2_|_6_|34_|38_|258|262|290|294|
  *  0.25 |_8_|12_0_40|44_|264|2684296|300|    0.0     0.25    0.5     0.75    1.0       x
  *       | 1 | 5 |33 |37 |257|261|289|293|
  *       |---0---|---4---|---32--|---36--|
  *       |_0_|_4_|32_|36_|256|260|288|292|______
  *      0.0     0.25    0.5     0.75    1.0       x
  *
 *
 *  In the diagram below are the target t = y[0], source s = x[0], and
 *  cell centers for the test.  And interaction list cell n = 45 at l = 3 will
 *  be used for the S|R translation test.  The cell contains the point o = x[1].
 *  Below are the coordinates for the points mentioned, their cell numbers, and
 *  the level we are testing.
 *    source - x[0] = [0.90, 0.925, 0.95] located in cell n = 511 at l = 3;
 *    point in interaction list cell - x[1] = [0.38, 0.06, 0.40] located in cell n = 45 at l = 3;
 *    target - y[0] = [0.15,0.03,0.20] located in cell n = 5 at l = 3;
 *
   *                 ________________________________________________________________
   *                |   |   |   |   |   |   |   |   ||   |   |   |   |   |   | s |   |
   *                |---|---|---|---|---|---|---|---||---|---|---|---|---|---|---z_b-|
   *                |___|___|___|___|___|___|___|___||___|___|___|___|___|___|___|___|
   *                |   |   |   |   |   |   |   |   ||   |   |   |   |   |   |   |   |
   *                |---|---|---|---|---|---|---|---||---|---|---|---|---|---|---|---|
   *                |___|___|___|___|___|___|___|___||___|___|___|___|___|___|___|___|
   *                |   |   |   |   |   |   |   |   ||   |   |   |   |   |   |   |   |
   *                |---|---|---|---|---|---|---|---||---|---|---|---|---|---|---|---|
   *                |___|___|___|___|___|___|___|___||___|___|___|___|___|___|___|___|
   *                |   |   |   |   |   |   |   |   ||   |   |   |   |   |   |   |   |
   *                |---|---|---|---|---|---|---|---||---|---|---|---|---|---|---|---|
   *                |   |   |   |   |   |   |   |   ||   |   |   |   |   |   |   |   |
   *                |===|===|===|===|===|===|===|===||===|===|===|===|===|===|===|===|
   *                |   |   |   |   #   |   | o |   ||   |   |   |   #   |   |   |   |
   *                |---i---|---i---#---i---|---i---||---|---|---|---#---|---|---|---|
   *                |___|__pn___|___#___|__pn___|___||___|___|___|___#___|___|___|___|
   *                |   |   |   |   #   |   |   |   ||   |   |   |   #   |   |   |   |                      t - target
   *                |---n---|---n---#---n---|---i---||---|---|---|---#---|---|---|---|                      + - x_a
   *                |   |   |   |   #   |   |   |   ||   |   |   |   #   |   |   |   |                      c - z_a,c - child center
   *                |###|###|###|#######|###|###|###||###|###|###|#######|###|###|###|                      p - z_a,p - parent center
   *                |   |   | t |   #   |   |   |   ||   |   |   |   #   |   |   |   |                      n - nearest neighbor
   *                |---n---|---c---#---n---|---i---||---|---|---|---#---|---|---|---|                      i - interaction list
   *                |___|__ p___|___#___|__pn___|___||___|___|___|___#___|___|___|___|                      o - interaction list
   *                |   |   |   |   #   |   |   |   ||   |   |   |   #   |   |   |   |                          cell for S|R
   *                |---n---|---n---#---n---|---i---||---|---|---|---#---|---|---|---|                     pn - parent neighbor
   *                |_ _|_ _|___|_ _#___|___|___|___||___|___|___|___#___|___|___|___|                      s - source
   *
 *
 *  1) test4.cc first performs a direct calculation of the potential function for the source s
 *     and the target t mentioned above.
 *  2) The code below then tests the member functions Potential::getSVector and Potential::getSCoeff
 *     of the class Potential in Potential.cc.  SVector returns the powers of the far-field series and
 *     getSCoeff returns the coefficients of the far-field approximation.  The two vectors are multiplied
 *     below to form the s-expansion and this far-field approximation is compared with the direct calculation
 *     for accuracy.
 *  3) The s-expansion calculated in (2) is then translated from the center of the cell (child cell)
 *     containing the source to the parent cell's center using the Potential::getSS member function which
 *     transforms a far-field series "from" one center "to" a far-field series having another center.
 *     The "from" and "to" arguments are the child-cell's center and parent-cell's center, respectively.
 *     We compare this to the direct calculation and truncated s-expansion approximation.
 *  4) The code then performs an S|R transformation of the far-field expansion centered at the parent cell
 *     of the source point from step (3) to the center of the parent cell containing the target point.
 *
 *  The code also calculates the value of the new series after the SS transformation.  All three
 *  results are output for comparison.
 *
 *
 *
 *
 *  The code below also checks (uses) that the following is working
 *    - FmmTree.cc member functions: constructor, printBoxInformation, getBox
 *    - Box.cc member functions: constructor, getCenter
 *    - Point.cc member function(s): constructor, getBoxIndex
 *    - Potential.cc member function(s): constructor, getSCoeff, getSVector, direct
 *
 */


int main()
{
  // note that the DEFAULT_NUM_LEVEL count starts on l = 1
  // Example:  If the default number of FMM levels is 5
  //           Then the levels are l = 0, l = 1, l = 2, l = 3, and l = 4
  int DEFAULT_NUM_LEVEL = 4;
  // abs_alpha is also the order_of_approximation - highest order of derivatives used in Taylor series approximation
  unsigned int abs_alpha = 4;
  // number of Taylor series terms used to approximate the potential function
  int p = (abs_alpha+1)*(abs_alpha+1);
  std::cout << "p = " << p << " and " << "order_of_approximation = " << abs_alpha << std::endl;

//  std::cout.precision(17);

  // p is the number of terms used in Taylor Series approximation
  // Setting p to 1 and adding number of derivatives at each order below
  // using (n+2)_C_2 combinations formula.
  //int p = 0;
  //int combination = 1; // number of derivatives for order |alpha| = 0
  //p = p + combination;
  //for (unsigned int i = 1; i <= abs_alpha; ++i)
  //{
  //  combination = 2*i+1;  // number derivatives for order |alpha| = i
  //  p = p + combination;
  //}

  //if (p == 16)
  //  std::cout << "p = " << p << " correct" << std::endl;
  //else // p \neq 20
  //  std::cout << "p != " << 20 << "incorrect" << std::endl;


//  Example1 example1(DEFAULT_NUM_LEVEL);


  std::vector<Point>  x, y;
  std::vector<double> u;


  // Below we use 8 source particles (2 in each coordinate direction)
  // per cell with pow(8,L-1) cells/boxes partitioning the domain
  unsigned int nSourceParticles = pow(2, 0);
  x.resize(nSourceParticles);
  //  y.resize(nTotalParticles);
  //  u.resize(nTotalParticles);
  //x.resize(1);

  // One target point
  y.resize(1);

  // u holds the charges of source particles
  u.resize(nSourceParticles);

  // The domain is the cube [0,1]^3
  // The cells per side is the number of cells along one side of the box
  // Dividing the length of a side by the number of cells gives the cell_length
  // We also calculate a quarter of the cell length
  double cells_per_side = pow(2.0,DEFAULT_NUM_LEVEL-1);
  std::cout << "Number of Boxes at the lowest level = " << cells_per_side * cells_per_side * cells_per_side << std::endl;
  double cell_length = (1.0-0.0)/cells_per_side;
  std::cout << "The cell length = " << cell_length << std::endl;
  double quarter_length = 0.25*cell_length;
  std::cout << "A Quarter of the length of a cell = " << quarter_length << std::endl;
  //  double three_quarter_length = 0.75*cell_length;

  // Let's run a check that the getBoxIndex function is working properly
  // and review where some points are located in what cells/boxes
  Point upper_right_corner(0.9375+quarter_length,quarter_length,0.9375+quarter_length);
  unsigned int cell_index_upper_right_corner = upper_right_corner.getBoxIndex(DEFAULT_NUM_LEVEL-1);
  std::cout << upper_right_corner.coordToString() << " is in box " << cell_index_upper_right_corner
		    << std::endl;
  // ********************************************************************
  // Creating the Source Points *****************************************
  // ********************************************************************
  //
  // Note: The points are the same as the 2D example.
  // Except that we are using xz-coordinates for xy-coordinates here.
  // Therefore, the extra y-coordinate fixed at y = quarter_length puts the points in the
  // 1st octant above the xz-plane.
  // Using the diagram above, we see the inner loop below increments along the x-axis
  // indicating that the order of the cells, w.r.t. the x and y vectors, containing the points will be
  //  0,  4,  32,  36, 256, 260, 288, 292,
  //  1,  5,  33,  37, 257, 261, 289, 293,    (working down (along) the x-axis and incrementing up the y-axis)
  //  8, 12,  40,  44, 264, 268, 296, 300,
  //  9, 13,  41,  45, 265, 269, 297, 301,
  // 64, 68,  96, 100, 320, 324, 352, 356,
  // 65, 69,  97, 101, 321, 325, 353, 357,
  // 72, 76, 104, 108, 328, 332, 360, 364,
  // 73, 77, 105, 109, 329, 333, 361, 365

  // index for points that will go into x (source) and y (targets) vectors
  unsigned int n_points = 0;
  //unsigned int flag_cell = 365;



  // level l = 2 refinement
  //                                                                                                                  y = 1.00
  //                                                                                                        1.00 ______ ______ ______ ______
  //                                                                                                            |      |      |      |      |
  //                                                                                y = 0.75                    |  27  |  31  |  59  |  63  |
  //                                                                                                        0.75|______3______|______7______|
  //                                                                      1.00 ______ ______ ______ ______      |      |      |      |      |
  //      |                                                                   |      |      |      |      |     |  26  |  30  |  58  |  62  |
  //                                              y = 0.50                    |  25  |  29  |o 57  |  61  | 0.50|______|______|______|______|
  //                                                                      0.75|______3______|______7______|     |      |      |      |      |
  //                                    1.00 ______ ______ ______ ______      |      |      |      |      |     |  19  |  23  |  51  |  55  |
  //     z                                  |      |      |      |      |     |  24  |  28  |  56  |  60  | 0.25|______2______|______6______|
  //      |      y = 0.25                   |  11  |  15  |  43  |  47  | 0.50|______|______|______|______|     |      |      |      |      |
  //      |                             0.75|______1______|______5______|     |      |      |      |      |     |  18  |  22  |  50  |  54  |
  // 1.00 |______ ______ ______ ______      |      |      |      |      |     |  17  |* 21  |  49  |  53  |     |______|______|______|______|
  //      |      |      |      |      |     |  10  |  14  |  42  |  46  | 0.25|______2______|______6______|    0.0    0.25   0.5    0.75   1.0       x
  //      |  9   |  13  |  41  |  45  | 0.50|______|______|______|______|     |      |      |      |      |
  // 0.75 |______1______|______5______|     |    x |      |      |      |     |  16  |  20  |  48  |  52  |
  //      |      |      |      |      |     |  3   |   7  |  35  |  39  |     |______|______|______|______|
  //      |  8   |  12  |  40  |  44  | 0.25|______0______|______4______|    0.0    0.25   0.5    0.75   1.0       x
  // 0.50 |______|______|______|______|     |      |      |      |      |
  //      |      |      |      |      |     |  2   |   6  |  34  |  38  |
  //      |  1   |  5   |  33  |  37  |     |______|______|______|______|
  // 0.25 |______0______|______4______|    0.0    0.25   0.5    0.75   1.0       x
  //      |      |      |      |      |
  //      |  0   |  4   |  32  |  36  |
  //      |______|______|______|______|______
  //     0.0    0.25   0.5    0.75   1.0       x
  //
  // looping over the lower front left corners of the Boxes making up the [0,1]^3 domain
  // Ex: Let the cell_length = 0.25 (level l = 2)
  //       for x_coord = 0.0 and y_coord = 0.0
  //         z_coord runs through 0.0, 0.0+0.25, 0.0+0.50, and 0.0+0.75
  //         - Seen above, these coordinates (0,0,0), (0,0,0.25), (0,0,0.50), (0,0,0.75)
  //           are the lower(z)-front(y)-left(x) corners of the front-left column of Boxes 0, 1, 8, and 9
  //       for x_coord = 0.0 and y_coord = 0.25
  //         z_coord runs through 0.0, 0.0+0.25, 0.0+0.50, and 0.0+0.75
  //         - Seen above, these coordinates (0,0.25,0), (0,0.25,0.25), (0,0.25,0.50), (0,0.25,0.75)
  //           are the lower(z)-front(y)-left(x) corners of the front-left column of Boxes 0, 1, 8, and 9
  //
  // Here we increment as if the cell length is 0.125 where the number of FMM levels would be 4
  // and the number of cells per side of the [0,1]^3 domain is 8
  // The quarter_length for this case would be 0.125/4 = 0.03125 and we set the source points based on
  // these values
/*
  for (double x_coord=0.0; x_coord<1.0; x_coord+=0.125)
    for (double y_coord=0.0; y_coord<1.0; y_coord+=0.125)
	  for (double z_coord=0.0; z_coord<1.0; z_coord+=0.125)
	  {
//	    Point flag_point(x_coord+quarter_length,y_coord+quarter_length,z_coord+quarter_length);
//	    unsigned int cell_index = flag_point.getBoxIndex(DEFAULT_NUM_LEVEL-1);
//        if(cell_index == 0 || cell_index == 5 || cell_index == 45 || cell_index == 365)
        if (x_coord == 0.0 && y_coord == 0.0 && z_coord == 0.0)  // in Box with Index = 0 for level = 3
        {
          std::cout << "x_coord = " << x_coord << std::endl;
          std::cout << "y_coord = " << y_coord << std::endl;
          std::cout << "z_coord = " << z_coord << std::endl;
		    y[0].setX(x_coord + 0.03125);       // source particle (point)
		    y[0].setY(y_coord + 0.03125);       // incrementing from lower left corner of cell
		    y[0].setZ(z_coord + 0.03125);       // using quarter_length for FMM with 4 levels
	        std::cout << "y[" << 0 << "] = " << y[0].coordToString() << std::endl;
	        std::cout << " and " << y[0].coordToString() << " is in box " << y[0].getBoxIndex(DEFAULT_NUM_LEVEL-1)
	    	  	      << std::endl;
        }
        if (x_coord == 0.875 && y_coord == 0.875 && z_coord == 0.875) // in Box with Index 511 for level l = 3
        {
          std::cout << "x_coord = " << x_coord << std::endl;
      	  std::cout << "y_coord = " << y_coord << std::endl;
      	  std::cout << "z_coord = " << z_coord << std::endl;
    	  for (unsigned int i = 0; i<2; ++i)      // stepping out a quarter length from lower front left corner of box
    	    for (unsigned int j = 0; j<2; ++j)
    	      for (unsigned int k = 0; k<2; ++k)
    	      {
    		    x[n_points].setX(x_coord + (2*i+1)*0.03125);       // source particle (point)
    		    x[n_points].setY(y_coord + (2*j+1)*0.03125);       // incrementing from lower left corner of cell
    		    x[n_points].setZ(z_coord + (2*k+1)*0.03125);       // using quarter_length for FMM with 4 levels
    	        std::cout << "x[" << n_points << "] = " << x[n_points].coordToString() << std::endl;
    	        std::cout << " and " << x[n_points].coordToString() << " is in box " << x[n_points].getBoxIndex(DEFAULT_NUM_LEVEL-1)
    	    	  	      << std::endl;
     		    ++n_points;
    	      }
        }
	  }
*/

// We hardcode the source x and target y particles to match with deal.ii quadrature points
   y[0].setX(0.026415608175648395);
   y[0].setY(0.026415608175648395);
   y[0].setZ(0.026415608175648395);

   x[0].setX(0.90141560817564848);
   x[0].setY(0.90141560817564836);
   x[0].setZ(0.90141560817564836);

/*
   x[1].setX(0.97358439182);
   x[1].setY(0.90141560818);
   x[1].setZ(0.90141560818);

   x[2].setX(0.90141560818);
   x[2].setY(0.97358439182);
   x[2].setZ(0.90141560818);

   x[3].setX(0.97358439182);
   x[3].setY(0.97358439182);
   x[3].setZ(0.90141560818);

   x[4].setX(0.90141560818);
   x[4].setY(0.90141560818);
   x[4].setZ(0.97358439182);

   x[5].setX(0.97358439182);
   x[5].setY(0.90141560818);
   x[5].setZ(0.97358439182);

   x[6].setX(0.901416);
   x[6].setY(0.973584);
   x[6].setZ(0.973584);

   x[7].setX(0.973584);
   x[7].setY(0.973584);
   x[7].setZ(0.973584);


   x[8].setX(0.026415608175648395);
   x[8].setY(0.52641560817564848);
   x[8].setZ(0.90141560817564836);

   x[9].setX(0.098584391824351608);
   x[9].setY(0.52641560817564848);
   x[9].setZ(0.90141560817564836);

   x[10].setX(0.026415608175648395);
   x[10].setY(0.59858439182435164);
   x[10].setZ(0.90141560817564836);

   x[11].setX(0.098584391824351608);
   x[11].setY(0.59858439182435164);
   x[11].setZ(0.90141560817564836);

   x[12].setX(0.026415608175648395);
   x[12].setY(0.52641560817564848);
   x[12].setZ(0.97358439182435164);

   x[13].setX(0.098584391824351608);
   x[13].setY(0.52641560817564848);
   x[13].setZ(0.97358439182435164);

   x[14].setX(0.026415608175648395);
   x[14].setY(0.59858439182435164);
   x[14].setZ(0.97358439182435164);

   x[15].setX(0.098584391824351608);
   x[15].setY(0.59858439182435164);
   x[15].setZ(0.97358439182435164);
*/

/*
   x[0].setX(0.026415608175648395);
   x[0].setY(0.52641560817564848);
   x[0].setZ(0.90141560817564836);

   x[1].setX(0.098584391824351608);
   x[1].setY(0.52641560817564848);
   x[1].setZ(0.90141560817564836);

   x[2].setX(0.026415608175648395);
   x[2].setY(0.59858439182435164);
   x[2].setZ(0.90141560817564836);

   x[3].setX(0.098584391824351608);
   x[3].setY(0.59858439182435164);
   x[3].setZ(0.90141560817564836);

   x[4].setX(0.026415608175648395);
   x[4].setY(0.52641560817564848);
   x[4].setZ(0.97358439182435164);

   x[5].setX(0.098584391824351608);
   x[5].setY(0.52641560817564848);
   x[5].setZ(0.97358439182435164);

   x[6].setX(0.026415608175648395);
   x[6].setY(0.59858439182435164);
   x[6].setZ(0.97358439182435164);

   x[7].setX(0.098584391824351608);
   x[7].setY(0.59858439182435164);
   x[7].setZ(0.97358439182435164);
*/
   //x_values[8] = 0.026415608175648395 0.52641560817564848 0.90141560817564836 is in box 201
   //x_values[9] = 0.098584391824351608 0.52641560817564836 0.90141560817564836 is in box 201
   //x_values[10] = 0.026415608175648392 0.59858439182435152 0.90141560817564848 is in box 201
   //x_values[11] = 0.098584391824351608 0.59858439182435164 0.90141560817564836 is in box 201
   //x_values[12] = 0.026415608175648395 0.52641560817564836 0.97358439182435164 is in box 201
   //x_values[13] = 0.098584391824351622 0.52641560817564836 0.97358439182435164 is in box 201
   //x_values[14] = 0.026415608175648395 0.59858439182435164 0.97358439182435175 is in box 201
   //x_values[15] = 0.098584391824351622 0.59858439182435164 0.97358439182435175 is in box 201



   // ********************************************************************
  // Creating the Target Points *****************************************
  // ********************************************************************
  // For this test only one target point
  // Note: In general the target points are same as source points
/*
  n_points = 0;

  for (double y_coord = 0.0; y_coord < 1.0; y_coord+=cell_length)
	for (double  x_coord = 0.0; x_coord < 1.0; x_coord+=cell_length)
	{
	  Point flag_point_2(x_coord+quarter_length,quarter_length,y_coord+quarter_length);
      unsigned int cell_index_2 = flag_point_2.getBoxIndex(DEFAULT_NUM_LEVEL-1);
	  if(cell_index_2 == 0)
	  {
	    y[n_points].setX(x_coord + quarter_length);       // target particle (point)
	    y[n_points].setZ(y_coord + quarter_length);       // lower left corner of cell
	    y[n_points].setY(quarter_length);
        std::cout << "y[" << n_points << "] = " << y[n_points].coordToString() << std::endl;
    	std::cout << " and " << y[n_points].coordToString() << " is in box " << y[n_points].getBoxIndex(DEFAULT_NUM_LEVEL-1)
    			    << std::endl;
//	    ++n_points;
//	    y[n_points].setX(x_coord + three_quarter_length); // source particle (point)
//	    y[n_points].setZ(y_coord + quarter_length);       // lower right corner of cell
//	    y[n_points].setY(quarter_length);
//    	std::cout << "x[" << n_points << "] = " << x[n_points].coordToString() << std::endl;
//	    ++n_points;
//	    y[n_points].setX(x_coord + quarter_length);       // source particle (point)
//	    y[n_points].setZ(y_coord + three_quarter_length); // upper left corner of cell
//	    y[n_points].setY(quarter_length);
//    	std::cout << "x[" << n_points << "] = " << x[n_points].coordToString() << std::endl;
//	    ++n_points;
//	    y[n_points].setX(x_coord + three_quarter_length); // source particle (point)
//	    y[n_points].setZ(y_coord + three_quarter_length); // upper right corner of cell
//	    y[n_points].setY(quarter_length);
//    	std::cout << "x[" << n_points << "] = " << x[n_points].coordToString() << std::endl;
//	    ++n_points;
	  }
      std::cout << "Looping through target points: "
                << " x_coord = " << x_coord << " y_coord = " << y_coord << std::endl;
	}
*/
  std::cout << "Target and Source points are set" << std::endl;

    int xBoxIndex = x[0].getBoxIndex(DEFAULT_NUM_LEVEL-1);
    int yBoxIndex = y[0].getBoxIndex(DEFAULT_NUM_LEVEL-1);
//    BoxGradient xBox(DEFAULT_NUM_LEVEL-1, xBoxIndex, p, abs_alpha);
//    BoxGradient yBox(DEFAULT_NUM_LEVEL-1, yBoxIndex, p, abs_alpha);
    Box xBox(DEFAULT_NUM_LEVEL-1, xBoxIndex, p, abs_alpha);
    Box yBox(DEFAULT_NUM_LEVEL-1, yBoxIndex, p, abs_alpha);
    std::cout << "Working on Boxes" << std::endl;

    Point center_xBoxChild = xBox.getCenter();
    Point center_yBoxChild = yBox.getCenter();

    std::cout << "yBoxChild cell n = " << yBoxIndex << " at level l = " << yBox.getLevel()
    		  << " has center = " << center_yBoxChild.coordToString() << std::endl;
    std::cout << "xBoxChild cell n = " << xBoxIndex << " at level l = " << xBox.getLevel()
    		  << " has center = " << center_xBoxChild.coordToString() << std::endl;

//    if (flag_cell == 365)
//    {
//      int xBoxParentIndex = x[0].getBoxIndex(DEFAULT_NUM_LEVEL-2);
//      int yBoxParentIndex = y[0].getBoxIndex(DEFAULT_NUM_LEVEL-2);
//      Box xBoxParent(DEFAULT_NUM_LEVEL-2, xBoxParentIndex, p, abs_alpha);
//      Box yBoxParent(DEFAULT_NUM_LEVEL-2, yBoxParentIndex, p, abs_alpha);
//      Point center_xBoxParent = xBoxParent.getCenter();
//      Point center_yBoxParent = yBoxParent.getCenter();
//      std::cout << "yBoxParent cell n = " << yBoxParentIndex << " at level l = " << yBoxParent.getLevel()
//    		    << " has center = " << center_yBoxParent.coordToString() << std::endl;
//      std::cout << "xBoxParent cell n = " << xBoxParentIndex << " at level l = " << xBoxParent.getLevel()
//    		    << " has center = " << center_xBoxParent.coordToString() << std::endl;
//    }


    std::cout << "We are here before setting the charges for the source points" << std::endl;
//  for (unsigned int i=0; i<u.size(); ++i)
//	  u[i] = 1.0;

/*
    u[0] = -0.48743173761521208;
    u[1] = -1.7922334787251739;
    u[2] = -0.48204223897677051;
    u[3] = -1.7724168784754299;
    u[4] = -0.53195571610036874;
    u[5] = -1.9559433045103198;
    u[6] = -0.52605747787099699;
    u[7] = -1.9342561241229774;

    u[8] = -1.4696897000236528;
    u[9] = -0.39970986488470323;
    u[10] = -1.3587323818075034;
    u[11] = -0.36953292707842128;
    u[12] = -1.6030130258631883;
    u[13] = -0.43596966078341737;
    u[14] = -1.481635518128485;
    u[15] = -0.40295875568154471;
*/

/*
    u[0] = -1.4696897000236528;
    u[1] = -0.39970986488470323;
    u[2] = -1.3587323818075034;
    u[3] = -0.36953292707842128;
    u[4] = -1.6030130258631883;
    u[5] = -0.43596966078341737;
    u[6] = -1.481635518128485;
    u[7] = -0.40295875568154471;

    u[8] = -0.48743173761521208;
    u[9] = -1.7922334787251739;
    u[10] = -0.48204223897677051;
    u[11] = -1.7724168784754299;
    u[12] = -0.53195571610036874;
    u[13] = -1.9559433045103198;
    u[14] = -0.52605747787099699;
    u[15] = -1.9342561241229774;
*/



    u[0] = -0.000028553269108821648;
/*
    u[0] = 1.4696897000236528;
    u[1] = 0.39970986488470323;
    u[2] = 1.6030130258631907;
    u[3] = 0.43596966078341742;
    u[4] = 1.3587323818075014;
    u[5] = 0.36953292707842134;
    u[6] = 1.481635518128485;
    u[7] = 0.40295875568154471;

    u[8] = 0.032075996874503299;
    u[9] = 0.11793995143449294;
   	u[10] = 0.11838563323222695;
   	u[11] = 0.43529109597367188;
   	u[12] = 0.029671446115807638;
   	u[13] = 0.10909867985026517;
   	u[14] = 0.10950752899770261;
    u[15] = 0.40264727242091125;
*/



/*
    		x_values[0] = 0.026415608175648395 0.52641560817564848 0.90141560817564836 and current at x-values[i] is J[0] = 0 -0.48743173761521208 0.032075996874503299 is in box 201
    		x_values[1] = 0.098584391824351608 0.52641560817564836 0.90141560817564836 and current at x-values[i] is J[1] = 0 -1.7922334787251739 0.11793995143449294 is in box 201
    		x_values[2] = 0.026415608175648392 0.59858439182435152 0.90141560817564848 and current at x-values[i] is J[2] = 0 -0.48204223897677051 0.11838563323222695 is in box 201
    		x_values[3] = 0.098584391824351608 0.59858439182435164 0.90141560817564836 and current at x-values[i] is J[3] = 0 -1.7724168784754299 0.43529109597367188 is in box 201
    		x_values[4] = 0.026415608175648395 0.52641560817564836 0.97358439182435164 and current at x-values[i] is J[4] = 0 -0.53195571610036874 0.029671446115807638 is in box 201
    		x_values[5] = 0.098584391824351622 0.52641560817564836 0.97358439182435164 and current at x-values[i] is J[5] = 0 -1.9559433045103198 0.10909867985026517 is in box 201
    		x_values[6] = 0.026415608175648395 0.59858439182435164 0.97358439182435175 and current at x-values[i] is J[6] = 0 -0.52605747787099699 0.10950752899770261 is in box 201
    		x_values[7] = 0.098584391824351622 0.59858439182435164 0.97358439182435175 and current at x-values[i] is J[7] = 0 -1.9342561241229774 0.40264727242091125 is in box 201
    		x_values[0] = 0.90141560817564848 0.90141560817564836 0.90141560817564836 and current at x-values[i] is J[0] = 0 -1.4696897000236528 1.4696897000236528 is in box 511
    		x_values[1] = 0.97358439182435164 0.90141560817564836 0.90141560817564836 and current at x-values[i] is J[1] = 0 -0.39970986488470323 0.39970986488470323 is in box 511
    		x_values[2] = 0.90141560817564836 0.97358439182435164 0.90141560817564848 and current at x-values[i] is J[2] = 0 -1.3587323818075034 1.6030130258631907 is in box 511
    		x_values[3] = 0.97358439182435164 0.97358439182435164 0.90141560817564836 and current at x-values[i] is J[3] = 0 -0.36953292707842128 0.43596966078341742 is in box 511
    		x_values[4] = 0.90141560817564848 0.90141560817564848 0.97358439182435164 and current at x-values[i] is J[4] = 0 -1.6030130258631883 1.3587323818075014 is in box 511
    		x_values[5] = 0.97358439182435164 0.90141560817564848 0.97358439182435164 and current at x-values[i] is J[5] = 0 -0.43596966078341737 0.36953292707842134 is in box 511
    		x_values[6] = 0.90141560817564848 0.97358439182435175 0.97358439182435175 and current at x-values[i] is J[6] = 0 -1.481635518128485 1.481635518128485 is in box 511
    		x_values[7] = 0.97358439182435164 0.97358439182435175 0.97358439182435175 and current at x-values[i] is J[7] = 0 -0.40295875568154471 0.40295875568154471 is in box 511
*/


  // Declaring the Potential Function needed by fmmtree data structure
  std::cout << "We are here before constructor for the potential function" << std::endl;
  Potential potential(p,abs_alpha);
//  PotentialGradient potential(p,abs_alpha);

  std::cout << "We are here before fmmtree construct" << std::endl;
  std::cout << "The number of source particles is " << x.size() << std::endl;
  std::cout << "The number of target particles is " << y.size() << std::endl;
  std::cout << "The potential function has p = " << potential.getP() << std::endl;
  std::cout << "The potential function has order_of_approximation = " << potential.getOrderApprox() << std::endl;

  // Declaring the FMM Tree data structure
  FmmTree fmmtree(DEFAULT_NUM_LEVEL,x,y,potential);
//  FmmTreeGradient fmmtree(DEFAULT_NUM_LEVEL,x,y,potential);

//  std::vector<std::vector<double> > direct = fmmtree.solveDirect(u);
//  std::vector<double> direct = fmmtree.solveDirect(u);
  std::vector<double> direct = fmmtree.solveDirect(u);
  for (unsigned int i=0; i<direct.size(); ++i)
    std::cout << "direct[" << i << "] = [" << direct[i] << "]" << std::endl;
//    std::cout << "direct[" << i << "] = [" << direct[i][0] << "," << direct[i][1] << "," << direct[i][2] << "]" << std::endl;

std::cin.ignore();
std::cin.get();

//  std::vector<std::vector<double> > indirect = fmmtree.solve(u);
//  std::vector<double> indirect = fmmtree.solve(u);
  std::vector<std::vector<std::complex<double> > > indirect = fmmtree.solve(u);

std::cin.ignore();
std::cin.get();

  for (unsigned int i=0; i<direct.size(); ++i)
  {
//    std::cout << "direct[" << i << "] = [" << direct[i][0] << "," << direct[i][1] << "," << direct[i][2] << "]" << " versus "
//    		  << "indirect[" << i << "] = [" << indirect[i][0] << "," << indirect[i][1] << "," << indirect[i][2] << "]" << "\n";
    std::cout << "direct[" << i << "] = [" << direct[i] << "]" << " versus "
    		  << "indirect[0][" << i << "] = [" << indirect[0][i] << "]" << "\n";
  }

std::cin.ignore();
std::cin.get();

  std::cout << "direct.size() = " << direct.size() << std::endl;
  double error = 0.0;
  for (unsigned int i = 0; i<direct.size(); ++i)
  {
//	  double tmp = std::abs(direct[i][0]-indirect[i][0])+std::abs(direct[i][1]-indirect[i][1])+std::abs(direct[i][2]-indirect[i][2]);
	  double tmp = std::abs(direct[i]-indirect[0][i].real());
	  if (tmp>error)
        error = tmp;
  }

  std::cout << "Error = " << error << "\n";



  std::cin.ignore();
  std::cin.get();

  std::vector<std::vector<double> > direct_grad = fmmtree.solveDirectGrad(u);
  for (unsigned int i=0; i<y.size(); ++i)
    std::cout << "direct_grad[" << i << "] = [" << direct_grad[0][i] << ", " << direct_grad[1][i] << ", "
                                                   << direct_grad[2][i] << "]" << std::endl;


  std::cin.ignore();
  std::cin.get();

    for (unsigned int i=0; i<y.size(); ++i)
    {
  //    std::cout << "direct[" << i << "] = [" << direct[i][0] << "," << direct[i][1] << "," << direct[i][2] << "]" << " versus "
  //    		  << "indirect[" << i << "] = [" << indirect[i][0] << "," << indirect[i][1] << "," << indirect[i][2] << "]" << "\n";
	    std::cout << "direct_grad[" << i << "] = [" << direct_grad[0][i] << ", " << direct_grad[1][i] << ", " << direct_grad[2][i] << "]" << " versus "
    		  << "indirect[" << i << "] = [" << indirect[1][i] << ", " << indirect[2][i] << ", " << indirect[3][i] << "]" << "\n";
        //std::cout << "direct_grad[" << i << "] = [" << direct_grad[i] << "]" << " versus "
      	//	  << "indirect[1][" << i << "] = [" << indirect[1][i] << "]" << "\n";
    }


    std::cout << "direct_grad.size() = " << direct_grad.size() << std::endl;
    double error_grad1 = 0.0;
    double error_grad2 = 0.0;
    for (unsigned int i = 0; i<y.size(); ++i)
    {
  //	  double tmp = std::abs(direct[i][0]-indirect[i][0])+std::abs(direct[i][1]-indirect[i][1])+std::abs(direct[i][2]-indirect[i][2]);
  	  double tmp1 = std::abs(direct_grad[0][i]-indirect[1][i].real());
  	  double tmp2 = std::abs(direct_grad[1][i]-indirect[2][i].real());
  	  if (tmp1>error_grad1)
          error_grad1 = tmp1;
  	  if (tmp2>error_grad2)
          error_grad2 = tmp2;
    }

    std::cout << "ErrorGrad = [" << error_grad1 << ", " << error_grad2 << "]" << "\n";


  std::cout << "Finished" << "\n";




  /*
   * The code below was used to check the FMM results above as well as troubleshoot
   * the FMM code.  It breaks down the FMM step by step
   *
    int xBoxIndex = x[0].getBoxIndex(DEFAULT_NUM_LEVEL-1);
    int yBoxIndex = y[0].getBoxIndex(DEFAULT_NUM_LEVEL-1);
    Box xBox(DEFAULT_NUM_LEVEL-1, xBoxIndex, p);
    Box yBox(DEFAULT_NUM_LEVEL-1, yBoxIndex, p);
    int xBoxParentIndex = x[0].getBoxIndex(DEFAULT_NUM_LEVEL-2);
    int yBoxParentIndex = y[0].getBoxIndex(DEFAULT_NUM_LEVEL-2);
    Box xBoxParent(DEFAULT_NUM_LEVEL-2, xBoxParentIndex, p);
    Box yBoxParent(DEFAULT_NUM_LEVEL-2, yBoxParentIndex, p);

    Point center_xBoxChild = xBox.getCenter();
    Point center_yBoxChild = yBox.getCenter();
    Point center_xBoxParent = xBoxParent.getCenter();
    Point center_yBoxParent = yBoxParent.getCenter();

    std::cout << "yBoxChild center cell n = 0 at level l = 3 is center_yBoxChild = " << center_yBoxChild.coordToString() << std::endl;
    std::cout << "xBoxChild center cell n = 365 at level l = 3 is center_xBoxChild =  = " << center_xBoxChild.coordToString() << std::endl;
    std::cout << "yBoxParent center cell n = 0 at level l = 2 is center_yBoxParent = " << center_yBoxParent.coordToString() << std::endl;
    std::cout << "xBoxParnet center cell n = 45 at level l = 2 is center_xBoxParent = " << center_xBoxParent.coordToString() << std::endl;

    std::vector<double> SCoeffChild = potential.getSCoeff(x[0].getCoord(), center_xBoxChild.getCoord());
    for (unsigned int i=0; i<SCoeffChild.size(); ++i)
      std::cout << "SCoeffChild[" << i << "] = " << SCoeffChild[i] << std::endl;
    std::cout << std::endl << std::endl;

    std::vector<double> SVecChild = potential.getSVector(y[0].getCoord(), center_xBoxChild.getCoord());
    std::vector<double> SCoeffParent = potential.getSS(center_xBoxChild.getCoord(), center_xBoxParent.getCoord(),
  		                                             SCoeffChild);
    for (unsigned int i=0; i<SCoeffParent.size(); ++i)
      std::cout << "SCoeffParent[" << i << "] = " << SCoeffParent[i] << std::endl;
    std::cout << std::endl << std::endl;


    std::vector<double> SVecParent = potential.getSVector(y[0].getCoord(), center_xBoxParent.getCoord());


    double s_expansion_child = 0.0;  double s_expansion_parent = 0.0;
    for (unsigned int i=0; i<SVecChild.size(); ++i)
    {
      s_expansion_child += SCoeffChild[i] * SVecChild[i];
      s_expansion_parent += SCoeffParent[i] * SVecParent[i];
    }

    std::cout << "s_expansion_child = " << s_expansion_child << std::endl;
    std::cout << "s_expansion_parent = " << s_expansion_parent << std::endl;



    std::vector<double> RCoeffParent = potential.getSR(center_xBoxParent.getCoord(), center_yBoxParent.getCoord(),
  		                                                     SCoeffParent);
    std::vector<double> RVecParent = potential.getRVector(y[0].getCoord(), center_yBoxParent.getCoord());
    std::vector<double> RCoeffChild = potential.getRR(center_yBoxParent.getCoord(), center_yBoxChild.getCoord(),
  		                                            RCoeffParent);
    std::vector<double> RVecChild = potential.getRVector(y[0].getCoord(), center_yBoxChild.getCoord());

    double r_expansion_parent = 0.0; double r_expansion_child = 0.0;
    for (unsigned int i=0; i<RVecParent.size(); ++i)
    {
      r_expansion_parent += RCoeffParent[i] * RVecParent[i];
      r_expansion_child += RCoeffChild[i] * RVecChild[i];
    }

    std::cout << "r_expansion_parent = " << r_expansion_parent << std::endl;
    std::cout << "r_expansion_child = " << r_expansion_child << std::endl;


    std::cout << "direct calculation = " << potential.direct(y[0].getCoord(), x[0].getCoord()) << std::endl;

  */


}











  /*******************************************************************************
   * far-field to far-field translation - S|S transformation
   *
   *
   *                 _____________________________
   *                |              |              |                         * x_b
   *                |              |              |
   *                |              |              |
   *                |       *      |              |
   *                |   *   z_a,c  |              |
   *                |   x_a        |              |
   *                |--------------*--------------|
   *                |              | z_a,p        |
   *                |              |              |
   *                |              |              |
   *                |              |              |
   *                |              |              |
   *                |______________|______________|
   *
   * Recall that for the far-field expansion the powers (and factorial) (1/alpha) (x_b - z_b)^{alpha}
   * are the coefficients of the series
   *
   *      f(x_b - x_a) = sum_{|alpha| >= 0} (1/alpha) (x_b - z_b)^{alpha} D^{alpha} f(z_b - x_a)
   *
   * and the derivatives  D^{alpha} f(z_b - x_a) are treated as the powers.
   * All the particles of a cell then have the same powers and we can combine the series of all
   * the particles into a single series.  The coefficients of this new series are a summation of
   * the coefficients of all the old series.
   *
   * Once we have a single series for each cell we perform the upward pass.  The upward pass consists
   * of an ss transformation that transforms a far-field series into another far-field series.
   * The old series is translated to a new location.  That is, the center of the old series is
   * changed to a new center.  In the upward pass the old center is the center z_a,c of the child cell.
   * The new center is the center z_a,p of the parent cell.
   *
   * Since we know where the new center will be and what the new powers will look like, we are really
   * only interested in what the new coefficients will be.  That is, the main goal of performing
   * the ss transformation is to take old coefficients to new coefficients.
   *
   * Recall from Main.cc that the new coefficients are determined by writing the powers
   * of the old series D^{alpha} f(y+t) centered at y + t = z_b,c the center of the child cell
   * as Taylor series expansions about the center y = z_b,p of the parent cell
   *
   *            D^{alpha} f(y+t) = sum_{|beta| = 0} b_{beta} D^{beta} f(y)
   *
   * Each of these series has the same powers and therefore can also be combined.  First we multiply
   * each equation above by the corresponding coefficient of the old series.  Then we sum all the equations.
   * The left hand side of the summation is the old series centered at y + t = z_b,c the center of the child cell.
   * The right-hand side will be our new series centered at y = z_b,p the center of the parent cell.  To write
   * our new series as a power series we collect like powers D^{beta} f(y) on the right-hand side.  From the like
   * terms have a power D^{beta} f(y) we factor off the power and the result is a coefficient for our new series
   * centered at y = z_b,p the center of the parent cell.
   *
   * Since we multiplied each equation by a coefficient of the old series, we can see that the coefficients of the
   * new series will be a linear combination of the coefficients of the old series.  We write each linear combination as
   * a row of a matrix.  This matrix is the ss-transformation matrix.  In Main.cc these linear combinations are written
   * as columns.
   *
   * In the code below, the source is x_b = x[0] and the target is x_a = y[0].
   * The parent center is y = z_b,p = center_xBoxParent and the child center
   * is y + t = z_b,c = center_xBox.  Therefore
   *
   *    t = y + t - y = z_b,c - z_b,p = center_xBox - center_xBoxParent = from - to
   *
   * and below we can see that in the call to potential.getSS the from = center_yBox and
   * to = center_yBoxParent
   *
   * The next call to potential.getSVector obtains the powers of the far-field expansion
   *
   *    D^{alpha} f(x_b - z_a,p)
   *
   * The powers are centered at the parent cell's center z_a,p = center_yBoxParent
   *******************************************************************************/


/*
  std::cout << "**************************************************" << std::endl;
  std::cout << "Far-Field to Far-Field Translation" << std::endl;
  std::cout << "**************************************************" << std::endl;
  std::cout << std::endl;


  // s-expansion shown above after translation to center of parent cell
  std::complex<double> s_expansion_translated;
  s_expansion_translated.real(0.0);
  s_expansion_translated.imag(0.0);

  Point center_xBoxParent = xBoxParent.getCenter();
  std::cout << "center_xBoxParent(x,y,z) = " << "center_xBoxParent(" << center_xBoxParent.getX() << ","
		    << center_xBoxParent.getY() << "," << center_xBoxParent.getZ() << ")" << std::endl;

  std::vector<double> SCoeffTranslated = potential.getSS(center_xBox.getCoord(), center_xBoxParent.getCoord(),
 		                                                   y_SCoeff);

  std::vector<double> SVecParentCenter = potential.getSVector(y[0].getCoord(), center_xBoxParent.getCoord());

//  for (unsigned int i=0; i < y_SVec_parent_center.size(); ++i)
//    std::cout << "y_SVec_parent_center[" << i << "] = " << y_SVec_parent_center[i] << std::endl;

  for (unsigned int i=0; i < SVecParentCenter.size(); ++i)
    s_expansion_translated += SCoeffTranslated[i] * SVecParentCenter[i];

  std::cout << "s_expansion_translated = " << s_expansion_translated << std::endl;

  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "**************************************************" << std::endl;
  std::cout << "**************************************************" << std::endl;

*/

  /*******************************************************************************
   * far-field to near-field translation - S|R transformation
   *
   *
   *                 _______________________________________________________________
   *                |   |   |   |   |   |   |   |   ||   |   |   |   #   |   | + |   |
   *                |---|---|---|---|---|---|---|---||---|---|---|---#---|---|--ffc--|
   *                |___|___|___|___|___|___|___|___||___|___|___|___#___|___p___|___|
   *                |   |   |   |   |   |   |   |   ||   |   |   |   #   |   |   |   |
   *                |---|---|---|---|---|---|---|---||---|---|---|---#---|---|---|---|
   *                |   |   |   |   |   |   |   |   ||   |   |   |   #   |   |   |   |
   *                |###|###|###|#######|###|###|###||###|###|###|###|###|###|###|###|
   *                |   |   |   |   |   |   |   |   ||   |   |   |   #   |   |   |   |
   *                |---|---|---|---|---|---|---|---||---|---|---|---#---|---|---|---|
   *                |___|___|___|___|___|___|___|___||___|___|___|___#___|___|___|___|
   *                |   |   |   |   |   |   |   |   ||   |   |   |   #   |   |   |   |
   *                |---|---|---|---|---|---|---|---||---|---|---|---#---|---|---|---|
   *                |   |   |   |   |   |   |   |   ||   |   |   |   #   |   |   |   |
   *                |===|===|===|===|===|===|===|===||===|===|===|===|===|===|===|===|
   *                |   |   |   |   #   |   | o |   ||   |   |   |   #   |   |   |   |
   *                |---i---|---i---#---i---|---i---||---|---|---|---#---|---|---|---|
   *                |___|__pn___|___#___|__pn___|___||___|___|___|___#___|___|___|___|
   *                |   |   |   |   #   |   |   |   ||   |   |   |   #   |   |   |   |                    ffc - far field cell
   *                |---n---|---n---#---n---|---i---||---|---|---|---#---|---|---|---|                      + - x_b
   *                |   |   |   |   #   |   |   |   ||   |   |   |   #   |   |   |   |                      c - z_a,c - parent cell
   *                |###|###|###|#######|###|###|###||###|###|###|###|###|###|###|###|                      p - z_a,p - parent cell
   *                |   |   | * |   #   |   |   |   ||   |   |   |   #   |   |   |   |                      n - nearest neighbor
   *                |---n---|---c---#---n---|---i---||---|---|---|---#---|---|---|---|                      i - interaction list
   *                |___|__ p___|___#___|__pn___|___||___|___|___|___#___|___|___|___|                      o - interaction list
   *                |   |   |   |   #   |   |   |   ||   |   |   |   #   |   |   |   |                          cell for S|R
   *                |---n---|---n---#---n---|---i---||---|---|---|---#---|---|---|---|                     pn - parent neighbor
   *                |_ _|_ _|___|_ _#___|___|___|___||___|___|___|___#___|___|___|___|                      * - x_a
   *
   * The S|R transformation is a shift of center and change of base
   * of a far-field approximation for a cell in the interaction list of a
   * cell of interest to a near-field approximation centered at the cell of
   * interest's center.
   *
   * Recall the far-field expansion coefficients are (1/alpha) (x_b - z_b)^{alpha}
   * and the powers are derivatives  D^{alpha} f(x_b - z_a) to make up the series
   *
   *      f(x_b - x_a) = sum_{|alpha| >= 0} (1/alpha) (x_b - z_b)^{alpha}  D^{alpha} f(z_b - x_a)
   *
   * The near-field expansion coefficients are (1/alpha) D^{alpha} f(x_b - z_a)
   * and the powers are (z_a - x_a)^{alpha} to make up the series
   *
   *      f(x_b - x_a) = sum_{|alpha| >= 0} (1/alpha) D^{alpha} f(x_b - z_a) (z_a - x_a)^{alpha}
   *
   * The s-expansion far-field approximation for the potential of the source x_b acting
   * on the target x_a is shifted and transformed into an r-expansion near-field approximation
   * for the potential of x_b acting on x_a.  The shift and transformation is part of the
   * downward pass of a potential function s-approximation from the peer cell of an interaction
   * list to the cell having the interaction list.
   *
   * In the code below the S|R transformation takes coefficients for a series centered
   * at the child cell's center to the new coefficients for the series centered at the parent
   * cell's center.
   * The s-expansion far-field approximation for the potential of the source x_b = x[0] acting
   * on the target x_a = y[0] is shifted and transformed into an r-expansion near-field approximation
   * for the potential of x_b = x[0] acting on x_a = y[0].  The shift is from the parent center of the
   * source to the parent center of the target.  The transformation from an s-expansion to an
   * r-expansion is necessary for convergence.
   *******************************************************************************/

/*
  std::cout << "**************************************************" << std::endl;
  std::cout << "Far-Field to Near-Field Translation" << std::endl;
  std::cout << "**************************************************" << std::endl;
  std::cout << std::endl;

  double r_expansion_s_transformed = 0.0;

  Point center_yBoxParent = yBoxParent.getCenter();

  // checking location of points
  // target - located in cell n = 5 at l = 3
  std::cout << "target y[0] = [" << y[0].getX() << "," << y[0].getY() << "," << y[0].getZ() << "]"
		    << std::endl;
  std::cout << "target y[0] is located in Box " << yBox.getIndex()
		    << " at level l = 3."<< std::endl;
  std::cout << "source x[0] = [" << x[0].getX() << "," << x[0].getY() << "," << x[0].getZ() << "]"
		    << std::endl;
  std::cout << "source x[0] is located in Box " << xBox.getIndex()
		    << " at level l = 3."<< std::endl;

  std::cout << "center_yBoxParent(x,y,z) = " << "center_yBoxParent(" << center_yBoxParent.getX() << ","
		    << center_yBoxParent.getY() << "," << center_yBoxParent.getZ() << ")" << std::endl;
  std::cout << "center_xBoxParent(x,y,z) = " << "center_xBoxParent(" << center_xBoxParent.getX() << ","
		    << center_xBoxParent.getY() << "," << center_xBoxParent.getZ() << ")" << std::endl;

//  for (unsigned int i=0; i<x2_1_SCoeff.size(); ++i)
//    std::cout << "x2_1_SCoeff[" << i << "] = " << x2_1_SCoeff[i] << std::endl;


  std::vector<double> SR_RCoeff = potential.getSR(center_xBoxParent.getCoord(), center_yBoxParent.getCoord(),
		                                                     SCoeffTranslated);
  for (unsigned int i=0; i<SR_RCoeff.size(); ++i)
    std::cout << "SR_RCoeff[" << i << "] = " << SR_RCoeff[i] << std::endl;

  std::vector<double> y_RVecSTransformed = potential.getRVector(y[0].getCoord(), center_yBoxParent.getCoord());

  for (unsigned int i=0; i<y_RVecSTransformed.size(); ++i)
    std::cout << "y_RVecSTransformed[" << i << "] = " << y_RVecSTransformed[i] << std::endl;

//  std::cout << "y_RCoeff_New.size() = " << y_RCoeff_New.size() << std::endl;
//  std::cout << "y_RVec_parent_center.size = " << y_RVec_parent_center.size() << std::endl;

  for (unsigned int i=0; i < y_RVecSTransformed.size(); ++i)
    r_expansion_s_transformed += SR_RCoeff[i] * y_RVecSTransformed[i];

  std::cout << "r_expansion_s_transformed = " << r_expansion_s_transformed << std::endl;

  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "**************************************************" << std::endl;
  std::cout << "**************************************************" << std::endl;

*/

  /*******************************************************************************
   * near-field expansion - r-expansion
   *
   * Performing r expansion
   * Similar to s-expansion, but the center of r-expansion is at center of box containing
   * the target point y.  The s-expansion center is at the center of box containing the
   * source point x.  Locations of center allow for convergence of the series.
   * In code below, the r-expansion center is at the center of the parent cell
   * containing the target point y.  In the next part of the code, the center is
   * shifted to the child cell center using an R|R translation
   *******************************************************************************/

/*

  double r_expansion = 0.0;

  std::vector<double> RCoeff = potential.getRCoeff(x[0].getCoord(), center_yBoxParent.getCoord());

  std::vector<double> RVec = potential.getRVector(y[0].getCoord(), center_yBoxParent.getCoord());

  for (unsigned int i=0; i<y_SVec.size(); ++i)
    r_expansion += RCoeff[i] * RVec[i];

  std::cout << "r_expansion = " << r_expansion << std::endl;


*/


  /*******************************************************************************
   * near-field to near-field translation
   *******************************************************************************/

/*

  // Performing R|R transformation that takes old coefficients for center at parent to new coefficients
  // having center at child.  The new approximation for potential of source x acting on target y
  // has its powers evaluated at the target minus the parent center (y - y_c^{(Parent(n),l-1)}.
  // where (Parent(n),l-1) is the parent cell of (n,l).
  double r_expansion_translated = 0.0;

  Point center_yBox = yBox.getCenter();

  std::vector<double> RCoeffTranslated = potential.getRR(center_yBoxParent.getCoord(), center_yBox.getCoord(),
		                                                   RCoeff);
  std::vector<double> RVecChildCenter = potential.getRVector(y[0].getCoord(), center_yBox.getCoord());


  for (unsigned int i=0; i < RVecChildCenter.size(); ++i)
    r_expansion_translated += RCoeffTranslated[i] * RVecChildCenter[i];

  std::cout << "r_expansion_translated = " << r_expansion_translated << std::endl;

*/




















































/*
#include<iostream>
#include<vector>
#include<cmath>

#include "Point.h"
#include "Util.h"
#include "Example1.h"
#include "Potential.h"
#include "FmmTree.h"
*/

/**
 *  test5.cc
 *  Created on: Dec 12, 2016
 *      Author: dbpc
 *
 *  The program test5.cc tests the Fast Multipole Method provided.
 *  The refinement level tested at is l = 3 ( l = 0 - no refinement, one cell).
 *  Below the code is a commented section that can be used to make sure that the code is
 *  working properly.  This section is based off of the workin test4.cc where the method has
 *  been broken its components - the S|S, S|R, and R|R transformations of Potential.cc
 *  The results of the commented section can be used to assemble the s and r-expansions for
 *  comparison with the results of this code.
 *
 *  Initially a single target point y[0] located in cell n = 0 was used to test the method.
 *  The points of three particular cells were used as source points for testing the code.
 *  The cells were n = 5, n = 45, and n = 365.  Each cell takes into account a different
 *  case for the Fast Multipole Method.  The points of cell n = 5 are near neighbors of the target
 *  point y[0] located in cell n = 0.  The source points of cell n = 45 are in the interactive list
 *  of y[0], and the source points of cell 365 are beyond the interactive list.
 *  The diagram of l = 3 below can be used to see the locations of these cells.
 *
 *  Using the cells n = 5, n = 45, and n = 365, it was possible to determine location of
 *  errors in the code.  For example, an error was determined to be occurring with the far-field
 *  calculations beyond the interaction list (for sources in cell 365 for example).  Tracing the
 *  code for this part of the FMM, it was foudn that the default p-value of 12 was being used by the class Box
 *  when its default constructor was used to create the fmmtree array of box objects.  The default
 *  constructor of class Box was then used to set the size of vectors c, d, and dtilde prior to the
 *  p-value of these instances being updated to the correct value used by this code - p = 20.  No segmentation
 *  fault run-time error was thrown for not having the correct vector sizes and the code compiled and ran
 *  producing erroneous results.  The error was fixed and a default p-value for box was changed to p = 20.
 *
 *  The results were also compared with maxima code test5.wxm.
 *
 *  The original 12 source points below in vector x that were used for debuggin were
 *
 *    x[0] = (0.03125, 0.03125, 0.03125)
 *    x[1] = (0.09375, 0.03125, 0.03125)
 *    x[2] = (0.03125, 0.03125, 0.09375)
 *    x[3] = (0.09375, 0.03125, 0.09375)
 *
 *    located in cell n = 0 of refinement level l = 3 (above)
 *
 *    x[4] = (0.15625, 0.03125, 0.15625)
 *    x[5] = (0.21875, 0.03125, 0.15625)
 *    x[6] = (0.15625, 0.03125, 0.21875)
 *    x[7] = (0.21875, 0.03125, 0.21875)
 *
 *    located in cell n = 5 of refinement level l = 3
 *
 *    x[8]  = (0.40625, 0.03125, 0.40625)
 *    x[9]  = (0.46875, 0.03125, 0.40625)
 *    x[10] = (0.40625, 0.03125, 0.46875)
 *    x[11] = (0.46875, 0.03125, 0.46875)
 *
 *    located in cell n = 45 of refinement level l = 3
 *
 *    x[12]  = (0.90625, 0.03125, 0.90625)
 *    x[13]  = (0.96875, 0.03125, 0.90625)
 *    x[14] = (0.90625, 0.03125, 0.96875)
 *    x[15] = (0.96875, 0.03125, 0.96875)
 *
 *    located in cell n = 365 of refinement level l = 3
 *
 *  A single target point below in vector y was tested
 *
 *    y[0] = (0.03125, 0.03125, 0.03125)
 *
 *  and all charge values in vector u corresponding to source points
 *  have value u[i] = 1.
 *
 *  At this point we can note where the source points x are with respect
 *  to the target y.
 *
 *  For this test program we use almost the same source and target particles
 *  as used in the 2D FMM test examples for refinement level l = 3 (below).
 *  Namely, we have 4 source points per cell and all the points are located in a plane
 *  parallel to the xz-coordinate plane.
 *
 *  The difference between this example and a 2D example is the addition of a
 *  third coordinate.  That is, the particles used in a 2D FMM example are
 *  the same except for the shift of the particles in the third coordinate
 *  direction.  We take the y-coordinate as the third coordinate direction.
 *  In contrast to the 2D examples, the particles do not lie in the xy-plane,
 *  but lie in a plane parallel to the xz-plane, shifted off the xz-plane into
 *  the first octant by the addition of the third coordinate direction - the y-coordinate .
 *  That is, the points for this example will have the same x-coordinate as the points in
 *  the 2D example, and the 2D points' y-coordinate will be the z-coordinate for this example.
 *  The y-coordinate is fixed for all points at y = 0.03125 - quarter_length of a cell.
 *
 *  The results of this code should be comparable to a 2D code since the particles are
 *  generate in the same order and have the same locations with respect to one another
 *  as in a 2D code.  We would like to make comparisons with the FMM 2D code, however, the
 *  2D code potential function is 1/(x-y) instead of 1/|x-y|.
 *  We can alter the 2D code if necessary, but new formulas for the translation operators
 *  would have to be determined.
 *
 *  Looking at the refinements below, we can now see where the sources are located with respect
 *  to the target.  The source point is located in cell n = 0.  The 8 particles are located in
 *  near neighbors of cell n = 0 in cell n = 0 and cell n = 5.  4 particles are located in cell
 *  n = 45 in the interactive list.  4 particles are located outside of the interactive list in
 *  the upper right front corner of the domain.
 *
 *  NEAR NEIGHBORS
 *
 *  The first part of the FMM code to check is the direct calculations of the potential for
 *  source points in near neighbors
 *  cells n = 0 and n = 5.
 *
 *  Recall the potential function is phi(t1,t2,t3) = 1 / sqrt(t1^2 + t2^2 + t3^2)
 *
 *  For target y[0] = (0.03125, 0.03125, 0.03125)
 *
 *      for source x[0] = (0.03125, 0.03125, 0.03125)
 *
 *          t = x[0] - y[0] = (0, 0, 0) and phi is undefined (singularity)
 *          and we do not include this point (do not include action of particle on itself)
 *
 *      for source x[1] = (0.09375, 0.03125, 0.03125)
 *
 *          t = x[1] - y[0] = (0.0625, 0.0000, 0.0000) and
 *          phi(t) = 1 / sqrt( (0.0625)^2 + (0.0000)^2 + (0.0000)^2 )
 *                 = 1 / sqrt( (1 / 16)^2 ) = 1 / (1 / 16) = 16
 *
 *      for source x[2] = (0.03125, 0.03125, 0.09375)
 *
 *          t = x[2] - y[0] = (0.0000, 0.0000, 0.0625)
 *          phi(t) = 1 / sqrt( (0.0000)^2 + (0.0000)^2 + (0.0625)^2 )
 *                 = 1 / sqrt( (1 / 16)^2 ) = 1 / (1 / 16) = 16
 *
 *      for source x[3] = (0.09375, 0.03125, 0.09375)
 *
 *          t = x[3] - y[0] = (0.0625, 0.0000, 0.0625)
 *          phi(t) = 1 / sqrt( (0.0625)^2 + (0.0000)^2 + (0.0625)^2 )
 *                 = 1 / sqrt( (1 / 16)^2 + (1 / 16)^2 ) = 1 / sqrt(2 / 256)
 *                 = 8 sqrt(2) ~ 11.313708
 *
 *      for source x[4] = (0.15625, 0.03125, 0.15625)
 *          phi(t) = 4 sqrt(2) ~ 5.656854
 *
 *      for source x[5] = (0.21875, 0.03125, 0.15625)
 *          phi(t) = 1 / sqrt( (3/16)^2 + 0^2 + (1/8)^2 )
 *                 = 1 / sqrt( 13 / 256 )
 *                 = 16 sqrt(13) / 13 ~ 4.437602
 *
 *      for source x[6] = (0.15625, 0.03125, 0.21875)
 *          phi(t) = 16 sqrt(13) / 13 ~ 4.437602
 *
 *      for source x[7] = (0.21875, 0.03125, 0.21875)
 *          phi(t) = 1 / sqrt( (3/16)^2 + 0^2 + (3/16)^2 )
 *                 = 1 / sqrt( 18 / 256 )
 *                 = 16 / (3 sqrt(2)) = 8 sqrt(2) / 3
 *                 ~ 3.771236
 *
 *  The sum of the potentials is
 *      Sum = 16 + 16 + 11.313708 + 5.656854 + 4.437602 + 4.437602 + 3.771236
 *          = 61.617002
 *
 *  Running the FMM with just these source points (cells n = 0 and n = 5 below) checks out
 *  and we get back 61.617.
 *
 *  We are now ready to check the interactive list.
 *
 *  INTERACTION LIST
 *
 *  For target y[0] = (0.03125, 0.03125, 0.03125) we have four source points
 *  in cell n = 45.
 *
 *    x[8]  = (0.40625, 0.03125, 0.40625)
 *    x[9]  = (0.46875, 0.03125, 0.40625)
 *    x[10] = (0.40625, 0.03125, 0.46875)
 *    x[11] = (0.46875, 0.03125, 0.46875)
 *
 *  We can see in the figure below that cell n = 45 is in the interaction list of cell n = 0.
 *
 *
 *
  * level l = 3 refinement
  *                                                                                                                                       y = 1.00
  *                                                                                                                             1.00 _______ _______ _______ _______
  *                                                                                                                                 |219|223|251|255|475|479|507|511|
  *                                                                                                 y = 0.875                       |---27--|---31--|---59--|---63--|
  *                                                                                                                             0.75|218|2221250|254|474|4785506|510|
  *                                                                                       1.00 _______ _______ _______ _______      |211|215|243|247|467|471|499|503|
  *               |                                                                           |217|221|249|253|473|477|505|509|     |---26--|---30--|---58--|---62--|
  *                                                           y = 0.75                        |---27--|---31--|---59--|---63--| 0.50|210|214|242|246|466|470|498|502|
  *                                                                                       0.75|216|2201248|252|472|4765504|508|     |155|159|187|191|411|415|443|447|
  *                                                 1.00 _______ _______ _______ _______      |209|213|241|245|465|469|497|501|     |---19--|---23--|---51--|---55--|
  *               z                                     |203|207|235|239|459|463|491|495|     |---26--|---30--|---58--|---62--| 0.25|154|1580186|190|410|4144442|446|
  *               |      y = 0.625                      |---25--|---29--|---57--|---61--| 0.50|208|212|240|244|464|468|496|500|     |147|151|179|183|403|407|435|439|
  *               |                                 0.75|202|2061234|238|458|4625490|494|     |153|157|185|189|409|413|441|445|     |---18--|---22--|---50--|---54--|
  *          1.00 |_______________________________      |195|199|227|231|451|455|483|487|     |---19--|---23--|---51--|---55--|     |146|150|178|182|402|406|434|438|
  *               |201|205|233|237|457|461|489|493|     |---24--|---28--|---56--|---60--| 0.25|152|1560184|188|408|4124440|444|    0.0     0.25    0.5    0 .75    1.0       x
  *               |---25--|---29--|---57--|---61--| 0.50|194|198|226|230|450|454|482|486|     |145|149|177|181|401|405|433|437|
  *          0.75 |200|2041232|236|456|4605488|492|     |139|143|171|175|395|399|427|431|     |---18--|---22--|---50--|---54--|
  *               |193|197|225|229|449|453|481|485|     |---17--|---21--|---49--|---53--|     |144|148|176|180|400|404|432|436|
  *               |---24--|---28--|---56--|---60--| 0.25|138|1420170|174|394|3984426|430|    0.0     0.25    0.5     0.75    1.0       x
  *          0.50 |192|196|224|228|448|452|480|484|     |131|135|163|167|387|391|419|423|
  *               |137|141|169|173|393|397|425|429|     |---16--|---20--|---48--|---52--|
  *               |---17--|---21--|---49--|---53--|     |130|134|162|166|386|390|418|422|
  *          0.25 |136|1400168|172|392|3964424|428|    0.0     0.25    0.5     0.75    1.0       x
  *               |129|133|161|165|385|389|417|421|
  *               |---16--|---20--|---48--|---52--|
  *               |128|132|160|164|384|388|416|420|______
  *              0.0     0.25    0.5     0.75    1.0       x
  *
  *
  *                                                                                                                               y = 0.50
  *                                                                                                                     1.00 _______ _______ _______ _______
  *                                                                                                                         |91 |95 |123|127|347|351|379|383|
  *                                                                                         y = 0.375                       |---11--|---15--|---43--|---47--|
  *                                                                                                                     0.75|90_|94_1122|126|346|3505378|382|
  *                                                                               1.00 _______ _______ _______ _______      |83 |87 |115|119|339|343|371|375|
  *       |                                                                           |89 |93 |121|125|345|349|377|381|     |---10--|---14--|---42--|---46--|
  *                                                   y = 0.25                        |---11--|---15--|---43--|---47--| 0.50|82_|86_|114|118|338|342|370|374|
  *                                                                               0.75|88_|92_1120|124|344|3485376|380|     |27 |31 | 59|63 |283|287|315|319|
  *                                         1.00 _______ _______ _______ _______      |81 |85 |113|117|337|341|369|373|     |---3---|---7---|---35--|---39--|
  *      z                                      |75 |79 |107|111|331|335|363|367|     |---10--|---14--|---42--|---46--| 0.25|26_|30_0_58|62_|282|2864314|318|
  *       |      y = 0.125                      |---9---|---13--|---41--|---45--| 0.50|80_|84_|112|116|336|340|368|372|     |19 |23 |51 |55 |275|279|307|311|
  *       |                                 0.75|74_|78_1106|110|330|3345362|366|     |25 |29 | 57| 61|281|285|313|317|     |---2---|---6---|---34--|---38--|
  *  1.00 |_______________________________      |67 |71 |99 |103|323|327|355|359|     |---3---|---7---|---35--|---39--|     |18_|22_|50_|54_|274|278|306|310|
  *       |73 |77 |105|109|329|333|361|365|     |---8---|---12--|---40--|---44--| 0.25|24_|28_0_56|_60|280|2844312|316|    0.0     0.25    0.5    0 .75    1.0       x
  *       |---9---|---13--|---41--|---45--| 0.50|66_|70_|98_|102|322|326|354|358|     |17 |21 |49 |53 |273|277|305|309|
  *  0.75 |72_|76_1104|108|328|3325360|364|     |11 |15 | 43|47 |267|271|299|303|     |---2---|---6---|---34--|---38--|
  *       |65 |69 |97 |101|321|325|353|357|     |---1---|---5---|---33--|---37--|     |16_|20_|48_|52_|272|276|304|308|
  *       |---8---|---12--|---40--|---44--| 0.25|10_|14_0_42|46_|266|2704298|302|    0.0     0.25    0.5     0.75    1.0       x
  *  0.50 |64_|68_|96_|100|320|324|352|356|     | 3 | 7 |35 |39 |259|263|291|295|
  *       | 9 |13 | 41|45 |265|269|297|301|     |---0---|---4---|---32--|---36--|
  *       |---1---|---5---|---33--|---37--|     |_2_|_6_|34_|38_|258|262|290|294|
  *  0.25 |_8_|12_0_40|44_|264|2684296|300|    0.0     0.25    0.5     0.75    1.0       x
  *       | 1 | 5 |33 |37 |257|261|289|293|
  *       |---0---|---4---|---32--|---36--|
  *       |_0_|_4_|32_|36_|256|260|288|292|______
  *      0.0     0.25    0.5     0.75    1.0       x
  *
 *
 *  For each source point in cell n = 45, we will i) obtain a far-field series (s-expansion)
 *  centered at the center of cell n = 45, ii) perform an SR translation of the series to the center
 *  of cell n = 0, iii) obtain the r-expansion for the potential of the source x[i] acting on target y[0],
 *  iv) compare the approximation for the potential of source acting on the target to the direct calculation
 *  of the potential function
 *
 *  The far field series (s-expansion)  psi(y-x) = sum (c_n S_n) = sum_{alpha} frac{ D^{alpha} f(y - x_c) }{ alpha! } (x_c - x)^{alpha}
 *  The coefficients are (x_c - x)^{alpha} / alpha!  and the powers are D^{alpha} f(y - x_c)
 *  The near field series (r-expansion)  psi(y-x) = sum (c_n S_n) = sum_{alpha} frac{ D^{alpha} f(x_c - x) }{ alpha! } (y - x_c)^{alpha}
 *  The coefficients are D^{alpha} f(x_c - x) / alpha!  and the powers are (y - x_c)^{alpha}
 *
 *  We show results for the first source point in n = 45 -  x[8]  = (0.40625, 0.03125, 0.40625)
 *
 *
 *  For target y[0] = (0.03125, 0.03125, 0.03125)
 *    and source x[8] = (0.40625, 0.03125, 0.40625)
 *
 *    The center of cell n = 0 containing y[0] is center_cell_0 = (0.0625,0.0625,0.0625)
 *    The center of cell n = 45 containing x[0] is center_cell_45 = (0.4375,0.0625,0.4375)
 *
 *    The first step is to collect the coefficients of the far-field expansion
 *
 *    		psi(x-y)  =  sum_{alpha >= 0} frac{ (x - x_*)^{alpha} }{ alpha! }  D^{alpha} psi(x_* - y)
 *
 *    where the coefficients frac{ (x - x_*)^{alpha} }{ alpha! } are collected
 *
 *    The next step is to perform the SR transformation from the center x_* of the interaction list cell n = 45 to the center of
 *    the cell n = 0 having the target y[0].  The far-field series
 *
 *			psi(x-y)  =  sum_{alpha >= 0} frac{ (x - x_*)^{alpha} }{ alpha! }  D^{alpha} psi(x_* - y)
 *
 *	  is transformed to the near-field series
 *
 *			psi(x-y)  =  sum_{alpha >= 0} frac{ D^{alpha} psi(x - y_*) }{ alpha! } (y_* - y)^{alpha}
 *
 *    Note that we want the near-field expansion in this form with powers (y_* - y)^{alpha} since we do not have to specify
 *    the target point y until the very end of the Fast Multipole Method.  The advantage of this strategy is only having to
 *    specify the source x and center's x_* and y_* as we start to approach the cell (possibly down a hierarchy of cells)
 *    containing the target.  Hence, the far-field expansions in the upward pass are general and for all targets far enough
 *    from the source cell.
 *
 *    The result of the SR transformation applied to the coefficients (SCoeff vector) of the far field series is the vector of
 *    near field coefficients RCoeff (approximately, since SR transformation is approximate) having the form
 *    a_n = frac{ D^{alpha} psi(x - y_*) }{ alpha! }.  Multiplying the coefficients by the near-field powers R_n = (y_* - y)^{alpha}
 *    we obtain the approximation
 *
 *                           psi(x-y) = 1.88562...
 *                                    = sum  a_n R_n
 *                                    ~  1.88571,
 *
 *    where 1.88562... was obtain from the direct calculation psi(x-y) = 1/|x-y|.
 *    Error was 0.0000913489
 *
 *    CELLS BEYOND THE INTERACTION LIST
 *
 *    Looking at the diagram above, the cell n = 365 is beyond the interaction list of the cell n = 0.
 *    Its parent, however, is in the interaction list of the parent of cell n = 0.  Therefore, we will treat
 *    the source in cell n = 365 at level l = 3 at the parent level l = 2, since the source falls into the
 *    interaction list of a parent cell containing the target in n = 0 at level l = 3.
 *
 *    For a source point in cell n = 365 at level l = 3, we will i) obtain a far-field series (s-expansion)
 *    centered at the center of cell n = 365 at l = 3, ii) perform an S|S translation from the child cell n = 365 at l = 3
 *    to the parent cell n = 45 at level l = 2, iii) perform an SR translation of the series from the parent center
 *    to the center of the parent of cell n = 0 at level l = 3, iv) translate the r-expansion for the potential of the
 *    source x[12] acting on target y[0] from the child center y_{*}^{child} to the parent center y_{*}^{parent}
 *    v) form the r-expansion approximating the potential function calculation 1/|x-y|
 *    vi) compare the approximation for the potential of source acting on the target to the direct calculation
 *    of the potential function
 *
 *    For target y[0] = (0.03125, 0.03125, 0.03125)
 *      and source x[12] = (0.90625, 0.03125, 0.90625)
 *
 *      The first step is the far-field expansion centered at center_cell_365
 *        The center of cell n = 0 containing y[0] is center_cell_0 = (0.0625,0.0625,0.0625)
 *        The center of cell n = 365 containing x[12] is center_cell_365 = (0.9375,0.0625,0.9375)
 *
 *      The second step is the S|S translation to the parent center (0.875, 0.125, 0.875) of n = 365
 *
 *      The third step is the S|R translation from the parent center of n = 365 to the parent center
 *      (0.125,0.125,0.125) of cell n = 0 ( at level l = 3).
 *
 *      The fourth step is the R|R translation of the near-field expansion of step 3 from the parent
 *      center of cell n = 0 to the (child) center of cell n = 0.  The result is
 *
 *                           psi(x-y) = 0.808122...
 *                                    = sum  a_n R_n
 *                                    ~  0.808235,
 *
 *      where 1.88562... was obtain from the direct calculation psi(x-y) = 1/|x-y|.
 *
 *      Error is 0.000112701
 *
 *
 */



/*
int main()
{
  // note that these refinement levels counts start on l = 1
  int DEFAULT_NUM_LEVEL = 4;                   // default refinement level
  // p represents the number of derivatives to approximate Taylor series to order |alpha|
  // The derivatives at order |alpha| = n are the combinations (n+2)_C_2 = ((n +2) \\ 2) = ((n+2)*(n+1))/(2*1)
  // We set p = 20 to approximate series to the third derivatives (1 + 3 + 6 + 10 = 20)
  //                                                               0th 1st 2nd 3rd  derivatives
  unsigned int abs_alpha = 3;
  // p is the number of terms used in Taylor Series approximation
  // Setting p to 1 and adding number of derivatives at each order below
  // using (n+2)_C_2 combinations formula.
  int p = 1;
  int combination = (2 * 1) / (2 * 1); // number of derivatives for order |alpha| = 0
  for (unsigned int i = 1; i <= abs_alpha; ++i)
  {
    combination = (i+2)*(i+1) / 2;  // number derivatives for order |alpha| = i
    p = p + combination;
  }

  if (p == 20)
    std::cout << "p = " << p << " correct" << std::endl;
  else // p \neq 20
    std::cout << "p != " << 20 << "incorrect" << std::endl;

//  Example1 example1(DEFAULT_NUM_LEVEL);


  double cells_per_side = pow(2.0,DEFAULT_NUM_LEVEL-1);
  double cell_length = (1.0-0.0)/cells_per_side;
  double quarter_length = 0.25*cell_length;
  double three_quarter_length = 0.75*cell_length;

  std::vector<Point>  x, y;
  std::vector<double> u;


  // The same number of particles as in the 2D FMM case are used and there are
  // 4 particles per cell with pow(4,L-1) cells partitioning the domain
//  unsigned int nTotalParticles = 4 * pow(4, DEFAULT_NUM_LEVEL - 1);

//  x.resize(nTotalParticles);
//  y.resize(nTotalParticles);
//  u.resize(nTotalParticles);
  x.resize(1);
  y.resize(1);
  u.resize(1);


  Point upper_right_corner(0.9375+quarter_length,quarter_length,0.9375+quarter_length);
  unsigned int cell_index_upper_right_corner = upper_right_corner.getBoxIndex(DEFAULT_NUM_LEVEL-1);
  std::cout << "point " << upper_right_corner.coordToString() << " is in box " << cell_index_upper_right_corner
		    << std::endl;
  // creating source points
  // note the points are the same as the 2D example, except
  // that we are using xz-coordinates for xy-coordinates here.  Therefore,
  // the extra y-coordinate fixed at y = quarter_length puts the points in the
  // 1st octant above the xz-plane.  Using the diagram above, we see the inner loop below
  // increments along the x-axis indicating that the order of the cells w.r.t. the x and y vectors
  // containing the points will be
  //  0,  4,  32,  36, 256, 260, 288, 292,
  //  1,  5,  33,  37, 257, 261, 289, 293,    (working doown the x-axis and incrementing up the y-axis)
  //  8, 12,  40,  44, 264, 268, 296, 300,
  //  9, 13,  41,  45, 265, 269, 297, 301,
  // 64, 68,  96, 100, 320, 324, 352, 356,
  // 65, 69,  97, 101, 321, 325, 353, 357,
  // 72, 76, 104, 108, 328, 332, 360, 364,
  // 73, 77, 105, 109, 329, 333, 361, 365

  // index for points that will go into x (source) and y (targets) vectors
  unsigned int n_points = 0;
  unsigned int flag_cell = 45;

  for (double y_coord=0.0; y_coord<1.0; y_coord+=cell_length)
	for (double  x_coord=0.0; x_coord<1.0; x_coord+=cell_length)
	{
	  Point flag_point(x_coord+quarter_length,quarter_length,y_coord+quarter_length);
	  unsigned int cell_index = flag_point.getBoxIndex(DEFAULT_NUM_LEVEL-1);
//      if(cell_index == 0 || cell_index == 5 || cell_index == 45 || cell_index == 365)
      if(cell_index == flag_cell)
      {
	    x[n_points].setX(x_coord + quarter_length);       // source particle (point)
	    x[n_points].setZ(y_coord + quarter_length);       // lower left corner of cell
	    x[n_points].setY(quarter_length);
    	std::cout << "x[" << n_points << "] = " << x[n_points].coordToString() << std::endl;
//	    ++n_points;
//	    x[n_points].setX(x_coord + three_quarter_length); // source particle (point)
//	    x[n_points].setZ(y_coord + quarter_length);       // lower right corner of cell
//	    x[n_points].setY(quarter_length);
//    	std::cout << "x[" << n_points << "] = " << x[n_points].coordToString() << std::endl;
//	    ++n_points;
//	    x[n_points].setX(x_coord + quarter_length);       // source particle (point)
//	    x[n_points].setZ(y_coord + three_quarter_length); // upper left corner of cell
//	    x[n_points].setY(quarter_length);
//    	std::cout << "x[" << n_points << "] = " << x[n_points].coordToString() << std::endl;
//	    ++n_points;
//	    x[n_points].setX(x_coord + three_quarter_length); // source particle (point)
//	    x[n_points].setZ(y_coord + three_quarter_length); // upper right corner of cell
//	    x[n_points].setY(quarter_length);
//    	std::cout << "x[" << n_points << "] = " << x[n_points].coordToString() << std::endl;
//	    ++n_points;
      }
	}

  // creating target points
  // for this test only one target point
  // note the target points are same as source points
  n_points = 0;

  for (double y_coord=0.0; y_coord<1.0; y_coord+=cell_length)
	for (double  x_coord=0.0; x_coord<1.0; x_coord+=cell_length)
	{
	  Point flag_point_2(x_coord+quarter_length,quarter_length,y_coord+quarter_length);
      unsigned int cell_index_2 = flag_point_2.getBoxIndex(DEFAULT_NUM_LEVEL-1);
	  if(cell_index_2 == 0)
	  {
	    y[n_points].setX(x_coord + quarter_length);       // target particle (point)
	    y[n_points].setZ(y_coord + quarter_length);       // lower left corner of cell
	    y[n_points].setY(quarter_length);
        std::cout << "y[" << n_points << "] = " << y[n_points].coordToString() << std::endl;
//	    ++n_points;
//	    y[n_points].setX(x_coord + three_quarter_length); // source particle (point)
//	    y[n_points].setZ(y_coord + quarter_length);       // lower right corner of cell
//	    y[n_points].setY(quarter_length);
//    	std::cout << "x[" << n_points << "] = " << x[n_points].coordToString() << std::endl;
//	    ++n_points;
//	    y[n_points].setX(x_coord + quarter_length);       // source particle (point)
//	    y[n_points].setZ(y_coord + three_quarter_length); // upper left corner of cell
//	    y[n_points].setY(quarter_length);
//    	std::cout << "x[" << n_points << "] = " << x[n_points].coordToString() << std::endl;
//	    ++n_points;
//	    y[n_points].setX(x_coord + three_quarter_length); // source particle (point)
//	    y[n_points].setZ(y_coord + three_quarter_length); // upper right corner of cell
//	    y[n_points].setY(quarter_length);
//    	std::cout << "x[" << n_points << "] = " << x[n_points].coordToString() << std::endl;
//	    ++n_points;
	  }
	}

    int xBoxIndex = x[0].getBoxIndex(DEFAULT_NUM_LEVEL-1);
    int yBoxIndex = y[0].getBoxIndex(DEFAULT_NUM_LEVEL-1);
    Box xBox(DEFAULT_NUM_LEVEL-1, xBoxIndex, p);
    Box yBox(DEFAULT_NUM_LEVEL-1, yBoxIndex, p);

    Point center_xBoxChild = xBox.getCenter();
    Point center_yBoxChild = yBox.getCenter();

    std::cout << "yBoxChild cell n = " << yBoxIndex << " at level l = " << yBox.getLevel()
    		  << " has center = " << center_yBoxChild.coordToString() << std::endl;
    std::cout << "xBoxChild cell n = " << xBoxIndex << " at level l = " << xBox.getLevel()
    		  << " has center = " << center_xBoxChild.coordToString() << std::endl;

    if (flag_cell == 365)
    {
      int xBoxParentIndex = x[0].getBoxIndex(DEFAULT_NUM_LEVEL-2);
      int yBoxParentIndex = y[0].getBoxIndex(DEFAULT_NUM_LEVEL-2);
      Box xBoxParent(DEFAULT_NUM_LEVEL-2, xBoxParentIndex, p);
      Box yBoxParent(DEFAULT_NUM_LEVEL-2, yBoxParentIndex, p);
      Point center_xBoxParent = xBoxParent.getCenter();
      Point center_yBoxParent = yBoxParent.getCenter();
      std::cout << "yBoxParent cell n = " << yBoxParentIndex << " at level l = " << yBoxParent.getLevel()
    		    << " has center = " << center_yBoxParent.coordToString() << std::endl;
      std::cout << "xBoxParent cell n = " << xBoxParentIndex << " at level l = " << xBoxParent.getLevel()
    		    << " has center = " << center_xBoxParent.coordToString() << std::endl;
    }




  for (unsigned int i=0; i<u.size(); ++i)
	  u[i] = 1.0;


  Potential potential(p);

  FmmTree fmmtree(DEFAULT_NUM_LEVEL,x,y,potential);


  std::vector<double> indirect = fmmtree.solve(u);
  std::vector<double> direct = fmmtree.solveDirect(u);

  for (unsigned int i=0; i<direct.size(); ++i)
  {
    std::cout << "direct[" << i << "] = " << direct[i] << " versus "
    		  << "indirect[" << i << "] = " << indirect[i] << "\n";
  }

  double error = 0.0;
  for (unsigned int i = 0; i<direct.size(); ++i)
  {
	  double tmp = std::abs(direct[i]-indirect[i]);
	  if (tmp>error)
        error = tmp;
  }

  std::cout << "Error = " << error << "\n";


  std::cout << "Finished" << "\n";


*/

  /*
   * The code below was used to check the FMM results above as well as troubleshoot
   * the FMM code.  It breaks down the FMM step by step
   *
    int xBoxIndex = x[0].getBoxIndex(DEFAULT_NUM_LEVEL-1);
    int yBoxIndex = y[0].getBoxIndex(DEFAULT_NUM_LEVEL-1);
    Box xBox(DEFAULT_NUM_LEVEL-1, xBoxIndex, p);
    Box yBox(DEFAULT_NUM_LEVEL-1, yBoxIndex, p);
    int xBoxParentIndex = x[0].getBoxIndex(DEFAULT_NUM_LEVEL-2);
    int yBoxParentIndex = y[0].getBoxIndex(DEFAULT_NUM_LEVEL-2);
    Box xBoxParent(DEFAULT_NUM_LEVEL-2, xBoxParentIndex, p);
    Box yBoxParent(DEFAULT_NUM_LEVEL-2, yBoxParentIndex, p);

    Point center_xBoxChild = xBox.getCenter();
    Point center_yBoxChild = yBox.getCenter();
    Point center_xBoxParent = xBoxParent.getCenter();
    Point center_yBoxParent = yBoxParent.getCenter();

    std::cout << "yBoxChild center cell n = 0 at level l = 3 is center_yBoxChild = " << center_yBoxChild.coordToString() << std::endl;
    std::cout << "xBoxChild center cell n = 365 at level l = 3 is center_xBoxChild =  = " << center_xBoxChild.coordToString() << std::endl;
    std::cout << "yBoxParent center cell n = 0 at level l = 2 is center_yBoxParent = " << center_yBoxParent.coordToString() << std::endl;
    std::cout << "xBoxParnet center cell n = 45 at level l = 2 is center_xBoxParent = " << center_xBoxParent.coordToString() << std::endl;

    std::vector<double> SCoeffChild = potential.getSCoeff(x[0].getCoord(), center_xBoxChild.getCoord());
    for (unsigned int i=0; i<SCoeffChild.size(); ++i)
      std::cout << "SCoeffChild[" << i << "] = " << SCoeffChild[i] << std::endl;
    std::cout << std::endl << std::endl;

    std::vector<double> SVecChild = potential.getSVector(y[0].getCoord(), center_xBoxChild.getCoord());
    std::vector<double> SCoeffParent = potential.getSS(center_xBoxChild.getCoord(), center_xBoxParent.getCoord(),
  		                                             SCoeffChild);
    for (unsigned int i=0; i<SCoeffParent.size(); ++i)
      std::cout << "SCoeffParent[" << i << "] = " << SCoeffParent[i] << std::endl;
    std::cout << std::endl << std::endl;


    std::vector<double> SVecParent = potential.getSVector(y[0].getCoord(), center_xBoxParent.getCoord());


    double s_expansion_child = 0.0;  double s_expansion_parent = 0.0;
    for (unsigned int i=0; i<SVecChild.size(); ++i)
    {
      s_expansion_child += SCoeffChild[i] * SVecChild[i];
      s_expansion_parent += SCoeffParent[i] * SVecParent[i];
    }

    std::cout << "s_expansion_child = " << s_expansion_child << std::endl;
    std::cout << "s_expansion_parent = " << s_expansion_parent << std::endl;



    std::vector<double> RCoeffParent = potential.getSR(center_xBoxParent.getCoord(), center_yBoxParent.getCoord(),
  		                                                     SCoeffParent);
    std::vector<double> RVecParent = potential.getRVector(y[0].getCoord(), center_yBoxParent.getCoord());
    std::vector<double> RCoeffChild = potential.getRR(center_yBoxParent.getCoord(), center_yBoxChild.getCoord(),
  		                                            RCoeffParent);
    std::vector<double> RVecChild = potential.getRVector(y[0].getCoord(), center_yBoxChild.getCoord());

    double r_expansion_parent = 0.0; double r_expansion_child = 0.0;
    for (unsigned int i=0; i<RVecParent.size(); ++i)
    {
      r_expansion_parent += RCoeffParent[i] * RVecParent[i];
      r_expansion_child += RCoeffChild[i] * RVecChild[i];
    }

    std::cout << "r_expansion_parent = " << r_expansion_parent << std::endl;
    std::cout << "r_expansion_child = " << r_expansion_child << std::endl;


    std::cout << "direct calculation = " << potential.direct(y[0].getCoord(), x[0].getCoord()) << std::endl;

  */

/*
}
*/

/*
 * Recall the FMM refinements 0 <= l <= 4
 * n - near neighbor
 * i - interaction list cell
 * c - cell
 *
 *     l = 0
 *      _______________________________________________________________________________
 *     |                                                                               |
 *     |                                                                               |
 *     |                                                                               |
 *     |                                                                               |
 *     |                                                                               |
 *     |                                                                               |
 *     |                                                                               |
 *     |                                                                               |
 *     |                                                                               |
 *     |                                                                               |
 *     |                                                                               |
 *     |                                                                               |
 *     |                                                                               |
 *     |                                                                               |
 *     |                                                                               |
 *     |                                                                               |
 *     |                                                                               |
 *     |                                                                               |
 *     |                                                                               |
 *     |                                                                               |
 *     |                                                                               |
 *     |                                                                               |
 *     |                                                                               |
 *     |                                                                               |
 *     |                                                                               |
 *     |                                                                               |
 *     |                                                                               |
 *     |                                                                               |
 *     |                                                                               |
 *     |                                                                               |
 *     |                                                                               |
 *     |                                                                               |
 *     |_______________________________________________________________________________|
 *
 *
 *     l = 1
 *      _______________________________________________________________________________
 *     |                                       =                                       |
 *     |                                       =                                       |
 *     |                                       =                                       |
 *     |                                       =                                       |
 *     |                                       =                                       |
 *     |                                       =                                       |
 *     |                                       =                                       |
 *     |                                       =                                       |
 *     |                                       =                                       |
 *     |                                       =                                       |
 *     |                                       =                                       |
 *     |                                       =                                       |
 *     |                                       =                                       |
 *     |                                       =                                       |
 *     |                                       =                                       |
 *     |                                       =                                       |
 *     |=========|===================|=================================================|
 *     |                                       =                                       |
 *     |                                       =                                       |
 *     |                                       =                                       |
 *     |                                       =                                       |
 *     |                                       =                                       |
 *     |                                       =                                       |
 *     |                                       =                                       |
 *     |                                       =                                       |
 *     |                                       =                                       |
 *     |                                       =                                       |
 *     |                                       =                                       |
 *     |                                       =                                       |
 *     |                                       =                                       |
 *     |                                       =                                       |
 *     |                                       =                                       |
 *     |_______________________________________________________________________________|
 *
 *     l = 2
 *      _______________________________________________________________________________
 *     |                   *                   =                   *                   |
 *     |                   *                   =                   *                   |
 *     |                   *                   =                   *                   |
 *     |        i          *         i         =         i         *         i         |
 *     |                   *                   =                   *                   |
 *     |                   *                   =                   *                   |
 *     |                   *                   =                   *                   |
 *     |***************************************=***************************************|
 *     |                   *                   =                   *                   |
 *     |                   *                   =                   *                   |
 *     |                   *                   =                   *                   |
 *     |         n         *         n         =         n         *         i         |
 *     |                   *                   =                   *                   |
 *     |                   *                   =                   *                   |
 *     |                   *                   =                   *                   |
 *     |                   *                   =                   *                   |
 *     |=========|===================|=================================================|
 *     |                   *                   =                   *                   |
 *     |                   *                   =                   *                   |
 *     |                   *                   =                   *                   |
 *     |         n         *         c         =         n         *         i         |
 *     |                   *                   =                   *                   |
 *     |                   *                   =                   *                   |
 *     |                   *                   =                   *                   |
 *     |***************************************=***************************************|
 *     |                   *                   =                   *                   |
 *     |                   *                   =                   *                   |
 *     |                   *                   =                   *                   |
 *     |         n         *         n         =         n         *         i         |
 *     |                   *                   =                   *                   |
 *     |                   *                   =                   *                   |
 *     |                   *                   =                   *                   |
 *     |_______________________________________________________________________________|
 *
 *
 *
 *     l = 3
 *      _______________________________________________________________________________
 *     |         ,         *         ,         =         ,         *         ,         |
 *     |         ,         *         ,         =         ,         *         ,         |
 *     |         ,         *         ,         =         ,         *         ,         |
 *     |,,,,,,,,,,,,,,,,,,,*,,,,,,,,,,,,,,,,,,,=,,,,,,,,,,,,,,,,,,,*,,,,|,,,,,,,,,,,,,,|
 *     |         ,         *         ,         =         ,         *         ,         |
 *     |         ,         *         ,         =         ,         *         ,         |
 *     |         ,         *         ,         =         ,         *         ,         |
 *     |***************************************=***************************************|
 *     |         ,         *         ,         =         ,         *         ,         |
 *     |    i    ,    i    *    i    ,    i    =    i    ,    i    *         ,         |
 *     |         ,         *         ,         =         ,         *         ,         |
 *     |,,,,,,,,,,,,,,,,,,,*,,,,,,,,,,,,,,,,,,,=,,,,,,,,,,,,,,,,,,,*,,,,,,,,,,,,,,,,,,,|
 *  -  |         ,         *         ,         =         ,         *         ,         |
 *  -  |    i    ,    n    *    n    ,    n    =    i    ,    i    *         ,         |
 *  -  |         ,         *         ,         =         ,         *         ,         |
 *  -  |         ,         *         ,         =         ,         *         ,         |
 *     |=========|===================|=================================================|
 *  -  |         ,         *         ,         =         ,         *         ,         |
 *  -  |    i    ,    n    *    c    ,    n    =    i    ,    i    *         ,         |
 *  -  |         ,         *         ,         =         ,         *         ,         |
 *  -  |,,,,,,,,,,,,,,,,,,,*,,,,,,,,,,,,,,,,,,,=,,,,,,,,,,,,,,,,,,,*,,,,,,,,,,,,,,,,,,,|
 *     |         ,         *         ,         =         ,         *         ,         |
 *     |    i    ,    n    *    n    ,    n    =    i    ,    i    *         ,         |
 *     |         ,         *         ,         =         ,         *         ,         |
 *     |***************************************=***************************************|
 *     |         ,         *         ,         =         ,         *         ,         |
 *     |    i    ,    i    *    i    ,    i    =    i    ,    i    *         ,         |
 *     |         ,         *         ,         =         ,         *         ,         |
 *     |,,,,,,,,,,,,,,,,,,,*,,,,,,,,,,,,,,,,,,,=,,,,,,,,,,,,,,,,,,,*,,,,,,,,,,,,,,,,,,,|
 *     |         ,         *         ,         =         ,         *         ,         |
 *     |    i    ,    i    *    i    ,    i    =    i    ,    i    *         ,         |
 *     |         ,         *         ,         =         ,         *         ,         |
 *     |_______________________________________________________________________________|
 *
 *
 *     l = 4
 *      _______________________________________________________________________________
 *     |    |    ,    |    *    |    ,    |    =    |    ,    |    *    |    ,    |    |
 *     |____|____,____|____*____|____,____|____=____|____,____|____*____|____,____|____|
 *     |    |    ,    |    *    |    ,    |    =    |    ,    |    *    |    ,    |    |
 *     |,,,,|,,,,,,,,,|,,,,*,,,,|,,,,,,,,,|,,,,=,,,,|,,,,,,,,,|,,,,*,,,,|,,,,,,,,,|,,,,|
 *     |    |    ,    |    *    |    ,    |    =    |    ,    |    *    |    ,    |    |
 *     |____|____,____|____*____|____,____|____=____|____,____|____*____|____,____|____|
 *     |    |    ,    |    *    |    ,    |    =    |    ,    |    *    |    ,    |    |
 *     |****|*********|************************=****|*********|************************|
 *     |    |    ,    |    *    |    ,    |    =    |    ,    |    *    |    ,    |    |
 *     |____|____,____|____*____|____,____|____=____|____,____|____*____|____,____|____|
 *     |    |    ,    |    *    |    ,    |    =    |    ,    |    *    |    ,    |    |
 *  -  |....|....,....|....*,,,,|,,,,,,,,,|,,,,=,,,,|,,,,,,,,,|,,,,*,,,,|,,,,,,,,,|,,,,|
 *  -  |    |    , i  | i  * i  | i  ,  i | i  =    |    ,    |    *    |    ,    |    |
 *  -  |____|____,____|____*____|____,____|____=____|____,____|____*____|____,____|____|
 *  -  |    |    , i  | i  * n  | n  ,  n | i  =    |    ,    |    *    |    ,    |    |
 *  -  |    |    ,    |    *    |    ,    |    =    |    ,    |    *    |    ,    |    |
 *     |=========|===================|===================|===================|=========|
 *  -  |    |    , i  | i  * n  | c  ,  n | i  =    |    ,    |    *    |    ,    |    |
 *  -  |____|____,____|____*____|____,____|____=____|____,____|____*____|____,____|____|
 *  -  |    |    , i  | i  * n  | n  ,  n | i  =    |    ,    |    *    |    ,    |    |
 *  -  |....|....,....|....*,,,,|,,,,,,,,,|,,,,=,,,,|,,,,,,,,,|,,,,*,,,,|,,,,,,,,,|,,,,|
 *     |    |    , i  | i  * i  | i  ,  i | i  =    |    ,    |    *    |    ,    |    |
 *     |____|____,____|____*____|____,____|____=____|____,____|____*____|____,____|____|
 *     |    |    , i  | i  * i  | i  ,  i | i  =    |    ,    |    *    |    ,    |    |
 *     |****|*********|************************=****|*********|************************|
 *     |    |    ,    |    *    |    ,    |    =    |    ,    |    *    |    ,    |    |
 *     |____|____,____|____*____|____,____|____=____|____,____|____*____|____,____|____|
 *     |    |    ,    |    *    |    ,    |    =    |    ,    |    *    |    ,    |    |
 *     |,,,,|,,,,,,,,,|,,,,*,,,,|,,,,,,,,,|,,,,=,,,,|,,,,,,,,,|,,,,*,,,,|,,,,,,,,,|,,,,|
 *     |    |    ,    |    *    |    ,    |    =    |    ,    |    *    |    ,    |    |
 *     |____|____,____|____*____|____,____|____=____|____,____|____*____|____,____|____|
 *     |    |    ,    |    *    |    ,    |    =    |    ,    |    *    |    ,    |    |
 *     |_______________________________________________________________________________|
 *
 *
 *
 *
 *
 *
 * Recall our work in two dimensions
 *        y - target
 *        x - source
 *        y_c^{(c)} - target cell center (child cell)
 *        x_c^{(c)} - source cell center (child cell)
 *        y_c^{(p)} - target cell center (parent cell)
 *        x_c^{(p)} - source cell center (parent cell)
 *   1) Upward Pass Part I
 *      a) Work at level l = L
 *      b) Far field expansions about for each source x and about each source cell's center x_c^{(c)}
 *
 *             ln(y-x)         = ln(y-x_c^{(c)}) + \sum_{m=1}^{\infty} \frac{-1}{m} \frac{ (x-x_c^{(c)})^m }{ (y-x_c^{(c)})^m }
 *                             = ln(y-x_c^{(c)}) + \sum_{m=1}^{\infty} \frac{-1}{m} (x-x_c^{(c)})^m    \frac{ (1 }{ (y-x_c^{(c)})^m }
 *                             = \sum_{m=0}^{\infty} b_m(x,x_c^{(c)}) S_m(y-x_c^{(c)}
 *        where
 *                b_m(x,x_c^{(c)}) = \frac{-1}{m} (x-x_c^{(c)})^m
 *                S_m(y-x_c^{(c)}) = \frac{ (1 }{ (y-x_c^{(c)})^m }
 *
 *        here need |x-x_c^{(c)}| < |y-x_c^{(c)}|
 *
 *   2) Upward Pass Part II
 *      a) form one power series for each cell at l = L
 *         i) combine far-field power series of all sources in each cell at l = L
 *      b) SS Translation of far field expansion for each cell's power series at l = L
 *         from child source center center x_c^{(c)} to parent source cell center x_c^{(p)}
 *         i)  only coefficients b_m(x,x_c^{(c)} are operated on in the translation process
 *             to obtain the new coefficients b_m(x,x_c^{(p)} for the new series that corresponds
 *             to the new powers S_m(y-x_c^{(p)}
 *         ii) power's are not determined in the translation process, but can always be calculated
 *             using new center x_c^{(p)}
 *      c) repeat steps (a) and (b) at level l = L-1
 *         i) continue until and including l = 2
 *
 *     l = 2
 *      _______________________________________________________________________________
 *     |                   *                   =                   *                   |
 *     |                   *                   =                   *                   |
 *     |                   *                   =                   *                   |
 *     |        i          *         i         =         i         *         i         |
 *     |                   *                   =                   *                   |
 *     |                   *                   =                   *                   |
 *     |                   *                   =                   *                   |
 *     |***************************************=***************************************|
 *     |                   *                   =                   *                   |
 *     |                   *                   =                   *                   |
 *     |                   *                   =                   *                   |
 *     |         n         *         n         =         n         *         i         |
 *     |                   *                   =                   *                   |
 *     |                   *                   =                   *                   |
 *     |                   *                   =                   *                   |
 *     |                   *                   =                   *                   |
 *     |=========|===================|=================================================|
 *     |                   *                   =                   *                   |
 *     |                   *                   =                   *                   |
 *     |                   *                   =                   *                   |
 *     |         n         *         c         =         n         *         i         |
 *     |                   *                   =                   *                   |
 *     |                   *                   =                   *                   |
 *     |                   *                   =                   *                   |
 *     |***************************************=***************************************|
 *     |                   *                   =                   *                   |
 *     |                   *                   =                   *                   |
 *     |                   *                   =                   *                   |
 *     |         n         *         n         =         n         *         i         |
 *     |                   *                   =                   *                   |
 *     |                   *                   =                   *                   |
 *     |                   *                   =                   *                   |
 *     |_______________________________________________________________________________|
 *
 *   3) Downward Pass Part I
 *      a) Work at level l = 2
 *      b) For each cell c at level l = 2
 *          i) determine the interaction list cells
 *         ii) SR translate the power series for each interaction list cell
 *             to the center y_c^{(p)} of the target cell c
 *        iii) the translation process operates on the coefficients b_m(x,x_c^{(p)} to
 *             determine the coefficients a_m(x,y_c^{(p)}) for the near field series
 *
 *               \sum_{m=0}^{\infty} a_m(x,y_c^{(p)}) R_m(y-y_c^{(p)}
 *
 *        Note
 *
 *             ln(y-x)         = ln(x-x_c^{(c)}) + \sum_{m=1}^{\infty} \frac{-1}{m} \frac{ (y-y_c^{(c)})^m }{ (x-y_c^{(c)})^m }
 *                             = ln(y-x_c^{(c)}) + \sum_{m=1}^{\infty} \frac{ -1 }{ m(x-y_c^{(c)})^m } (y-y_c^{(c)})^m
 *                             = \sum_{m=0}^{\infty} a_m(x,y_c^{(c)}) R_m(y-y_c^{(c)}
 *        where
 *                a_m(x,y_c^{(c)}) = \frac{ -1 }{ m(x-y_c^{(c)})^m }
 *                R_m(y-y_c^{(c)}) = (y-y_c^{(c)})^m
 *
 *        here need |y-y_c^{(c)}| < |x-y_c^{(c)}|
 *
 *
 *         iv) form one near-field power series for cell by combining all near-field power series
 *             translated to cell (can combine since each series has same powers)
 *       Notes: * At l = 2 there are only near neighbors and an interaction list for each cell
 *                (nothing more)
 *              * The near neighbors n are too close to approximate with a far-field series
 *              * If we were at the lowest FMM level, we would calculate the interaction with the near neighbor
 *                source points and targets in cell c directly (and do not approximate with a far-field series)
 *
 *     l = 3
 *      _______________________________________________________________________________
 *     |         ,         *         ,         =         ,         *         ,         |
 *     |         ,         *         ,         =         ,         *         ,         |
 *     |         ,         *         ,         =         ,         *         ,         |
 *     |,,,,,,,,,,,,,,,,,,,*,,,,,,,,,,,,,,,,,,,=,,,,,,,,,,,,,,,,,,,*,,,,|,,,,,,,,,,,,,,|
 *     |         ,         *         ,         =         ,         *         ,         |
 *     |         ,         *         ,         =         ,         *         ,         |
 *     |         ,         *         ,         =         ,         *         ,         |
 *     |***************************************=***************************************|
 *     |         ,         *         ,         =         ,         *         ,         |
 *     |    i    ,    i    *    i    ,    i    =    i    ,    i    *         ,         |
 *     |         ,         *         ,         =         ,         *         ,         |
 *     |,,,,,,,,,,,,,,,,,,,*,,,,,,,,,,,,,,,,,,,=,,,,,,,,,,,,,,,,,,,*,,,,,,,,,,,,,,,,,,,|
 *  -  |         ,         *         ,         =         ,         *         ,         |
 *  -  |    i    ,    n    *    n    ,    n    =    i    ,    i    *         ,         |
 *  -  |         ,         *         ,         =         ,         *         ,         |
 *  -  |         ,         *         ,         =         ,         *         ,         |
 *     |=========|===================|=================================================|
 *  -  |         ,         *         ,         =         ,         *         ,         |
 *  -  |    i    ,    n    *    c    ,    n    =    i    ,    i    *         ,         |
 *  -  |         ,         *         ,         =         ,         *         ,         |
 *  -  |,,,,,,,,,,,,,,,,,,,*,,,,,,,,,,,,,,,,,,,=,,,,,,,,,,,,,,,,,,,*,,,,,,,,,,,,,,,,,,,|
 *     |         ,         *         ,         =         ,         *         ,         |
 *     |    i    ,    n    *    n    ,    n    =    i    ,    i    *         ,         |
 *     |         ,         *         ,         =         ,         *         ,         |
 *     |***************************************=***************************************|
 *     |         ,         *         ,         =         ,         *         ,         |
 *     |    i    ,    i    *    i    ,    i    =    i    ,    i    *         ,         |
 *     |         ,         *         ,         =         ,         *         ,         |
 *     |,,,,,,,,,,,,,,,,,,,*,,,,,,,,,,,,,,,,,,,=,,,,,,,,,,,,,,,,,,,*,,,,,,,,,,,,,,,,,,,|
 *     |         ,         *         ,         =         ,         *         ,         |
 *     |    i    ,    i    *    i    ,    i    =    i    ,    i    *         ,         |
 *     |         ,         *         ,         =         ,         *         ,         |
 *     |_______________________________________________________________________________|
 *
 *
 *    4) Downward Pass Part II
 *       a) Step down to the next level (For each level 3 <= l <= L )
 *       b) for each cell at level (first level for Part II is l = 3)
 *           i) perform an RR translation of the near-field expansion determined
 *              in Downward Pass Part I for the parent cell (from parent to child)
 *
 *                \sum_{m=0}^{\infty} a_m(x,y_c^{(c)}) R_m(y-y_c^{(c)}
 *
 *          ii) determine the interaction list for the cell
 *         iii) SR translate the far-field power series for each interaction list cell
 *          iv) form a single series for the RR and SR translated series to the cell
 *              (since powers for each series are the same)
 *    5) FMM Calculation
 *       a) Work at the lowest level l = L
 *       b) For each cell at l = L
 *          i) For each target in the cell
 *             - calculate the interaction between all sources beyond the near neighbors
 *               using the single power series that has been translated and formed for the cell
 *             - calculate the interaction between the target and all the sources in the near
 *               neighbors by direct calculation (since sources too close for approximation)
 *
 */
