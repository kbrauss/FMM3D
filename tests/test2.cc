/*
 * test2.cc
 *
 *  Created on: Mar 29, 2017
 *      Author: dbpc
 */



#include<iostream>
#include<vector>
#include<cmath>

#include "Point.h"
#include "Util.h"
#include "Potential.h"
#include "FmmTree.h"


/**
 *  test2.cc
 *  Created on: March 29, 2017
 *      Author: dbpc
 *
 *  The program below tests approximation of \psi(\mathbf{r}) = 1 / \|\mathbf{r}\| by a truncated near-field
 *  series.  The series is truncated to order |alpha| = 3 resulting in p = 1 + 3 + 5 + 7 = 16 terms.
 *
 *  The far-field expansion has the form sum_{i=1}^{p} ( b_m (x_i, x^s_*) S_m(y,x^s_*) ) where x_i is the source
 *  point, x^s_* is the source cell center, and y is the target cell.  Using Cartesian coordinates the powers
 *  and factorial (1/alpha!) (x_b - z_b)^{alpha} are the coefficients b_m(x_b, z_b) and the derivatives
 *  D^{alpha} f(z_b - x_a) are the powers S_m(x_a,z_b).  The notation x_i = x_a represents the source point,
 *  x^s_* = z_a is the source center, and y = x_b is the target point.  The far-field expansion in the Taylor
 *  series is
 *
 *      f(x_b - x_a) = sum_{|alpha| >= 0} [ ( (1 / (alpha !)) (x_a - z_a)^{alpha} ) D^{alpha} f(z_a - x_b) ]
 *
 *  The series converges when || x_a - z_a || < || z_a - x_b || (or || x_i - x^s_* || < || x^s_* - y ||=|| y - x_* ||)
 *  In Dehnen's coordinate system the same series has the form
 *
 *      \frac{1}{|x_a - x_b|} = \frac{1}{|(x_a - z_A) - (x_b - z_A)|}
 *                            = \sum_{n=0}^{\infty} \sum_{m=-n}^{n} \gamma_n^{m*} (x_a - z_A) \theta_n^{m} (x_b - z_A)
 *
 *  The coefficients of the series are \gamma_n^{m*} (x_a - z_A) and the powers are \theta_n^{m} (x_b - z_A)
 *
 *  For a near-field series sum_{i=1}^{p} ( a_m (y, x^t_*) R_m(y,x^t_*) ) we replace the variables in the far-field
 *  expansion z_a with z_b (source cell center for the target cell center) and swap x_a and x_b.  The coefficients
 *  a_m (y, x^t_*) are (1/alpha!) D^{alpha} f(z_b - x_a) and the powers R_m(y,x^t_*) are (x_b - z_b)^{alpha}.
 *  The near-field expansion as a Taylor series in Cartesian coordinates is
 *
 *      f(x_b - x_a) = sum_{|alpha| >= 0} [ ( (1 / (alpha !)) (x_b - z_b)^{alpha} ) D^{alpha} f(z_b - x_a) ]
 *
 *  The series converges when || x_b - z_b || < || z_b - x_a || (or || y - x^t_* || < || x^t_* - x_i ||=|| y - x_* ||)
 *  In Dehnen's coordinate system the same series has the form
 *
 *      \frac{1}{|x_a - x_b|} = \frac{1}{|(x_b - z_B) - (x_a - z_B)|}
 *                            = \sum_{n=0}^{\infty} \sum_{m=-n}^{n} \gamma_n^{m*} (x_b - z_B) \theta_n^{m} (x_a - z_B)
 *
 *  and the coefficients of the series are \theta_n^{m} (x_a - z_B) and the powers are \gamma_n^{m*} (x_b - z_B)
 *
 *  level l = 2 refinement
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
 *  cell centers for the test.  An interaction list cell n = 45 at l =
 *  is also given, and contains the point o = x[1].
 *  Below are the coordinates for the points mentioned, their cell numbers, and
 *  the cell number for level l = 3.
 *
 *    source - x[0] = [0.90, 0.925, 0.95] located in cell n = 511 at l = 3;
 *    point in interaction list cell - x[1] = [0.38, 0.06, 0.40] located in cell n = 45 at l = 3;
 *    target - y[0] = [0.15,0.03,0.20] located in cell n = 5 at l = 3;
 *
 *  The locations for the source, interaction list point, and target are shown below
 *  The y-axis is into the screen.  For ease of view, depth of the points in the y-direction is not shown
 *
 *                                 z-axis |
 *                                        |
 *                                        |
 *
 * "    "  "    "  "    "  "    "  1.0000  ________________________________________________________________
 *                                        |   |   |   |   |   |   |   |   ||   |   |   |   |   |   | s |   |
 * 0.9375                                 |---|---|---|---|---|---|---|---||---|---|---|---|---|---|---z_a-|
 * "    "  0.8750                         |___|___|___|___|___|___|___|___||___|___|___|___|___|___|___|___|
 * 0.8125                                 |   |   |   |   |   |   |   |   ||   |   |   |   |   |   |   |   |
 *                                        |---|---|---|---|---|---|---|---||---|---|---|---|---|---|---|---|
 * "    "  "    "  0.7500                 |___|___|___|___|___|___|___|___||___|___|___|___|___|___|___|___|
 *                                        |   |   |   |   |   |   |   |   ||   |   |   |   |   |   |   |   |
 * 0.6875                                 |---|---|---|---|---|---|---|---||---|---|---|---|---|---|---|---|
 * "    "  0.6250                         |___|___|___|___|___|___|___|___||___|___|___|___|___|___|___|___|
 * 0.5625                                 |   |   |   |   |   |   |   |   ||   |   |   |   |   |   |   |   |
 *                                        |---|---|---|---|---|---|---|---||---|---|---|---|---|---|---|---|
 *                                        |   |   |   |   |   |   |   |   ||   |   |   |   |   |   |   |   |
 * "    "  "    "  "    "  0.5000         |===|===|===|===|===|===|===|===||===|===|===|===|===|===|===|===|
 *                                        |   |   |   |   #   |   |   |   ||   |   |   |   #   |   |   |   |
 * 0.4375                                 |---i---|---i---#---i---|---i---||---|---|---|---#---|---|---|---|
 * "    "  0.3750                         |___|__pn___|___#___|__pn_o_|___||___|___|___|___#___|___|___|___|
 * 0.3125                                 |   |   |   |   #   |   |   |   ||   |   |   |   #   |   |   |   |
 *                                        |---n---|---n---#---n---|---i---||---|---|---|---#---|---|---|---|
 *                                        |   |   |   |   #   |   |   |   ||   |   |   |   #   |   |   |   |
 * "    "  "    "  0.2500                 |###|###|###|#######|###|###|###||###|###|###|#######|###|###|###|
 *                                        |   |   | t |   #   |   |   |   ||   |   |   |   #   |   |   |   |
 * 0.1875                                 |---n---|---c---#---n---|---i---||---|---|---|---#---|---|---|---|
 * "    "  0.1250                         |___|__ p___|___#___|__pn___|___||___|___|___|___#___|___|___|___|
 * 0.0625                                 |   |   |   |   #   |   |   |   ||   |   |   |   #   |   |   |   |
 *                                        |---n---|---n---#---n---|---i---||---|---|---|---#---|---|---|---|
 * "    "  "    "  "    "  "    "  0.0000 |_ _|_ _|___|_ _#___|___|___|___||___|___|___|___#___|___|___|___|    _ _ _
 *                                      0.00     0.125   0.25    0.375    0.5             0.75    0.875   1.00        x-axis
 *
 * l = 4   l = 3   l = 2   l = 1   l = 0                                                   t - target = y[0] = z_b
 * -- levels for refinements above                                                         c - z_a,c - target cell center
 *                                                                                         s - source = x[0] = z_a
 *                                                                                         c - z_a,c - source cell center
 *                                                                                         p - z_a,p - parent center
 *                                                                                         n - nearest neighbor center
 *                                                                                         i - interaction list cell center
 *                                                                                         o - interaction list point x[1]
 *                                                                                         pn - parent neighbor center
 *
 *
 *  1) The code below performs a direct calculation of the potential function for the source s
 *     and the target t mentioned above.
 *  2) The code below then tests the member functions Potential::getRVector and Potential::getRCoeff
 *     of the class Potential in Potential.cc by performing an r-expansion approximation to be compared
 *     to the direct calculation.  The Potential.cc member functions RVector returns the powers of the
 *     near-field series as a matrix and the member function getRCoeff returns the coefficients of the
 *     near-field approximation as a matrix.  The corresponding elements of the two matrices are multiplied
 *     together to form the r-expansion and this near-field approximation is compared with the direct calculation
 *     for accuracy.
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
  std::cout << "center_yBox(x,y,z) = " << "center_yBox(" << center_yBox.getX() << ","
		    << center_yBox.getY() << "," << center_yBox.getZ() << ")" << std::endl;

  // obtaining far-field coefficients (s-expansion coefficients)
  // far field coefficients are (z_a - x_a)^{\alpha} = (center_yBox - y[0])^{\alpha}
  // Note: y[0] is the target ( in Main.cc it is x )
  std::cout << "We are here" << std::endl;
  std::vector<std::vector<std::complex<double> > > x_RCoeff = potential.getRCoeff(x[0].getCoord(), center_yBox.getCoord());

  // obtaining far-field powers (s-expansion powers)
  // far field powers are D^{\alpha} f (x_b - z_a) = D^{\alpha} f (x[0] - center_yBox)
  std::vector<std::vector<std::complex<double> > > x_RVec = potential.getRVector(y[0].getCoord(), center_yBox.getCoord());

  for (unsigned int i=0; i<=abs_alpha; ++i)
    for (unsigned int j=0; j<=abs_alpha; ++j)
    {
//      std::cout << "We are here" << std::endl;
      r_expansion = r_expansion + x_RCoeff[j][i] * x_RVec[j][i];
    }

  std::cout << "r_expansion = " << r_expansion.real() << " + i" << r_expansion.imag() << std::endl;

  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "**************************************************" << std::endl;
  std::cout << "**************************************************" << std::endl;

}



