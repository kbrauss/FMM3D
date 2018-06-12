/*
 * test6.cc
 *
 *  Created on: Apr 23, 2017
 *      Author: dbpc
 */

/*
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
*/
int main()
{
  // note that these refinement levels counts start on l = 1
  int DEFAULT_NUM_LEVEL = 4;                   // default refinement level
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

  // ********************************************************************
  // Creating the Target Points *****************************************
  // ********************************************************************
  // For this test only one target point
  // Note: In general the target points are same as source points
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


