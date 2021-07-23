scale_f = 1.0e-3;
cl__1 = 0.02*scale_f;

Mesh.CharacteristicLengthExtendFromBoundary = 0;
// point definition section
// coordinates extracted from
// HFSS in mm, scaled as per scale_f

// Point(1) = {-0.05*scale_f, -0.015*scale_f, 0.018*scale_f, cl__1};
Point(2) = {-0.05*scale_f, -0.015*scale_f, 0.025*scale_f, cl__1};
Point(3) = {-0.06*scale_f, -0.015*scale_f, 0.025*scale_f, cl__1};
// Point(4) = {-0.06*scale_f, -0.015*scale_f, 0.022*scale_f, cl__1};
// Point(5) = {-0.095*scale_f, -0.015*scale_f, 0.022*scale_f, cl__1};
// Point(6) = {-0.095*scale_f, -0.015*scale_f, 0.029*scale_f, cl__1};
// Point(7) = {-0.2*scale_f, -0.015*scale_f, 0.029*scale_f, cl__1};
// Point(8) = {-0.2*scale_f, -0.015*scale_f, 0.025*scale_f, cl__1};
// Point(9) = {-0.105*scale_f, -0.015*scale_f, 0.025*scale_f, cl__1};
// Point(10) = {-0.105*scale_f, -0.015*scale_f, 0.018*scale_f, cl__1};
// Point(11) = {-0.105*scale_f, -0.005*scale_f, 0.025*scale_f, cl__1};
// Point(14) = {-0.2*scale_f, -0.005000000000000001*scale_f, 0.025*scale_f, cl__1};
// Point(15) = {-0.2*scale_f, -0.005*scale_f, 0.029*scale_f, cl__1};
// Point(18) = {-0.095*scale_f, -0.005000000000000001*scale_f, 0.029*scale_f, cl__1};
// Point(23) = {-0.095*scale_f, -0.005000000000000001*scale_f, 0.022*scale_f, cl__1};
// Point(28) = {-0.06*scale_f, -0.005000000000000001*scale_f, 0.022*scale_f, cl__1};
// Point(32) = {-0.105*scale_f, -0.005000000000000001*scale_f, 0.018*scale_f, cl__1};
// Point(33) = {-0.05*scale_f, -0.005*scale_f, 0.018*scale_f, cl__1};
Point(39) = {-0.05*scale_f, -0.065*scale_f, 0.025*scale_f, cl__1};
Point(43) = {-0.05*scale_f, -0.015*scale_f, 0.029*scale_f, cl__1};
Point(44) = {-0.05*scale_f, -0.065*scale_f, 0.029*scale_f, cl__1};

// total 21 points

Point(54) = {-0.06*scale_f, -0.015*scale_f, 0.029*scale_f, cl__1};
Point(59) = {-0.06*scale_f, -0.075*scale_f, 0.025*scale_f, cl__1};
Point(60) = {-0.06*scale_f, -0.075*scale_f, 0.029*scale_f, cl__1};

//total point 24

Point(63) = {0.065*scale_f, -0.065*scale_f, 0.029*scale_f, cl__1};
Point(64) = {0.065*scale_f, -0.065*scale_f, 0.025*scale_f, cl__1};
Point(67) = {0.07500000000000001*scale_f, -0.075*scale_f, 0.025*scale_f, cl__1};
Point(68) = {0.075*scale_f, -0.075*scale_f, 0.029*scale_f, cl__1};
Point(71) = {0.065*scale_f, 0.065*scale_f, 0.029*scale_f, cl__1};
Point(72) = {0.065*scale_f, 0.065*scale_f, 0.025*scale_f, cl__1};
Point(75) = {0.075*scale_f, 0.075*scale_f, 0.025*scale_f, cl__1};
Point(76) = {0.075*scale_f, 0.075*scale_f, 0.029*scale_f, cl__1};
Point(79) = {-0.065*scale_f, 0.065*scale_f, 0.029*scale_f, cl__1};
Point(80) = {-0.065*scale_f, 0.065*scale_f, 0.025*scale_f, cl__1};
Point(83) = {-0.075*scale_f, 0.075*scale_f, 0.025*scale_f, cl__1};
Point(84) = {-0.075*scale_f, 0.075*scale_f, 0.029*scale_f, cl__1};
Point(87) = {-0.065*scale_f, -0.08000000000000002*scale_f, 0.029*scale_f, cl__1};
Point(88) = {-0.065*scale_f, -0.08*scale_f, 0.025*scale_f, cl__1};
Point(91) = {-0.075*scale_f, -0.09*scale_f, 0.025*scale_f, cl__1};
Point(92) = {-0.075*scale_f, -0.09*scale_f, 0.029*scale_f, cl__1};
Point(95) = {0.08000000000000002*scale_f, -0.08*scale_f, 0.029*scale_f, cl__1};
Point(96) = {0.08*scale_f, -0.08*scale_f, 0.025*scale_f, cl__1};
Point(99) = {0.09*scale_f, -0.09*scale_f, 0.025*scale_f, cl__1};
Point(100) = {0.09*scale_f, -0.09*scale_f, 0.029*scale_f, cl__1};

// this section 20 points


Point(103) = {0.08*scale_f, 0.08*scale_f, 0.029*scale_f, cl__1};
Point(104) = {0.08*scale_f, 0.08*scale_f, 0.025*scale_f, cl__1};
Point(107) = {0.09*scale_f, 0.09*scale_f, 0.025*scale_f, cl__1};
Point(108) = {0.09*scale_f, 0.09*scale_f, 0.029*scale_f, cl__1};
Point(111) = {-0.08*scale_f, 0.08*scale_f, 0.029*scale_f, cl__1};
Point(112) = {-0.08*scale_f, 0.08*scale_f, 0.025*scale_f, cl__1};
Point(115) = {-0.09*scale_f, 0.09*scale_f, 0.025*scale_f, cl__1};
Point(116) = {-0.09*scale_f, 0.09*scale_f, 0.029*scale_f, cl__1};
Point(119) = {-0.08*scale_f, -0.095*scale_f, 0.029*scale_f, cl__1};
Point(120) = {-0.08*scale_f, -0.095*scale_f, 0.025*scale_f, cl__1};
Point(123) = {-0.09*scale_f, -0.105*scale_f, 0.025*scale_f, cl__1};
Point(124) = {-0.09*scale_f, -0.105*scale_f, 0.029*scale_f, cl__1};
Point(127) = {0.095*scale_f, -0.095*scale_f, 0.029*scale_f, cl__1};
Point(128) = {0.095*scale_f, -0.095*scale_f, 0.025*scale_f, cl__1};
Point(131) = {0.105*scale_f, -0.105*scale_f, 0.025*scale_f, cl__1};
Point(132) = {0.105*scale_f, -0.105*scale_f, 0.029*scale_f, cl__1};
Point(135) = {0.095*scale_f, 0.095*scale_f, 0.029*scale_f, cl__1};
Point(136) = {0.095*scale_f, 0.095*scale_f, 0.025*scale_f, cl__1};
Point(139) = {0.105*scale_f, 0.105*scale_f, 0.025*scale_f, cl__1};
Point(140) = {0.105*scale_f, 0.105*scale_f, 0.029*scale_f, cl__1};
Point(143) = {-0.095*scale_f, 0.095*scale_f, 0.029*scale_f, cl__1};
Point(144) = {-0.095*scale_f, 0.095*scale_f, 0.025*scale_f, cl__1};
Point(147) = {-0.105*scale_f, 0.105*scale_f, 0.025*scale_f, cl__1};
Point(148) = {-0.105*scale_f, 0.105*scale_f, 0.029*scale_f, cl__1};
Point(151) = {-0.095*scale_f, 0.005000000000000003*scale_f, 0.029*scale_f, cl__1};
Point(152) = {-0.095*scale_f, 0.005*scale_f, 0.025*scale_f, cl__1};
Point(155) = {-0.105*scale_f, 0.015*scale_f, 0.025*scale_f, cl__1};
Point(156) = {-0.105*scale_f, 0.015*scale_f, 0.029*scale_f, cl__1};
Point(157) = {-0.2*scale_f, 0.005000000000000001*scale_f, 0.029*scale_f, cl__1};

// total 29 points

Point(186) = {-0.2*scale_f, 0.015*scale_f, 0.029*scale_f, cl__1};
Point(189) = {-0.2*scale_f, 0.015*scale_f, 0.025*scale_f, cl__1};


// total 2 points

Point(206) = {-0.2*scale_f, 0.005*scale_f, 0.025*scale_f, cl__1};

Point(229) = {-0.245*scale_f, 0.015*scale_f, 0.025*scale_f, cl__1};
Point(230) = {-0.245*scale_f, 0.015*scale_f, 0.029*scale_f, cl__1};
Point(231) = {-0.242*scale_f, 0.015*scale_f, 0.029*scale_f, cl__1};
Point(232) = {-0.242*scale_f, 0.015*scale_f, 0.025*scale_f, cl__1};
Point(233) = {-0.255*scale_f, 0.005*scale_f, 0.029*scale_f, cl__1};
Point(234) = {-0.255*scale_f, 0.005*scale_f, 0.025*scale_f, cl__1};
Point(235) = {-0.242*scale_f, 0.005*scale_f, 0.025*scale_f, cl__1};
Point(236) = {-0.242*scale_f, 0.005*scale_f, 0.029*scale_f, cl__1};

// total 9 points


Point(243) = {-0.245*scale_f, 0.145*scale_f, 0.025*scale_f, cl__1};
Point(244) = {-0.245*scale_f, 0.145*scale_f, 0.029*scale_f, cl__1};
Point(247) = {-0.255*scale_f, 0.155*scale_f, 0.029*scale_f, cl__1};
Point(248) = {-0.255*scale_f, 0.155*scale_f, 0.025*scale_f, cl__1};
Point(251) = {0.145*scale_f, 0.145*scale_f, 0.025*scale_f, cl__1};
Point(252) = {0.145*scale_f, 0.145*scale_f, 0.029*scale_f, cl__1};
Point(255) = {0.155*scale_f, 0.155*scale_f, 0.029*scale_f, cl__1};
Point(256) = {0.155*scale_f, 0.155*scale_f, 0.025*scale_f, cl__1};
Point(259) = {0.145*scale_f, -0.145*scale_f, 0.025*scale_f, cl__1};
Point(260) = {0.145*scale_f, -0.145*scale_f, 0.029*scale_f, cl__1};
Point(263) = {0.155*scale_f, -0.155*scale_f, 0.029*scale_f, cl__1};
Point(264) = {0.155*scale_f, -0.155*scale_f, 0.025*scale_f, cl__1};
Point(267) = {-0.245*scale_f, -0.145*scale_f, 0.025*scale_f, cl__1};
Point(268) = {-0.245*scale_f, -0.145*scale_f, 0.029*scale_f, cl__1};
Point(271) = {-0.255*scale_f, -0.155*scale_f, 0.029*scale_f, cl__1};
Point(272) = {-0.255*scale_f, -0.155*scale_f, 0.025*scale_f, cl__1};
Point(275) = {-0.245*scale_f, -0.01499999999999999*scale_f, 0.025*scale_f, cl__1};
Point(276) = {-0.245*scale_f, -0.015*scale_f, 0.029*scale_f, cl__1};
Point(279) = {-0.255*scale_f, -0.005000000000000003*scale_f, 0.029*scale_f, cl__1};
Point(280) = {-0.255*scale_f, -0.005*scale_f, 0.025*scale_f, cl__1};
Point(281) = {-0.242*scale_f, -0.015*scale_f, 0.025*scale_f, cl__1};
Point(296) = {-0.242*scale_f, -0.005*scale_f, 0.025*scale_f, cl__1};
Point(299) = {-0.242*scale_f, -0.005*scale_f, 0.029*scale_f, cl__1};

// total 23 points

Point(316) = {-0.242*scale_f, -0.015*scale_f, 0.029*scale_f, cl__1};

// total 1 point









Line(1) = {43, 2};
Line(2) = {2, 3};
Line(3) = {54, 3};
// Line(4) = {4, 5};
// Line(5) = {5, 6};
// Line(6) = {6, 7};
// Line(7) = {7, 8};
// Line(8) = {8, 9};
// Line(9) = {9, 10};
// Line(10) = {10, 1};
//Line(11) = {11, 12};
// Line(11) = {11, 9};

//Line(12) = {12, 8};   // no line 12
//Line(12) = {12, 13};

// Line(13) = {8, 14};
//Line(13) = {13, 14};

// Line(14) = {14, 11};

// wrong line
//Line(15) = {15, 16};
// Line(15) = {15, 7};

//Line(16) = {16, 17};  // no line 16
//Line(16) = {7, 17};

//Line(17) = {17, 18};
// Line(17) = {6, 18};
// total 15 lines
// Line(18) = {18, 15};
//Line(19) = {19, 20};
//Line(20) = {20, 21};
//Line(20) = {7, 21};  // this line is wrong
//Line(21) = {14, 8};
//Line(22) = {8, 16};
//Line(22) = {22, 19};
// Line(23) = {23, 18};

//Line(25) = {25, 26};
//Line(24) = {18, 6};
//Line(24) = {24, 25};
//Line(25) = {25, 26};
//Line(26) = {26, 23};
// Line(26) = {5, 23};

//Line(27) = {27, 28};
//Line(28) = {28, 29};
//Line(29) = {29, 30};
//Line(30) = {30, 27};
// Line(27) = {4, 28};
// Line(28) = {28, 23};
//Line(29) = {23, 5}; // same as line 26
//Line(30) = {5, 4};
//Line(31) = {31, 32};
//Line(32) = {32, 33};
//Line(33) = {33, 34};
//Line(34) = {34, 31};
// Line(31) = {10, 32};
// Line(32) = {32, 33};
// Line(33) = {33, 1};
//Line(34) = {1, 10};

//Line(35) = {35, 36};
//Line(36) = {36, 37};
//Line(37) = {37, 38};
//Line(38) = {38, 35};
//Line(35) = {9, 11};
// Line(36) = {11, 32};
//Line(37) = {32, 10};
//Line(38) = {10, 9};

//Line(39) = {39, 40};
//Line(40) = {40, 41};
//Line(41) = {41, 42};
//Line(42) = {42, 43};
//Line(43) = {43, 44};

Line(39) = {39, 2};
//Line(40) = {2, 1};
//Line(41) = {1, 33};
// Line(42) = {33, 43};
Line(43) = {43, 44};
Line(44) = {44, 39};


//Line(45) = {43, 33};
//Line(46) = {33, 32};
//Line(47) = {32, 11};
//Line(48) = {11, 14};
 // Line(49) = {14, 15};

//Line(50) = {15, 18};

//Line(51) = {18, 23};   // line 23
//Line(52) = {23, 28};  // line 28

////////////////////////////////////////

  // Line(53) = {28, 54};
  Line(54) = {54, 43};

  Line(58) = {3, 59};
  Line(59) = {59, 60};
  Line(60) = {60, 54};

// total 34 lines
////////////////////////////////
//Line(55) = {54, 28};
//Line(56) = {28, 4};
//Line(57) = {4, 3};

/////////////////////////////////////////
Line(61) = {44, 63};
Line(62) = {39, 64};
Line(63) = {60, 68};
Line(64) = {59, 67};
Line(65) = {64, 63};
Line(66) = {67, 68};



Line(68) = {71, 63};
Line(69) = {72, 64};
Line(70) = {76, 68};
Line(71) = {75, 67};
Line(72) = {71, 72};
Line(73) = {75, 76};


Line(74) = {71, 79};
Line(75) = {72, 80};
Line(76) = {76, 84};
Line(77) = {75, 83};
Line(78) = {79, 80};
Line(79) = {84, 83};


Line(80) = {87, 79};
Line(81) = {88, 80};
Line(82) = {92, 84};
Line(83) = {91, 83};
Line(84) = {87, 88};
Line(85) = {92, 91};

Line(86) = {87, 95};
Line(87) = {88, 96};
Line(88) = {92, 100};
Line(89) = {91, 99};
Line(90) = {95, 96};
Line(91) = {100, 99};


Line(92) = {103, 95};
Line(93) = {104, 96};
Line(94) = {108, 100};
Line(95) = {107, 99};
Line(96) = {103, 104};
Line(97) = {108, 107};

Line(98) = {103, 111};
Line(99) = {104, 112};
Line(100) = {108, 116};
Line(101) = {107, 115};
Line(102) = {111, 112};
Line(103) = {116, 115};

Line(104) = {119, 111};
Line(105) = {120, 112};
Line(106) = {124, 116};
Line(107) = {123, 115};
Line(108) = {119, 120};
Line(109) = {124, 123};

Line(110) = {119, 127};
Line(111) = {120, 128};
Line(112) = {124, 132};
Line(113) = {123, 131};
Line(114) = {127, 128};
Line(115) = {132, 131};


Line(116) = {135, 127};
Line(117) = {136, 128};
Line(118) = {140, 132};
Line(119) = {139, 131};
Line(120) = {135, 136};
Line(121) = {140, 139};

Line(122) = {135, 143};
Line(123) = {136, 144};
Line(124) = {140, 148};
Line(125) = {139, 147};
Line(126) = {143, 144};
Line(127) = {148, 147};

Line(128) = {151, 143};
Line(129) = {152, 144};
Line(130) = {156, 148};
Line(131) = {155, 147};
Line(132) = {151, 152};
Line(133) = {156, 155};

Line(134) = {151, 157};
Line(135) = {152, 206};
Line(136) = {156, 186};
Line(137) = {155, 189};
Line(138) = {157, 206};
Line(139) = {186, 189};

///////////////////////////////////////
// total lines 78



Line(140) = {189, 206};
Line(141) = {186, 157};

Line(142) = {231, 232};
Line(143) = {230, 229};
Line(144) = {231, 230};
Line(145) = {232, 229};


Line(146) = {231, 236};
Line(147) = {232, 235};
Line(148) = {235, 236};

Line(149) = {233, 236};
Line(150) = {234, 235};
Line(151) = {234, 233};


Line(152) = {230, 244};
Line(153) = {229, 243};
Line(154) = {233, 247};
Line(155) = {234, 248};
Line(156) = {244, 243};
Line(157) = {247, 248};

Line(158) = {252, 244};
Line(159) = {251, 243};
Line(160) = {255, 247};
Line(161) = {256, 248};
Line(162) = {252, 251};
Line(163) = {255, 256};

Line(164) = {252, 260};
Line(165) = {251, 259};
Line(166) = {255, 263};
Line(167) = {256, 264};
Line(168) = {260, 259};
Line(169) = {263, 264};

Line(170) = {268, 260};
Line(171) = {267, 259};
Line(172) = {271, 263};
Line(173) = {272, 264};
Line(174) = {268, 267};
Line(175) = {271, 272};


Line(176) = {268, 276};
Line(177) = {267, 275};
Line(178) = {271, 279};
Line(179) = {272, 280};
Line(180) = {279, 280};
Line(181) = {276, 275};

Line(182) = {316, 276};
Line(183) = {281, 275};
Line(184) = {299, 279};
Line(185) = {296, 280};
Line(186) = {299, 296};
Line(187) = {316, 281};


Line(188) = {299, 316};
Line(189) = {296, 281};

////////////////////////
// total lines 50
// total lines 162


// surface defination


// Line Loop(2) = {7, 15, 49, 13};
// Plane Surface(2) = {2};  // port 3 Surface

Line Loop(3) = {138, 141, -139, -140};
Plane Surface(3) = {3};  // port 4 Surface

///////////////////////////////////////////////
// Line Loop(4) = {-7, -8, -9, -10, -1, -2, -3, -4, -5, -6};
// Plane Surface(4) = {4};  // port 4 Surface

// Line Loop(5) = {14, 36, 32, 42, -54, -53, 28, 23, 18, -49};
// Plane Surface(5) = {5};  // port 4 Surface

// Line Loop(6) = {-14, -11, 8, -13};
// Plane Surface(6) = {6};  //

// Line Loop(7) = {31, 9, 11, -36};
// Plane Surface(7) = {7};

// Line Loop(8) = {-32, -33, 10, -31};
// Plane Surface(8) = {8};

////////////////
// Line Loop(9) = {6, -17, -18, -15};
// Plane Surface(9) = {9};

Line Loop(10) = {-60, -59, -58, -3};
Plane Surface(10) = {10};

Line Loop(11) = {-54, 3, -2, -1};
Plane Surface(11) = {11};

// Line Loop(12) = {-39, 1, 33, -42, -43, -44};
// Plane Surface(12) = {12};

Line Loop(13) = {-43, 1, -39, -44};
Plane Surface(13) = {13};

Line Loop(14) = {-64, -66, 63, 59};
Plane Surface(14) = {14};

Line Loop(15) = {-61, 65, 62, 44};
Plane Surface(15) = {15};

Line Loop(16) = {68, -72, -69, -65};
Plane Surface(16) = {16};

Line Loop(17) = {71, -73, -70, 66};
Plane Surface(17) = {17};

Line Loop(18) = {75, 72, -74, -78};
Plane Surface(18) = {18};

Line Loop(19) = {73, -77, 79, 76};
Plane Surface(19) = {19};

Line Loop(20) = {-81, 78, 80, -84};
Plane Surface(20) = {20};

Line Loop(21) = {83, 85, -82, -79};
Plane Surface(21) = {21};

Line Loop(22) = {87, 84, -86, -90};
Plane Surface(22) = {22};

Line Loop(23) = {88, -85, -89, 91};
Plane Surface(23) = {23};

Line Loop(24) = {92, -96, -93, 90};
Plane Surface(24) = {24};

Line Loop(25) = {-94, -91, 95, 97};
Plane Surface(25) = {25};

Line Loop(26) = {99, 96, -98, -102};
Plane Surface(26) = {26};

Line Loop(27) = {-101, 103, 100, -97};
Plane Surface(27) = {27};

Line Loop(28) = {102, 104, -108, -105};
Plane Surface(28) = {28};

Line Loop(29) = {-103, 107, 109, -106};
Plane Surface(29) = {29};

Line Loop(30) = {-110, -114, 111, 108};
Plane Surface(30) = {30};

Line Loop(31) = {-113, 115, 112, -109};
Plane Surface(31) = {31};

Line Loop(32) = {-117, 114, 116, -120};
Plane Surface(32) = {32};

Line Loop(33) = {119, 121, -118, -115};
Plane Surface(33) = {33};

Line Loop(34) = {-125, 127, 124, -121};
Plane Surface(34) = {34};

Line Loop(35) = {-122, -126, 123, 120};
Plane Surface(35) = {35};

Line Loop(36) = {128, -132, -129, 126};
Plane Surface(36) = {36};

Line Loop(37) = {-127, 131, 133, -130};
Plane Surface(37) = {37};

Line Loop(38) = {135, 132, -134, -138};
Plane Surface(38) = {38};

Line Loop(39) = {-137, 139, 136, -133};
Plane Surface(39) = {39};

// two large surface for the inductor
//Line Loop(40) = {-54, -43, -61, 68, -74, 80, -86, 92, -98, 104, -110, 116, -122, 128, -134, 141, 136, -130, 124, -118, 112, -106, 100, -94, 88, -82, 76, -70, 63, -60};

Line Loop(40) = {60, -63, 70, -76, 82, -88, 94, -100, 106, -112, 118, -124, 130, -136, -141, 134, -128, 122, -116, 110, -104, 98, -92, 86, -80, 74, -68, 61, 43, 54};
Plane Surface(40) = {40};

//Line Loop(41) = {-2, -58, -64, 71, -77, 83, -89, 95, -101, 107, -113, 119, -125, 131, -137, -140, 135, -129, 123, -117, 111, -105, 99, -93, 87, -81, 75, -69, 62, -39};

Line Loop(41) = {39, -62, 69, -75, 81, -87, 93, -99, 105, -111, 117, -123, 129, -135, 140, 137, -131, 125, -119, 113, -107, 101, -95, 89, -83, 77, -71, 64, 58, 2};
Plane Surface(41) = {41};



Line Loop(42) = {-187, 189, 186, -188};
Plane Surface(42) = {42};  // port 1 Surface

Line Loop(43) = {148, 147, 142, -146};
Plane Surface(43) = {43};  // port 2 Surface

Line Loop(44) = {183, 187, -182, -181};
Plane Surface(44) = {44};

Line Loop(45) = {-177, 181, 176, -174};
Plane Surface(45) = {45};

Line Loop(46) = {175, -178, -180, 179};
Plane Surface(46) = {46};

Line Loop(47) = {171, 174, -170, -168};
Plane Surface(47) = {47};

Line Loop(48) = {172, -175, -173, 169};
Plane Surface(48) = {48};

Line Loop(49) = {-165, 168, 164, -162};
Plane Surface(49) = {49};

Line Loop(50) = {-166, -169, 167, 163};
Plane Surface(50) = {50};

Line Loop(51) = {159, 162, -158, -156};
Plane Surface(51) = {51};

Line Loop(52) = {157, 160, -163, -161};
Plane Surface(52) = {52};

Line Loop(53) = {-153, 156, 152, -143};
Plane Surface(53) = {53};

Line Loop(54) = {-157, 155, -151, -154};
Plane Surface(54) = {54};

Line Loop(55) = {144, 143, -145, -142};
Plane Surface(55) = {55};

Line Loop(56) = {149, 151, -150, -148};
Plane Surface(56) = {56};

Line Loop(57) = {184, 180, -185, -186};
Plane Surface(57) = {57};

Line Loop(58) = {182, 188, -184, 178, -172, 166, -160, 154, -149, 146, -144, -152, 158, -164, 170, -176};
Plane Surface(58) = {58};

Line Loop(59) = {185, -189, -183, 177, -171, 165, -159, 153, 145, -147, 150, -155, 161, -167, 173, -179};
Plane Surface(59) = {59};



// port mapping file should be
// 2
// 1 3
// 2 4
// port defination
Physical Surface("Port1") = {2};
Physical Surface("Port2") = {3};
Physical Surface("Surface1") = {10, 11, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41};
Physical Surface("Port3") = {42};
Physical Surface("Port4") = {43};
Physical Surface("Surface2") = {44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59};
