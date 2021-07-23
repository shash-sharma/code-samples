
// Multiple conductor crossover bus inspired by:
// Zhu, Song, White - Algorithms in FastImp - A Fast and Wide-Band Impedance Extraction Program for Complicated 3-D Geometries - 2005
// and
// Lee et. al - Design and Signal Integrity Analysis of HBM Interposer in 2.5D TB Bandwidth Graphics Module - 2015

// Shashwat Sharma
// Mar 06, 2018

scale = 1.0e-1;
lc = scale*0.25e-3;
l_gnd = scale*5e-3; // length of each conductor
l_line = scale*6e-3; // length of each conductor
w_gnd = scale*0.20e-3; // width of each conductor
h_gnd = scale*0.06e-3; // height of each conductor
w_line = scale*0.09e-3; // width of each conductor
h_line = scale*0.08e-3; // height of each conductor
bend_shift = 25*w_line; // lateral shift between the two ends of the bends
bend_point = 0.2*l_line; // how far along the line the bend begins and ends

z1 = scale*1.0e-3; // z coordinate of the centre of the first (top-most) layer
z2 = scale*0.6e-3; // z coordinate of the centre of the second layer
z3 = scale*0.2e-3; // z coordinate of the centre of the third layer
z4 = scale*1.6e-3; // z coordinate of the top part of the z-varying layer
z5 = scale*1.6e-3; // z coordinate of the bottom part of the z-varying layer

// z5 = scale*0.9e-3; // z coordinate of the bottom part of the z-varying layer

separation = scale*0.33e-3; // centre-to-centre lateral separation between conductors

nConductors1 = 9; // number of conductors in layer 1
nConductors2 = 13; // number of conductors in layer 2
nConductors3 = 9; // number of conductors in layer 3
// nConductors1 = 0; // number of conductors in layer 1
// nConductors2 = 0; // number of conductors in layer 2
// nConductors3 = 0; // number of conductors in layer 3
nConductors4 = 0; // number of conductors in z-varying layer

// specify mesh size
Mesh.CharacteristicLengthExtendFromBoundary = 0;
//Mesh.CharacteristicLengthMax = 4*lc;

// Array properties
nBuses = 2;
array_sep = 0.04*l_line;

kk = 1;
mm = 1;
nn = 1;
jj = 1;
qq = 1;

For ww In {0:nBuses-1}

    // First layer - with bends
    x_bl = 0.0 - ww*(array_sep+l_gnd-bend_shift-array_sep); y_bl = 0.0; z_bl = z1;
    x_br = x_bl+w_line; y_br = y_bl; z_br = z_bl;
    x_tr = x_br; y_tr = y_bl; z_tr = z_br+h_line;
    x_tl = x_bl; y_tl = y_bl; z_tl = z_tr;

    x2_bl = x_bl; y2_bl = y_bl+bend_point; z2_bl = z_bl;
    x2_br = x2_bl+w_line; y2_br = y2_bl; z2_br = z2_bl;
    x2_tr = x2_br; y2_tr = y2_bl; z2_tr = z2_br+h_line;
    x2_tl = x2_bl; y2_tl = y2_bl; z2_tl = z2_tr;

    x3_bl = x_bl+bend_shift*(-1)^ww; y3_bl = y_bl+(l_line-bend_point); z3_bl = z_bl;
    x3_br = x3_bl+w_line; y3_br = y3_bl; z3_br = z3_bl;
    x3_tr = x3_br; y3_tr = y3_bl; z3_tr = z3_br+h_line;
    x3_tl = x3_bl; y3_tl = y3_bl; z3_tl = z3_tr;

    x4_bl = x_bl+bend_shift*(-1)^ww; y4_bl = y_bl+l_line; z4_bl = z_bl;
    x4_br = x4_bl+w_line; y4_br = y4_bl; z4_br = z4_bl;
    x4_tr = x4_br; y4_tr = y4_bl; z4_tr = z4_br+h_line;
    x4_tl = x4_bl; y4_tl = y4_bl; z4_tl = z4_tr;


    For ii In {0:nConductors1-1}

        Point(kk) = {x_bl+separation*ii, y_bl, z_bl, lc}; kk++;
        Point(kk) = {x_br+separation*ii, y_br, z_br, lc}; kk++;
        Point(kk) = {x_tr+separation*ii, y_tr, z_tr, lc}; kk++;
        Point(kk) = {x_tl+separation*ii, y_tl, z_tl, lc}; kk++;

        Line(mm) = {kk-4,kk-3}; mm++;
        Line(mm) = {kk-3,kk-2}; mm++;
        Line(mm) = {kk-2,kk-1}; mm++;
        Line(mm) = {kk-1,kk-4}; mm++;

        Line Loop(nn) = {-(mm-1), -(mm-2), -(mm-3), -(mm-4)}; nn++;


        Point(kk) = {x2_bl+separation*ii, y2_bl, z2_bl, lc}; kk++;
        Point(kk) = {x2_br+separation*ii, y2_br, z2_br, lc}; kk++;
        Point(kk) = {x2_tr+separation*ii, y2_tr, z2_tr, lc}; kk++;
        Point(kk) = {x2_tl+separation*ii, y2_tl, z2_tl, lc}; kk++;

        Line(mm) = {kk-4,kk-3}; mm++;
        Line(mm) = {kk-3,kk-2}; mm++;
        Line(mm) = {kk-2,kk-1}; mm++;
        Line(mm) = {kk-1,kk-4}; mm++;

        Line(mm) = {kk-8,kk-4}; mm++;
        Line(mm) = {kk-7,kk-3}; mm++;
        Line(mm) = {kk-6,kk-2}; mm++;
        Line(mm) = {kk-5,kk-1}; mm++;

        Line Loop(nn) = {mm-1, -(mm-6), -(mm-2), mm-10}; nn++;
        Line Loop(nn) = {-(mm-4), mm-12, mm-3, -(mm-8)}; nn++;
        Line Loop(nn) = {mm-2, -(mm-7), -(mm-3), mm-11}; nn++;
        Line Loop(nn) = {-(mm-1), mm-9, mm-4, -(mm-5)}; nn++;


        Point(kk) = {x3_bl+separation*ii, y3_bl, z3_bl, lc}; kk++;
        Point(kk) = {x3_br+separation*ii, y3_br, z3_br, lc}; kk++;
        Point(kk) = {x3_tr+separation*ii, y3_tr, z3_tr, lc}; kk++;
        Point(kk) = {x3_tl+separation*ii, y3_tl, z3_tl, lc}; kk++;

        Line(mm) = {kk-4,kk-3}; mm++;
        Line(mm) = {kk-3,kk-2}; mm++;
        Line(mm) = {kk-2,kk-1}; mm++;
        Line(mm) = {kk-1,kk-4}; mm++;

        Line(mm) = {kk-8,kk-4}; mm++;
        Line(mm) = {kk-7,kk-3}; mm++;
        Line(mm) = {kk-6,kk-2}; mm++;
        Line(mm) = {kk-5,kk-1}; mm++;

        Line Loop(nn) = {mm-1, -(mm-6), -(mm-2), mm-14}; nn++;
        Line Loop(nn) = {-(mm-4), mm-16, mm-3, -(mm-8)}; nn++;
        Line Loop(nn) = {mm-2, -(mm-7), -(mm-3), mm-15}; nn++;
        Line Loop(nn) = {-(mm-1), mm-13, mm-4, -(mm-5)}; nn++;

        Point(kk) = {x4_bl+separation*ii, y4_bl, z4_bl, lc}; kk++;
        Point(kk) = {x4_br+separation*ii, y4_br, z4_br, lc}; kk++;
        Point(kk) = {x4_tr+separation*ii, y4_tr, z4_tr, lc}; kk++;
        Point(kk) = {x4_tl+separation*ii, y4_tl, z4_tl, lc}; kk++;

        Line(mm) = {kk-4,kk-3}; mm++;
        Line(mm) = {kk-3,kk-2}; mm++;
        Line(mm) = {kk-2,kk-1}; mm++;
        Line(mm) = {kk-1,kk-4}; mm++;

        Line Loop(nn) = {mm-4, mm-3, mm-2, mm-1}; nn++;

        Line(mm) = {kk-8,kk-4}; mm++;
        Line(mm) = {kk-7,kk-3}; mm++;
        Line(mm) = {kk-6,kk-2}; mm++;
        Line(mm) = {kk-5,kk-1}; mm++;

        Line Loop(nn) = {mm-1, -(mm-6), -(mm-2), mm-14}; nn++;
        Line Loop(nn) = {-(mm-4), mm-16, mm-3, -(mm-8)}; nn++;
        Line Loop(nn) = {mm-2, -(mm-7), -(mm-3), mm-15}; nn++;
        Line Loop(nn) = {-(mm-1), mm-13, mm-4, -(mm-5)}; nn++;

        Plane Surface(jj) = {nn-14}; jj++;
        Plane Surface(jj) = {nn-5}; jj++;
        Plane Surface(jj) = {nn-13}; jj++;
        Plane Surface(jj) = {nn-12}; jj++;
        Plane Surface(jj) = {nn-11}; jj++;
        Plane Surface(jj) = {nn-10}; jj++;
        Plane Surface(jj) = {nn-9}; jj++;
        Plane Surface(jj) = {nn-8}; jj++;
        Plane Surface(jj) = {nn-7}; jj++;
        Plane Surface(jj) = {nn-6}; jj++;
        Plane Surface(jj) = {nn-4}; jj++;
        Plane Surface(jj) = {nn-3}; jj++;
        Plane Surface(jj) = {nn-2}; jj++;
        Plane Surface(jj) = {nn-1}; jj++;

        If (qq == 7 || qq == 10 || qq == 1 || qq == 118)
            Physical Surface(Sprintf("Port%g",qq)) = {jj-14}; qq++;
        Else
            Physical Surface(Sprintf("Surface%g",qq)) = {jj-14}; qq++;
        EndIf

        If (qq == 8 || qq == 11 || qq == 2 || qq == 119)
            Physical Surface(Sprintf("Port%g",qq)) = {jj-13}; qq++;
        Else
            Physical Surface(Sprintf("Surface%g",qq)) = {jj-13}; qq++;
        EndIf

        // If (qq != 1 && qq != 79 && qq != 2 && qq != 80 && qq != 7 && qq != 10 && qq != 70 && qq != 73 && qq != 8 && qq != 11 && qq != 71 && qq != 74)
            Physical Surface(Sprintf("Surface%g",qq)) = {jj-12, jj-11, jj-10, jj-9, jj-8, jj-7, jj-6, jj-5, jj-4, jj-3, jj-2, jj-1}; qq++;
        // EndIf

    EndFor


    // Second layer
    x_bl = -(l_gnd - (2*nConductors2-1)*w_gnd)/2 - ww*(array_sep+l_gnd); y_bl = w_gnd+(l_gnd - (2*nConductors2-1)*w_gnd)/2; z_bl = z2;
    x_br = x_bl; y_br = y_bl-w_gnd; z_br = z_bl;
    x_tr = x_br; y_tr = y_br; z_tr = z_br+h_gnd;
    x_tl = x_bl; y_tl = y_bl; z_tl = z_tr;

    x2_bl = x_bl+l_gnd; y2_bl = y_bl; z2_bl = z_bl;
    x2_br = x2_bl; y2_br = y2_bl-w_gnd; z2_br = z2_bl;
    x2_tr = x2_br; y2_tr = y2_br; z2_tr = z2_br+h_gnd;
    x2_tl = x2_bl; y2_tl = y2_bl; z2_tl = z2_tr;

    For ii In {0:nConductors2-1}

        Point(kk) = {x_bl, y_bl+separation*ii, z_bl, lc}; kk++;
        Point(kk) = {x_br, y_br+separation*ii, z_br, lc}; kk++;
        Point(kk) = {x_tr, y_tr+separation*ii, z_tr, lc}; kk++;
        Point(kk) = {x_tl, y_tl+separation*ii, z_tl, lc}; kk++;

        Line(mm) = {kk-4,kk-3}; mm++;
        Line(mm) = {kk-3,kk-2}; mm++;
        Line(mm) = {kk-2,kk-1}; mm++;
        Line(mm) = {kk-1,kk-4}; mm++;

        Line Loop(nn) = {-(mm-1), -(mm-2), -(mm-3), -(mm-4)}; nn++;

        Point(kk) = {x2_bl, y2_bl+separation*ii, z2_bl, lc}; kk++;
        Point(kk) = {x2_br, y2_br+separation*ii, z2_br, lc}; kk++;
        Point(kk) = {x2_tr, y2_tr+separation*ii, z2_tr, lc}; kk++;
        Point(kk) = {x2_tl, y2_tl+separation*ii, z2_tl, lc}; kk++;

        Line(mm) = {kk-4,kk-3}; mm++;
        Line(mm) = {kk-3,kk-2}; mm++;
        Line(mm) = {kk-2,kk-1}; mm++;
        Line(mm) = {kk-1,kk-4}; mm++;

        Line Loop(nn) = {mm-4, mm-3, mm-2, mm-1}; nn++;

        Line(mm) = {kk-8,kk-4}; mm++;
        Line(mm) = {kk-7,kk-3}; mm++;
        Line(mm) = {kk-6,kk-2}; mm++;
        Line(mm) = {kk-5,kk-1}; mm++;

        Line Loop(nn) = {mm-1, -(mm-6), -(mm-2), mm-10}; nn++;
        Line Loop(nn) = {-(mm-4), mm-12, mm-3, -(mm-8)}; nn++;
        Line Loop(nn) = {mm-2, -(mm-7), -(mm-3), mm-11}; nn++;
        Line Loop(nn) = {-(mm-1), mm-9, mm-4, -(mm-5)}; nn++;

        Plane Surface(jj) = {nn-6}; jj++;
        Plane Surface(jj) = {nn-5}; jj++;
        Plane Surface(jj) = {nn-4}; jj++;
        Plane Surface(jj) = {nn-3}; jj++;
        Plane Surface(jj) = {nn-2}; jj++;
        Plane Surface(jj) = {nn-1}; jj++;

        Physical Surface(Sprintf("Surface%g",qq)) = {jj-6}; qq++;
        Physical Surface(Sprintf("Surface%g",qq)) = {jj-5}; qq++;
        Physical Surface(Sprintf("Surface%g",qq)) = {jj-4, jj-3, jj-2, jj-1}; qq++;

    EndFor


    // Third layer
    x_bl = 0.0 - ww*(array_sep+l_gnd - nConductors3*w_line - 4*separation); y_bl = 0.0; z_bl = z3;
    x_br = x_bl+w_line; y_br = y_bl; z_br = z_bl;
    x_tr = x_br; y_tr = y_bl; z_tr = z_br+h_line;
    x_tl = x_bl; y_tl = y_bl; z_tl = z_tr;

    x2_bl = x_bl; y2_bl = y_bl+l_line; z2_bl = z_bl;
    x2_br = x2_bl+w_line; y2_br = y2_bl; z2_br = z2_bl;
    x2_tr = x2_br; y2_tr = y2_bl; z2_tr = z2_br+h_line;
    x2_tl = x2_bl; y2_tl = y2_bl; z2_tl = z2_tr;

    For ii In {0:nConductors3-1}

        Point(kk) = {x_bl+separation*ii, y_bl, z_bl, lc}; kk++;
        Point(kk) = {x_br+separation*ii, y_br, z_br, lc}; kk++;
        Point(kk) = {x_tr+separation*ii, y_tr, z_tr, lc}; kk++;
        Point(kk) = {x_tl+separation*ii, y_tl, z_tl, lc}; kk++;

        Line(mm) = {kk-4,kk-3}; mm++;
        Line(mm) = {kk-3,kk-2}; mm++;
        Line(mm) = {kk-2,kk-1}; mm++;
        Line(mm) = {kk-1,kk-4}; mm++;

        Line Loop(nn) = {-(mm-1), -(mm-2), -(mm-3), -(mm-4)}; nn++;

        Point(kk) = {x2_bl+separation*ii, y2_bl, z2_bl, lc}; kk++;
        Point(kk) = {x2_br+separation*ii, y2_br, z2_br, lc}; kk++;
        Point(kk) = {x2_tr+separation*ii, y2_tr, z2_tr, lc}; kk++;
        Point(kk) = {x2_tl+separation*ii, y2_tl, z2_tl, lc}; kk++;

        Line(mm) = {kk-4,kk-3}; mm++;
        Line(mm) = {kk-3,kk-2}; mm++;
        Line(mm) = {kk-2,kk-1}; mm++;
        Line(mm) = {kk-1,kk-4}; mm++;

        Line Loop(nn) = {mm-4, mm-3, mm-2, mm-1}; nn++;

        Line(mm) = {kk-8,kk-4}; mm++;
        Line(mm) = {kk-7,kk-3}; mm++;
        Line(mm) = {kk-6,kk-2}; mm++;
        Line(mm) = {kk-5,kk-1}; mm++;

        Line Loop(nn) = {mm-1, -(mm-6), -(mm-2), mm-10}; nn++;
        Line Loop(nn) = {-(mm-4), mm-12, mm-3, -(mm-8)}; nn++;
        Line Loop(nn) = {mm-2, -(mm-7), -(mm-3), mm-11}; nn++;
        Line Loop(nn) = {-(mm-1), mm-9, mm-4, -(mm-5)}; nn++;

        Plane Surface(jj) = {nn-6}; jj++;
        Plane Surface(jj) = {nn-5}; jj++;
        Plane Surface(jj) = {nn-4}; jj++;
        Plane Surface(jj) = {nn-3}; jj++;
        Plane Surface(jj) = {nn-2}; jj++;
        Plane Surface(jj) = {nn-1}; jj++;

        If (qq == 71 || qq == 74)
            Physical Surface(Sprintf("Port%g",qq)) = {jj-6}; qq++;
        Else
            Physical Surface(Sprintf("Surface%g",qq)) = {jj-6}; qq++;
        EndIf

        If (qq == 70 || qq == 73)    
            Physical Surface(Sprintf("Port%g",qq)) = {jj-5}; qq++;
        Else
            Physical Surface(Sprintf("Surface%g",qq)) = {jj-5}; qq++;
        EndIf
        
        Physical Surface(Sprintf("Surface%g",qq)) = {jj-4, jj-3, jj-2, jj-1}; qq++;

    EndFor

EndFor

bend_point = 3*bend_point;

// z-varying layer - with bends in the z-direction
x_bl = 0.0 - (l_gnd-bend_shift) + 5*separation + (nConductors4-1/2)*w_line; y_bl = 0.0; z_bl = z4;
x_br = x_bl+w_line; y_br = y_bl; z_br = z_bl;
x_tr = x_br; y_tr = y_bl; z_tr = z_br+h_line;
x_tl = x_bl; y_tl = y_bl; z_tl = z_tr;

x2_bl = x_bl; y2_bl = y_bl+bend_point-w_line; z2_bl = z_bl;
x2_br = x2_bl+w_line; y2_br = y2_bl; z2_br = z2_bl;
x2_tr = x2_br; y2_tr = y2_bl; z2_tr = z2_br+h_line;
x2_tl = x2_bl; y2_tl = y2_bl; z2_tl = z2_tr;

x3_bl = x_bl; y3_bl = y_bl+bend_point; z3_bl = z_bl;
x3_br = x3_bl+w_line; y3_br = y3_bl; z3_br = z3_bl;
x3_tr = x3_br; y3_tr = y3_bl; z3_tr = z3_br+h_line;
x3_tl = x3_bl; y3_tl = y3_bl; z3_tl = z3_tr;

x4_bl = x_bl; y4_bl = y3_bl-w_line; z4_bl = z5;
x4_br = x4_bl+w_line; y4_br = y4_bl; z4_br = z5;
x4_tr = x4_br; y4_tr = y4_bl; z4_tr = z4_br+h_line;
x4_tl = x4_bl; y4_tl = y4_bl; z4_tl = z4_tr;

x5_bl = x_bl; y5_bl = y3_bl; z5_bl = z5;
x5_br = x5_bl+w_line; y5_br = y5_bl; z5_br = z5;
x5_tr = x5_br; y5_tr = y5_bl; z5_tr = z5_br+h_line;
x5_tl = x5_bl; y5_tl = y5_bl; z5_tl = z5_tr;

x6_bl = x_bl; y6_bl = y2_bl+(l_line-bend_point+w_line); z6_bl = z5;
x6_br = x6_bl+w_line; y6_br = y6_bl; z6_br = z5;
x6_tr = x6_br; y6_tr = y6_bl; z6_tr = z6_br+h_line;
x6_tl = x6_bl; y6_tl = y6_bl; z6_tl = z6_tr;

For ii In {0:nConductors4-1}

    Point(kk) = {x_bl+separation*ii, y_bl, z_bl, lc}; kk++;
    Point(kk) = {x_br+separation*ii, y_br, z_br, lc}; kk++;
    Point(kk) = {x_tr+separation*ii, y_tr, z_tr, lc}; kk++;
    Point(kk) = {x_tl+separation*ii, y_tl, z_tl, lc}; kk++;
    Line(mm) = {kk-4,kk-3}; mm++;
    Line(mm) = {kk-3,kk-2}; mm++;
    Line(mm) = {kk-2,kk-1}; mm++;
    Line(mm) = {kk-1,kk-4}; mm++;

    Line Loop(nn) = {-(mm-1), -(mm-2), -(mm-3), -(mm-4)}; nn++;

    Point(kk) = {x2_bl+separation*ii, y2_bl, z2_bl, lc}; kk++;
    Point(kk) = {x2_br+separation*ii, y2_br, z2_br, lc}; kk++;
    Point(kk) = {x2_tr+separation*ii, y2_tr, z2_tr, lc}; kk++;
    Point(kk) = {x2_tl+separation*ii, y2_tl, z2_tl, lc}; kk++;
    Line(mm) = {kk-4,kk-3}; mm++;
    Line(mm) = {kk-3,kk-2}; mm++;
    Line(mm) = {kk-2,kk-1}; mm++;
    Line(mm) = {kk-1,kk-4}; mm++;

    Line(mm) = {kk-8,kk-4}; mm++;
    Line(mm) = {kk-7,kk-3}; mm++;
    Line(mm) = {kk-6,kk-2}; mm++;
    Line(mm) = {kk-5,kk-1}; mm++;

    Line Loop(nn) = {mm-1, -(mm-6), -(mm-2), mm-10}; nn++;
    Line Loop(nn) = {-(mm-4), mm-12, mm-3, -(mm-8)}; nn++;
    Line Loop(nn) = {mm-2, -(mm-7), -(mm-3), mm-11}; nn++;
    Line Loop(nn) = {-(mm-1), mm-9, mm-4, -(mm-5)}; nn++;

    Point(kk) = {x3_bl+separation*ii, y3_bl, z3_bl, lc}; kk++;
    Point(kk) = {x3_br+separation*ii, y3_br, z3_br, lc}; kk++;
    Point(kk) = {x3_tr+separation*ii, y3_tr, z3_tr, lc}; kk++;
    Point(kk) = {x3_tl+separation*ii, y3_tl, z3_tl, lc}; kk++;
    Line(mm) = {kk-4,kk-3}; mm++;
    Line(mm) = {kk-3,kk-2}; mm++;
    Line(mm) = {kk-2,kk-1}; mm++;
    Line(mm) = {kk-1,kk-4}; mm++;

    Line(mm) = {kk-8,kk-4}; mm++;
    Line(mm) = {kk-7,kk-3}; mm++;
    Line(mm) = {kk-6,kk-2}; mm++;
    Line(mm) = {kk-5,kk-1}; mm++;

    Line Loop(nn) = {-(mm-6), -(mm-2), mm-14, mm-1}; nn++;
    Line Loop(nn) = {mm-2, -(mm-7), -(mm-3), mm-15}; nn++;
    Line Loop(nn) = {mm-7, mm-6, mm-5, mm-8}; nn++;
    Line Loop(nn) = {-(mm-5), -(mm-1), mm-13, mm-4}; nn++;


    Point(kk) = {x4_bl+separation*ii, y4_bl, z4_bl, lc}; kk++;
    Point(kk) = {x4_br+separation*ii, y4_br, z4_br, lc}; kk++;
    Point(kk) = {x4_tr+separation*ii, y4_tr, z4_tr, lc}; kk++;
    Point(kk) = {x4_tl+separation*ii, y4_tl, z4_tl, lc}; kk++;
    Line(mm) = {kk-4,kk-3}; mm++;
    Line(mm) = {kk-3,kk-2}; mm++;
    Line(mm) = {kk-2,kk-1}; mm++;
    Line(mm) = {kk-1,kk-4}; mm++;


    Point(kk) = {x5_bl+separation*ii, y5_bl, z5_bl, lc}; kk++;
    Point(kk) = {x5_br+separation*ii, y5_br, z5_br, lc}; kk++;
    Point(kk) = {x5_tr+separation*ii, y5_tr, z5_tr, lc}; kk++;
    Point(kk) = {x5_tl+separation*ii, y5_tl, z5_tl, lc}; kk++;
    Line(mm) = {kk-4,kk-3}; mm++;
    Line(mm) = {kk-3,kk-2}; mm++;
    Line(mm) = {kk-2,kk-1}; mm++;
    Line(mm) = {kk-1,kk-4}; mm++;

    Line(mm) = {kk-8,kk-4}; mm++;
    Line(mm) = {kk-7,kk-3}; mm++;
    Line(mm) = {kk-6,kk-2}; mm++;
    Line(mm) = {kk-5,kk-1}; mm++;

    Line Loop(nn) = {-(mm-5), -(mm-1), mm-9, (mm-4)}; nn++;
    Line Loop(nn) = {-(mm-10), -(mm-11), -(mm-12), -(mm-9)}; nn++;
    Line Loop(nn) = {mm-11, mm-2, -(mm-7), -(mm-3)}; nn++;
    Line Loop(nn) = {-(mm-4), mm-12, mm-3, -(mm-8)}; nn++;

    Line(mm) = {kk-2,kk-11}; mm++;
    Line(mm) = {kk-1,kk-12}; mm++;
    Line(mm) = {kk-5,kk-16}; mm++;
    Line(mm) = {kk-6,kk-15}; mm++;

    Line Loop(nn) = {(mm-3), -(mm-20), -(mm-2), (mm-5)}; nn++;
    Line Loop(nn) = {(mm-2), (mm-32), -(mm-1), (mm-14)}; nn++;
    Line Loop(nn) = {mm-1, mm-19, -(mm-4), -(mm-6)}; nn++;
    Line Loop(nn) = {(mm-4), -(mm-24), -(mm-3), -(mm-10)}; nn++;

    Point(kk) = {x6_bl+separation*ii, y6_bl, z6_bl, lc}; kk++;
    Point(kk) = {x6_br+separation*ii, y6_br, z6_br, lc}; kk++;
    Point(kk) = {x6_tr+separation*ii, y6_tr, z6_tr, lc}; kk++;
    Point(kk) = {x6_tl+separation*ii, y6_tl, z6_tl, lc}; kk++;
    Line(mm) = {kk-4,kk-3}; mm++;
    Line(mm) = {kk-3,kk-2}; mm++;
    Line(mm) = {kk-2,kk-1}; mm++;
    Line(mm) = {kk-1,kk-4}; mm++;

    Line Loop(nn) = {(mm-1), (mm-4), (mm-3), (mm-2)}; nn++;

    Line(mm) = {kk-8,kk-4}; mm++;
    Line(mm) = {kk-7,kk-3}; mm++;
    Line(mm) = {kk-6,kk-2}; mm++;
    Line(mm) = {kk-5,kk-1}; mm++;

    Line Loop(nn) = {(mm-1), -(mm-6), -(mm-2), mm-18}; nn++;
    Line Loop(nn) = {-(mm-1), mm-17, mm-4, -(mm-5)}; nn++;
    Line Loop(nn) = {-(mm-4), (mm-20), (mm-3), -(mm-8)}; nn++;
    Line Loop(nn) = {(mm-2), -(mm-7), -(mm-3), (mm-19)}; nn++;

    Plane Surface(jj) = {nn-22}; jj++;
    Plane Surface(jj) = {nn-5}; jj++;
    Plane Surface(jj) = {nn-21}; jj++;
    Plane Surface(jj) = {nn-20}; jj++;
    Plane Surface(jj) = {nn-19}; jj++;
    Plane Surface(jj) = {nn-18}; jj++;
    Plane Surface(jj) = {nn-17}; jj++;
    Plane Surface(jj) = {nn-16}; jj++;
    Plane Surface(jj) = {nn-15}; jj++;
    Plane Surface(jj) = {nn-14}; jj++;
    Plane Surface(jj) = {nn-13}; jj++;
    Plane Surface(jj) = {nn-12}; jj++;
    Plane Surface(jj) = {nn-11}; jj++;
    Plane Surface(jj) = {nn-10}; jj++;
    Plane Surface(jj) = {nn-9}; jj++;
    Plane Surface(jj) = {nn-8}; jj++;
    Plane Surface(jj) = {nn-7}; jj++;
    Plane Surface(jj) = {nn-6}; jj++;
    Plane Surface(jj) = {nn-4}; jj++;
    Plane Surface(jj) = {nn-3}; jj++;
    Plane Surface(jj) = {nn-2}; jj++;
    Plane Surface(jj) = {nn-1}; jj++;

    If (qq == 193 || qq == 196)
        Physical Surface(Sprintf("Port%g",qq)) = {jj-22}; qq++;
    Else
        Physical Surface(Sprintf("Surface%g",qq)) = {jj-22}; qq++;
    EndIf

    If (qq == 194 || qq == 197)
        Physical Surface(Sprintf("Port%g",qq)) = {jj-21}; qq++;
    Else
        Physical Surface(Sprintf("Surface%g",qq)) = {jj-21}; qq++;
    EndIf

    Physical Surface(Sprintf("Surface%g",qq)) = {jj-20, jj-19, jj-18, jj-17, jj-16, jj-15, jj-14, jj-13,jj-12, jj-11, jj-10, jj-9, jj-8, jj-7, jj-6, jj-5, jj-4, jj-3, jj-2, jj-1}; qq++;

EndFor

Coherence;
