/* 
 * Copyright (C) 2007, 2008, 2009, 2010, 2011 M. T. Homer Reid
 * This file is part of libAmosBessel. 
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/*
 * libAmosBessel.h  -- a C++ wrapper around the original fortran routines 
 *                  -- of T. E. Amos for evaluating bessel-type functions
 *                  -- at complex arguments 
 *
 * homer reid       -- 2/2011
 *
 *
 * --------------------------------------------------------------
 *
 * this library exports two routines: 
 *  
 * (1) int AmosBessel(char WhichFunction, cdouble z,
 *                    double MinOrder, int NumOrders,
 *                    int Scale, cdouble *f);
 *
 * (2) int AmosAiry(char WhichFunction, cdouble z,
 *                  int Scale, cdouble *f);
 * 
 * --------------------------------------------------------------
 * 
 * the function parameters are:
 *
 *  WhichFunction: a character that selects which function you 
 *                 want to compute (see below) 
 *
 *              z: the complex argument of the bessel/airy function
 *
 *       MinOrder: 
 *      NumOrders: bessel functions of orders 
 *                 MinOrder, MinOrder+1, ... MinOrder+NumOrders-1
 *                 will be computed. (Note that MinOrder need not be 
 *                 an integer or even a half-integer, but may be  
 *                 any arbitrary real number >=0 ).
 *                 
 *          Scale: if Scale==1, scaled versions of the bessel functions 
 *                 are computed. 
 *                 the definition of 'scaled' depends on the function 
 *                 in question (see below). 
 *
 *              f: on entry, must point to a buffer with room for 
 *                 NumOrders cdoubles (for AmosBessel) or 1 cdouble
 *                 (for AmosAiry).
 *            
 *                 on return, f is filled in with the values of the 
 *                 function requested. for example, if WhichFunction='J',
 *                 MinOrder=3.25, NumOrders=3, then on return we have         
 *            
 *                 f[0]=J_{3.25}(z), f[1]=J_{4.25}(z), f[2]=J_{5.25}(z)
 *            
 *                 where J_{\nu}(z) is the regular cylindrical bessel 
 *                 function of order \nu.
 *            
 * --------------------------------------------------------------
 *            
 *  the return value from both functions is the IERR parameter in
 *  the fortran routine (=0 for a successful computation) 
 *            
 * --------------------------------------------------------------
 *            
 * Values of the 'WhichFunction' field:
 *
 * == for AmosBessel:
 *
 *  'J': regular cylindrical bessel J_n
 *  'Y': irregular cylindrical bessel Y_n
 *  'I': modified regular cylindrical bessel I_n
 *  'K': modified irregular cylindrical bessel K_n
 *  'O': type-1 cylindrical hankel function H^(1)_n
 *  'T': type-2 cylindricak hankel function H^(2)_n
 *
 *  'j': regular cylindrical bessel j_n
 *  'y': irregular cylindrical bessel y_n
 *  'i': modified regular cylindrical bessel i_n
 *  'k': modified irregular cylindrical bessel k_n
 *  'o': type-1 spherical hankel function h^(1)_n
 *  't': type-2 spherical hankel function h^(2)_n
 *
 * == for AmosAiry:
 *
 *  'A': value of Airy function Ai(z)
 *  'a': derivative of Airy function dAi(z) /dz
 *  'B': value of Airy function Bi(z)
 *  'b': derivative of Airy function dBi(z) /dz
 *
 * --------------------------------------------------------------
 *
 * as noted above, the meaning of 'scaling' differs for different
 * functions. if you set Scale=1 on entry to AmosBessel or AmosAiry,  
 * the function returned will be: 
 *
 * WhichFunction | Quantity returned 
 * ---------------------------------
 *
 *     'J'       | e^{ -|im z| } * J_\nu(z) 
 *     'j'       | e^{ -|im z| } * j_\nu(z) 
 *
 *     'Y'       | e^{ -|im z| } * Y_\nu(z) 
 *     'y'       | e^{ -|im z| } * y_\nu(z) 
 *
 *     'I'       | e^{ -|re z| } * I_\nu(z) 
 *     'i'       | e^{ -|re z| } * i_\nu(z) 
 *
 *     'K'       | e^{ +z } * K_\nu(z) 
 *     'k'       | e^{ +z } * k_\nu(z) 
 *
 *     'O'       | e^{ -iz } * H^(1)_\nu(z) 
 *     'o'       | e^{ -iz } * h^(1)_\nu(z) 
 *
 *     'T'       | e^{ +iz } * H^(2)_\nu(z) 
 *     't'       | e^{ +iz } * h^(2)_\nu(z) 
 *
 *
 *     'A'       | e^{ 2/3 * z^(3/2) }* Ai(z)
 *     'a'       | e^{ 2/3 * z^(3/2) }* dAi(z)/dz
 *
 *     'B'       | e^{ -2/3 * z^(3/2) }* Bi(z)
 *     'b'       | e^{ -2/3 * z^(3/2) }* dBi(z)/dz
 */

#ifndef LIBAMOSBESSEL_H
#define LIBAMOSBESSEL_H

#include<complex>

#ifndef cdouble
  typedef _Complex double cdouble;
#endif

int AmosBesselCPP(char WhichFunction, std::complex<double> z, double MinOrder, int NumOrders, 
               int Scale, std::complex<double> *f);

int AmosBessel(char WhichFunction, cdouble z, 
                double MinOrder, int NumOrders, 
                int Scale, cdouble *f);

int AmosAiry(char WhichFunction, cdouble Z, int Scale, cdouble *f);

#endif
