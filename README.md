************************************************************************

#GREAT = General Relativistic Eigenmode Analysis Tool

************************************************************************

  This code computes the oscillation modes (eigenmodes) of a spherically
  symmetric self-gravitating object. It has been developed to study the 
  oscillation modes (p, g and f-modes) of proto-neutron stars formed in the 
  collapse of the cores of massive stars. A detailed description of the equations 
  implemented can be found at:

  A. Torres-Forne, P. Cerda-Duran, A. Passamonti & J.A. Font, MNRAS, 474, 5272 (2018)
  
  A. Torres-Forne, P. Cerda-Duran, A. Passamonti & J.A. Font, in prep. (2018)

  The code comes as a library, which should be called from a main code providing
  the data to me analyzed. We provide a test code (src/test.f90) as an example of
  how to call the routines in the library.


##REQUIREMENTS:

1. Linux of unix-like operating system.

2. Modern fortran compiler (e.g. latest version of gfortran).

3. LAPACK (http://www.netlib.org/lapack/).

##INSTALLATION:

1. Create a "local_settings" file necessary to setup the Makefile. A working
    example can be found at "local_settings_example".

2. Execute "make". Libraries and modules will be created at lib and include,
    respectively

3. To compile the test code (src/test.f90) execute "make test". The executable
    bin/test will be created.
 
##USAGE:

Create a "parameters" file. An example can be found at "parameters_test".

To use it as a library, compile and link your code adding:

    /<PATH TO GREAT>/lib/libgreat.a
    
    -I/<PATH TO GREAT/include/

To use the test code, execute /bin/test.


##CREDITS

Copyright (C) 2018 Alejandro Torres-Forne and Pablo Cerda-Duran

Copyright (C) 2009 Anthony SCEMAMA for the content of the file src/module_fast_inv.f90

##LICENSE

GREAT is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

GREAT is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with GREAT.  If not, see <http://www.gnu.org/licenses/>.


##CONTACT:

Alejandro Torres-Forne and Pablo Cerda-Duran

Universitat de AVlencia

c/ Dr. Moliner, 50

E46100 Burjassot (Valencia), Spain

alejandro.torres(at)uv.es

pablo.cerda(at)uv.es
