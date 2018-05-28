!**********************************************************************
! Test program:
! Example of the use of the module analysis module with a simple example.
! The test is described in the appendix of 
! Torres-Forne et al MNRAS 474, 6272 (2018)
!**********************************************************************
!
! GPL  License:
! -------------
!
!    GREAT = General Relativistic Eigenmode Analysis Tool
!    Copyright (C) 2018 Alejandro Torres-Forne and Pablo Cerda-Duran
!
!    This file is part of GREAT.
!
!    GREAT is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    GREAT is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with GREAT.  If not, see <http://www.gnu.org/licenses/>.
!
!************************************************************************

program main

  use module_param
  use module_background
  use module_eigen
  use module_mode_analysis


  implicit none

  integer :: i, nt

  !> input/output directory
  character(len=255) :: output_directory,input_directory, parfile

  !> parameters 
  type (parameters_t) :: param  
  !> backgroun data (1D profiles)
  type(BGData_t) :: data

  !> error flag
  integer :: ierr

  !> Size of the grid (without ghost cells)
  integer :: m
  !> index location of the shock (last cell inside the shock)
  integer :: iR

  !> variables needed for this example
  real*8 :: delta_r
  real*8 :: gamma_one


  ! ------------------------------
  input_directory = "./Data/"
  output_directory = "./output/"
  parfile = "./parameters"
  m = 300
  ! ------------------------------


  ! ---- Read parameter file  -------
  call Read_parameters(param, parfile, ierr)
  if (ierr.ne.0) stop
  

  ! ------ nt loops over different times of the simulation. In this example ---
  ! ------ we compute only one time -------------------------------------------
  do nt = 1, 1
          
     ! ...... Allocate arrays ................
     call Init_background_data (data, m, ierr)
     

     ! ..... not used. Just to keep track of where the data comes from 
     data%input_directory = input_directory
     ! ..... Output directory. 
     data%output_directory = output_directory
    
     ! time index
     data%nt = nt
     ! simulation time (since this example only has one time we set it to zero)
     data%time = 0.0d0
     ! index location of the shock (last cell inside the shock)
     ! for this example we choose it to be equal to the last cell
     ! in the grid, but does not have to be like this.
     data%iR = m
     

     delta_r = 1.0d0 / dble(data%m)
     gamma_one = 4.0d0 / 3.0d0
     
     ! ..... Fill arrays with background data ............
     ! We use units of G=c=1 and the remaining length scale in cm.
     do i = 1, m
        ! ... radius (in cm)
        data%r (i) =  (dble(i) - 0.5d0) * delta_r
        ! ... rest mass density (cm^-2)
        data%rho (i) = 1.0d0
        ! ... sound speed squared (dimensionless)
        data%c_sound_squared (i) = 1.0d0
        ! ... pressure (cm^-2)
        data%p (i) = data%c_sound_squared (i) * data%rho (i) / gamma_one
        ! ... specific internal energy (dimensionless)
        data%eps (i) = data%p (i) / data%rho (i) / (gamma_one - 1.0d0)
        
        ! ---- Metric: set to mikowsky -----
        ! conformal factor (dimensionless)
        data%phi (i) = 1.0d0
        ! lapse (dimensionless)
        data%alpha (i) = 1.0d0
     enddo

     ! If this flag is set to true, it indicates that the data does not contain the value
     ! of N^2 (Brunt-Vaisala frequency) and will be computed internally.
     data%calculate_n2 = .true.

     ! ...... Fill ghost cells ...........
     call Impose_boundary_conditions_bg (data, ierr)

     ! ....... Compute eigenmodes ........
     ! It computes eigenfrequencies and eigenmodes.
     ! Eigenfrequencies are appended to the freqs.dat file.
     ! Eigenfunctions for each time/nt are stored into  
     ! eigen_<nt>.dat. All output in the output_directory.
     call Perform_mode_analysis (data, param, ierr)

     ! ......... Prints an error message if ierr != 0 .....
     call Print_error_message (ierr)


     ! ........ Output background quantities .........
     ! It includes not only the values set above but also
     ! internally calculated quantities such as the Brunt-Vaisala
     ! frequency and the Lamb frequency
     call Output_backgroud_data(data, param, 2, ierr)

     ! ......... Deallocate arrays in data ..........
     call Free_background_data (data, ierr)

  enddo

end program main



