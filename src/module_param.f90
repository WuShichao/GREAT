!************************************************************************
! Module to manage the parameters used in the eigenmode computation
!************************************************************************
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

module module_param

  implicit none

  private

  type, public :: parameters_t
     !> l = degree of the spherical harmonic decomposition
     integer :: l
     !> Compute with Cowling approximation
     logical :: cowling
     !> Compute with delta_phi 
     logical :: cal_dphi
     !>  Mode approximation:
     !>      1.All modes
     !>      2.Compute P mode approximation
     !>      3.Compute G mode approximation 
     integer :: mode
     !> Data format
     !> 1. Janka
     !> 2. Coconut (silo)
     !> 3. Martin
     !> 4. Haakon
     integer :: input_mode
     !> Minimun/maximum value of the frequency for the eigenvalue search
     real*8 :: fmin, fmax
     !> Units for the frequency
     !> 0: code units (G=C=1 + cm), i.e. cm^-1
     !> 1: Hz
     integer :: funits
     !> Number of frequencies for the coarse eigenfrequency search
     integer :: nwmax
     !> logf_mode
     !>  .true.: Coarse search equally spaced in logarithm of f
     !>  .false.: equally spaced in f
     logical :: logf_mode
     !> Tolerance for eigenvalue calculation
     real*8 :: tol
     !> Maximum number of iterations in eigenvalue calculation
     integer :: maxiter
     !> Interpolate to a new grid
     logical :: interpolate_grid
     !> Interpolation factor ???? (only 2 or 4 allowed)
     integer :: int_factor
     !> Enable restart
     logical :: restart
     !> Index of the first file
     integer :: t_ini
     !> Time-step to read files
     integer :: t_step
     !> Select how to integrate equations
     !>     1. Explicit integration for all equations (6x6)
     !>     2. Explicit integration & elliptic solver
     integer :: integration_mode
     !> Use N^2 from files
     logical :: use_precomputed_n2
     !Select the ggrav computation formula
     !> 1.     ggrav = -dalpha_dr
     !> 2.     ggrav = dp_dr/(rho h)
     !> 3.     Mean value of 1 and 2 
     integer :: ggrav_mode
     !> Enable newtonian calculation of ggrav
     logical :: newtonian
     !> regularize first ip_coeff points
     integer :: ip_coeffs

     !> Parameter to control the transition between full calculation (pfrac=1)
     !> and the p-mode limit (pfrac=0)
     real*8 :: pfrac = 1.0d0
     !> Parameter to control the transition between full calculation (gfrac=1)
     !> and the g-mode limit (gfrac=0)
     real*8 :: gfrac = 1.0d0
     !> Generate files outputing eigenvalue frequencies and eigenfunctions
     logical :: output_eigen = .true.
     !> Verbose mode
     logical :: verb = .true.

  end type parameters_t

  public :: Read_parameters

contains
  !=================================================================
  !   This subroutine reads in parameters
  !=================================================================
  subroutine Read_parameters(param, parfile, ierr)

    implicit none

    ! -- arguments ----
    type (parameters_t), intent (inout) :: param
    character (len=255) :: parfile
    integer, intent (out) :: ierr
    ! ----------------

    ierr = 0

    call get_integer_parameter ("l=", param%l, parfile, ierr)

    call get_logical_parameter ("cowling=", param%cowling, parfile, ierr)

    call get_logical_parameter ("cal_dphi=", param%cal_dphi, parfile, ierr)

    call get_integer_parameter ("mode=", param%mode, parfile, ierr)

    call get_integer_parameter ("input_mode=", param%input_mode, parfile, ierr)

    call get_real_parameter ("fmin=", param%fmin, parfile, ierr)
    call get_real_parameter ("fmax=", param%fmax, parfile, ierr)

    call get_integer_parameter ("funits=", param%funits, parfile, ierr)  

    call get_integer_parameter ("nwmax=", param%nwmax, parfile, ierr)  

    call get_logical_parameter ("logf_mode=", param%logf_mode, parfile, ierr)

    call get_real_parameter ("tol=", param%tol, parfile, ierr)
    call get_integer_parameter ("maxiter=", param%maxiter, parfile, ierr)  


    call get_logical_parameter ("interpolate_grid=", param%interpolate_grid, parfile, ierr)
    call get_integer_parameter ("int_factor=", param%int_factor, parfile, ierr)

    call get_logical_parameter ("restart=", param%restart, parfile, ierr)
    call get_integer_parameter ("t_ini=", param%t_ini, parfile, ierr)
    call get_integer_parameter ("time_step=", param%t_step, parfile, ierr)


    call get_integer_parameter ("integration_mode", param%integration_mode, parfile, ierr)

    call get_logical_parameter ("use_precomputed_n2=", param%use_precomputed_n2, parfile, ierr)

    call get_integer_parameter ("ggrav_mode", param%ggrav_mode, parfile, ierr)

    call get_logical_parameter ("newtonian", param%newtonian, parfile, ierr)

    call get_integer_parameter ("ip_coeffs", param%ip_coeffs, parfile, ierr)

    return

  contains ! --- the next subroutines are not visible from outside 

    !     ===================================================================================
    subroutine get_string_parameter (name_string, value, parfile, ierr)

      implicit none


      character*(*) name_string, value
      character (len=255) :: parfile
      integer, intent (out) :: ierr

      character*(100) line_string
      character*(50) name
      integer i, j, l, ll

      ierr = 0


      i = index (name_string, '=')
!!!   index of equal sign
      if(i.gt.0)then
!!!   name of parameter
         name = name_string (1 : i - 1)
      else
         name = name_string
      endif

!!!   store unoccupied unit handle in unit
      open (unit = 9, file = parfile, status = 'unknown')
!!!   open paramter file
10    format (a)
!!!   character string format

15    continue
      read (9, 10, end = 19, err = 15) line_string
!!!   read line with format 10, goto 19 at end, goto 15 (redo) at error, store in line
      i = index (line_string, '=')
!!!   index of equal sign
      j = index (line_string, '#')
!!!   index of comment char

      if (i.eq.0.or.j.eq.1) goto 15
!!!   if the whole line is a comment or there is no
!!!   equal sign, then go on to the next line    

      if(j.gt.0.and.j.lt.i) goto 15
!!!   if there is an equal sign, but it is in a comment
!!!   then go on to the next line

      l = len (line_string)
!!!   if there is a comment in the line, make sure
!!!   we do not include it as part of the value

      if (j.gt.0) l = j - 1
!!!   strip out leading blanks

      ll = 0
17    continue
      ll = ll + 1
      if (line_string (ll : ll).eq.' ') goto 17
      if (name.eq.line_string (ll : i - 1)) then
         value = line_string (i + 1 : l)
         goto 20
      endif
      goto 15

19    continue

      write (*, *) "ERROR: could not read parameter", name
      ierr = 125

20    continue

      close (9)

      return
    end subroutine get_string_parameter



    !     ===================================================================================
    subroutine get_real_parameter (name, value, parfile, ierr)

      implicit none

      character (len=255) :: parfile
      integer, intent (out) :: ierr

      character*(*) name
      character*(50) name_tmp
      character*(100) value_string
      real*8 value
      integer i

      ierr = 0

      call get_string_parameter (name, value_string, parfile, ierr)

!!!   make sure that there is a decimal point in the string
      if (index (value_string, '.').eq.0) then
         i = index (name, '=')
         if (i.eq.0) i = len (name)
         name_tmp = name (1 : i - 1)
         print*
         print*, 'Fatal problem in:'
         print*, '   Parameter reader'
         print*
         print*, 'Diagnosis:'
         print*, '   Could not find decimal point in real parameter ' // &
              name_tmp
         print*
         print*, 'Action:'
         print*, '   Exiting'
         print*
         stop
      endif

!!!   use internal file to convert string to real
      read (value_string,12) value
12    format (e20.15)

      return
    end subroutine get_real_parameter



    !     ===================================================================================
    subroutine get_integer_parameter (name, value, parfile, ierr)

      implicit none

      character (len=255) :: parfile
      integer, intent (out) :: ierr

      character*(*) name
      character*(50) value_string
      integer value

      ierr = 0


      call get_string_parameter (name, value_string, parfile, ierr)

!!!   use internal file to convert string to integer
      read (value_string, 13) value
13    format (i10)
      return
    end subroutine get_integer_parameter



    !     ===================================================================================
    subroutine get_logical_parameter (name, value, parfile, ierr)

      implicit none

      character (len=255) :: parfile
      integer, intent (out) :: ierr


      character*(*) name
      character*(50) value_string
      logical value

      ierr = 0


      call get_string_parameter (name, value_string, parfile, ierr)
      read (value_string, 14) value
14    format (l10)

      return
    end subroutine get_logical_parameter

  end subroutine read_parameters

end module module_param
