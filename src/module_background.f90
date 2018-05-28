!************************************************************************
! Module to manage the background data use to perform the mode analysis
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



module module_background

  use module_param

  implicit none

  private

  real*8, public, parameter :: c_light=2.99792458d10
  real*8, public, parameter :: pi = 3.1415926535897932384626433832d0 
  real*8, public, parameter :: m_solar=1.989d33
  real*8, public, parameter :: g_grav=6.673d-8
  real*8, public, parameter :: p_geom_factor = g_grav / c_light ** 4
  real*8, public, parameter :: rho_geom_factor = g_grav / c_light ** 2

  logical, save :: bg_file_exists = .false.
  logical, save :: glo_file_exists = .false.

  !> This type holds the background data (1D profiles) for a given time
  type, public :: BGData_t

     ! -------------------------------------------------
     ! ----- variables needed for the computation ------

     !> Size of the grid
     integer :: m
     !> Number of ghost cells
     integer :: g
     !> Time of the profile
     real*8 :: time
     !> number of the file
     integer :: nt
     !> Shock location (index of the last point INSIDE the shock)
     integer :: iR
     !> directory with the data
     character (len=255) :: input_directory
     !> directory for the output
     character (len=255) :: output_directory

     !> radius
     real*8, dimension (:), allocatable :: r
     !> rest mass density
     real*8, dimension (:), allocatable :: rho
     !> specific internal energy
     real*8, dimension (:), allocatable :: eps
     !> pressure
     real*8, dimension (:), allocatable :: p
     !> speed of sound square
     real*8, dimension (:), allocatable :: c_sound_squared
     !> conformal factor
     real*8, dimension (:), allocatable :: phi
     !> lapse
     real*8, dimension (:), allocatable :: alpha


     ! -------------------------------------------------
     ! -------------------------------------------------

     ! ------ derivatives coefficients (computed in the module) ----
     real*8, dimension (:), allocatable :: coeff_dm
     real*8, dimension (:), allocatable :: coeff_d0
     real*8, dimension (:), allocatable :: coeff_dp

     ! ------- derivatives of variables (computed in the module) ---
     real*8, dimension (:), allocatable :: dphi_dr
     real*8, dimension (:), allocatable :: dalpha_dr
     real*8, dimension (:), allocatable :: de_dr
     real*8, dimension (:), allocatable :: dq_dr
     real*8, dimension (:), allocatable :: dp_dr
     real*8, dimension (:), allocatable :: drho_dr
     real*8, dimension (:), allocatable :: drhoh_dr
     real*8, dimension (:), allocatable :: dlnq_dr
     real*8, dimension (:), allocatable :: dlnrhoh_dr

     ! ------- extra variables (computed in the module) ------------
     !> relativistic enthalpy (h = 1 + eps + P / rho)
     real*8, dimension (:), allocatable :: h
     !> adiabatic index
     real*8, dimension (:), allocatable :: gamma_one
     !> G = - grad log alpha or G = grap P / rho h
     real*8, dimension (:), allocatable :: ggrav
     !> Brunt Vaisala frequency square
     real*8, dimension (:), allocatable :: n2
     !> Lamb frequency
     real*8, dimension (:), allocatable :: lamb2
     !> Buoyancy term 
     real*8, dimension (:), allocatable :: Bstar
     !> Q = alpha * phi
     real*8, dimension (:), allocatable :: Q
     !> 1 / cs^2
     real*8, dimension (:), allocatable :: inv_cs2
     !> integrated mass
     real*8, dimension (:), allocatable :: mass_v


     ! ----- extra variables (not computed) -------
     !> radial velocity
     real*8, dimension (:), allocatable :: v_1
     !> electron number fraction
     real*8, dimension (:), allocatable :: y_e
     !> temperature ????
     real*8, dimension (:), allocatable :: t
     !> entropy per baryon
     real*8, dimension (:), allocatable :: entropy
     !> Newtonian potential
     real*8, dimension (:), allocatable :: U


     !> flag to compute n2
     logical :: calculate_n2 = .true.


     !> total mass
     real*8 :: mass


  end type BGData_t
  
  public :: Init_background_data
  public :: Free_background_data
  public :: Output_backgroud_data
  public :: Impose_boundary_conditions_bg
  public :: Compute_shock_location_bg
  public :: Check_for_errors_bg
  public :: Compute_integrated_quantities_bg

contains
  
  !****************************************************************************
  ! Allocate and initialize 1D data
  !****************************************************************************  
  subroutine Init_background_data (data, m, ierr)
    implicit none

    ! -----  Arguments  --------
    type(BGData_t), intent (inout) :: data
    integer, intent (in) :: m
    integer, intent (out) :: ierr
    ! --------------------------

    integer :: g

    ierr = 0

    write (*, *) "Allocating memory: grid", m

    g = 1
    data%m = m
    data%g = g

    allocate (data%r (1-g : m+g))
    allocate (data%rho (1-g : m+g))
    allocate (data%eps (1-g : m+g))
    allocate (data%p (1-g : m+g))
    allocate (data%c_sound_squared (1-g : m+g))
    allocate (data%phi  (1-g : m+g))
    allocate (data%alpha  (1-g : m+g))

    allocate (data%coeff_dm(1:m))
    allocate (data%coeff_d0(1:m))
    allocate (data%coeff_dp(1:m))

    allocate (data%dphi_dr(1-g:m+g))
    allocate (data%dalpha_dr(1-g:m+g))
    allocate (data%de_dr(1-g:m+g))
    allocate (data%dq_dr(1-g:m+g))
    allocate (data%dp_dr(1-g:m+g))
    allocate (data%drho_dr(1-g:m+g))
    allocate (data%drhoh_dr(1-g:m+g))
    allocate (data%dlnq_dr(1-g:m+g))
    allocate (data%dlnrhoh_dr(1-g:m+g))

    allocate (data%h (1-g : m+g))
    allocate (data%gamma_one  (1-g : m+g))
    allocate (data%ggrav (1-g : m+g))
    allocate (data%n2(1-g : m+g))
    allocate (data%lamb2(1-g : m+g))
    allocate (data%Bstar(1-g : m+g))
    allocate (data%Q(1-g : m+g))
    allocate (data%inv_cs2(1-g : m+g))
    allocate (data%mass_v(1-g : m+g))

    allocate (data%v_1 (1-g : m+g))
    allocate (data%y_e(1-g : m+g))
    allocate (data%t  (1-g : m+g))
    allocate (data%entropy(1-g : m+g))
    allocate (data%U(1-g : m+g))

    data%r = 0.0d0
    data%rho = 0.0d0
    data%eps = 0.0d0
    data%p = 0.0d0
    data%c_sound_squared = 0.0d0
    data%phi = 1.0d0
    data%alpha = 1.0d0

    data%coeff_dm = 0.0d0
    data%coeff_d0 = 0.0d0
    data%coeff_dp = 0.0d0

    data%dphi_dr = 0.0d0
    data%dalpha_dr = 0.0d0
    data%de_dr = 0.0d0
    data%dq_dr = 0.0d0
    data%dp_dr = 0.0d0
    data%drho_dr = 0.0d0
    data%drhoh_dr = 0.0d0
    data%dlnq_dr = 0.0d0
    data%dlnrhoh_dr = 0.0d0

    data%h = 1.0d0
    data%gamma_one = 0.0d0
    data%ggrav = 0.0d0
    data%n2 = 0.0d0
    data%lamb2 = 0.0d0
    data%Bstar = 0.0d0
    data%Q = 0.0d0
    data%inv_cs2 = 0.0d0
    data%mass_v = 0.0d0

    data%v_1=0.0d0
    data%y_e=0.0d0
    data%t = 0.0d0
    data%entropy = 0.0d0
    data%U = 0.0d0

  end subroutine Init_background_data

  !****************************************************************************
  ! Deallocate 1D data
  !****************************************************************************  
  subroutine Free_background_data (data, ierr)
    implicit none
    ! -----  Arguments  --------
    type(BGData_t), intent (inout) :: data
    integer, intent (out) :: ierr
    ! --------------------------

    ierr = 0

    deallocate (data%r)
    deallocate (data%rho)
    deallocate (data%eps)
    deallocate (data%p)
    deallocate (data%c_sound_squared)
    deallocate (data%phi)
    deallocate (data%alpha)

    deallocate (data%coeff_dm)
    deallocate (data%coeff_d0)
    deallocate (data%coeff_dp)

    deallocate (data%dphi_dr)
    deallocate (data%dalpha_dr)
    deallocate (data%de_dr)
    deallocate (data%dq_dr)
    deallocate (data%dp_dr)
    deallocate (data%drho_dr)
    deallocate (data%drhoh_dr)
    deallocate (data%dlnq_dr)
    deallocate (data%dlnrhoh_dr)

    deallocate (data%h)
    deallocate (data%gamma_one)
    deallocate (data%ggrav)
    deallocate (data%n2)
    deallocate (data%lamb2)
    deallocate (data%Bstar)
    deallocate (data%Q)
    deallocate (data%inv_cs2)
    deallocate (data%mass_v)

    deallocate (data%v_1)
    deallocate (data%y_e)
    deallocate (data%t)
    deallocate (data%entropy)
    deallocate (data%U)
 

  end subroutine Free_background_data



! ==============================================================================================================
! ==============================================================================================================
  subroutine Output_backgroud_data(data, param, output_format, ierr)
    
    
    implicit none

    ! ------ Arguments --------- 
    ! Stores 2D profile data at each time
    type(BGData_t), intent (inout) :: data
    type(parameters_t), intent (inout) :: param
    !> Format of the outupt file
    !>     1: csv file
    !>     2: ASCII for gnuplot
    integer, intent (in) :: output_format
    integer, intent (out) :: ierr
    ! --------------------------

    character (len=255) :: fileout
    integer i , m, nt
    real*8 :: delta_r

    ierr = 0

    ! local copy
    m = data%m
    nt = data%nt

    select case (output_format)
    case (1) ! -------------- CSV format ---------------

       write(fileout,"(a11,i0.4,a4)") "quantities_",nt,'.csv'
       fileout = trim(data%output_directory)//"/"//trim(fileout)
       write (*, *) "Output background data :", fileout   
       open (22, file = fileout,STATUS='UNKNOWN', ACTION='WRITE')

       write(22,*) "time,","r,","rho,","h,","p,","c_s2,","alpha,", &
            "lamb2,","n2,","phi,","ggrav",",","v_1", &
            ",","entropy",",","t",",","delta_r",",","l",",", &
            "total_mass",",","mass",",","dpi_phi",",","eps",",","dp_dr",&
            ",","dalpha_dr"

       do i = 1, data%iR
          delta_r = data%r (i+1) - data%r (i)
          write(22,*) data%time,",",data%r(i),",",data%rho(i),",", &
               data%h(i),",",data%p(i),",",data%c_sound_squared(i),",", &
               data%alpha(i),",",data%lamb2(i),",",data%n2(i),",", &
               data%phi(i),",",data%ggrav(i),"," &
               ,data%v_1(i),",",data%entropy(i), &
               ",",data%t(i),",",delta_r,",",param%l,&
               ",",data%mass,",",data%mass_v(i),",",data%dphi_dr(i),",",data%eps(i),",",data%dp_dr(i),&
               ",",data%dalpha_dr(i)
       enddo
       close(22)

    case (2)

       if (param%restart) then
          bg_file_exists = .true.
          glo_file_exists = .true.
       endif


       fileout = trim(data%output_directory)//"/"//"background.dat"
       if (.not.bg_file_exists) then
          write (*, *) "Output background data :", fileout
          open (22, file = fileout,status='unknown', action='write')
          write (22, "(25A25)") &
                 "time [s]" &
               , "r [cm]" &
               , "rho [cm^-2]" &
               , "eps " &
               , "p [cm^-1]" &
               , "cs^2" &
               , "phi" &
               , "alpha" &
               , "h" &
               , "Gamma_1" &
               , "G [cm^-1]" &
               , "N^2 [cm^-2]" &
               , "L^2 [cm^-2]" &
               , "B [cm^-1]" &
               , "Q" &
               , "1/cs^2" &
               , "v_1" &
               , "Y_e" &
               , "T" &
               , "s" &
               , "U" &
               , "M [cm]" 

          bg_file_exists = .true.
       else
          open (22, file = fileout,status='old', position='append', action='write')
       endif

       do i = 1 - data%g, data%m + data%g 
          write (22, "(25E25.15)") &
               data%time &
               , data%r (i) &
               , data%rho (i) &
               , data%eps (i) &
               , data%p (i) &
               , data%c_sound_squared (i) &
               , data%phi (i) &
               , data%alpha (i) &
               , data%h (i) &
               , data%gamma_one (i) &
               , data%ggrav (i) &
               , data%n2 (i) &
               , data%lamb2 (i) &
               , data%Bstar (i) &
               , data%Q (i) &
               , data%inv_cs2 (i) &
               , data%v_1 (i)&
               , data%y_e (i)&
               , data%t (i)&
               , data%entropy (i)&
               , data%U (i)&
               , data%mass_v (i)

       enddo
       write (22, *) " "
       close (22)

       fileout = trim(data%output_directory)//"/"//"global.dat"
       if (.not.glo_file_exists) then
          write (*, *) "Output global quantities :", fileout
          open (22, file =fileout,status='unknown', action='write')
          write (22, "(A25,3A10,25A25)") &
               "time [s]", &
               "nt", &
               "m", &
               "iR", &
               "r [cm]", &
               "Mass [cm]"
          
          glo_file_exists = .true.
       else
          open (22, file = fileout,status='old', position='append', action='write')
       endif
       write (22, "(E25.15,3I10,25E25.15)") &
            data%time &
            , data%nt &
            , data%m &
            , data%iR &
            , data%r(data%iR) &
            , data%mass 
       close (22)


    case default
       write (*, *) "ERROR in Output_background_data"
       write (*, *) "invalid output_format:", output_format
       ierr = 333
    end select

  end subroutine Output_backgroud_data



  !****************************************************************************
  ! Impose boundary conditions
  !****************************************************************************  
  subroutine Impose_boundary_conditions_bg (data, ierr)
    implicit none
    ! -----  Arguments  --------
    type(BGData_t), intent (inout) :: data
    integer, intent (out) :: ierr
    ! --------------------------

    integer :: i, m, g

    ierr = 0

    m = data%m
    g = data%g

    !   ---- Compute inner ghost cells for mean quantities ----
    do i = 1, g
       data%r (1-i) = -data%r(i)
       data%rho (1-i) = data%rho (i)
       data%eps (1-i) = data%eps (i)
       data%p (1-i) = data%p (i)
       data%c_sound_squared (1-i) = data%c_sound_squared (i)
       data%phi (1-i) = data%phi (i)
       data%alpha (1-i) = data%alpha (i)

       data%h (1-i) = data%h (i)
       data%gamma_one (1-i) = data%gamma_one (i)
       data%ggrav (1-i) = - data%ggrav (i)
       data%n2 (1-i) = - data%n2 (i)
       data%lamb2 (1-i) = - data%lamb2 (i)
       data%Bstar (1-i) = - data%Bstar (i)
       data%Q (1-i) = data%Q (i)
       data%inv_cs2 (1-i) = data%inv_cs2 (i)
       data%mass_v (1-i) = - data%mass_v (i)

       data%v_1 (1-i) = - data%v_1 (i)
       data%y_e (1-i) = data%y_e (i)
       data%t (1-i) = data%t (i)
       data%entropy (1-i) = data%entropy (i)
       data%U (1-i) = data%U (i)
    enddo

    do i = 1, g
       data%r (m+i) = data%r (m) + i*(data%r (m) - data%r (m-1))
       data%rho (m+i) = data%rho (m)
       data%eps (m+i) = data%eps (m)
       data%p (m+i) = data%p (m)
       data%c_sound_squared (m+i) = data%c_sound_squared (m)
       data%phi (m+i) = data%phi (m)
       data%alpha (m+i) = data%alpha (m)

       data%h (m+i) = data%h (m)
       data%gamma_one (m+i) = data%gamma_one (m)
       data%ggrav (m+i) = data%ggrav (m)
       data%n2 (m+i) = data%n2 (m)
       data%lamb2 (m+i) = data%lamb2 (m)
       data%Bstar (m+i) = data%Bstar (m)
       data%Q (m+i) = data%Q (m)
       data%inv_cs2 (m+i) = data%inv_cs2 (m)
       data%mass_v (m+i) = data%mass_v (m)

       data%v_1 (m+i) = data%v_1 (m)
       data%y_e (m+i) = data%y_e (m)
       data%t (m+i) = data%t (m)
       data%entropy (m+i) = data%entropy (m)
       data%U (m+i) = data%U (m)
    enddo


  end subroutine Impose_boundary_conditions_bg


  !****************************************************************************
  ! Compute shock location
  !****************************************************************************  
  subroutine Compute_shock_location_bg (data, ierr)
    implicit none
    ! -----  Arguments  --------
    type(BGData_t), intent (inout) :: data
    integer, intent (out) :: ierr
    ! --------------------------

    integer :: i

    ierr = 0

    ! Compute shock location (iR is the last point inside the shock)
    data%iR = 0 ! if not detected iR=0
    do i = 1, data%m
       if (data%c_sound_squared(i).lt.data%v_1(i)**2) then
          data%iR = i-1
          exit
       endif
    enddo
    
  end subroutine Compute_shock_location_bg


  !****************************************************************************
  ! Check for errors in data
  !****************************************************************************  
  subroutine Check_for_errors_bg (data, ierr)
    implicit none
    ! -----  Arguments  --------
    type(BGData_t), intent (inout) :: data
    integer, intent (out)  :: ierr
    ! --------------------------

    ierr = 0

    if (any(isnan(data%r))) ierr = 11
    if (any(isnan(data%rho))) ierr = 12
    if (any(isnan(data%eps))) ierr = 13
    if (any(isnan(data%p))) ierr = 14
    if (any(isnan(data%c_sound_squared))) ierr = 15
    if (any(isnan(data%phi))) ierr = 16
    if (any(isnan(data%alpha))) ierr = 17

    if (any(isnan(data%coeff_dm))) ierr = 21
    if (any(isnan(data%coeff_d0))) ierr = 22
    if (any(isnan(data%coeff_dp))) ierr = 23

    if (any(isnan(data%dphi_dr))) ierr = 31
    if (any(isnan(data%dalpha_dr))) ierr = 32
    if (any(isnan(data%de_dr))) ierr = 33
    if (any(isnan(data%dq_dr))) ierr = 34
    if (any(isnan(data%dp_dr))) ierr = 35
    if (any(isnan(data%drho_dr))) ierr = 36
    if (any(isnan(data%drhoh_dr))) ierr = 37
    if (any(isnan(data%dlnq_dr))) ierr = 38
    if (any(isnan(data%dlnrhoh_dr))) ierr = 39

    if (any(isnan(data%h))) ierr = 41
    if (any(isnan(data%gamma_one))) ierr = 42
    if (any(isnan(data%ggrav))) ierr = 43
    if (any(isnan(data%n2))) ierr = 44
    if (any(isnan(data%lamb2))) ierr = 45
    if (any(isnan(data%Bstar))) ierr = 46
    if (any(isnan(data%Q))) ierr = 47
    if (any(isnan(data%inv_cs2))) ierr = 48
    if (any(isnan(data%mass_v))) ierr = 49

    if (any(isnan(data%v_1))) ierr = 51
    if (any(isnan(data%y_e))) ierr = 52
    if (any(isnan(data%t))) ierr = 53
    if (any(isnan(data%entropy))) ierr = 54
    if (any(isnan(data%U))) ierr = 55


    return

  end subroutine Check_for_errors_bg

  
  !****************************************************************************
  ! Compute integrated quantities
  !****************************************************************************  
  subroutine Compute_integrated_quantities_bg (data, ierr)
    implicit none
    ! -----  Arguments  --------
    type(BGData_t), intent (inout) :: data
    integer, intent (out) :: ierr
    ! --------------------------

    integer :: i
    real*8 :: r_if_p, r_if_m
    real*8 :: dvol

    ierr = 0

    data%mass = 0.0d0
    do i = 1, data%iR
       
       r_if_p = 0.5d0 * (data%r (i+1) + data%r (i))
       r_if_m = 0.5d0 * (data%r (i) + data%r (i-1))

       dvol = 4.0d0 / 3.0d0 * pi * ( r_if_p ** 3 - r_if_m ** 3 )
       ! We neglect the Lorentz factor since we are assuming
       ! small velocities anyway
       data%mass = data%mass + data%rho (i) * data%phi (i) ** 6 * dvol
       data%mass_v (i) = data%mass
       
    enddo

  end subroutine Compute_integrated_quantities_bg
  
end module module_background
