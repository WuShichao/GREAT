!************************************************************************
! Module to integrate the eigenmode equations
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

module module_mode_integrator

  use module_param
  use module_background
  use module_eigen
  use module_fast_inv

  implicit none
  private

  ! Coefficients matrix for the integrator
  type, public :: mode_coeff_t
     !> size of the grid
     integer :: iR
     !> Size of the coefficients matrix
     integer :: n
     !> Coefficients matrix
     real*8, dimension (:,:,:), allocatable :: matrix
     !> Solution vector
     real*8, allocatable, dimension (:,:) :: u
     !> r
     real*8, allocatable, dimension (:) :: r     
     !> needed for the outer bc
     real*8, allocatable, dimension (:,:) :: lim
     integer, allocatable, dimension (:,:) :: rpower
     logical, allocatable, dimension (:,:) :: modify
     

  end type mode_coeff_t

  public :: Initialize_coefficients_matrix
  public :: Free_coefficients_matrix
  public :: Compute_coefficients_full
  public :: Set_inner_boundary_conditions
  public :: Compute_outer_boundary_conditions
  public :: Copy_solution
  public :: Integrate_equation
  public :: Invert_matrix
!  public :: Compute_coefficients_Cowling

contains

  !=========================================================================
  ! Initializes coefficients matrix
  !=========================================================================

  subroutine Initialize_coefficients_matrix (data, param, coeff, ierr)

    ! -----  Arguments  --------
    ! input
    type(BGData_t), intent (inout) :: data
    type (parameters_t), intent (inout) ::param    
    ! output
    type (mode_coeff_t), intent (inout) :: coeff
    integer, intent (out) :: ierr
    ! --------------------------

    ierr = 0

    coeff%iR = data%iR
    if (param%cowling) then ! Cowling: 2x2 matrix
       coeff%n = 2
    elseif(.not.param%cal_dphi) then ! Non-cowling but delta psi = 0 : 4x4 matrix stored in 6x6 matrix
       coeff%n = 4
    else ! Non-Cowling (general): 6x6 matrix
       coeff%n = 6
    endif
    allocate (coeff%matrix (1:coeff%n, 1:coeff%n, 1:coeff%iR)) 
    allocate (coeff%u (1:coeff%n, 1:coeff%iR)) 
    allocate (coeff%r (1:coeff%iR)) 
    allocate (coeff%lim (1:coeff%n, 1:coeff%n))
    allocate (coeff%rpower (1:coeff%n, 1:coeff%n))
    allocate (coeff%modify (1:coeff%n, 1:coeff%n))

    coeff%matrix = 0.0d0
    coeff%u = 0.0d0
    coeff%r (1:coeff%iR) = data%r (1:data%iR)
    coeff%lim = 0.0d0
    coeff%rpower = 0
    coeff%modify = .false.

  end subroutine Initialize_coefficients_matrix

  !=========================================================================
  ! Free coefficients matrix
  !=========================================================================

  subroutine Free_coefficients_matrix (coeff, ierr)

    ! -----  Arguments  --------
    type (mode_coeff_t), intent (inout) :: coeff
    ! output
    integer, intent (out) :: ierr
    ! --------------------------

    ierr = 0
    coeff%iR = 0
    coeff%n = 0
    deallocate (coeff%matrix) 
    deallocate (coeff%u) 
    deallocate (coeff%r)
    deallocate (coeff%lim)
    deallocate (coeff%rpower)
    deallocate (coeff%modify)

  end subroutine Free_coefficients_matrix


  !=========================================================================
  ! Check_for_errors
  !=========================================================================

  subroutine Check_for_coefficients_errors (coeff, ierr)

    ! -----  Arguments  --------
    ! input
    type (mode_coeff_t), intent (inout) :: coeff
    ! output
    integer, intent (out) :: ierr
    ! --------------------------

    integer :: n1, n2

    ierr = 0

    do n1 = 1, coeff%n
       do n2 = 1, coeff%n
          if (any(isnan(coeff%matrix(n1, n2, :))))  then
             Write (*, *) "NaN found in coefficients matrix:", n1, n2
          endif
       enddo
    enddo
    ierr = 200

  end subroutine Check_for_coefficients_errors



  ! ==============================================================================================================
  ! Computes the coefficient matrix corresponding to the non-Cowling approximation (6x6)
  ! ==============================================================================================================  
  subroutine Compute_coefficients_full (data, param, w2, coeff, ierr)
    
    implicit none
    
    ! ------- Arguments -------
    ! input
    type(BGData_t), intent (inout) :: data
    type(parameters_t), intent (inout) :: param
    real*8 :: w2
    ! output
    type (mode_coeff_t), intent (inout) :: coeff
    integer, intent (out) :: ierr
    ! ----------------------------

    integer :: i
    
    ! Auxiliary variables
    real*8 :: aa, bb, cc , dd

    ierr = 0
    
    if (param%l.ne.0) then !---- l != 0 ---------
      
       do i = 1, data%iR
          
          ! Temporary variables
          aa  = data%alpha(i) / data%phi(i)**5.0 / w2 * (data%dlnrhoh_dr(i) - (1.0d0 + data%inv_cs2(i))*data%ggrav(i))
          
          bb = 2.0d0 * pi * (data%rho(i) * data%h(i) + 5.0d0*data%p(i)) * data%phi(i)**4.0
          
          cc = 2.0d0 * pi * data%rho(i)  * data%h(i) * data%alpha(i) * data%phi(i)**5.0
          
          dd = 6.0d0 + data%inv_cs2(i)

          !_____ Equation for eta_r ______
          coeff%matrix (1,1,i) = -2.0d0/data%r(i) - data%ggrav(i)* data%inv_cs2(i) - 6.0 * data%dphi_dr(i)/data%phi(i)
          !coeff%matrix (1,1,i) = -1.0d0 * (2.0d0/r(i) + 1.0/gamma_one(i)*data%dp_dr(i)/p(i) + 6.0 * data%dphi_dr(i)/phi(i))
          coeff%matrix (1, 2, i) = &
               param%l*(param%l+1.0)/data%r(i)**2 &
               - data%phi(i)**4 * data%inv_cs2(i) / data%alpha(i)**2 * w2
          if (.not.param%cowling) then
             coeff%matrix (1,3,i) = 0.0d0
             coeff%matrix (1,4,i) = data%inv_cs2(i)/data%Q(i)
             if (param%cal_dphi) then
                coeff%matrix (1,5,i) = 0.0d0
                coeff%matrix (1,6,i) = -dd/data%phi(i)
             endif
          endif
          !_____ Equation for eta_p ______
          coeff%matrix (2,1,i) = 1.0d0 - data%n2(i)/w2
          coeff%matrix (2,2,i) = - data%dlnq_dr(i)  + data%ggrav(i) * (1.0d0+data%inv_cs2(i))
          if (.not.param%cowling) then
             coeff%matrix (2,3,i) = 0.0d0
             coeff%matrix (2,4,i) = aa
             if (param%cal_dphi) then
                coeff%matrix (2,5,i) = 0.0d0
                coeff%matrix (2,6,i) = -aa*data%alpha(i)
             endif
          endif

          if (.not.param%cowling) then

             !______ Equation for K ______
             coeff%matrix (3,1,i) = -cc * data%Bstar(i)
             coeff%matrix (3,2,i) = cc * dd * data%phi(i)**4 * w2 / data%alpha(i)**2 
             coeff%matrix (3,3,i) = -2.0/data%r(i)
             coeff%matrix (3,4,i) = param%l*(param%l+1)/data%r(i)**2 &
                  - 2.0d0*pi*data%phi(i)**4 * ( 5.0d0*(data%rho(i)*data%h(i)-data%p(i))&
                  +data%rho(i)*data%h(i)*data%inv_cs2(i) )
             if (param%cal_dphi) then
                coeff%matrix (3,5,i) = 0.0d0
                coeff%matrix (3,6,i) = 2.0d0*pi*data%phi(i)**4*data%alpha(i) &
                     * (10.0d0*data%rho(i)*data%h(i) + 20.0d0*data%p(i) &
                     + data%rho(i)*data%h(i)*data%inv_cs2(i))
             endif

             !______ Equation for delta Q ______
             coeff%matrix (4,1,i) = 0.0d0
             coeff%matrix (4,2,i) = 0.0d0
             coeff%matrix (4,3,i) = 1.0d0
             coeff%matrix (4,4,i) = 0.0d0
             if (param%cal_dphi) then
                coeff%matrix (4,5,i) = 0.0d0
                coeff%matrix (4,6,i) = 0.0d0
             endif

             if (param%cal_dphi) then

                !______ Equation for Psi _________
                coeff%matrix (5,1,i) = 2.0d0*pi*data%rho(i)*data%h(i)*data%phi(i)**5 *data%Bstar(i)
                coeff%matrix (5,2,i) = -2.0d0*pi*data%rho(i)*data%h(i)*data%phi(i)**5 &
                     * data%phi(i)**4*w2*data%inv_cs2(i)/data%alpha(i)**2.0
                coeff%matrix (5,3,i) = 0.0d0
                coeff%matrix (5,4,i) = 2.0d0*pi*data%phi(i)**4.0 * data%rho(i)*data%h(i) &
                     * data%inv_cs2(i)/data%alpha(i)
                coeff%matrix (5,5,i) = -2.0d0/data%r(i)
                coeff%matrix (5,6,i) = param%l*(param%l+1)/data%r(i)**2 &
                     - 2.0d0*pi * data%phi(i)**4.0 * ( 5.0d0*data%rho(i)*(1.0d0+data%eps(i)) &
                     + data%rho(i)*data%h(i)*data%inv_cs2(i) )

                !_______ Equation for delta psi __________
                coeff%matrix (6,1,i) = 0.0d0
                coeff%matrix (6,2,i) = 0.0d0
                coeff%matrix (6,3,i) = 0.0d0
                coeff%matrix (6,4,i) = 0.0d0
                coeff%matrix (6,5,i) = 1.0d0
                coeff%matrix (6,6,i) = 0.0d0

             endif

          endif
          
       enddo

    else ! ---- l = 0 ---------

       do i = 1, data%iR
          !__________ Equation for eta_r ____________
          coeff%matrix (1,1,i) = - 2.0d0/data%r(i) - data%ggrav(i)* data%inv_cs2(i) &
               - 6.0 * data%dphi_dr(i)/ data%phi(i)
          coeff%matrix (1,2,i) = - data%inv_cs2(i) /(data%rho(i)*data%h(i))
          if (.not.param%cowling) then
             coeff%matrix (1,3,i) = 0.0d0
             coeff%matrix (1,4,i) = 0.0d0
             if (param%cal_dphi) then
                coeff%matrix (1,5,i) = 0.0d0
                coeff%matrix (1,6,i) = -6.0d0/data%phi(i)
             endif
          endif
          !__________ Equation for delta P ___________
          coeff%matrix (2,1,i) = -data%rho(i)*data%h(i)/data%alpha(i)**2.0*data%phi(i)**4.0 &
               *(data%n2(i) - w2)
          coeff%matrix (2,2,i) = data%ggrav(i) * (1.0d0 + data%inv_cs2(i))
          if (.not.param%cowling) then
             coeff%matrix (2,3,i) = -data%rho(i)*data%h(i)/(data%alpha(i)*data%phi(i))
             coeff%matrix (2,4,i) = -data%rho(i)*data%h(i)/(data%alpha(i)*data%phi(i)) &
                  * (data%ggrav(i) -  data%dphi_dr(i)/data%phi(i))
             if (param%cal_dphi) then
                coeff%matrix (2,5,i) = data%rho(i)*data%h(i)/data%phi(i)
                coeff%matrix (2,6,i) = -data%rho(i)*data%h(i)/data%phi(i)*data%dphi_dr(i)/data%phi(i)
             endif
          endif

          if (.not.param%cowling) then
             !____________ Equation for K ___________
             coeff%matrix (3,1,i) = -2.0d0* pi * data%rho(i) * data%h(i) &
                  * data%alpha(i) * data%phi(i)**5.0 * data%Bstar(i)
             coeff%matrix (3,2,i) = 2.0d0* pi * data%alpha(i) * data%phi(i)**5.0 &
                  * (6.0d0 + 2.0d0*data%inv_cs2(i))
             coeff%matrix (3,3,i) = -2.0d0/data%r(i)
             coeff%matrix (3,4,i) = 2.0d0* pi * data%phi(i)**4.0 &
                  * (data%rho(i) * data%h(i) + 5.0d0*data%p(i))
             if (param%cal_dphi) then
                coeff%matrix (3,5,i) = 0.0d0
                coeff%matrix (3,6,i) = 8.0d0* pi *data%alpha(i)* data%phi(i)**4.0 &
                     * (data%rho(i) * data%h(i) + 5.0d0*data%p(i))
             endif

             !__________ Equation for delta Q _________
             coeff%matrix (4,1,i) = 0.0d0
             coeff%matrix (4,2,i) = 0.0d0
             coeff%matrix (4,3,i) = 1.0d0
             coeff%matrix (4,4,i) = 0.0d0
             if (param%cal_dphi) then
                coeff%matrix (4,5,i) = 0.0d0
                coeff%matrix (4,6,i) = 0.0d0
             endif

             if (param%cal_dphi) then

                !____________ Equation for Psi ___________
                coeff%matrix (5,1,i) = 2.0d0 * pi * data%rho(i)  *data%h(i) * data%phi(i)**5.0 * data%Bstar(i)
                coeff%matrix (5,2,i) = - 2.0d0 * pi * data%phi(i)**5.0 * data%inv_cs2(i)
                coeff%matrix (5,3,i) = 0.0d0
                coeff%matrix (5,4,i) = 0.0d0
                coeff%matrix (5,5,i) = -2.0d0 / data%r(i)
                coeff%matrix (5,6,i) = -10.0d0 * pi * data%phi(i)**4.0 * (data%rho(i)*(1.0d0+data%eps(i)))   


                !__________ Equation for delta psi _________
                coeff%matrix (6,1,i) = 0.0d0
                coeff%matrix (6,2,i) = 0.0d0
                coeff%matrix (6,3,i) = 0.0d0
                coeff%matrix (6,4,i) = 0.0d0
                coeff%matrix (6,5,i) = 1.0d0
                coeff%matrix (6,6,i) = 0.0d0
             endif
          endif
       enddo
    endif

    call Interpolate_coefficients_full (data, param, coeff, ierr)
    
    call Check_for_coefficients_errors (coeff, ierr)


  end subroutine Compute_coefficients_full
  
 


  ! ==============================================================================================================
  ! Regularizes first points by interpolating to the center (non-Cowling)
  ! ==============================================================================================================
  subroutine Interpolate_coefficients_full (data, param, coeff, ierr)
    
    implicit none
    
    ! ------- Arguments -------    
    ! input 
    type(BGData_t), intent (inout) :: data
    type(parameters_t), intent (inout) :: param
    ! output
    type (mode_coeff_t), intent (inout) :: coeff
    integer, intent (out) :: ierr
    ! --------------------------

    integer :: i, n1, n2, ip
    real*8 :: tmp


    ierr = 0

    if (param%ip_coeffs.eq.0) return

    ip = param%ip_coeffs+1

    coeff%modify = .false.

    ! previous version: 11, 12, 21, 22, 24

    !   Define coeffient limits r=0
    if (param%l.ne.0) then !---- l != 0 ---------
   
       coeff%lim (1,1) = -2.0d0
       coeff%rpower (1,1) = 1
       coeff%modify (1, 1) = .true.

       coeff%lim (1,2) = param%l*(param%l+1.0d0)
       coeff%rpower (1,2) = 2
       coeff%modify (1, 2)= .true.

       coeff%lim (2,1) = 1.0d0
       coeff%rpower (2,1) = 0
       coeff%modify (2, 1)= .true.

       coeff%lim (2,2) = 0.0d0
       coeff%rpower (2,2) = 0
       coeff%modify (2, 2)= .true.

       if (.not.param%cowling) then

          coeff%lim (2,4) = 0.0d0
          coeff%rpower (2,4) = 0
          coeff%modify (2, 4) = .true.


          coeff%lim (3,1) = 0.0d0
          coeff%rpower (3,1) = 0
          coeff%modify (3, 1) = .true.
          !    coeff%lim (3,4) = l*(l+1.0d0)
          !    coeff%rpower (3,4) = 2
          !    coeff%modify (3, 4) = .true.

          if (param%cal_dphi) then

             coeff%lim (2,6) = 0.0d0
             coeff%rpower (2,6) = 0
             coeff%modify (2, 6) = .true.

             coeff%lim (5,1) = 0.0d0
             coeff%rpower (5,1) = 0
             coeff%modify (5, 1) = .true.

             !    coeff%lim (5,6) = l*(l+1.0d0)
             !    coeff%rpower (5, 6) = 2
             !    coeff%modify (5, 1) = .true.
          endif

       endif

    else !--------------- l = 0 ----------

       coeff%lim (1,1) = -2.0d0
       coeff%rpower (1,1) = 1
       coeff%modify (1, 1) = .true.

       coeff%lim (2,2) = 0.0d0
       coeff%rpower (2,2) = 0
       coeff%modify (2, 2) = .true.

       if (.not.param%cowling) then

          coeff%lim (2,4) = 0.0d0
          coeff%rpower (2,4) = 0
          coeff%modify (2, 4) = .true.

          coeff%lim (3,1) = 0.0d0
          coeff%rpower (3,1) = 0
          coeff%modify (3, 1) = .true.

          if (param%cal_dphi) then
             coeff%lim (2,6) = 0.0d0
             coeff%rpower (2,6) = 0
             coeff%modify (2, 6) = .true.

             coeff%lim (5,1) = 0.0d0
             coeff%rpower (5,1) = 0
             coeff%modify (5, 1) = .true.
          endif

       endif

    endif


    do n1 = 1, coeff%n
       do n2 = 1, coeff%n
          if (coeff%modify (n1, n2)) then
             do i=1,ip-1 
                tmp = coeff%lim (n1, n2) &
                     + ( coeff%matrix (n1, n2, ip) * data%r(ip) ** coeff%rpower (n1, n2) - coeff%lim (n1, n2)) &
                     / data%r (ip) * data%r (i)
                coeff%matrix (n1, n2, i) = tmp / data%r (i) ** coeff%rpower (n1, n2)
             enddo
          endif
       enddo
    enddo

  end subroutine Interpolate_coefficients_full


! ==============================================================================================================
! Set inner boundary conditions
! ==============================================================================================================
  subroutine Set_inner_boundary_conditions (ufree, data, param, coeff, ierr)

    implicit none 

    ! ------- Arguments -------    
    ! input     
    real*8, dimension (1:6) :: ufree
    type (BGData_t), intent (inout) :: data    
    type (parameters_t), intent (inout) :: param    
    type (mode_coeff_t), intent (inout) :: coeff
    ! output
    integer, intent (out) :: ierr
    ! --------------------------

    ierr = 0

    coeff%u (:,1) = 0.0d0

    ! ---------------- Cowling --------------------------
    ! one free parameters: eta_r
    if (param%cowling) then 
       
       if (param%l.eq.0) then ! ----- l = 0 -------
          coeff%u (1, 1) = ufree (1) ! eta_r
          coeff%u (2, 1) = 0.0d0 ! delta_p
       else                   ! ----- l |= 0 ------
          coeff%u (1, 1) = ufree (1) ! eta_r
          coeff%u (2, 1) = coeff%u (1, 1) * coeff%r(1) / param%l  ! eta_p
       endif

    !---------  Non-cowling and delta psi = 0 -----------
    ! two free parameters: eta_r and delta Q
    elseif(.not.param%cal_dphi) then 

       if (param%l.eq.0) then ! ----- l = 0 -------
          
          coeff%u (1, 1) = ufree (1) ! eta_r
          coeff%u (4, 1) = ufree (4) ! delta_Q

          coeff%u (2, 1) =  & ! delta_p
               -data%p (1) * data%gamma_one (1) &
               * coeff%u (1, 1) / data%r (1) 
          coeff%u (3, 1) = 0.0d0 ! K

       else                   ! ----- l |= 0 ------

          coeff%u (1, 1) = ufree (1) ! eta_r
          coeff%u (4, 1) = ufree (4) ! delta_Q

          coeff%u (2, 1) = coeff%u (1, 1) * coeff%r(1) / param%l  ! eta_p
          coeff%u (3, 1) =  coeff%u (4, 1) * param%l / coeff%r(1)! K
       endif



    ! -----------  Non-Cowling (general) ----------------
    ! three free parameters: eta_r, delta Q and delta psi
    else 

       if (param%l.eq.0) then ! ----- l = 0 -------

          coeff%u (1, 1) = ufree (1) ! eta_r
          coeff%u (4, 1) = ufree (4) ! delta_Q
          coeff%u (6, 1) = ufree (6) ! delta psi

          coeff%u (2, 1) =  & ! delta_p
               -data%p (1) * data%gamma_one (1) &
               * (6.0d0 * coeff%u(6, 1)/data%phi(1) + coeff%u (1, 1) / data%r (1) )
          coeff%u (3, 1) = 0.0d0 ! K
          coeff%u (5, 1) = 0.0d0 ! Psi

       else                   ! ----- l |= 0 ------

          coeff%u (1, 1) = ufree (1) ! eta_r
          coeff%u (4, 1) = ufree (4) ! delta_Q
          coeff%u (6, 1) = ufree (6) ! delta psi

          coeff%u (2, 1) = coeff%u (1, 1) * coeff%r(1) / param%l  ! eta_p
          coeff%u (3, 1) =  coeff%u (4, 1) * param%l / coeff%r(1)! K
          coeff%u (5, 1) = coeff%u (6, 1) * param%l / coeff%r(1) ! Psi
       endif


    endif


  end subroutine Set_inner_boundary_conditions


! ==============================================================================================================
! Compute outer boundary conditions
! ==============================================================================================================
  subroutine Compute_outer_boundary_conditions (func, param, coeff, ierr)
    
    implicit none 
    
    ! ------- Arguments -------    
    ! input     
    real*8, dimension (1:6), intent (inout) :: func
    type (parameters_t), intent (inout) :: param    
    type (mode_coeff_t), intent (inout) :: coeff
    ! output
    integer, intent (out) :: ierr
    ! --------------------------

    integer :: iR

    ierr = 0
    
    iR = coeff%iR

    func = 0.0d0
    
    ! ---------------- Cowling --------------------------
    if (param%cowling) then 

       ! eta_r = 0
       func (1) = coeff%u (1,iR)
       
    !---------  Non-cowling and delta psi = 0 -----------
    ! two free parameters: eta_r and delta Q
    elseif(.not.param%cal_dphi) then 
       
       func (1) = coeff%u (1, iR)
       func (2) = coeff%u (3, iR) + coeff%u (4, iR) * (param%l+1)/coeff%r(iR)
       
    ! -----------  Non-Cowling (general) ----------------
    ! three free parameters: eta_r, delta Q and delta psi
    else 
       
       func (1) = coeff%u (1,iR)
       func (2) = coeff%u (3, iR) + coeff%u (4, iR) * (param%l+1)/coeff%r(iR)
       func (3) = coeff%u (5, iR) + coeff%u (6, iR) * (param%l+1)/coeff%r(iR)

    endif


  end subroutine Compute_outer_boundary_conditions

! ==============================================================================================================
! Copy solution vector to variables
! ==============================================================================================================
  subroutine Copy_solution (data, param, coeff, eigen, ierr)

    implicit none 

    ! ------- Arguments -------    
    ! input     
    type(BGData_t), intent (inout) :: data
    type (parameters_t), intent (inout) ::param    
    type (mode_coeff_t), intent (inout) :: coeff 
    ! output
    type(Eigen_t), intent (inout) :: eigen
    integer, intent (out) :: ierr
    ! --------------------------

    integer :: iR

    ierr = 0

    iR = coeff%iR

    ! ---------------- Cowling --------------------------
    if (param%cowling) then 
       
       if (param%l.eq.0) then 

          eigen%nr = coeff%u (1,:) 
          eigen%delta_p = coeff%u (2,:) 

          eigen%np = 0.0d0
          eigen%k_cap = 0.0d0
          eigen%delta_Q = 0.0d0
          eigen%phi_cap = 0.0d0
          eigen%delta_phi = 0.0d0

       else

          eigen%nr = coeff%u (1,:) 
          eigen%np = coeff%u (2,:) 

          eigen%delta_p (1:iR) = -data%p (1:iR) * data%gamma_one (1:iR) &
               * coeff%u (1, 1:iR) / data%r (1:iR)

          eigen%k_cap = 0.0d0
          eigen%delta_Q = 0.0d0
          eigen%phi_cap = 0.0d0
          eigen%delta_phi = 0.0d0


       endif

    !---------  Non-cowling and delta psi = 0 -----------
    elseif(.not.param%cal_dphi) then 


       if (param%l.eq.0) then 

          eigen%nr (:) = coeff%u (1,:) 
          eigen%delta_p (:) = coeff%u (2,:) 
          
          eigen%np = 0.0d0
          eigen%k_cap = coeff%u (3,:) 
          eigen%delta_Q (:) = coeff%u (4,:) 
          eigen%phi_cap = 0.0d0
          eigen%delta_phi (:) = 0.0d0
          
       else
          
          eigen%nr (:) = coeff%u (1,:) 
          eigen%np (:) = coeff%u (2,:) 
          
          eigen%delta_p (1:iR) = -data%p (1:iR) * data%gamma_one (1:iR) &
               * coeff%u (1, 1:iR) / data%r (1:iR)
          eigen%k_cap = coeff%u (3,:) 
          eigen%delta_Q (:) = coeff%u (4,:) 
          eigen%phi_cap = 0.0d0
          eigen%delta_phi (:) = 0.0d0

       endif


    ! -----------  Non-Cowling (general) ----------------
    else 


       if (param%l.eq.0) then 

          eigen%nr (:) = coeff%u (1,:) 
          eigen%delta_p (:) = coeff%u (2,:) 
          
          eigen%np = 0.0d0
          eigen%k_cap = coeff%u (3,:) 
          eigen%delta_Q (:) = coeff%u (4,:) 
          eigen%phi_cap = coeff%u (5,:) 
          eigen%delta_phi (:) = coeff%u (6,:)

          
       else

          eigen%nr (:) = coeff%u (1,:) 
          eigen%np (:) = coeff%u (2,:) 
          
          eigen%delta_p (1:iR) = -data%p (1:iR) * data%gamma_one (1:iR) &
               * (6.0d0 * coeff%u(6, 1:iR)/data%phi(1:iR) + coeff%u (1, 1:iR) / data%r (1:iR) )
          eigen%k_cap = coeff%u (3,:) 
          eigen%delta_Q (:) = coeff%u (4,:) 
          eigen%phi_cap = coeff%u (5,:) 
          eigen%delta_phi (:) = coeff%u (6, :)          

       endif


    endif


  end subroutine Copy_solution

! ==============================================================================================================
! Integrates the equations given by its coefficients form
! ==============================================================================================================
  subroutine Integrate_equation (coeff, ierr)

    implicit none 

    ! ------- Arguments -------    
    ! input 
    type (mode_coeff_t), intent (inout) :: coeff
    ! output
    integer, intent (out) :: ierr
    ! --------------------------


    ! ------- Arguments -------

    integer :: i

    real*8, dimension(1:coeff%n,1:coeff%n) :: Mn,Nn,Mn_inv,aux
    real*8, dimension (1:coeff%n) :: Un, Un1

    real*8 :: delta_r
    
    do i = 1, coeff%iR-1
       
       Un (:) = coeff%u (:,i)

       delta_r = coeff%r (i+1) - coeff%r (i)

       !1. Fill matrix
       call fill_matrix(coeff%matrix(:,:,i+1),-delta_r,Mn, ierr)
       
       call fill_matrix(coeff%matrix(:,:,i),delta_r,Nn, ierr)
       
       !2. Calculate inverse matrix
       call Invert_matrix (Mn, Mn_inv, ierr)
       !3. Update values
       aux = MATMUL(Mn_inv, Nn)
       Un1 = MATMUL(aux,Un)
       
       coeff%u (:,i+1) = Un1
       
    enddo
   
  end subroutine Integrate_equation

  ! ==============================================================================================================
  ! ==============================================================================================================
  subroutine fill_matrix(coeff_matrix, delta_r, Mn, ierr)
    
    implicit none

    ! --- Arguments -------
    ! input 
    real*8, dimension(:,:), intent(in) :: coeff_matrix
    real*8, intent (in) :: delta_r
    ! output
    real*8, dimension(size(coeff_matrix,1),size(coeff_matrix,2)),intent(out) :: Mn
    integer, intent (out) :: ierr
    !---------------------

    integer n1

    ierr = 0

    Mn = 0.5d0 * delta_r * coeff_matrix 
    do n1 = 1,size(coeff_matrix,1)
       Mn (n1, n1) = Mn (n1, n1) + 1.0d0
    enddo

    
  end subroutine fill_matrix
 

  ! ============================================================================================================== 
  ! Returns the inverse of a matrix calculated by finding the LU
  ! decomposition.  Depends on LAPACK.
  ! ==============================================================================================================

  subroutine Invert_matrix(A, Ainv, ierr)

    implicit none

    !------ Arguments ----------
    real*8, dimension(:,:), intent(inout) :: A
    real*8, dimension(size(A,1),size(A,2)), intent(inout) :: Ainv
    integer :: ierr
    !--------------------
    real*8, dimension(size(A,1)) :: work  ! work array for LAPACK
    integer, dimension(size(A,1)) :: ipiv   ! pivot indices
    integer :: n, info
    real*8 :: detA

    ! External procedures defined in LAPACK
    external DGETRF
    external DGETRI

    ierr = 0

    n = size(A,1)

    if (n.le.5) then ! Small matrices (5x5 or smaller)

       Ainv = A
       select case (n)
       case (2) ! 2x2
          call invert2 (Ainv, n, n, detA) 
       case (3) ! 3x3
          call invert3 (Ainv, n, n, detA) 
       case (4) ! 4x4
          call invert4 (Ainv, n, n, detA) 
       case (5) ! 5x5
          call invert5 (Ainv, n, n, detA) 
       end select
       
       if (abs(detA).lt.tiny(1.0d0)) then ! Singular matrix
          ierr = 1111
          return
       else
          Ainv = Ainv / detA
       endif

    else ! General case 

       ! Store A in Ainv to prevent it from being overwritten by LAPACK
       Ainv = A

       ! DGETRF computes an LU factorization of a general M-by-N matrix A
       ! using partial pivoting with row interchanges.
       call DGETRF(n, n, Ainv, n, ipiv, info)

       if (info /= 0) then
!          print*,info
!          write (*, *) 'Matrix is numerically singular!'
          ierr = info
          return
       end if

       ! DGETRI computes the inverse of a matrix using the LU factorization
       ! computed by DGETRF.
       call DGETRI(n, Ainv, n, ipiv, work, n, info)

       if (info /= 0) then
!          write (*, *) 'Matrix inversion failed!', n,"x",n
          ierr = info
          return
       end if

    endif


  end subroutine Invert_matrix



!!$  ! ==============================================================================================================
!!$  ! Computes the coefficient matrix corresponding to the Cowling approximation (2x2)
!!$  ! ==============================================================================================================
!!$  subroutine Compute_coefficients_Cowling (data, param, w2, coeff, ierr)  
!!$    implicit none
!!$    
!!$    ! -----  Arguments  --------
!!$    ! input
!!$    type(BGData_t), intent (inout) :: data
!!$    type(parameters_t), intent (inout) :: param
!!$    real*8 :: w2
!!$    ! output
!!$    type (mode_coeff_t), intent (inout) :: coeff
!!$    integer, intent (out) :: ierr
!!$    ! ----------------------------
!!$
!!$    integer :: i
!!$
!!$    ierr = 0
!!$    
!!$    coeff%matrix = 0.0d0
!!$
!!$    if (param%l.ne.0) then !---- l != 0 ---------
!!$       do i = 1, data%iR
!!$
!!$          ! _____ Equation for eta_r______
!!$          ! .. eta_r
!!$          coeff%matrix (1, 1, i)= &
!!$               - 2.0d0/data%r(i) &
!!$               - data%ggrav(i) * data%inv_cs2(i)  &
!!$               - 6.0 *data%dphi_dr(i)/data%phi(i)
!!$          ! .. eta_p
!!$          coeff%matrix (1, 2, i) = &
!!$               param%l*(param%l+1.0)/data%r(i)**2 &
!!$               - data%phi(i)**4 * data%inv_cs2(i) / data%alpha(i)**2 * w2
!!$
!!$          ! ______ Equation for eta_p ______
!!$          ! .. eta_r
!!$          coeff%matrix (2, 1, i) = 1.0d0 - data%n2(i)/w2
!!$          ! .. eta_p
!!$          coeff%matrix (2, 2, i) = &
!!$               - data%dlnq_dr(i) &
!!$               + data%ggrav(i) * (1.0d0+data%inv_cs2(i))
!!$       enddo
!!$
!!$       call Interpolate_coefficients_Cowling (data, param, coeff, ierr)
!!$       
!!$    else ! ---- l = 0 ---------
!!$       
!!$       ierr = 101
!!$       print*, "ERROR: l=0 for Cowling not implemented"
!!$
!!$    endif
!!$
!!$    call Check_for_coefficients_errors (coeff, ierr)
!!$
!!$  end subroutine Compute_coefficients_Cowling
!!$ 
!!$  ! ==============================================================================================================
!!$  ! Regularizes first points by interpolating to the center (Cowling)
!!$  ! ==============================================================================================================
!!$  subroutine Interpolate_coefficients_Cowling (data, param, coeff, ierr)
!!$    
!!$    implicit none
!!$    
!!$    ! ------- Arguments -------    
!!$    ! input 
!!$    type(BGData_t), intent (inout) :: data
!!$    type(parameters_t), intent (inout) :: param
!!$    ! output
!!$    type (mode_coeff_t), intent (inout) :: coeff
!!$    integer, intent (out) :: ierr
!!$    ! --------------------------
!!$
!!$    integer :: i, n1, n2, ip
!!$    real*8 :: tmp
!!$
!!$    real*8 :: lim (1:coeff%n, 1:coeff%n)
!!$    integer :: rpower (1:coeff%n, 1:coeff%n)
!!$
!!$    ierr = 0
!!$
!!$    ip = param%ip_coeffs+1
!!$
!!$    !   Define coeffient limits r=0
!!$    lim (1, 1) = -2.0 ! /r
!!$    rpower (1, 1) = 1
!!$    lim (1, 2) = param%l*(param%l+1) ! /r^2
!!$    rpower (1, 2) = 2
!!$    lim (2, 1) = 1.0
!!$    rpower (2, 1) = 0
!!$    lim (2, 2) = 0.0
!!$    rpower (2, 2) = 0
!!$
!!$
!!$    do n1= 1, coeff%n
!!$       do n2= 1, coeff%n
!!$          do i=1,ip-1 
!!$             tmp = lim (n1, n2) &
!!$                  + ( coeff%matrix (n1, n2, ip) * data%r(ip) ** rpower (n1, n2) - lim (n1, n2)) &
!!$                  / data%r (ip) * data%r (i)
!!$             coeff%matrix (n1, n2, i) = tmp / data%r (i) ** rpower (n1, n2)
!!$          enddo
!!$       enddo
!!$    enddo
!!$
!!$  end subroutine Interpolate_coefficients_Cowling



end module module_mode_integrator
