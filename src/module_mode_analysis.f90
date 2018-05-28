!************************************************************************
! Module to manage to compute the eigenmodes of a 1D equilibrium model
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

module  module_mode_analysis
  
  use module_param
  use module_background
  use module_eigen
  use module_mode_integrator

  implicit none

  private
    
  public :: Perform_mode_analysis
  public :: Print_error_message


contains

 ! ==============================================================================================================
 ! ==============================================================================================================
  subroutine Perform_mode_analysis (data, param, ierr)

    implicit none

    ! ------- Arguments -------
    type(BGData_t), intent (inout) :: data
    type(parameters_t), intent (inout) :: param
    integer, intent (out) :: ierr
    ! ----------------------------

    integer ::  nw,i

    real*8 :: logf
    
    ! Frequencies in the coarse search grid
    real*8, allocatable, dimension (:) :: freqs
    ! Value of the function that has to be zero at the eigenvalue
    real*8, allocatable, dimension (:) :: bc_func


    real*8 :: xa, xb, ya, yb
    real*8 :: ym,xm
    integer :: type

    type (Eigen_t) :: eigen
    integer :: ne, nemax

    type (mode_coeff_t) :: coeff


    ierr = 0
    

    write (*, *) "-----------------------------------------------------------------"
    write (*, *) data%nt,data%time," s"
    write (*, *) "Shock location: ",data%iR,data%r(data%iR)


    allocate (bc_func(1:param%nwmax))
    allocate (freqs(1:param%nwmax))
 
    ! ----Print Parameters -----------
    write (*,*)  "Parameters: "
    
    ! ----------------------------

    if (param%force_compute_n2) data%calculate_n2 = .true.

    if (param%cowling) then
       write (*,*)  "   Cowling aprox."
    else
       write (*,*)  "   No cowling"
    endif
    write (*,*) "   ggrav_mode =",param%ggrav_mode
    if(param%cal_dphi) then
       write (*,*) "   Calculating delta_phi"
    else
       write (*,*) "   Not calculating delta_phi"
    endif
    write (*,*) "   iR, m= ",data%iR, data%m

    if (param%newtonian) then
       write (*,*) "   Newtonian calculation"
    else
       write (*,*) "   GR calculation"
    endif
    
 
    !------ Calculate coefficients for derivatives -------
    data%coeff_dp = 0.0d0
    data%coeff_d0 = 0.0d0
    data%coeff_dm = 0.0d0
    do i = 1, data%iR -1
       data%coeff_dp (i) = (data%r(i) - data%r(i-1)) &
            / ( (data%r(i+1) - data%r (i-1)) * (data%r(i+1) - data%r(i)) )
       data%coeff_d0 (i) = - (2.0d0 * data%r(i) - data%r (i+1) - data%r(i-1)) &
            / ( (data%r(i) - data%r (i-1)) * (data%r(i+1) - data%r(i)) )
       data%coeff_dm (i) = - (data%r(i+1) - data%r(i)) &
            / ( (data%r(i) - data%r (i-1)) * (data%r(i+1) - data%r(i-1)) )
    enddo
    ! one sided derivatives for the last point
    data%coeff_dp (data%iR) = 0.0d0
    data%coeff_d0 (data%iR) =  1.0d0 /(data%r(data%m)-data%r(data%m-1))
    data%coeff_dm (data%iR) = -1.0d0 /(data%r(data%m)-data%r(data%m-1))

    !--- Impose boundary conditions on input variables to ensure that 
    !    derivatives are possible -----------------------------------
    call Impose_boundary_conditions_bg (data, ierr)

    
    !------- 2nd order derivatives of background quantities -------
    
    do i = 1, data%m

       data%dphi_dr(i) = &
            data%phi (i + 1) *  data%coeff_dp (i) &
            + data%phi (i) *  data%coeff_d0 (i)&
            + data%phi (i - 1) *  data%coeff_dm (i)
       data%dalpha_dr(i) = &
            data%alpha (i + 1) *  data%coeff_dp (i) &
            + data%alpha (i) *  data%coeff_d0 (i) &
            + data%alpha (i - 1) *  data%coeff_dm (i)
       data%de_dr(i) = &
            (data%rho (i+1) * (1.0d0 + data%eps (i+1))) *  data%coeff_dp (i) &
            +  (data%rho (i) * (1.0d0 + data%eps (i))) *  data%coeff_d0 (i)&
            +  (data%rho (i-1) * (1.0d0 + data%eps (i-1))) *  data%coeff_dm (i)
       data%dq_dr(i) = &
            ((data%rho (i+1) *data%h (i+1) *  data%phi (i+1) ** 4 / data%alpha (i+1) ** 2) &
            *  data%coeff_dp (i) &
            +  (data%rho (i) *data%h (i) *  data%phi (i) ** 4 / data%alpha (i) ** 2) &
            *  data%coeff_d0 (i) &
            + (data%rho (i-1) *data%h (i-1) *  data%phi (i-1) ** 4 / data%alpha (i-1) ** 2) &
            * data%coeff_dm (i))
       data%dp_dr(i) = &
            data%p (i + 1) *  data%coeff_dp (i) &
            + data%p (i) *  data%coeff_d0 (i)&
            + data%p (i - 1) *  data%coeff_dm (i)
       data%drho_dr(i)= data%rho (i + 1)*  data%coeff_dp (i) &
            + data%rho (i) *  data%coeff_d0 (i ) &
            + data%rho (i - 1) *  data%coeff_dm (i)
       data%drhoh_dr(i)= &
            data%rho (i + 1)*data%h(i+1) *  data%coeff_dp (i) &
            + data%rho (i)*data%h(i) *  data%coeff_d0 (i ) &
            + data%rho (i - 1)*data%h(i-1) *  data%coeff_dm (i)
       data%dlnq_dr(i) = data%dq_dr(i) &
            / (data%rho(i)*data%h(i)*data%phi(i)**4/data%alpha(i)**2)
       data%dlnrhoh_dr(i)=data%drhoh_dr(i)/(data%rho(i)*data%h(i))
       
    enddo

    !--------- Compute extra quantities -----------------------------

    ! ... relativistic enthalpy
    where (abs(data%rho)>0.0d0) & ! non-zero rho
         data%h = 1.0d0 + data%eps + data%p / data%rho
    ! ... adiabatic index 
    where (abs(data%p)>0.0d0) & ! non-zero P
         data%gamma_one = data%rho*data%h*data%c_sound_squared / data%p
    ! ... gravitational acceleration
    select case (param%ggrav_mode)
    case (1)
       if (param%newtonian) then
          data%ggrav = -data%dalpha_dr
       else
          data%ggrav = -data%dalpha_dr/ data%alpha
       endif
    case (2)
       if (param%newtonian) then
          where (abs(data%rho)>0.0d0) data%ggrav = data%dp_dr / data%rho
       else
          where (abs(data%rho*data%h)>0.0d0) data%ggrav = data%dp_dr / data%rho / data%h
       endif
    case (3)
       if (param%newtonian) then
          where (abs(data%rho)>0.0d0) &
               data%ggrav = 0.5*(data%dp_dr/data%rho-data%dalpha_dr)
       else
          where (abs(data%rho*data%h)>0.0d0) &
               data%ggrav = 0.5*(data%dp_dr/data%rho/data%h-data%dalpha_dr/ data%alpha)
       endif
    case default
       write(*,*) "ERROR: perform_mode_analysis"
       write (*, *) "invalid ggrav_mode = ", param%ggrav_mode
       ierr = 10
       return
    end select
    !..... Inverse of sound speed square
    where (abs(data%c_sound_squared)>0.0d0)
       data%inv_cs2 = 1.0d0 / data%c_sound_squared
    elsewhere
       data%inv_cs2 = huge(1.0d0)
    end where
    !.... Buoyancy 
    where (abs(data%rho*data%h)>0.0d0) &
         data%Bstar = data%de_dr /data%rho/data%h &
         - data%dp_dr/data%rho/data%h*data%inv_cs2
    !.... Brunt Vaisala frequency
    if (data%calculate_n2)  then
       write (*,*)  "Calculating Brunt-Vaisala"
       data%n2 = data%alpha**2.0/data%phi**4.0 * data%ggrav * data%Bstar
    else
       write (*,*)  "Using Brunt-Vaisala from file"
    endif
    !.... Lamb frequency squared
    data%lamb2 = data%alpha**2/data%phi**4.0*data%c_sound_squared &
         * param%l * ( param%l + 1.0d0 )/ data%r ** 2
    !.... Q = phi*alpha
    data%Q = data%alpha * data%phi


    ! -------- Select calculation mode ---------------
    select case (param%mode)
    case (1) ! Full analysis
       write (*,*) "Analyzing all modes"
    case (2) ! Only p-modes (set buoyancy to zero)
       write (*,*)  "Analyzing p-modes"
       data%Bstar = 0.0d0
       data%n2 = 0.0d0
    case (3) ! Only g-modes (set cs2 to infinity and recompute Bstar)
       write (*,*)  "Analyzing g-modes"
       data%inv_cs2 = 0.0d0
       where (abs(data%rho*data%h)>0.0d0) &
            data%Bstar = data%de_dr /data%rho/data%h
    end select

    !----- Compute integrated quantities (mass ) -------------------
    call Compute_integrated_quantities_bg (data, ierr)

    !----- Impose boundary conditions to all quantities -------------
    call Impose_boundary_conditions_bg (data, ierr)

    !---- Check for errors in data-----
    call Check_for_errors_bg (data, ierr)
    if (ierr.ne.0) then
       write (*,*) "ERROR in perform_data_analysis"
       write (*, *) "data contains nans. ierr=", ierr
       return
    endif


    call Initialize_eigenvectors (data, eigen, ierr)

    call Initialize_coefficients_matrix (data, param, coeff, ierr)
    print*, "coeff%n =", coeff%n



    ! ------------------- Coarse search for eigenvalues ----------------------
    write (*,*)  "Coarse frequency search"
    loop_frequency: do nw = 1, param%nwmax !Frequency loop

       if (param%logf_mode) then
          logf = log10(param%fmin) &
               + dble(nw -1) / dble(param%nwmax-1) &
               * (log10(param%fmax) - log10(param%fmin)) 
          freqs(nw) = 10.0d0**logf
       else
          freqs(nw) = param%fmin + dble(nw -1) / dble(param%nwmax-1) &
               * (param%fmax - param%fmin)
       endif

       call Integrate_eigenfunction_equation (freqs (nw), data, param, eigen, coeff, bc_func(nw), ierr)

    end do loop_frequency

    !-------------------- Fine search for eigenvalues ----------- ----------------

    write (*,*)  "Fine frequency search"
    ne = 0
    loop_frequency2: do nw=1,param%nwmax-1 
       if (bc_func(nw)*bc_func(nw+1).lt.0.0) then

          xa = freqs(nw)
          xb = freqs(nw+1)
          ya = bc_func(nw)
          yb = bc_func(nw+1)

          call  bisection (data, param, eigen, coeff,&
               xa, xb, ya, yb, xm, ym, type, ierr)
          
          if (type.eq.0) then
             ne = ne + 1
             eigen%nw = nw
             eigen%freq = xm
             eigen%ne = ne
             call Normalize_eigenvectors (eigen, ierr)
             call Output_eigenvectors (data, param, eigen, ierr)
             write (*,*)  ne, nw, xm, xb-xa, "(",(xb-xa)/xm,")", ym
          elseif (type.eq.1) then
             write (*,*)  ne, nw, xm, xb-xa, "(",(xb-xa)/xm,")", ym, " <----- asymptote"
          elseif (type.eq.2) then
             write (*,*)  ne, nw, xm, xb-xa, "(",(xb-xa)/xm,")", ym, " <----- not converged"
          endif
          
       endif
    enddo loop_frequency2

    nemax = ne 
    
    print*, "Eigenvalues found: ", nemax
    
    call Free_eigenvectors (eigen, ierr)

    call Free_coefficients_matrix (coeff, ierr)

    
    deallocate (bc_func)
    deallocate (freqs)

  end subroutine perform_mode_analysis



  
! ==============================================================================================================
! ==============================================================================================================
  subroutine Integrate_eigenfunction_equation (freq, data, param, eigen, coeff, bc_func, ierr)
    implicit none

    !------ Arguments ----!
    real*8, intent (in) :: freq
    type(BGData_t), intent (inout) :: data
    type(parameters_t), intent (inout) :: param
    type(Eigen_t), intent (inout) :: eigen
    type (mode_coeff_t), intent (inout) :: coeff
    real*8, intent (out) :: bc_func
    integer, intent (out) :: ierr
    !---------------------

    integer :: j
    integer :: iR


    real*8, dimension (:,:), allocatable :: S,S_inv

    real*8 :: w2

    real*8 :: ufree (1:6), func (1:6)


    real*8 :: x,x1,x0
    real*8 :: z,z1,z0
    real*8 :: f0,f1,f2,dF_x,dF_z,f01,f10
    real*8 :: g0,g1,g2,dG_x,dG_z,g01,g10
    real*8 :: aux(2)
    real*8 :: errorx, errorz

    real*8, dimension (:), allocatable :: FG


    logical :: converged

    iR = data%iR

    allocate (S (1:2,1:2))
    allocate (S_inv (1:2,1:2))
    allocate(FG(1:2))

    !Initilize variables
    select case (param%funits)
    case (0) ! frequency in cm^-1
       w2 = ((2.0d0*pi*freq))**2
    case (1) ! frequency in Hz
       w2 = ((2.0d0*pi*freq)/c_light)**2    
    case (2) ! angular frequency in cm^-1
       w2 = freq**2   
    case default
       write (*, *) "ERROR in integrate_eigenfucntion_equation"
       write (*, *) "funits = ", param%funits
       return
       ierr = 1234
    end select
    ! ----------------- Cowling approximation --------------------------------------------
    if (param%cowling) then


       call Compute_coefficients_full (data, param, w2,coeff, ierr)
!       call  Compute_coefficients_Cowling (data, param, w2, coeff, ierr)

       ufree (1) = 0.1d0

       call Set_inner_boundary_conditions (ufree, data, param, coeff, ierr)
       
       call Integrate_equation (coeff, ierr)

       call Compute_outer_boundary_conditions (func, param, coeff, ierr)
       
       call Copy_solution (data, param, coeff, eigen, ierr)
       
    ! ----------------- non-Cowling approximation with delta Q (but delta psi = 0) ------------
    elseif(.not.param%cal_dphi) then

       call Compute_coefficients_full (data, param, w2,coeff, ierr)

       ufree (1) = 0.1d0
       x0 = 1.0d-8
       x1 = 1.0d-7
       

       ufree (4) = x0
       call Set_inner_boundary_conditions (ufree, data, param, coeff, ierr)       
       call Integrate_equation (coeff, ierr)       
       call Compute_outer_boundary_conditions (func, param, coeff, ierr)
       f0 = func (2)

       ufree (4) = x1
       call Set_inner_boundary_conditions (ufree, data, param, coeff, ierr)       
       call Integrate_equation (coeff, ierr)       
       call Compute_outer_boundary_conditions (func, param, coeff, ierr)
       f1 = func (2)

       converged = .false.
       do j = 1, param%maxiter
          
          x = x1 - (x1-x0)/(f1-f0)*f1

          ufree (4) = x
          call Set_inner_boundary_conditions (ufree, data, param, coeff, ierr)       
          call Integrate_equation (coeff, ierr)       
          call Compute_outer_boundary_conditions (func, param, coeff, ierr)
          f2 = func (2)
             
          if (abs(x-x1)/abs(x1).lt.param%tol) then
             converged = .true.
             exit
          endif
          f0 = f1
          f1 = f2
          x0 = x1
          x1 = x      
       enddo       

       if (.not.converged) then
!          write (*,*) "ERROR: shooting method not converged",i,freq,abs(x-x1)/abs(x1),x,x1
          ierr = 500
       endif

       call Copy_solution (data, param, coeff, eigen, ierr)

    ! ----------------- non-Cowling approximation with delta Q and delta psi ------------
    else
       

       call Compute_coefficients_full (data, param, w2,coeff, ierr)
       
       ufree (1) = 0.1d0
       if (param%l.ne.0) then 
          x0 = 1.0d-8
          z0 = 1.0d-8
          x1 = 1.0d-7
          z1 = 1.0d-7
       else
          x0 = 1.0d-8
          z0 = 1.0d-1
          x1 = 1.0d-7
          z1 = 1.0d-2
       endif

       ufree (4) = x0
       ufree (6) = z0
       call Set_inner_boundary_conditions (ufree, data, param, coeff, ierr)       
       call Integrate_equation (coeff, ierr)       
       call Compute_outer_boundary_conditions (func, param, coeff, ierr)
       f0 = func (2)
       g0 = func (3)

       ufree (4) = x1
       ufree (6) = z1
       call Set_inner_boundary_conditions (ufree, data, param, coeff, ierr)       
       call Integrate_equation (coeff, ierr)       
       call Compute_outer_boundary_conditions (func, param, coeff, ierr)
       f1 = func (2)
       g1 = func (3)

       converged = .false.
       do j = 1, param%maxiter
       
          ufree (4) = x1
          ufree (6) = z0
          call Set_inner_boundary_conditions (ufree, data, param, coeff, ierr)       
          call Integrate_equation (coeff, ierr)       
          call Compute_outer_boundary_conditions (func, param, coeff, ierr)
          f10 = func (2)
          g10 = func (3)
          
          ufree (4) = x0
          ufree (6) = z1
          call Set_inner_boundary_conditions (ufree, data, param, coeff, ierr)       
          call Integrate_equation (coeff, ierr)       
          call Compute_outer_boundary_conditions (func, param, coeff, ierr)
          f01 = func (2)
          g01 = func (3)

          
          !Derivatives
          dF_x = (f1-f01)/(x1-x0)          
          dG_x = (g1-g01)/(x1-x0)
          dF_z = (f1-f10)/(z1-z0)          
          dG_z = (g1-g10)/(z1-z0)
          
          !Secant matrix
          S(1,1) = dF_x
          S(1,2) = dF_z
          S(2,1) = dG_x
          S(2,2) = dG_z
          !Compute the inverse
          call Invert_matrix (S, S_inv, ierr)
          ! Singular matrix 
          if (ierr.ne.0) then
             S_inv (1, 1) = 1.0d0 / S(1, 1)
             S_inv (2, 2) = 1.0d0 / S(2, 2)
             S_inv (1, 2) = 0.0d0
             S_inv (2, 1) = 0.0d0
             ierr = 0
          endif
          
          !Function
          FG(1) = f1
          FG(2) = g1
          
          aux =  MATMUL(S_inv,FG)
             
          x = x1 - aux(1)
          z = z1 - aux(2)
          
          ufree (4) = x
          ufree (6) = z
          call Set_inner_boundary_conditions (ufree, data, param, coeff, ierr)       
          call Integrate_equation (coeff, ierr)       
          call Compute_outer_boundary_conditions (func, param, coeff, ierr)
          f2 = func (2)
          g2 = func (3)
          
          errorx = abs(x-x1)/abs(x)
          errorz = abs(z-z1)/abs(z)
    
          if (sqrt(errorx**2+errorz**2).lt.param%tol) then
             converged = .true.
             exit
          endif

          if (abs(aux(2)).gt.tiny(0.0d0)) then
             g0 = g1
             g1 = g2
             z0 = z1
             z1 = z
          endif

          if (abs(aux(1)).gt.tiny(0.0d0)) then
             f0 = f1
             f1 = f2         
             x0 = x1
             x1 = x
          endif
       enddo
       if (.not.converged) then
!          write (*,*) "ERROR: shooting method not converged", freq,sqrt(errorx**2+errorz**2)
          ierr = 500
       endif
       
       call Copy_solution (data, param, coeff, eigen, ierr)
   
    endif
    !------------------ Integration completed ...........
    

    deallocate (S)
    deallocate (S_inv)
    deallocate(FG)

    ! return one function (eta_r at the shock) that has to be zero for the eigenvalue
    bc_func = func (1)


  end subroutine Integrate_eigenfunction_equation


  ! ==============================================================================================================
  ! ==============================================================================================================  
  subroutine bisection (data, param, eigen, coeff, xa, xb, ya, yb, xm, ym, type, ierr)
    
    implicit none
    !------ Arguments ----!
    type(BGData_t), intent (inout) :: data
    type(parameters_t), intent (inout) :: param
    type(Eigen_t), intent (inout) :: eigen
    type (mode_coeff_t) :: coeff
    !> Extreme points of the interval
    real*8, intent (inout) :: xa,xb
    !> Values of the fuction at extreme points 
    real*8, intent (inout) :: ya,yb
    !  Output
    !> type of change of sign
    !>             0: zero
    !>             1: asymptote
    integer, intent (out) :: type    
    real*8, intent (out) :: xm,ym
    integer, intent (out) :: ierr
    !----------------------
    
    integer :: n
    integer, parameter  :: nmax = 100 ! Ensures decrease of 30 orders of magnitude in x

    real*8 :: ymean_old,ymean
    real*8 :: xtol, ytol


    logical error
    
    type = 0
    xtol = param%tol
!    ytol = 1.0d-5
    ytol = param%tol

    error = .true.
    
    do n = 1, nmax
       
       xm = (xa+xb)*0.5d0
       
       call Integrate_eigenfunction_equation (xm, data, param, eigen, coeff, ym, ierr)

       ! compute mean value to check if is a zero or an asymptota
       ymean_old = 0.5*(abs(ya)+abs(yb))
       ! bisection
       if (ya*ym.lt.0.0d0) then
          xb=xm
          yb=ym
       else
          xa=xm
          ya=ym
       endif
       ymean = 0.5*(abs(ya)+abs(yb))

       ! if ymean is growing it is an asymptote
       if (n.gt.1.and.ymean.gt.ymean_old) then 
          ! if shooting did not converge but it is an asymptote, reset ierr=0
          if (ierr.eq.500) ierr = 0
          type = 1
          error = .false.
          exit
       endif
       ! stop criterium
       if (n.gt.1.and.(abs(ym).lt.ytol.and.abs((xa-xb)/xa).lt.xtol)) then
          error = .false.
          exit
       endif

    enddo

    if (error) then
!       write (*,*)  "Bisection not converged", xm, xa, xb, ym, ya, yb &
!            , ymean, ymean_old
!       write (*,*)  abs(xb-xa)/abs(xb)
!       read*
       type = 2
       ierr = 700
    endif

    return

  end subroutine bisection


  subroutine Print_error_message (ierr)

    integer, intent (in) :: ierr

    if (ierr.eq.0) return

    write (*, *) "***********************************************"
    write (*, *) "ERROR:"

    select case (ierr)
    case (500)
       write (*,*) "ERROR: shooting method not converged"
    case (700)
       write (*,*)  "Bisection not converged"
    case default
       write (*, *) " ierr =", ierr
    end select
    write (*, *) "***********************************************"
       


  end subroutine Print_error_message


end module module_mode_analysis
