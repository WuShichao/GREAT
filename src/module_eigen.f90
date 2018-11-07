!************************************************************************
! Module to manage the eigenvalues and eigenfunctions
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

module module_eigen

  use module_param
  use module_background

  implicit none
  
  private

  logical, save :: freqs_file_exists = .false.

  type, public :: Eigen_t
     !> Size of the grid
     integer :: iR
     !> radial grid
     real*8, allocatable, dimension (:) :: r
     !> eta_r: radial perturbation
     real*8, allocatable, dimension (:) :: nr
     !> eta_perp: perpendicular perturbation
     real*8, allocatable, dimension (:) :: np
     !> delta P: Pressure perturbation
     real*8, allocatable, dimension (:) :: delta_p
     !> metric perturbation delta Q
     real*8, allocatable, dimension (:) :: delta_Q
     !> metric perturbation delta Psi
     real*8, allocatable, dimension (:) :: delta_phi
     !> radial derivative of delta Q
     real*8, allocatable, dimension (:) :: k_cap
     !> radial derivative of delta psi
     real*8, allocatable, dimension (:) :: phi_cap
     !> frequency
     real*8 :: freq
     !> nw
     integer :: nw
     !> index
     integer :: ne
  end type Eigen_t

  public :: Initialize_eigenvectors
  public :: Free_eigenvectors
  public :: Normalize_eigenvectors
  public :: Output_eigenvectors


contains



  !****************************************************************************
  ! Initialize eigenvectors
  !****************************************************************************  
  subroutine Initialize_eigenvectors (data, eigen, ierr)
    implicit none
    ! -----  Arguments  --------
    type(BGData_t), intent (inout) :: data
    type(Eigen_t), intent (inout) :: eigen
    integer, intent (out) :: ierr
    ! --------------------------
    
    ierr = 0

    allocate (eigen%r (1:data%iR))
    allocate (eigen%nr (1:data%iR))
    allocate (eigen%np (1:data%iR))
    allocate (eigen%delta_p (1:data%iR))
    allocate (eigen%delta_Q (1:data%iR))
    allocate (eigen%delta_phi (1:data%iR))
    allocate (eigen%k_cap (1:data%iR))
    allocate (eigen%phi_cap (1:data%iR))
    
    eigen%iR = data%iR
    eigen%r (1:data%iR) = data%r (1:data%iR)
    eigen%nr = 0.0d0
    eigen%np = 0.0d0
    eigen%delta_p = 0.0d0
    eigen%delta_Q = 0.0d0
    eigen%delta_phi = 0.0d0
    eigen%k_cap = 0.0d0
    eigen%phi_cap = 0.0d0

  end subroutine Initialize_eigenvectors

  !****************************************************************************
  ! Free eigenvectors
  !****************************************************************************  
  subroutine Free_eigenvectors (eigen, ierr)
    implicit none
    ! -----  Arguments  --------
    type(Eigen_t), intent (inout) :: eigen
    integer, intent (out) :: ierr
    ! --------------------------
    
    ierr = 0

    deallocate (eigen%r)
    deallocate (eigen%nr)
    deallocate (eigen%np)
    deallocate (eigen%delta_p)
    deallocate (eigen%delta_Q)
    deallocate (eigen%delta_phi)
    deallocate (eigen%k_cap)
    deallocate (eigen%phi_cap)
    
  end subroutine Free_eigenvectors


  !****************************************************************************
  ! Normalize eigenvectors
  !****************************************************************************  
  subroutine Normalize_eigenvectors (eigen, ierr)
    implicit none
    ! -----  Arguments  --------
    type(Eigen_t), intent (inout) :: eigen
    integer, intent (out) :: ierr
    ! --------------------------

    real*8 :: factor

    ierr = 0

    factor = max(maxval(abs(eigen%nr)), maxval(abs(eigen%np/eigen%r)))

    eigen%nr = eigen%nr / factor
    eigen%np = eigen%np / factor
    eigen%delta_p = eigen%delta_p / factor
    eigen%delta_Q = eigen%delta_Q / factor
    eigen%delta_phi = eigen%delta_phi / factor
    eigen%k_cap = eigen%k_cap / factor
    eigen%phi_cap = eigen%phi_cap / factor


  end subroutine Normalize_eigenvectors

  !****************************************************************************
  ! Output eigenvectors
  !****************************************************************************  
  subroutine Output_eigenvectors (data, param, eigen, ierr)
    implicit none
    ! -----  Arguments  --------
    type(BGData_t), intent (inout) :: data
    type(parameters_t), intent (inout) :: param
    type(Eigen_t), intent (inout) :: eigen
    integer, intent (out) :: ierr
    ! --------------------------

    character (len=255) :: fileout
    character (len=25) :: fstring

    integer :: i

    ierr = 0

    select case (param%funits)
    case (0) 
       fstring = "freq [cm^-1]"
    case (1) 
       fstring = "freq [Hz]"
    case (2) 
       fstring = "omega [cm^-1]"
    end select


    write (fileout,"(a6,i0.4,a4)") "eigen_",data%nt,".dat"
    fileout = trim(data%output_directory)//"/"//trim(fileout)

    if (eigen%ne.eq.1) then
       write (*, *) "Output eigenfunctions: ", fileout
       open (22, file=fileout,status='unknown', action='write')
       write (22, "(A2,I8,A20)") "# ", eigen%iR, "cells"
       write (22, "(A2,A23,25A25)") "# ",&
            "time [s]", &
            trim(fstring), &
            "r [cm]", &
            "eta_r", &
            "eta_perp/r", &
            "delta P", &
            "delta Q", &
            "delta psi", &
            "K", &
            "Psi"
       
    else
       open (22, file =fileout,status='old', position='append', action='write')
    endif

    do i = 1, eigen%iR
       write (22, "(25E25.16)") &
            data%time, &
            eigen%freq, &
            eigen%r(i), &
            eigen%nr(i), &
            eigen%np(i)/data%r(i), &
            eigen%delta_p(i), &
            eigen%delta_Q(i), &
            eigen%delta_phi(i), &
            eigen%k_cap(i), &
            eigen%phi_cap(i)
       
    enddo
    write (22, *) " "

    close (22)

    if (param%restart) freqs_file_exists = .true.

    fileout = trim(data%output_directory)//"/"//"freqs.dat"

    if (.not.freqs_file_exists) then
       write (*, *) "Output eigenvalues: ", fileout
       open (22, file=fileout,status='unknown', action='write')
       write (22, "(A2,A23,25A25)") "# ", "time [s]", "nt", trim(fstring)
       freqs_file_exists = .true.
    else
       open (22, file =fileout,status='old', position='append', action='write')
    endif

    write (22, "(E25.15,I25,E25.15)") data%time, eigen%ne, eigen%freq

    close (22)



  end subroutine Output_eigenvectors


end module module_eigen
