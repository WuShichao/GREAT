!****************************************************************************************
! Module with fast inversion routines for matrices smaller or equal than 5x5
! Subroutines from QMC=Chem (https://github.com/scemama/qmcchem) with GPL License 
! downloaded from: 
! https://raw.githubusercontent.com/scemama/qmcchem/master/src/TOOLS/invert.irp.f
! Wrapped in a module for easy use with minor modifications
!*****************************************************************************************
!
! GPL  License:
! -------------
!
!    GREAT = General Relativistic Eigenmode Analysis Tool
!
!    Copyright (C) 2009 Anthony SCEMAMA 
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
module module_fast_inv

  implicit none
 
  private

  public :: invert2
  public :: invert3
  public :: invert4
  public :: invert5


contains

!***********************************************************************************************************************************
subroutine invert2(a,LDA,na,det_l)
 implicit none
 double precision :: a (LDA,na)
 integer          :: LDA
 integer          :: na
 double precision :: det_l
 double precision :: b(2,2)

 b(1,1) = a(1,1)
 b(2,1) = a(2,1)
 b(1,2) = a(1,2)
 b(2,2) = a(2,2)
 det_l = a(1,1)*a(2,2) - a(1,2)*a(2,1)
 a(1,1) = b(2,2)
 a(2,1) = -b(2,1)
 a(1,2) = -b(1,2)
 a(2,2) = b(1,1)
end subroutine invert2

!***********************************************************************************************************************************
subroutine invert3(a,LDA,na,det_l)
 implicit none
 double precision, intent(inout) :: a (LDA,na)
 integer, intent(in)             :: LDA
 integer, intent(in)             :: na
 double precision, intent(inout) :: det_l
 double precision :: b(4,3)
 !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: b
 integer :: i
 det_l = a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2)) &
        -a(1,2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1)) &
        +a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))
 do i=1,4
  b(i,1) = a(i,1)
  b(i,2) = a(i,2)
  b(i,3) = a(i,3)
 enddo
 a(1,1) =  b(2,2)*b(3,3) - b(2,3)*b(3,2)
 a(2,1) =  b(2,3)*b(3,1) - b(2,1)*b(3,3)
 a(3,1) =  b(2,1)*b(3,2) - b(2,2)*b(3,1)

 a(1,2) =  b(1,3)*b(3,2) - b(1,2)*b(3,3)
 a(2,2) =  b(1,1)*b(3,3) - b(1,3)*b(3,1)
 a(3,2) =  b(1,2)*b(3,1) - b(1,1)*b(3,2)

 a(1,3) =  b(1,2)*b(2,3) - b(1,3)*b(2,2)
 a(2,3) =  b(1,3)*b(2,1) - b(1,1)*b(2,3)
 a(3,3) =  b(1,1)*b(2,2) - b(1,2)*b(2,1)

end subroutine invert3

!***********************************************************************************************************************************
subroutine invert4(a,LDA,na,det_l)
 implicit none
 double precision, intent(inout) :: a (LDA,na)
 integer, intent(in)             :: LDA
 integer, intent(in)             :: na
 double precision, intent(inout) :: det_l
 double precision :: b(4,4)
 !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: b
 integer :: i
 det_l =  a(1,1)*(a(2,2)*(a(3,3)*a(4,4)-a(3,4)*a(4,3))  &
                 -a(2,3)*(a(3,2)*a(4,4)-a(3,4)*a(4,2))  &
                 +a(2,4)*(a(3,2)*a(4,3)-a(3,3)*a(4,2))) &
         -a(1,2)*(a(2,1)*(a(3,3)*a(4,4)-a(3,4)*a(4,3))  &
                 -a(2,3)*(a(3,1)*a(4,4)-a(3,4)*a(4,1))  &
                 +a(2,4)*(a(3,1)*a(4,3)-a(3,3)*a(4,1))) &
         +a(1,3)*(a(2,1)*(a(3,2)*a(4,4)-a(3,4)*a(4,2))  &
                 -a(2,2)*(a(3,1)*a(4,4)-a(3,4)*a(4,1))  &
                 +a(2,4)*(a(3,1)*a(4,2)-a(3,2)*a(4,1))) &
         -a(1,4)*(a(2,1)*(a(3,2)*a(4,3)-a(3,3)*a(4,2))  &
                 -a(2,2)*(a(3,1)*a(4,3)-a(3,3)*a(4,1))  &
                 +a(2,3)*(a(3,1)*a(4,2)-a(3,2)*a(4,1)))
 do i=1,4
  b(1,i) = a(1,i)
  b(2,i) = a(2,i)
  b(3,i) = a(3,i)
  b(4,i) = a(4,i)
 enddo
 
 a(1,1) =  b(2,2)*(b(3,3)*b(4,4)-b(3,4)*b(4,3))-b(2,3)*(b(3,2)*b(4,4)-b(3,4)*b(4,2))+b(2,4)*(b(3,2)*b(4,3)-b(3,3)*b(4,2))
 a(2,1) = -b(2,1)*(b(3,3)*b(4,4)-b(3,4)*b(4,3))+b(2,3)*(b(3,1)*b(4,4)-b(3,4)*b(4,1))-b(2,4)*(b(3,1)*b(4,3)-b(3,3)*b(4,1))
 a(3,1) =  b(2,1)*(b(3,2)*b(4,4)-b(3,4)*b(4,2))-b(2,2)*(b(3,1)*b(4,4)-b(3,4)*b(4,1))+b(2,4)*(b(3,1)*b(4,2)-b(3,2)*b(4,1))
 a(4,1) = -b(2,1)*(b(3,2)*b(4,3)-b(3,3)*b(4,2))+b(2,2)*(b(3,1)*b(4,3)-b(3,3)*b(4,1))-b(2,3)*(b(3,1)*b(4,2)-b(3,2)*b(4,1))

 a(1,2) = -b(1,2)*(b(3,3)*b(4,4)-b(3,4)*b(4,3))+b(1,3)*(b(3,2)*b(4,4)-b(3,4)*b(4,2))-b(1,4)*(b(3,2)*b(4,3)-b(3,3)*b(4,2))
 a(2,2) =  b(1,1)*(b(3,3)*b(4,4)-b(3,4)*b(4,3))-b(1,3)*(b(3,1)*b(4,4)-b(3,4)*b(4,1))+b(1,4)*(b(3,1)*b(4,3)-b(3,3)*b(4,1))
 a(3,2) = -b(1,1)*(b(3,2)*b(4,4)-b(3,4)*b(4,2))+b(1,2)*(b(3,1)*b(4,4)-b(3,4)*b(4,1))-b(1,4)*(b(3,1)*b(4,2)-b(3,2)*b(4,1))
 a(4,2) =  b(1,1)*(b(3,2)*b(4,3)-b(3,3)*b(4,2))-b(1,2)*(b(3,1)*b(4,3)-b(3,3)*b(4,1))+b(1,3)*(b(3,1)*b(4,2)-b(3,2)*b(4,1))

 a(1,3) =  b(1,2)*(b(2,3)*b(4,4)-b(2,4)*b(4,3))-b(1,3)*(b(2,2)*b(4,4)-b(2,4)*b(4,2))+b(1,4)*(b(2,2)*b(4,3)-b(2,3)*b(4,2))
 a(2,3) = -b(1,1)*(b(2,3)*b(4,4)-b(2,4)*b(4,3))+b(1,3)*(b(2,1)*b(4,4)-b(2,4)*b(4,1))-b(1,4)*(b(2,1)*b(4,3)-b(2,3)*b(4,1))
 a(3,3) =  b(1,1)*(b(2,2)*b(4,4)-b(2,4)*b(4,2))-b(1,2)*(b(2,1)*b(4,4)-b(2,4)*b(4,1))+b(1,4)*(b(2,1)*b(4,2)-b(2,2)*b(4,1))
 a(4,3) = -b(1,1)*(b(2,2)*b(4,3)-b(2,3)*b(4,2))+b(1,2)*(b(2,1)*b(4,3)-b(2,3)*b(4,1))-b(1,3)*(b(2,1)*b(4,2)-b(2,2)*b(4,1))

 a(1,4) = -b(1,2)*(b(2,3)*b(3,4)-b(2,4)*b(3,3))+b(1,3)*(b(2,2)*b(3,4)-b(2,4)*b(3,2))-b(1,4)*(b(2,2)*b(3,3)-b(2,3)*b(3,2))
 a(2,4) =  b(1,1)*(b(2,3)*b(3,4)-b(2,4)*b(3,3))-b(1,3)*(b(2,1)*b(3,4)-b(2,4)*b(3,1))+b(1,4)*(b(2,1)*b(3,3)-b(2,3)*b(3,1))
 a(3,4) = -b(1,1)*(b(2,2)*b(3,4)-b(2,4)*b(3,2))+b(1,2)*(b(2,1)*b(3,4)-b(2,4)*b(3,1))-b(1,4)*(b(2,1)*b(3,2)-b(2,2)*b(3,1))
 a(4,4) =  b(1,1)*(b(2,2)*b(3,3)-b(2,3)*b(3,2))-b(1,2)*(b(2,1)*b(3,3)-b(2,3)*b(3,1))+b(1,3)*(b(2,1)*b(3,2)-b(2,2)*b(3,1))

end subroutine invert4

!***********************************************************************************************************************************
subroutine invert5(a,LDA,na,det_l)
 implicit none
 double precision, intent(inout) :: a (LDA,na)
 integer, intent(in)             :: LDA
 integer, intent(in)             :: na
 double precision, intent(inout) :: det_l
 double precision :: b(5,5)
 !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: b
 integer :: i
 det_l = a(1,1)*(a(2,2)*(a(3,3)*(a(4,4)*a(5,5)-a(4,5)*a(5,4))-a(3,4)*( &
      a(4,3)*a(5,5)-a(4,5)*a(5,3))+a(3,5)*(a(4,3)*a(5,4)-a(4,4)*a(5,3)))- &
 a(2,3)*(a(3,2)*(a(4,4)*a(5,5)-a(4,5)*a(5,4))-a(3,4)*(a(4,2)*a(5,5)- &
 a(4,5)*a(5,2))+a(3,5)*(a(4,2)*a(5,4)-a(4,4)*a(5,2)))+a(2,4)*(a(3,2)*( &
 a(4,3)*a(5,5)-a(4,5)*a(5,3))-a(3,3)*(a(4,2)*a(5,5)-a(4,5)*a(5,2))+ &
 a(3,5)*(a(4,2)*a(5,3)-a(4,3)*a(5,2)))-a(2,5)*(a(3,2)*(a(4,3)*a(5,4)- &
 a(4,4)*a(5,3))-a(3,3)*(a(4,2)*a(5,4)-a(4,4)*a(5,2))+a(3,4)*(a(4,2)* &
 a(5,3)-a(4,3)*a(5,2))))-a(1,2)*(a(2,1)*(a(3,3)*(a(4,4)*a(5,5)-a(4,5)* &
 a(5,4))-a(3,4)*(a(4,3)*a(5,5)-a(4,5)*a(5,3))+a(3,5)*(a(4,3)*a(5,4)- &
 a(4,4)*a(5,3)))-a(2,3)*(a(3,1)*(a(4,4)*a(5,5)-a(4,5)*a(5,4))-a(3,4)*( &
 a(4,1)*a(5,5)-a(4,5)*a(5,1))+a(3,5)*(a(4,1)*a(5,4)-a(4,4)*a(5,1)))+ &
 a(2,4)*(a(3,1)*(a(4,3)*a(5,5)-a(4,5)*a(5,3))-a(3,3)*(a(4,1)*a(5,5)- &
 a(4,5)*a(5,1))+a(3,5)*(a(4,1)*a(5,3)-a(4,3)*a(5,1)))-a(2,5)*(a(3,1)*( &
 a(4,3)*a(5,4)-a(4,4)*a(5,3))-a(3,3)*(a(4,1)*a(5,4)-a(4,4)*a(5,1))+ &
 a(3,4)*(a(4,1)*a(5,3)-a(4,3)*a(5,1))))+a(1,3)*(a(2,1)*(a(3,2)*(a(4,4)* &
 a(5,5)-a(4,5)*a(5,4))-a(3,4)*(a(4,2)*a(5,5)-a(4,5)*a(5,2))+a(3,5)*( &
 a(4,2)*a(5,4)-a(4,4)*a(5,2)))-a(2,2)*(a(3,1)*(a(4,4)*a(5,5)-a(4,5)* &
 a(5,4))-a(3,4)*(a(4,1)*a(5,5)-a(4,5)*a(5,1))+a(3,5)*(a(4,1)*a(5,4)- &
 a(4,4)*a(5,1)))+a(2,4)*(a(3,1)*(a(4,2)*a(5,5)-a(4,5)*a(5,2))-a(3,2)*( &
 a(4,1)*a(5,5)-a(4,5)*a(5,1))+a(3,5)*(a(4,1)*a(5,2)-a(4,2)*a(5,1)))- &
 a(2,5)*(a(3,1)*(a(4,2)*a(5,4)-a(4,4)*a(5,2))-a(3,2)*(a(4,1)*a(5,4)- &
 a(4,4)*a(5,1))+a(3,4)*(a(4,1)*a(5,2)-a(4,2)*a(5,1))))-a(1,4)*(a(2,1)*( &
 a(3,2)*(a(4,3)*a(5,5)-a(4,5)*a(5,3))-a(3,3)*(a(4,2)*a(5,5)-a(4,5)* &
 a(5,2))+a(3,5)*(a(4,2)*a(5,3)-a(4,3)*a(5,2)))-a(2,2)*(a(3,1)*(a(4,3)* &
 a(5,5)-a(4,5)*a(5,3))-a(3,3)*(a(4,1)*a(5,5)-a(4,5)*a(5,1))+a(3,5)*( &
 a(4,1)*a(5,3)-a(4,3)*a(5,1)))+a(2,3)*(a(3,1)*(a(4,2)*a(5,5)-a(4,5)* &
 a(5,2))-a(3,2)*(a(4,1)*a(5,5)-a(4,5)*a(5,1))+a(3,5)*(a(4,1)*a(5,2)- &
 a(4,2)*a(5,1)))-a(2,5)*(a(3,1)*(a(4,2)*a(5,3)-a(4,3)*a(5,2))-a(3,2)*( &
 a(4,1)*a(5,3)-a(4,3)*a(5,1))+a(3,3)*(a(4,1)*a(5,2)-a(4,2)*a(5,1))))+ &
 a(1,5)*(a(2,1)*(a(3,2)*(a(4,3)*a(5,4)-a(4,4)*a(5,3))-a(3,3)*(a(4,2)* &
 a(5,4)-a(4,4)*a(5,2))+a(3,4)*(a(4,2)*a(5,3)-a(4,3)*a(5,2)))-a(2,2)*( &
 a(3,1)*(a(4,3)*a(5,4)-a(4,4)*a(5,3))-a(3,3)*(a(4,1)*a(5,4)-a(4,4)* &
 a(5,1))+a(3,4)*(a(4,1)*a(5,3)-a(4,3)*a(5,1)))+a(2,3)*(a(3,1)*(a(4,2)* &
 a(5,4)-a(4,4)*a(5,2))-a(3,2)*(a(4,1)*a(5,4)-a(4,4)*a(5,1))+a(3,4)*( &
 a(4,1)*a(5,2)-a(4,2)*a(5,1)))-a(2,4)*(a(3,1)*(a(4,2)*a(5,3)-a(4,3)* &
 a(5,2))-a(3,2)*(a(4,1)*a(5,3)-a(4,3)*a(5,1))+a(3,3)*(a(4,1)*a(5,2)- &
 a(4,2)*a(5,1))))

 do i=1,5
  b(1,i) = a(1,i)
  b(2,i) = a(2,i)
  b(3,i) = a(3,i)
  b(4,i) = a(4,i)
  b(5,i) = a(5,i)
 enddo
 
 a(1,1) = &
 (b(2,2)*(b(3,3)*(b(4,4)*b(5,5)-b(4,5)*b(5,4))-b(3,4)*(b(4,3)*b(5,5)-b(4,5)*b(5,3))+b(3,5)*(b(4,3)*b(5,4)-b(4,4)*b(5,3)))-b(2,3)* &
 (b(3,2)*(b(4,4)*b(5,5)-b(4,5)*b(5,4))-b(3,4)*(b(4,2)*b(5,5)-b(4,5)*b(5,2))+b(3,5)*(b(4,2)*b(5,4)-b(4,4)*b(5,2)))+b(2,4)* &
 (b(3,2)*(b(4,3)*b(5,5)-b(4,5)*b(5,3))-b(3,3)*(b(4,2)*b(5,5)-b(4,5)*b(5,2))+b(3,5)*(b(4,2)*b(5,3)-b(4,3)*b(5,2)))-b(2,5)* &
 (b(3,2)*(b(4,3)*b(5,4)-b(4,4)*b(5,3))-b(3,3)*(b(4,2)*b(5,4)-b(4,4)*b(5,2))+b(3,4)*(b(4,2)*b(5,3)-b(4,3)*b(5,2))))
 a(2,1) = &
 (-b(2,1)*(b(3,3)*(b(4,4)*b(5,5)-b(4,5)*b(5,4))-b(3,4)*(b(4,3)*b(5,5)-b(4,5)*b(5,3))+b(3,5)*(b(4,3)*b(5,4)-b(4,4)*b(5,3)))+b(2,3)* &
 (b(3,1)*(b(4,4)*b(5,5)-b(4,5)*b(5,4))-b(3,4)*(b(4,1)*b(5,5)-b(4,5)*b(5,1))+b(3,5)*(b(4,1)*b(5,4)-b(4,4)*b(5,1)))-b(2,4)* &
 (b(3,1)*(b(4,3)*b(5,5)-b(4,5)*b(5,3))-b(3,3)*(b(4,1)*b(5,5)-b(4,5)*b(5,1))+b(3,5)*(b(4,1)*b(5,3)-b(4,3)*b(5,1)))+b(2,5)* &
 (b(3,1)*(b(4,3)*b(5,4)-b(4,4)*b(5,3))-b(3,3)*(b(4,1)*b(5,4)-b(4,4)*b(5,1))+b(3,4)*(b(4,1)*b(5,3)-b(4,3)*b(5,1))))
 a(3,1) = &
 (b(2,1)*(b(3,2)*(b(4,4)*b(5,5)-b(4,5)*b(5,4))-b(3,4)*(b(4,2)*b(5,5)-b(4,5)*b(5,2))+b(3,5)*(b(4,2)*b(5,4)-b(4,4)*b(5,2)))-b(2,2)* &
 (b(3,1)*(b(4,4)*b(5,5)-b(4,5)*b(5,4))-b(3,4)*(b(4,1)*b(5,5)-b(4,5)*b(5,1))+b(3,5)*(b(4,1)*b(5,4)-b(4,4)*b(5,1)))+b(2,4)* &
 (b(3,1)*(b(4,2)*b(5,5)-b(4,5)*b(5,2))-b(3,2)*(b(4,1)*b(5,5)-b(4,5)*b(5,1))+b(3,5)*(b(4,1)*b(5,2)-b(4,2)*b(5,1)))-b(2,5)* &
 (b(3,1)*(b(4,2)*b(5,4)-b(4,4)*b(5,2))-b(3,2)*(b(4,1)*b(5,4)-b(4,4)*b(5,1))+b(3,4)*(b(4,1)*b(5,2)-b(4,2)*b(5,1))))
 a(4,1) = &
 (-b(2,1)*(b(3,2)*(b(4,3)*b(5,5)-b(4,5)*b(5,3))-b(3,3)*(b(4,2)*b(5,5)-b(4,5)*b(5,2))+b(3,5)*(b(4,2)*b(5,3)-b(4,3)*b(5,2)))+b(2,2)* &
 (b(3,1)*(b(4,3)*b(5,5)-b(4,5)*b(5,3))-b(3,3)*(b(4,1)*b(5,5)-b(4,5)*b(5,1))+b(3,5)*(b(4,1)*b(5,3)-b(4,3)*b(5,1)))-b(2,3)* &
 (b(3,1)*(b(4,2)*b(5,5)-b(4,5)*b(5,2))-b(3,2)*(b(4,1)*b(5,5)-b(4,5)*b(5,1))+b(3,5)*(b(4,1)*b(5,2)-b(4,2)*b(5,1)))+b(2,5)* &
 (b(3,1)*(b(4,2)*b(5,3)-b(4,3)*b(5,2))-b(3,2)*(b(4,1)*b(5,3)-b(4,3)*b(5,1))+b(3,3)*(b(4,1)*b(5,2)-b(4,2)*b(5,1))))
 a(5,1) = &
 (b(2,1)*(b(3,2)*(b(4,3)*b(5,4)-b(4,4)*b(5,3))-b(3,3)*(b(4,2)*b(5,4)-b(4,4)*b(5,2))+b(3,4)*(b(4,2)*b(5,3)-b(4,3)*b(5,2)))-b(2,2)* &
 (b(3,1)*(b(4,3)*b(5,4)-b(4,4)*b(5,3))-b(3,3)*(b(4,1)*b(5,4)-b(4,4)*b(5,1))+b(3,4)*(b(4,1)*b(5,3)-b(4,3)*b(5,1)))+b(2,3)* &
 (b(3,1)*(b(4,2)*b(5,4)-b(4,4)*b(5,2))-b(3,2)*(b(4,1)*b(5,4)-b(4,4)*b(5,1))+b(3,4)*(b(4,1)*b(5,2)-b(4,2)*b(5,1)))-b(2,4)* &
 (b(3,1)*(b(4,2)*b(5,3)-b(4,3)*b(5,2))-b(3,2)*(b(4,1)*b(5,3)-b(4,3)*b(5,1))+b(3,3)*(b(4,1)*b(5,2)-b(4,2)*b(5,1))))

 a(1,2) = &
 (-b(1,2)*(b(3,3)*(b(4,4)*b(5,5)-b(4,5)*b(5,4))-b(3,4)*(b(4,3)*b(5,5)-b(4,5)*b(5,3))+b(3,5)*(b(4,3)*b(5,4)-b(4,4)*b(5,3)))+b(1,3)* &
 (b(3,2)*(b(4,4)*b(5,5)-b(4,5)*b(5,4))-b(3,4)*(b(4,2)*b(5,5)-b(4,5)*b(5,2))+b(3,5)*(b(4,2)*b(5,4)-b(4,4)*b(5,2)))-b(1,4)* &
 (b(3,2)*(b(4,3)*b(5,5)-b(4,5)*b(5,3))-b(3,3)*(b(4,2)*b(5,5)-b(4,5)*b(5,2))+b(3,5)*(b(4,2)*b(5,3)-b(4,3)*b(5,2)))+b(1,5)* &
 (b(3,2)*(b(4,3)*b(5,4)-b(4,4)*b(5,3))-b(3,3)*(b(4,2)*b(5,4)-b(4,4)*b(5,2))+b(3,4)*(b(4,2)*b(5,3)-b(4,3)*b(5,2))))
 a(2,2) = &
 (b(1,1)*(b(3,3)*(b(4,4)*b(5,5)-b(4,5)*b(5,4))-b(3,4)*(b(4,3)*b(5,5)-b(4,5)*b(5,3))+b(3,5)*(b(4,3)*b(5,4)-b(4,4)*b(5,3)))-b(1,3)* &
 (b(3,1)*(b(4,4)*b(5,5)-b(4,5)*b(5,4))-b(3,4)*(b(4,1)*b(5,5)-b(4,5)*b(5,1))+b(3,5)*(b(4,1)*b(5,4)-b(4,4)*b(5,1)))+b(1,4)* &
 (b(3,1)*(b(4,3)*b(5,5)-b(4,5)*b(5,3))-b(3,3)*(b(4,1)*b(5,5)-b(4,5)*b(5,1))+b(3,5)*(b(4,1)*b(5,3)-b(4,3)*b(5,1)))-b(1,5)* &
 (b(3,1)*(b(4,3)*b(5,4)-b(4,4)*b(5,3))-b(3,3)*(b(4,1)*b(5,4)-b(4,4)*b(5,1))+b(3,4)*(b(4,1)*b(5,3)-b(4,3)*b(5,1))))
 a(3,2) = &
 (-b(1,1)*(b(3,2)*(b(4,4)*b(5,5)-b(4,5)*b(5,4))-b(3,4)*(b(4,2)*b(5,5)-b(4,5)*b(5,2))+b(3,5)*(b(4,2)*b(5,4)-b(4,4)*b(5,2)))+b(1,2)* &
 (b(3,1)*(b(4,4)*b(5,5)-b(4,5)*b(5,4))-b(3,4)*(b(4,1)*b(5,5)-b(4,5)*b(5,1))+b(3,5)*(b(4,1)*b(5,4)-b(4,4)*b(5,1)))-b(1,4)* &
 (b(3,1)*(b(4,2)*b(5,5)-b(4,5)*b(5,2))-b(3,2)*(b(4,1)*b(5,5)-b(4,5)*b(5,1))+b(3,5)*(b(4,1)*b(5,2)-b(4,2)*b(5,1)))+b(1,5)* &
 (b(3,1)*(b(4,2)*b(5,4)-b(4,4)*b(5,2))-b(3,2)*(b(4,1)*b(5,4)-b(4,4)*b(5,1))+b(3,4)*(b(4,1)*b(5,2)-b(4,2)*b(5,1))))
 a(4,2) = &
 (b(1,1)*(b(3,2)*(b(4,3)*b(5,5)-b(4,5)*b(5,3))-b(3,3)*(b(4,2)*b(5,5)-b(4,5)*b(5,2))+b(3,5)*(b(4,2)*b(5,3)-b(4,3)*b(5,2)))-b(1,2)* &
 (b(3,1)*(b(4,3)*b(5,5)-b(4,5)*b(5,3))-b(3,3)*(b(4,1)*b(5,5)-b(4,5)*b(5,1))+b(3,5)*(b(4,1)*b(5,3)-b(4,3)*b(5,1)))+b(1,3)* &
 (b(3,1)*(b(4,2)*b(5,5)-b(4,5)*b(5,2))-b(3,2)*(b(4,1)*b(5,5)-b(4,5)*b(5,1))+b(3,5)*(b(4,1)*b(5,2)-b(4,2)*b(5,1)))-b(1,5)* &
 (b(3,1)*(b(4,2)*b(5,3)-b(4,3)*b(5,2))-b(3,2)*(b(4,1)*b(5,3)-b(4,3)*b(5,1))+b(3,3)*(b(4,1)*b(5,2)-b(4,2)*b(5,1))))
 a(5,2) = &
 (-b(1,1)*(b(3,2)*(b(4,3)*b(5,4)-b(4,4)*b(5,3))-b(3,3)*(b(4,2)*b(5,4)-b(4,4)*b(5,2))+b(3,4)*(b(4,2)*b(5,3)-b(4,3)*b(5,2)))+b(1,2)* &
 (b(3,1)*(b(4,3)*b(5,4)-b(4,4)*b(5,3))-b(3,3)*(b(4,1)*b(5,4)-b(4,4)*b(5,1))+b(3,4)*(b(4,1)*b(5,3)-b(4,3)*b(5,1)))-b(1,3)* &
 (b(3,1)*(b(4,2)*b(5,4)-b(4,4)*b(5,2))-b(3,2)*(b(4,1)*b(5,4)-b(4,4)*b(5,1))+b(3,4)*(b(4,1)*b(5,2)-b(4,2)*b(5,1)))+b(1,4)* &
 (b(3,1)*(b(4,2)*b(5,3)-b(4,3)*b(5,2))-b(3,2)*(b(4,1)*b(5,3)-b(4,3)*b(5,1))+b(3,3)*(b(4,1)*b(5,2)-b(4,2)*b(5,1))))

 a(1,3) = &
 (b(1,2)*(b(2,3)*(b(4,4)*b(5,5)-b(4,5)*b(5,4))-b(2,4)*(b(4,3)*b(5,5)-b(4,5)*b(5,3))+b(2,5)*(b(4,3)*b(5,4)-b(4,4)*b(5,3)))-b(1,3)* &
 (b(2,2)*(b(4,4)*b(5,5)-b(4,5)*b(5,4))-b(2,4)*(b(4,2)*b(5,5)-b(4,5)*b(5,2))+b(2,5)*(b(4,2)*b(5,4)-b(4,4)*b(5,2)))+b(1,4)* &
 (b(2,2)*(b(4,3)*b(5,5)-b(4,5)*b(5,3))-b(2,3)*(b(4,2)*b(5,5)-b(4,5)*b(5,2))+b(2,5)*(b(4,2)*b(5,3)-b(4,3)*b(5,2)))-b(1,5)* &
 (b(2,2)*(b(4,3)*b(5,4)-b(4,4)*b(5,3))-b(2,3)*(b(4,2)*b(5,4)-b(4,4)*b(5,2))+b(2,4)*(b(4,2)*b(5,3)-b(4,3)*b(5,2))))
 a(2,3) = &
 (-b(1,1)*(b(2,3)*(b(4,4)*b(5,5)-b(4,5)*b(5,4))-b(2,4)*(b(4,3)*b(5,5)-b(4,5)*b(5,3))+b(2,5)*(b(4,3)*b(5,4)-b(4,4)*b(5,3)))+b(1,3)* &
 (b(2,1)*(b(4,4)*b(5,5)-b(4,5)*b(5,4))-b(2,4)*(b(4,1)*b(5,5)-b(4,5)*b(5,1))+b(2,5)*(b(4,1)*b(5,4)-b(4,4)*b(5,1)))-b(1,4)* &
 (b(2,1)*(b(4,3)*b(5,5)-b(4,5)*b(5,3))-b(2,3)*(b(4,1)*b(5,5)-b(4,5)*b(5,1))+b(2,5)*(b(4,1)*b(5,3)-b(4,3)*b(5,1)))+b(1,5)* &
 (b(2,1)*(b(4,3)*b(5,4)-b(4,4)*b(5,3))-b(2,3)*(b(4,1)*b(5,4)-b(4,4)*b(5,1))+b(2,4)*(b(4,1)*b(5,3)-b(4,3)*b(5,1))))
 a(3,3) = &
 (b(1,1)*(b(2,2)*(b(4,4)*b(5,5)-b(4,5)*b(5,4))-b(2,4)*(b(4,2)*b(5,5)-b(4,5)*b(5,2))+b(2,5)*(b(4,2)*b(5,4)-b(4,4)*b(5,2)))-b(1,2)* &
 (b(2,1)*(b(4,4)*b(5,5)-b(4,5)*b(5,4))-b(2,4)*(b(4,1)*b(5,5)-b(4,5)*b(5,1))+b(2,5)*(b(4,1)*b(5,4)-b(4,4)*b(5,1)))+b(1,4)* &
 (b(2,1)*(b(4,2)*b(5,5)-b(4,5)*b(5,2))-b(2,2)*(b(4,1)*b(5,5)-b(4,5)*b(5,1))+b(2,5)*(b(4,1)*b(5,2)-b(4,2)*b(5,1)))-b(1,5)* &
 (b(2,1)*(b(4,2)*b(5,4)-b(4,4)*b(5,2))-b(2,2)*(b(4,1)*b(5,4)-b(4,4)*b(5,1))+b(2,4)*(b(4,1)*b(5,2)-b(4,2)*b(5,1))))
 a(4,3) = &
 (-b(1,1)*(b(2,2)*(b(4,3)*b(5,5)-b(4,5)*b(5,3))-b(2,3)*(b(4,2)*b(5,5)-b(4,5)*b(5,2))+b(2,5)*(b(4,2)*b(5,3)-b(4,3)*b(5,2)))+b(1,2)* &
 (b(2,1)*(b(4,3)*b(5,5)-b(4,5)*b(5,3))-b(2,3)*(b(4,1)*b(5,5)-b(4,5)*b(5,1))+b(2,5)*(b(4,1)*b(5,3)-b(4,3)*b(5,1)))-b(1,3)* &
 (b(2,1)*(b(4,2)*b(5,5)-b(4,5)*b(5,2))-b(2,2)*(b(4,1)*b(5,5)-b(4,5)*b(5,1))+b(2,5)*(b(4,1)*b(5,2)-b(4,2)*b(5,1)))+b(1,5)* &
 (b(2,1)*(b(4,2)*b(5,3)-b(4,3)*b(5,2))-b(2,2)*(b(4,1)*b(5,3)-b(4,3)*b(5,1))+b(2,3)*(b(4,1)*b(5,2)-b(4,2)*b(5,1))))
 a(5,3) = &
 (b(1,1)*(b(2,2)*(b(4,3)*b(5,4)-b(4,4)*b(5,3))-b(2,3)*(b(4,2)*b(5,4)-b(4,4)*b(5,2))+b(2,4)*(b(4,2)*b(5,3)-b(4,3)*b(5,2)))-b(1,2)* &
 (b(2,1)*(b(4,3)*b(5,4)-b(4,4)*b(5,3))-b(2,3)*(b(4,1)*b(5,4)-b(4,4)*b(5,1))+b(2,4)*(b(4,1)*b(5,3)-b(4,3)*b(5,1)))+b(1,3)* &
 (b(2,1)*(b(4,2)*b(5,4)-b(4,4)*b(5,2))-b(2,2)*(b(4,1)*b(5,4)-b(4,4)*b(5,1))+b(2,4)*(b(4,1)*b(5,2)-b(4,2)*b(5,1)))-b(1,4)* &
 (b(2,1)*(b(4,2)*b(5,3)-b(4,3)*b(5,2))-b(2,2)*(b(4,1)*b(5,3)-b(4,3)*b(5,1))+b(2,3)*(b(4,1)*b(5,2)-b(4,2)*b(5,1))))

 a(1,4) = &
 (-b(1,2)*(b(2,3)*(b(3,4)*b(5,5)-b(3,5)*b(5,4))-b(2,4)*(b(3,3)*b(5,5)-b(3,5)*b(5,3))+b(2,5)*(b(3,3)*b(5,4)-b(3,4)*b(5,3)))+b(1,3)* &
 (b(2,2)*(b(3,4)*b(5,5)-b(3,5)*b(5,4))-b(2,4)*(b(3,2)*b(5,5)-b(3,5)*b(5,2))+b(2,5)*(b(3,2)*b(5,4)-b(3,4)*b(5,2)))-b(1,4)* &
 (b(2,2)*(b(3,3)*b(5,5)-b(3,5)*b(5,3))-b(2,3)*(b(3,2)*b(5,5)-b(3,5)*b(5,2))+b(2,5)*(b(3,2)*b(5,3)-b(3,3)*b(5,2)))+b(1,5)* &
 (b(2,2)*(b(3,3)*b(5,4)-b(3,4)*b(5,3))-b(2,3)*(b(3,2)*b(5,4)-b(3,4)*b(5,2))+b(2,4)*(b(3,2)*b(5,3)-b(3,3)*b(5,2))))
 a(2,4) = &
 (b(1,1)*(b(2,3)*(b(3,4)*b(5,5)-b(3,5)*b(5,4))-b(2,4)*(b(3,3)*b(5,5)-b(3,5)*b(5,3))+b(2,5)*(b(3,3)*b(5,4)-b(3,4)*b(5,3)))-b(1,3)* &
 (b(2,1)*(b(3,4)*b(5,5)-b(3,5)*b(5,4))-b(2,4)*(b(3,1)*b(5,5)-b(3,5)*b(5,1))+b(2,5)*(b(3,1)*b(5,4)-b(3,4)*b(5,1)))+b(1,4)* &
 (b(2,1)*(b(3,3)*b(5,5)-b(3,5)*b(5,3))-b(2,3)*(b(3,1)*b(5,5)-b(3,5)*b(5,1))+b(2,5)*(b(3,1)*b(5,3)-b(3,3)*b(5,1)))-b(1,5)* &
 (b(2,1)*(b(3,3)*b(5,4)-b(3,4)*b(5,3))-b(2,3)*(b(3,1)*b(5,4)-b(3,4)*b(5,1))+b(2,4)*(b(3,1)*b(5,3)-b(3,3)*b(5,1))))
 a(3,4) = &
 (-b(1,1)*(b(2,2)*(b(3,4)*b(5,5)-b(3,5)*b(5,4))-b(2,4)*(b(3,2)*b(5,5)-b(3,5)*b(5,2))+b(2,5)*(b(3,2)*b(5,4)-b(3,4)*b(5,2)))+b(1,2)* &
 (b(2,1)*(b(3,4)*b(5,5)-b(3,5)*b(5,4))-b(2,4)*(b(3,1)*b(5,5)-b(3,5)*b(5,1))+b(2,5)*(b(3,1)*b(5,4)-b(3,4)*b(5,1)))-b(1,4)* &
 (b(2,1)*(b(3,2)*b(5,5)-b(3,5)*b(5,2))-b(2,2)*(b(3,1)*b(5,5)-b(3,5)*b(5,1))+b(2,5)*(b(3,1)*b(5,2)-b(3,2)*b(5,1)))+b(1,5)* &
 (b(2,1)*(b(3,2)*b(5,4)-b(3,4)*b(5,2))-b(2,2)*(b(3,1)*b(5,4)-b(3,4)*b(5,1))+b(2,4)*(b(3,1)*b(5,2)-b(3,2)*b(5,1))))
 a(4,4) = &
 (b(1,1)*(b(2,2)*(b(3,3)*b(5,5)-b(3,5)*b(5,3))-b(2,3)*(b(3,2)*b(5,5)-b(3,5)*b(5,2))+b(2,5)*(b(3,2)*b(5,3)-b(3,3)*b(5,2)))-b(1,2)* &
 (b(2,1)*(b(3,3)*b(5,5)-b(3,5)*b(5,3))-b(2,3)*(b(3,1)*b(5,5)-b(3,5)*b(5,1))+b(2,5)*(b(3,1)*b(5,3)-b(3,3)*b(5,1)))+b(1,3)* &
 (b(2,1)*(b(3,2)*b(5,5)-b(3,5)*b(5,2))-b(2,2)*(b(3,1)*b(5,5)-b(3,5)*b(5,1))+b(2,5)*(b(3,1)*b(5,2)-b(3,2)*b(5,1)))-b(1,5)* &
 (b(2,1)*(b(3,2)*b(5,3)-b(3,3)*b(5,2))-b(2,2)*(b(3,1)*b(5,3)-b(3,3)*b(5,1))+b(2,3)*(b(3,1)*b(5,2)-b(3,2)*b(5,1))))
 a(5,4) = &
 (-b(1,1)*(b(2,2)*(b(3,3)*b(5,4)-b(3,4)*b(5,3))-b(2,3)*(b(3,2)*b(5,4)-b(3,4)*b(5,2))+b(2,4)*(b(3,2)*b(5,3)-b(3,3)*b(5,2)))+b(1,2)* &
 (b(2,1)*(b(3,3)*b(5,4)-b(3,4)*b(5,3))-b(2,3)*(b(3,1)*b(5,4)-b(3,4)*b(5,1))+b(2,4)*(b(3,1)*b(5,3)-b(3,3)*b(5,1)))-b(1,3)* &
 (b(2,1)*(b(3,2)*b(5,4)-b(3,4)*b(5,2))-b(2,2)*(b(3,1)*b(5,4)-b(3,4)*b(5,1))+b(2,4)*(b(3,1)*b(5,2)-b(3,2)*b(5,1)))+b(1,4)* &
 (b(2,1)*(b(3,2)*b(5,3)-b(3,3)*b(5,2))-b(2,2)*(b(3,1)*b(5,3)-b(3,3)*b(5,1))+b(2,3)*(b(3,1)*b(5,2)-b(3,2)*b(5,1))))

 a(1,5) = &
 (b(1,2)*(b(2,3)*(b(3,4)*b(4,5)-b(3,5)*b(4,4))-b(2,4)*(b(3,3)*b(4,5)-b(3,5)*b(4,3))+b(2,5)*(b(3,3)*b(4,4)-b(3,4)*b(4,3)))-b(1,3)* &
 (b(2,2)*(b(3,4)*b(4,5)-b(3,5)*b(4,4))-b(2,4)*(b(3,2)*b(4,5)-b(3,5)*b(4,2))+b(2,5)*(b(3,2)*b(4,4)-b(3,4)*b(4,2)))+b(1,4)* &
 (b(2,2)*(b(3,3)*b(4,5)-b(3,5)*b(4,3))-b(2,3)*(b(3,2)*b(4,5)-b(3,5)*b(4,2))+b(2,5)*(b(3,2)*b(4,3)-b(3,3)*b(4,2)))-b(1,5)* &
 (b(2,2)*(b(3,3)*b(4,4)-b(3,4)*b(4,3))-b(2,3)*(b(3,2)*b(4,4)-b(3,4)*b(4,2))+b(2,4)*(b(3,2)*b(4,3)-b(3,3)*b(4,2))))
 a(2,5) = &
 (-b(1,1)*(b(2,3)*(b(3,4)*b(4,5)-b(3,5)*b(4,4))-b(2,4)*(b(3,3)*b(4,5)-b(3,5)*b(4,3))+b(2,5)*(b(3,3)*b(4,4)-b(3,4)*b(4,3)))+b(1,3)* &
 (b(2,1)*(b(3,4)*b(4,5)-b(3,5)*b(4,4))-b(2,4)*(b(3,1)*b(4,5)-b(3,5)*b(4,1))+b(2,5)*(b(3,1)*b(4,4)-b(3,4)*b(4,1)))-b(1,4)* &
 (b(2,1)*(b(3,3)*b(4,5)-b(3,5)*b(4,3))-b(2,3)*(b(3,1)*b(4,5)-b(3,5)*b(4,1))+b(2,5)*(b(3,1)*b(4,3)-b(3,3)*b(4,1)))+b(1,5)* &
 (b(2,1)*(b(3,3)*b(4,4)-b(3,4)*b(4,3))-b(2,3)*(b(3,1)*b(4,4)-b(3,4)*b(4,1))+b(2,4)*(b(3,1)*b(4,3)-b(3,3)*b(4,1))))
 a(3,5) = &
 (b(1,1)*(b(2,2)*(b(3,4)*b(4,5)-b(3,5)*b(4,4))-b(2,4)*(b(3,2)*b(4,5)-b(3,5)*b(4,2))+b(2,5)*(b(3,2)*b(4,4)-b(3,4)*b(4,2)))-b(1,2)* &
 (b(2,1)*(b(3,4)*b(4,5)-b(3,5)*b(4,4))-b(2,4)*(b(3,1)*b(4,5)-b(3,5)*b(4,1))+b(2,5)*(b(3,1)*b(4,4)-b(3,4)*b(4,1)))+b(1,4)* &
 (b(2,1)*(b(3,2)*b(4,5)-b(3,5)*b(4,2))-b(2,2)*(b(3,1)*b(4,5)-b(3,5)*b(4,1))+b(2,5)*(b(3,1)*b(4,2)-b(3,2)*b(4,1)))-b(1,5)* &
 (b(2,1)*(b(3,2)*b(4,4)-b(3,4)*b(4,2))-b(2,2)*(b(3,1)*b(4,4)-b(3,4)*b(4,1))+b(2,4)*(b(3,1)*b(4,2)-b(3,2)*b(4,1))))
 a(4,5) = &
 (-b(1,1)*(b(2,2)*(b(3,3)*b(4,5)-b(3,5)*b(4,3))-b(2,3)*(b(3,2)*b(4,5)-b(3,5)*b(4,2))+b(2,5)*(b(3,2)*b(4,3)-b(3,3)*b(4,2)))+b(1,2)* &
 (b(2,1)*(b(3,3)*b(4,5)-b(3,5)*b(4,3))-b(2,3)*(b(3,1)*b(4,5)-b(3,5)*b(4,1))+b(2,5)*(b(3,1)*b(4,3)-b(3,3)*b(4,1)))-b(1,3)* &
 (b(2,1)*(b(3,2)*b(4,5)-b(3,5)*b(4,2))-b(2,2)*(b(3,1)*b(4,5)-b(3,5)*b(4,1))+b(2,5)*(b(3,1)*b(4,2)-b(3,2)*b(4,1)))+b(1,5)* &
 (b(2,1)*(b(3,2)*b(4,3)-b(3,3)*b(4,2))-b(2,2)*(b(3,1)*b(4,3)-b(3,3)*b(4,1))+b(2,3)*(b(3,1)*b(4,2)-b(3,2)*b(4,1))))
 a(5,5) = &
 (b(1,1)*(b(2,2)*(b(3,3)*b(4,4)-b(3,4)*b(4,3))-b(2,3)*(b(3,2)*b(4,4)-b(3,4)*b(4,2))+b(2,4)*(b(3,2)*b(4,3)-b(3,3)*b(4,2)))-b(1,2)* &
 (b(2,1)*(b(3,3)*b(4,4)-b(3,4)*b(4,3))-b(2,3)*(b(3,1)*b(4,4)-b(3,4)*b(4,1))+b(2,4)*(b(3,1)*b(4,3)-b(3,3)*b(4,1)))+b(1,3)* &
 (b(2,1)*(b(3,2)*b(4,4)-b(3,4)*b(4,2))-b(2,2)*(b(3,1)*b(4,4)-b(3,4)*b(4,1))+b(2,4)*(b(3,1)*b(4,2)-b(3,2)*b(4,1)))-b(1,4)* &
 (b(2,1)*(b(3,2)*b(4,3)-b(3,3)*b(4,2))-b(2,2)*(b(3,1)*b(4,3)-b(3,3)*b(4,1))+b(2,3)*(b(3,1)*b(4,2)-b(3,2)*b(4,1))))

end subroutine invert5



end module module_fast_inv


