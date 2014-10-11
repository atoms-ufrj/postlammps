module mConstants

!   This file is part of Postlammps.
!
!    Postlammps is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    Postlammps is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with Postlammps.  If not, see <http://www.gnu.org/licenses/>.

integer,      parameter :: rb = 8
integer,      parameter :: sl = 256
character(3), parameter :: csl = "256"
real(rb), parameter :: zero  = 0.0_rb, &
                       half  = 0.5_rb, &
                       one   = 1.0_rb, &
                       two   = 2.0_rb, &
                       three = 3.0_rb
end module mConstants
