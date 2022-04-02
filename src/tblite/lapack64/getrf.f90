! This file is part of tblite.
! SPDX-Identifier: LGPL-3.0-or-later
!
! tblite is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! tblite is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with tblite.  If not, see <https://www.gnu.org/licenses/>.

module tblite_lapack_getrf
   use mctc_env, only : sp, dp, i8
   implicit none
   private

   public :: wrap_getrf


   !> Computes an LU factorization of a general M-by-N matrix A
   !> using partial pivoting with row interchanges.
   !>
   !> The factorization has the form
   !>    A = P * L * U
   !> where P is a permutation matrix, L is lower triangular with unit
   !> diagonal elements (lower trapezoidal if m > n), and U is upper
   !> triangular (upper trapezoidal if m < n).
   interface wrap_getrf
      module procedure :: wrap_sgetrf
      module procedure :: wrap_dgetrf
   end interface wrap_getrf


   !> Computes an LU factorization of a general M-by-N matrix A
   !> using partial pivoting with row interchanges.
   !>
   !> The factorization has the form
   !>    A = P * L * U
   !> where P is a permutation matrix, L is lower triangular with unit
   !> diagonal elements (lower trapezoidal if m > n), and U is upper
   !> triangular (upper trapezoidal if m < n).
   interface lapack_getrf
      pure subroutine sgetrf(m, n, a, lda, ipiv, info)
         import :: sp, i8
         integer(i8), intent(in) :: lda
         real(sp), intent(inout) :: a(lda, *)
         integer(i8), intent(out) :: ipiv(*)
         integer(i8), intent(out) :: info
         integer(i8), intent(in) :: m
         integer(i8), intent(in) :: n
      end subroutine sgetrf
      pure subroutine dgetrf(m, n, a, lda, ipiv, info)
         import :: dp, i8
         integer(i8), intent(in) :: lda
         real(dp), intent(inout) :: a(lda, *)
         integer(i8), intent(out) :: ipiv(*)
         integer(i8), intent(out) :: info
         integer(i8), intent(in) :: m
         integer(i8), intent(in) :: n
      end subroutine dgetrf
   end interface lapack_getrf

contains

subroutine wrap_sgetrf(amat, ipiv, info)
   real(sp), intent(inout) :: amat(:, :)
   integer(i8), intent(out) :: ipiv(:)
   integer(i8), intent(out) :: info
   integer(i8) :: m, n, lda
   lda = max(1, size(amat, 1))
   m = size(amat, 1, kind=i8)
   n = size(amat, 2, kind=i8)
   call lapack_getrf(m, n, amat, lda, ipiv, info)
end subroutine wrap_sgetrf


subroutine wrap_dgetrf(amat, ipiv, info)
   real(dp), intent(inout) :: amat(:, :)
   integer(i8), intent(out) :: ipiv(:)
   integer(i8), intent(out) :: info
   integer(i8) :: m, n, lda
   lda = max(1, size(amat, 1))
   m = size(amat, 1, kind=i8)
   n = size(amat, 2, kind=i8)
   call lapack_getrf(m, n, amat, lda, ipiv, info)
end subroutine wrap_dgetrf

end module tblite_lapack_getrf
