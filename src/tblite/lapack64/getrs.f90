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

module tblite_lapack_getrs
   use mctc_env, only : sp, dp, i8
   implicit none
   private

   public :: wrap_getrs


   !> Solves a system of linear equations
   !>    A * X = B  or  A**T * X = B
   !> with a general N-by-N matrix A using the LU factorization computed
   !> by ?GETRF.
   interface wrap_getrs
      module procedure :: wrap_sgetrs
      module procedure :: wrap_dgetrs
   end interface wrap_getrs


   !> Solves a system of linear equations
   !>    A * X = B  or  A**T * X = B
   !> with a general N-by-N matrix A using the LU factorization computed
   !> by ?GETRF.
   interface lapack_getrs
      pure subroutine sgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info)
         import :: sp, i8
         integer(i8), intent(in) :: lda
         integer(i8), intent(in) :: ldb
         real(sp), intent(in) :: a(lda, *)
         integer(i8), intent(in) :: ipiv(*)
         real(sp), intent(inout) :: b(ldb, *)
         character(len=1), intent(in) :: trans
         integer(i8), intent(out) :: info
         integer(i8), intent(in) :: n
         integer(i8), intent(in) :: nrhs
      end subroutine sgetrs
      pure subroutine dgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info)
         import :: dp, i8
         integer(i8), intent(in) :: lda
         integer(i8), intent(in) :: ldb
         real(dp), intent(in) :: a(lda, *)
         integer(i8), intent(in) :: ipiv(*)
         real(dp), intent(inout) :: b(ldb, *)
         character(len=1), intent(in) :: trans
         integer(i8), intent(out) :: info
         integer(i8), intent(in) :: n
         integer(i8), intent(in) :: nrhs
      end subroutine dgetrs
   end interface lapack_getrs

contains

subroutine wrap_sgetrs(amat, bmat, ipiv, info, trans)
   real(sp), intent(in) :: amat(:, :)
   real(sp), intent(inout) :: bmat(:, :)
   integer(i8), intent(in) :: ipiv(:)
   integer(i8), intent(out) :: info
   character(len=1), intent(in), optional :: trans
   character(len=1) :: tra
   integer(i8) :: n, nrhs, lda, ldb
   if (present(trans)) then
      tra = trans
   else
      tra = 'n'
   end if
   lda = max(1, size(amat, 1))
   ldb = max(1, size(bmat, 1))
   n = size(amat, 2, kind=i8)
   nrhs = size(bmat, 2, kind=i8)
   call lapack_getrs(tra, n, nrhs, amat, lda, ipiv, bmat, ldb, info)
end subroutine wrap_sgetrs


subroutine wrap_dgetrs(amat, bmat, ipiv, info, trans)
   real(dp), intent(in) :: amat(:, :)
   real(dp), intent(inout) :: bmat(:, :)
   integer(i8), intent(in) :: ipiv(:)
   integer(i8), intent(out) :: info
   character(len=1), intent(in), optional :: trans
   character(len=1) :: tra
   integer(i8) :: n, nrhs, lda, ldb
   if (present(trans)) then
      tra = trans
   else
      tra = 'n'
   end if
   lda = max(1, size(amat, 1))
   ldb = max(1, size(bmat, 1))
   n = size(amat, 2, kind=i8)
   nrhs = size(bmat, 2, kind=i8)
   call lapack_getrs(tra, n, nrhs, amat, lda, ipiv, bmat, ldb, info)
end subroutine wrap_dgetrs

end module tblite_lapack_getrs
