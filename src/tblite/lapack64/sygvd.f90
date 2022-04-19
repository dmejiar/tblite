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

module tblite_lapack_sygvd
   use mctc_env, only : sp, dp, error_type, fatal_error, i8, i4
   use tblite_output_format, only : format_string
   use tblite_scf_solver, only : solver_type
   implicit none
   private

   public :: sygvd_solver


   interface lapack_sygvd
      pure subroutine ssygvd(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, &
            & iwork, liwork, info)
         import :: sp, i8
         integer(i8), intent(in) :: lda
         integer(i8), intent(in) :: ldb
         real(sp), intent(inout) :: a(lda, *)
         real(sp), intent(inout) :: b(ldb, *)
         real(sp), intent(out) :: w(*)
         integer(i8), intent(in) :: itype
         character(len=1), intent(in) :: jobz
         character(len=1), intent(in) :: uplo
         integer(i8), intent(out) :: info
         integer(i8), intent(in) :: n
         real(sp), intent(inout) :: work(*)
         integer(i8), intent(in) :: lwork
         integer(i8), intent(inout) :: iwork(*)
         integer(i8), intent(in) :: liwork
      end subroutine ssygvd
      pure subroutine dsygvd(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, &
            & iwork, liwork, info)
         import :: dp, i8
         integer(i8), intent(in) :: lda
         integer(i8), intent(in) :: ldb
         real(dp), intent(inout) :: a(lda, *)
         real(dp), intent(inout) :: b(ldb, *)
         real(dp), intent(out) :: w(*)
         integer(i8), intent(in) :: itype
         character(len=1), intent(in) :: jobz
         character(len=1), intent(in) :: uplo
         integer(i8), intent(out) :: info
         integer(i8), intent(in) :: n
         real(dp), intent(inout) :: work(*)
         integer(i8), intent(in) :: lwork
         integer(i8), intent(inout) :: iwork(*)
         integer(i8), intent(in) :: liwork
      end subroutine dsygvd
   end interface lapack_sygvd


   type, extends(solver_type) :: sygvd_solver
      private
      integer(i8) :: n = 0
      integer(i8), allocatable :: iwork(:)
      real(sp), allocatable :: swork(:)
      real(sp), allocatable :: sbmat(:, :)
      real(dp), allocatable :: dwork(:)
      real(dp), allocatable :: dbmat(:, :)
   contains
      procedure :: solve_sp
      procedure :: solve_dp
   end type sygvd_solver

contains

subroutine solve_sp(self, hmat, smat, eval, error)
   class(sygvd_solver), intent(inout) :: self
   real(sp), contiguous, intent(inout) :: hmat(:, :)
   real(sp), contiguous, intent(in) :: smat(:, :)
   real(sp), contiguous, intent(inout) :: eval(:)
   type(error_type), allocatable, intent(out) :: error
   integer(i8) :: info, lswork, liwork

   if (self%n == 0) then
      self%n = size(hmat, 1)
   end if
   if (.not.allocated(self%swork)) then
      allocate(self%swork(1 + 6*self%n + 2*self%n**2))
   end if
   if (.not.allocated(self%iwork)) then
      allocate(self%iwork(3 + 5*self%n))
   end if
   self%sbmat = smat
   lswork = size(self%swork)
   liwork = size(self%iwork)

   call lapack_sygvd(1_i8, 'v', 'u', self%n, hmat, self%n, self%sbmat, self%n, eval, &
      & self%swork, lswork, self%iwork, liwork, info)

   call handle_info(error, info)

end subroutine solve_sp

subroutine solve_dp(self, hmat, smat, eval, error)
   class(sygvd_solver), intent(inout) :: self
   real(dp), contiguous, intent(inout) :: hmat(:, :)
   real(dp), contiguous, intent(in) :: smat(:, :)
   real(dp), contiguous, intent(inout) :: eval(:)
   type(error_type), allocatable, intent(out) :: error
   integer(i8) :: info, ldwork, liwork

   if (self%n == 0) then
      self%n = size(hmat, 1, kind=i8)
   end if
   if (.not.allocated(self%dwork)) then
      allocate(self%dwork(1 + 6*self%n + 2*self%n**2))
   end if
   if (.not.allocated(self%iwork)) then
      allocate(self%iwork(3 + 5*self%n))
   end if
   self%dbmat = smat
   ldwork = size(self%dwork)
   liwork = size(self%iwork)

   call lapack_sygvd(1_i8, 'v', 'u', self%n, hmat, self%n, self%dbmat, self%n, eval, &
      & self%dwork, ldwork, self%iwork, liwork, info)

   call handle_info(error, info)

end subroutine solve_dp

subroutine handle_info(error, info)
   type(error_type), allocatable, intent(out) :: error
   integer(i8), intent(in) :: info

   if (info /= 0) then
      call fatal_error(error, "(sygvd) failed to solve eigenvalue problem.&
         & info="//format_string(int(info), '(i0)'))
   end if
end subroutine handle_info

end module tblite_lapack_sygvd
