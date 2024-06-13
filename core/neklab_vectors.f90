      module neklab_vectors
      !---------------------------------------
      !-----     LightKrylov Imports     -----
      !---------------------------------------
      ! Default real kind.
      use LightKrylov, only: dp
      ! Abstract types for real- and complex-valued vectors.
      use LightKrylov, only: abstract_vector_rdp, abstract_vector_cdp
      
      implicit none
      include "SIZE"
      include "TOTAL"
      private
      
      public :: nek2vec, vec2nek
      public :: outpost_vec
      public :: nopcopy
      
      integer, parameter :: lv = lx1*ly1*lz1*lelv
      !! Local number of grid points for the velocity mesh.
      integer, parameter :: lp = lx2*ly2*lz2*lelv
      !! Local number of grid points for the pressure mesh.
      
      interface nek2vec
      module procedure nek2vec_prt
      module procedure nek2vec_std
      end interface
      
      interface vec2nek
      module procedure vec2nek_std
      module procedure vec2nek_prt
      end interface
      
      interface outpost_vec
      module procedure outpost_dvector
      end interface
      
      !----------------------------------------
      !-----     NEK REAL VECTOR TYPE     -----
      !----------------------------------------
      
      type, extends(abstract_vector_rdp), public :: nek_dvector
      !! Type definition for Nek5000 state-vector (real).
      real(kind=dp), dimension(lv) :: vx, vy, vz
      !! Velocity components of the state vector.
      real(kind=dp), dimension(lp) :: pr
      !! Pressure components of the state vector.
      real(kind=dp), dimension(lv, ldimt) :: theta
      !! Temperature + passive scalars components.
      contains
      private
      procedure, pass(self), public :: zero => nek_dzero
      !! Sets a vector to zero.
      procedure, pass(self), public :: rand => nek_drand
      !! Create a random vector.
      procedure, pass(self), public :: scal => nek_dscal
      !! Scale a vector such that \( \mathbf{x} = \alpha \mathbf{x}$ with     $                        \alpha \in \mathbb{R} \).
      procedure, pass(self), public :: axpby => nek_daxpby
      !! Add (in-place) two vectors such that \( \mathbf{x} = \alpha \mathbf{x} + \beta \mathbf{y} \) with \( \alpha \) and \( \beta \in \mathbb{R} \).
      procedure, pass(self), public :: dot => nek_ddot
      !! Compute the \( \ell_2 \) inner-product between two vectors.
      end type nek_dvector
      
      !-------------------------------------------
      !-----     NEK COMPLEX VECTOR TYPE     -----
      !-------------------------------------------
      
      type, extends(abstract_vector_cdp), public :: nek_zvector
      type(nek_dvector) :: re
      type(nek_dvector) :: im
      contains
      private
      procedure, pass(self), public :: zero => nek_zzero
      !! Sets a vector to zero.
      procedure, pass(self), public :: rand => nek_zrand
      !! Create a random vector.
      procedure, pass(self), public :: scal => nek_zscal
      !! Scale a vector such that \( \mathbf{x} = \alpha \mathbf{x}$ with     $                        \alpha \in \mathbb{R} \).
      procedure, pass(self), public :: axpby => nek_zaxpby
      !! Add (in-place) two vectors such that \( \mathbf{x} = \alpha \mathbf{x} + \beta \mathbf{y} \) with \( \alpha \) and \( \beta \in \mathbb{R} \).
      procedure, pass(self), public :: dot => nek_zdot
      !! Compute the \( \ell_2 \) inner-product between two vectors.
      end type nek_zvector
      
      contains
      
      !----------------------------------------------------------------------------
      !-----                                                                  -----
      !-----     DEFINITION OF THE TYPE BOUND PROCEDURES FOR REAL VECTORS     -----
      !-----                                                                  -----
      !----------------------------------------------------------------------------
      
      subroutine nek_dzero(self)
      class(nek_dvector), intent(inout) :: self
      !! Vector to be zeroed-out.
      call self%scal(0.0_dp)
      return
      end subroutine nek_dzero
      
      subroutine nek_drand(self, ifnorm)
      class(nek_dvector), intent(inout) :: self
      logical, optional, intent(in) :: ifnorm
      logical :: normalize
      integer :: i, n, ieg, iel
      real(kind=dp) :: xl(ldim), fcoeff(3), alpha
      
      n = nx1*ny1*nz1*nelv
      if (present(ifnorm)) then
      normalize = ifnorm
      else
      normalize = .false.
      end if
      
      do i = 1, n
      ieg = lglel(iel)
      xl(1) = xm1(i, 1, 1, 1)
      xl(2) = ym1(i, 1, 1, 1)
      if (if3d) xl(3) = zm1(i, 1, 1, 1)
      
      call random_number(fcoeff); fcoeff = fcoeff*1.0e4_dp
      self%vx(i) = self%vx(i) + mth_rand(i, 1, 1, 1, xl, fcoeff)
      
      call random_number(fcoeff); fcoeff = fcoeff*1.0e4_dp
      self%vy(i) = self%vy(i) + mth_rand(i, 1, 1, 1, xl, fcoeff)
      end do
      
      ! Face averaging.
      call opdssum(self%vx, self%vy, self%vz)
      call opcolv(self%vx, self%vy, self%vz, vmult)
      call dsavg(self%vx)
      call dsavg(self%vy)
      if (if3d) call dsavg(self%vz)
      call bcdirvc(self%vx, self%vy, self%vz, v1mask, v2mask, v3mask)
      
      if (normalize) then
      alpha = self%norm()
      call self%scal(1.0_dp/alpha)
      end if
      
      return
      end subroutine nek_drand
      
      subroutine nek_dscal(self, alpha)
      class(nek_dvector), intent(inout) :: self
      real(kind=dp), intent(in) :: alpha
      call dscal(lv, alpha, self%vx, 1)
      call dscal(lv, alpha, self%vy, 1)
      if (if3d) call dscal(lv, alpha, self%vz, 1)
      call dscal(lp, alpha, self%pr, 1)
      if (ifto) call dscal(lv, alpha, self%theta(:, 1), 1)
      return
      end subroutine nek_dscal
      
      subroutine nek_daxpby(self, alpha, vec, beta)
      class(nek_dvector), intent(inout) :: self
      real(kind=dp), intent(in) :: alpha
      class(abstract_vector_rdp), intent(in) :: vec
      real(kind=dp), intent(in) :: beta
      select type (vec)
      type is (nek_dvector)
      call daxpby(lv, beta, vec%vx, 1, alpha, self%vx, 1)
      call daxpby(lv, beta, vec%vy, 1, alpha, self%vy, 1)
      if (if3d) call daxpby(lv, beta, vec%vz, 1, alpha, self%vz, 1)
      call daxpby(lp, beta, vec%pr, 1, alpha, self%pr, 1)
      if (ifto) call daxpby(lv, beta, vec%theta(:, 1), 1, alpha, self%theta(:, 1), 1)
      end select
      return
      end subroutine nek_daxpby
      
      real(kind=dp) function nek_ddot(self, vec) result(alpha)
      class(nek_dvector), intent(in) :: self
      class(abstract_vector_rdp), intent(in) :: vec
      real(kind=dp), external :: glsc3
      integer :: i, n
      
      n = nx1*ny1*nz1*nelv
      select type (vec)
      type is (nek_dvector)
      ! Kinetic energy contribution.
      alpha = glsc3(self%vx, bm1, vec%vx, n) + glsc3(self%vy, bm1, vec%vy, n)
      if (if3d) then
      alpha = alpha + glsc3(self%vz, bm1, vec%vz)
      end if
      
      ! Thermal energy contribution.
      if (ifto) then
      alpha = alpha + glsc3(self%theta(:, 1), bm1, vec%theta(:, 1), n)
      end if
      
      ! Whatever contribution from additional scalars.
      if (ldimt > 1) then
      do i = 2, ldimt
      if (ifpsco(i - 1)) alpha = alpha + glsc3(self%theta(:, i), bm1, vec%theta(:, i), n)
      end do
      end if
      end select
      
      return
      end function nek_ddot
      !-------------------------------------------------------------------------------
      !-----                                                                     -----
      !-----     DEFINITION OF THE TYPE BOUND PROCEDURES FOR COMPLEX VECTORS     -----
      !-----                                                                     -----
      !-------------------------------------------------------------------------------
      
      subroutine nek_zzero(self)
      class(nek_zvector), intent(inout) :: self
      !! Vector to be zeroed-out.
      call self%scal(cmplx(0.0_dp, 0.0_dp, kind=dp))
      return
      end subroutine nek_zzero
      
      subroutine nek_zrand(self, ifnorm)
      class(nek_zvector), intent(inout) :: self
      logical, optional, intent(in) :: ifnorm
      logical :: normalize
      real(kind=dp) :: alpha
      
      ! Optional args.
      if (present(ifnorm)) then
      normalize = ifnorm
      else
      normalize = .false.
      end if
      
      ! Generate random vector.
      call self%re%rand(); call self%im%rand()
      
      ! Normalize if needed.
      if (normalize) then
      alpha = self%norm()
      call self%scal(cmplx(1.0_dp, 0.0_dp, kind=dp)/alpha)
      end if
      
      return
      end subroutine nek_zrand
      
      subroutine nek_zscal(self, alpha)
      class(nek_zvector), intent(inout) :: self
      complex(kind=dp), intent(in) :: alpha
      type(nek_zvector), allocatable :: wrk
      
      ! Scratch array.
      wrk = self
      ! Scale complex vector.
      call self%re%axpby(alpha%re, wrk%im, -alpha%im)
      call self%im%axpby(alpha%re, wrk%re, alpha%im)
      
      return
      end subroutine nek_zscal
      
      subroutine nek_zaxpby(self, alpha, vec, beta)
      class(nek_zvector), intent(inout) :: self
      complex(kind=dp), intent(in) :: alpha
      class(abstract_vector_cdp), intent(in) :: vec
      complex(kind=dp), intent(in) :: beta
      type(nek_zvector), allocatable :: wrk
      
      select type (vec)
      type is (nek_zvector)
      ! Scratch array.
      wrk = vec
      ! Scale vectors before addition.
      call self%scal(alpha); call wrk%scal(beta)
      ! Vector addition.
      call self%re%add(wrk%re); call self%im%add(wrk%im)
      end select
      return
      end subroutine nek_zaxpby
      
      complex(kind=dp) function nek_zdot(self, vec) result(alpha)
      class(nek_zvector), intent(in) :: self
      class(abstract_vector_cdp), intent(in) :: vec
      real(kind=dp) :: alpha_r, alpha_i
      select type (vec)
      type is (nek_zvector)
      alpha_r = self%re%dot(vec%re) + self%im%dot(vec%im)
      alpha_i = self%re%dot(vec%im) - self%im%dot(vec%re)
      alpha = cmplx(alpha_r, alpha_i, kind=dp)
      end select
      return
      end function nek_zdot
      
      !-----------------------------------------
      !-----                               -----
      !-----     NEK-RELATED UTILITIES     -----
      !-----                               -----
      !-----------------------------------------
      
      subroutine nek2vec_prt(vec, vx_, vy_, vz_, pr_, t_)
      include "SIZE"
      type(nek_dvector), intent(out) :: vec
      real(kind=dp), dimension(lx1*ly1*lz1*lelv, lpert), intent(in) :: vx_
      real(kind=dp), dimension(lx1*ly1*lz1*lelv, lpert), intent(in) :: vy_
      real(kind=dp), dimension(lx1*ly1*lz1*lelv, lpert), intent(in) :: vz_
      real(kind=dp), dimension(lx2*ly2*lz2*lelv, lpert), intent(in) :: pr_
      real(kind=dp), dimension(lx1*ly1*lz1*lelt, ldimt, lpert), intent(in) :: t_
      
      call nopcopy(vec%vx, vec%vy, vec%vz, vec%pr, vec%theta, vx_(:, 1), vy_(:, 1), vz_(:, 1), pr_(:, 1), t_(:, :, 1))
      
      return
      end subroutine nek2vec_prt
      
      subroutine nek2vec_std(vec, vx_, vy_, vz_, pr_, t_)
      include "SIZE"
      type(nek_dvector), intent(out) :: vec
      real(kind=dp), dimension(lx1, ly1, lz1, lelv), intent(in) :: vx_
      real(kind=dp), dimension(lx1, ly1, lz1, lelv), intent(in) :: vy_
      real(kind=dp), dimension(lx1, ly1, lz1, lelv), intent(in) :: vz_
      real(kind=dp), dimension(lx2, ly2, lz2, lelv), intent(in) :: pr_
      real(kind=dp), dimension(lx1, ly1, lz1, lelt, ldimt), intent(in) :: t_
      
      call nopcopy(vec%vx, vec%vy, vec%vz, vec%pr, vec%theta, vx_, vy_, vz_, pr_, t_)
      
      return
      end subroutine nek2vec_std
      
      subroutine vec2nek_std(vx_, vy_, vz_, pr_, t_, vec)
      include "SIZE"
      type(nek_dvector), intent(in) :: vec
      real(kind=dp), dimension(lx1, ly1, lz1, lelv), intent(out) :: vx_
      real(kind=dp), dimension(lx1, ly1, lz1, lelv), intent(out) :: vy_
      real(kind=dp), dimension(lx1, ly1, lz1, lelv), intent(out) :: vz_
      real(kind=dp), dimension(lx2, ly2, lz2, lelv), intent(out) :: pr_
      real(kind=dp), dimension(lx1, ly1, lz1, lelt, ldimt), intent(out) :: t_
      
      call nopcopy(vx_, vy_, vz_, pr_, t_, vec%vx, vec%vy, vec%vz, vec%pr, vec%theta)
      
      return
      end subroutine vec2nek_std
      
      subroutine vec2nek_prt(vx_, vy_, vz_, pr_, t_, vec)
      include "SIZE"
      type(nek_dvector), intent(in) :: vec
      real(kind=dp), dimension(lx1*ly1*lz1*lelv, lpert), intent(out) :: vx_
      real(kind=dp), dimension(lx1*ly1*lz1*lelv, lpert), intent(out) :: vy_
      real(kind=dp), dimension(lx1*ly1*lz1*lelv, lpert), intent(out) :: vz_
      real(kind=dp), dimension(lx2*ly2*lz2*lelv, lpert), intent(out) :: pr_
      real(kind=dp), dimension(lx1*ly1*lz1*lelt, ldimt, lpert), intent(out) :: t_
      
      call nopcopy(vx_(:, 1), vy_(:, 1), vz_(:, 1), pr_(:, 1), t_(:, :, 1), vec%vx, vec%vy, vec%vz, vec%pr, vec%theta)
      
      return
      end subroutine vec2nek_prt
      
      subroutine nopcopy(a1, a2, a3, a4, a5, b1, b2, b3, b4, b5)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      integer :: n, k
      real(kind=dp), intent(inout) :: a1(1), a2(1), a3(1), a4(1), a5(lx1*ly1*lz1*lelt, 1)
      real(kind=dp), intent(in) :: b1(1), b2(1), b3(1), b4(1), b5(lx1*ly1*lz1*lelt, 1)
      n = nx1*ny1*nz1*nelv
      call copy(a1, b1, n)
      call copy(a2, b2, n)
      if (if3D) call copy(a3, b3, n)
      if (ifpo) call copy(a4, b4, nx2*ny2*nz2*nelv)
      if (ifto) call copy(a5(1, 1), b5(1, 1), lx1*ly1*lz1*nelfld(2))
      if (ldimt > 1) then
      do k = 1, npscal
      if (ifpsco(k)) call copy(a5(1, k + 1), b5(1, k + 1), lx1*ly1*lz1*nelfld(k + 2))
      end do
      end if
      return
      end subroutine nopcopy
      
      real(kind=dp) function mth_rand(ix, iy, iz, ieg, xl, fcoeff) !generate random number
      include 'INPUT'           ! IF3D
      integer ix, iy, iz, ieg
      real(kind=dp) xl(LDIM), fcoeff(3)
      mth_rand = fcoeff(1)*(ieg + xl(1)*sin(xl(2))) + fcoeff(2)*ix*iy + fcoeff(3)*ix
      if (IF3D) mth_rand = fcoeff(1)*(ieg + xl(NDIM)*sin(mth_rand)) + fcoeff(2)*iz*ix + fcoeff(3)*iz
      mth_rand = 1.0e3_dp*sin(mth_rand)
      mth_rand = 1.0e3_dp*sin(mth_rand)
      mth_rand = cos(mth_rand)
      return
      end function mth_rand
      
      subroutine outpost_dvector(vec, prefix)
      type(nek_dvector), intent(in) :: vec
      character(len=3), intent(in) :: prefix
      call outpost(vec%vx, vec%vy, vec%vz, vec%pr, vec%theta, prefix)
      return
      end subroutine outpost_dvector
      
      end module neklab_vectors
