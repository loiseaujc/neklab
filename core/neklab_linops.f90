      module neklab_linops
      !---------------------------------------
      !-----     LightKrylov Imports     -----
      !---------------------------------------
      ! Default real kind.
         use LightKrylov, only: dp
      ! Abstract types for real-valued linear operators and vectors.
         use LightKrylov, only: abstract_linop_rdp, abstract_vector_rdp
         use LightKrylov, only: Id_rdp, axpby_linop_rdp
         ! Abstract types for complex-valued linear operator and vectors.
         use LightKrylov, only: abstract_linop_cdp, abstract_vector_cdp
         use LightKrylov, only: Id_cdp, axpby_linop_cdp
         ! Iterative linear solver for the resolvent computation.
         use LightKrylov, only: gmres, gmres_dp_opts
         use LightKrylov, only: rtol_dp, atol_dp
      ! Extensions of the abstract vector types to nek data format.
         use neklab_vectors
         implicit none
         include "SIZE"
         include "TOTAL"
         include "ADJOINT"
         private

         complex(kind=dp), parameter :: one_cdp = cmplx(1.0_dp, 0.0_dp, kind=dp)
         complex(kind=dp), parameter :: im_dp   = cmplx(0.0_dp, 1.0_dp, kind=dp)
         real(kind=dp), parameter :: pi_ = 4.0_dp * atan(1.0_dp)
         procedure(neklab_forcing_interface), pointer, public :: neklab_forcing => dummy_neklab_forcing

         interface
            subroutine neklab_forcing_interface(ffx, ffy, ffz, ix, iy, iz, ieg)
               import dp
               real(kind=dp), intent(inout) :: ffx, ffy, ffz
               integer      , intent(in)    :: ix, iy, iz, ieg
            end subroutine neklab_forcing_interface
         end interface

      !------------------------------------------
      !-----     EXPONENTIAL PROPAGATOR     -----
      !------------------------------------------
      
         type, extends(abstract_linop_rdp), public :: exptA_linop
            real(kind=dp) :: tau
            type(nek_dvector) :: baseflow
         contains
            private
            procedure, pass(self), public :: matvec => exptA_matvec
            procedure, pass(self), public :: rmatvec => exptA_rmatvec
            procedure, pass(self), public :: init => init_exptA
         end type exptA_linop

      !--------------------------------------
      !-----     RESOLVENT OPERATOR     -----
      !--------------------------------------

         type, extends(abstract_linop_cdp) :: zexptA_linop
            real(kind=dp) :: tau
            type(nek_dvector) :: baseflow
         contains
            private
            procedure, pass(self), public :: matvec => zexptA_matvec
            procedure, pass(self), public :: rmatvec => zexptA_rmatvec
            procedure, pass(self), public :: init => init_zexptA
         end type zexptA_linop

         type, extends(abstract_linop_cdp), public :: resolvent_linop
            real(kind=dp) :: omega
            type(nek_dvector) :: baseflow
         contains
            private
            procedure, pass(self), public :: matvec => resolvent_matvec
            procedure, pass(self), public :: rmatvec => resolvent_rmatvec
         end type resolvent_linop

         type(nek_zvector) :: resolvent_forcing

      contains
      
      !--------------------------------------------------
      !-----     TYPE-BOUND PROCEDURE FOR exptA     -----
      !--------------------------------------------------
      
         subroutine init_exptA(self)
            class(exptA_linop), intent(in) :: self
            call nekgsync()
      
      ! Setup Nek5000 logical flags for perturbative solver.
            ifpert = .true.; ifbase = .false.
            call bcast(ifpert, lsize); call bcast(ifbase, lsize)
      
      ! Force single perturbation mode.
            if (param(31) > 1) then
               if (nid == 0) write (6, *) "neklab does not support (yet) npert > 1."
               call nek_end()
            else
               param(31) = 1; npert = 1
            end if
      
      ! Deactivate OIFS.
            if (ifchar) then
            if (nid == 0) then
               write (6, *) "WARNING : OIFS is not available for linearized solver."
               write (6, *) "          Turning it off."
            end if
            ifchar = .false.
            end if
      
      ! Force CFL to 0.5
            if (param(26) > 0.5) then
            if (nid == 0) then
               write (6, *) "WARNING : Target CFL is larger than 0.5"
               write (6, *) "          Forcing it to 0.5"
            end if
            param(26) = 0.5_dp
            end if
      
      ! Compute appropriate step size.
            param(10) = self%tau
            call compute_cfl(ctarg, self%baseflow%vx, self%baseflow%vy, self%baseflow%vz, 1.0_dp)
            dt = param(26)/ctarg; nsteps = ceiling(self%tau/dt)
            dt = self%tau/nsteps; param(12) = dt
            call compute_cfl(ctarg, self%baseflow%vx, self%baseflow%vy, self%baseflow%vz, dt)
            if (nid == 0) write (6, *) "INFO : Current CFL and target CFL :", ctarg, param(26)
            lastep = 0; fintim = nsteps*dt
      
      ! Force contant time step.
            param(12) = -abs(param(12))
      
      ! Broadcast parameters.
            call bcast(param, 200*wdsize)
      
            return
         end subroutine init_exptA
      
         subroutine exptA_matvec(self, vec_in, vec_out)
            class(exptA_linop), intent(in) :: self
            class(abstract_vector_rdp), intent(in) :: vec_in
            class(abstract_vector_rdp), intent(out) :: vec_out
      
      ! Nek-related setup.
            ifadj = .false.; lastep = 0; fintim = param(10)
            call bcast(ifadj, lsize)
      
      ! Force the baseflow field.
            call vec2nek(vx, vy, vz, pr, t, self%baseflow)
      
      ! Sets the initial condition for Nek5000's linearized solver.
            select type (vec_in)
            type is (nek_dvector)
               call vec2nek(vxp, vyp, vzp, prp, tp, vec_in)
            end select
      
      ! Integrate the equations forward in time.
            time = 0.0_dp ; lastep = 0 ; fintim = param(10)
            do istep = 1, nsteps
               call nek_advance()
            end do
      
      ! Copy the final solution to vector.
            select type (vec_out)
            type is (nek_dvector)
               call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)
            end select
      
            return
         end subroutine exptA_matvec
      
         subroutine exptA_rmatvec(self, vec_in, vec_out)
            class(exptA_linop), intent(in) :: self
            class(abstract_vector_rdp), intent(in) :: vec_in
            class(abstract_vector_rdp), intent(out) :: vec_out
      ! Nek-related setup.
            ifadj = .true.; lastep = 0; fintim = param(10)
            call bcast(ifadj, lsize)
      
      ! Force the baseflow field.
            call vec2nek(vx, vy, vz, pr, t, self%baseflow)
      
      ! Sets the initial condition for Nek5000's linearized solver.
            select type (vec_in)
            type is (nek_dvector)
               call vec2nek(vxp, vyp, vzp, prp, tp, vec_in)
            end select
      
      ! Integrate the equations forward in time.
            time = 0.0_dp ; lastep = 0 ; fintim = param(10)
            do istep = 1, nsteps
               call nek_advance()
            end do
      
      ! Copy the final solution to vector.
            select type (vec_out)
            type is (nek_dvector)
               call nek2vec(vec_out, vxp, vyp, vzp, prp, tp)
            end select
      
            ifadj = .false.
            return
         end subroutine exptA_rmatvec
 
      !---------------------------------------------------
      !-----     TYPE-BOUND PROCEDURE FOR zexptA     -----
      !---------------------------------------------------
      
         subroutine init_zexptA(self)
            class(zexptA_linop), intent(in) :: self
            call nekgsync()
      
      ! Setup Nek5000 logical flags for perturbative solver.
            ifpert = .true.; ifbase = .false.
            call bcast(ifpert, lsize); call bcast(ifbase, lsize)
      
      ! Force single perturbation mode.
            if (param(31) > 1) then
               if (nid == 0) write (6, *) "neklab does not support (yet) npert > 1."
               call nek_end()
            else
               param(31) = 1; npert = 1
            end if
      
      ! Deactivate OIFS.
            if (ifchar) then
            if (nid == 0) then
               write (6, *) "WARNING : OIFS is not available for linearized solver."
               write (6, *) "          Turning it off."
            end if
            ifchar = .false.
            end if
      
      ! Force CFL to 0.5
            if (param(26) > 0.5) then
            if (nid == 0) then
               write (6, *) "WARNING : Target CFL is larger than 0.5"
               write (6, *) "          Forcing it to 0.5"
            end if
            param(26) = 0.5_dp
            end if
      
      ! Compute appropriate step size.
            param(10) = self%tau
            call compute_cfl(ctarg, self%baseflow%vx, self%baseflow%vy, self%baseflow%vz, 1.0_dp)
            dt = param(26)/ctarg; nsteps = ceiling(self%tau/dt)
            dt = self%tau/nsteps; param(12) = dt
            call compute_cfl(ctarg, self%baseflow%vx, self%baseflow%vy, self%baseflow%vz, dt)
            if (nid == 0) write (6, *) "INFO : Current CFL and target CFL :", ctarg, param(26)
            lastep = 0; fintim = nsteps*dt
      
      ! Force contant time step.
            param(12) = -abs(param(12))
      
      ! Broadcast parameters.
            call bcast(param, 200*wdsize)
      
            return
         end subroutine init_zexptA
      
         subroutine zexptA_matvec(self, vec_in, vec_out)
            class(zexptA_linop), intent(in) :: self
            class(abstract_vector_cdp), intent(in) :: vec_in
            class(abstract_vector_cdp), intent(out) :: vec_out
      
      ! Nek-related setup.
            ifadj = .false.; lastep = 0; fintim = param(10)
            call bcast(ifadj, lsize)
      
      ! Force the baseflow field.
            call vec2nek(vx, vy, vz, pr, t, self%baseflow)
     
            !-----------------------------
            !-----     REAL PART     -----
            !-----------------------------

      ! Sets the initial condition for Nek5000's linearized solver.
            select type (vec_in)
            type is (nek_zvector)
               call vec2nek(vxp, vyp, vzp, prp, tp, vec_in%re)
            end select
      
      ! Integrate the equations forward in time.
            time = 0.0_dp ; lastep = 0 ; fintim = param(10)
            do istep = 1, nsteps
               call nek_advance()
            end do
            call outpost(vxp, vyp, vzp, prp, tp, "prt")
      
      ! Copy the final solution to vector.
            select type (vec_out)
            type is (nek_zvector)
               call nek2vec(vec_out%re, vxp, vyp, vzp, prp, tp)
            end select

            !-----------------------------
            !-----     IMAG PART     -----
            !-----------------------------

      ! Sets the initial condition for Nek5000's linearized solver.
            select type (vec_in)
            type is (nek_zvector)
               call vec2nek(vxp, vyp, vzp, prp, tp, vec_in%im)
            end select
      
      ! Integrate the equations forward in time.
            time = 0.0_dp ; lastep = 0 ; fintim = param(10)
            do istep = 1, nsteps
               call nek_advance()
            end do
            call outpost(vxp, vyp, vzp, prp, tp, "adj")
      
      ! Copy the final solution to vector.
            select type (vec_out)
            type is (nek_zvector)
               call nek2vec(vec_out%im, vxp, vyp, vzp, prp, tp)
            end select
      
            return
         end subroutine zexptA_matvec
      
         subroutine zexptA_rmatvec(self, vec_in, vec_out)
            class(zexptA_linop), intent(in) :: self
            class(abstract_vector_cdp), intent(in) :: vec_in
            class(abstract_vector_cdp), intent(out) :: vec_out
      ! Nek-related setup.
            ifadj = .true.; lastep = 0; fintim = param(10)
            call bcast(ifadj, lsize)
      
      ! Force the baseflow field.
            call vec2nek(vx, vy, vz, pr, t, self%baseflow)
     
            !------------------------------
            !-----     REAL PART      -----
            !------------------------------

      ! Sets the initial condition for Nek5000's linearized solver.
            select type (vec_in)
            type is (nek_zvector)
               call vec2nek(vxp, vyp, vzp, prp, tp, vec_in%re)
            end select
      
      ! Integrate the equations forward in time.
            time = 0.0_dp ; lastep = 0 ; fintim = param(10)
            do istep = 1, nsteps
               call nek_advance()
            end do
      
      ! Copy the final solution to vector.
            select type (vec_out)
            type is (nek_zvector)
               call nek2vec(vec_out%re, vxp, vyp, vzp, prp, tp)
            end select

            !------------------------------
            !-----     IMAG PART      -----
            !------------------------------

      ! Sets the initial condition for Nek5000's linearized solver.
            select type (vec_in)
            type is (nek_zvector)
               call vec2nek(vxp, vyp, vzp, prp, tp, vec_in%im)
            end select
      
      ! Integrate the equations forward in time.
            time = 0.0_dp ; lastep = 0 ; fintim = param(10)
            do istep = 1, nsteps
               call nek_advance()
            end do
      
      ! Copy the final solution to vector.
            select type (vec_out)
            type is (nek_zvector)
               call nek2vec(vec_out%im, vxp, vyp, vzp, prp, tp)
            end select
      
            ifadj = .false.
            return
         end subroutine zexptA_rmatvec

         !------------------------------------------------------------
         !-----     TYPE-BOUND PROCEDURE FOR RESOLVENT LINOP     -----
         !------------------------------------------------------------

          subroutine resolvent_matvec(self, vec_in, vec_out)
            class(resolvent_linop), intent(in) :: self
            class(abstract_vector_cdp), intent(in) :: vec_in
            class(abstract_vector_cdp), intent(out) :: vec_out

            !----- Internal variables -----
            ! Complex-valued exponential propagator.
            class(exptA_linop), allocatable :: exptA
            class(Id_rdp), allocatable :: Id
            type(nek_zvector) :: b
            real(kind=dp) :: tau

            !------------------------------------------
            !-----     Setting-up the problem     -----
            !------------------------------------------

               ! Integration time.
               if (self%omega == 0.0_dp) then
                  tau = 1.0_dp
               else
                  tau = 2*pi_ / abs(self%omega)
               endif

               ! Initialize required operators.
               exptA = exptA_linop(tau, self%baseflow) ; call exptA%init()
               Id = Id_rdp()

            !-------------------------------------
            !-----     Computing the rhs     -----
            !-------------------------------------

            ! Initialize rhs vector.
            call b%zero()

            ! Set the forcing function.
            neklab_forcing => neklab_resolvent_forcing

            ! ----- Real part -----

            ! Set initial condition to zero.
            call opzero(vxp, vyp, vzp)

            ! Time integration.
            time = 0.0_dp ; lastep = 0 ; fintim = param(10)
            do istep = 1, nsteps
               ! Zero-out the resolvent forcing.
               call resolvent_forcing%zero()
               ! Update the resolvent forcing for the corresponding time.
               call resolvent_forcing%add(vec_in)
               call resolvent_forcing%scal(exp(im_dp*self%omega*time))
               ! Nek5000 time-step.
               call nek_advance()
            enddo

            ! Copy the real-part of the response.
            call nek2vec(b%re, vxp, vyp, vzp, prp, tp)
            call outpost_vec(b%re, "bre")

            ! ! ----- Imaginary part -----
            !
            ! ! Set initial condition to zero.
            ! call opzero(vxp, vyp, vzp)
            !
            ! ! Time integration.
            ! time = 0.0_dp ; lastep = 0 ; fintim = param(10)
            ! do istep = 1, nsteps
            !    ! Zero-out the resolvent forcing.
            !    call resolvent_forcing%zero()
            !    ! Update the resolvent forcing for the corresponding time.
            !    call resolvent_forcing%add(vec_in)
            !    call resolvent_forcing%scal(exp(im_dp*self%omega*time))
            !    resolvent_forcing%re = resolvent_forcing%im
            !    ! Nek5000 time-step.
            !    call nek_advance()
            ! enddo
            !
            ! ! Copy the imaginary-part of the response.
            ! call nek2vec(b%im, vxp, vyp, vzp, prp, tp)
            ! call outpost_vec(b%im, "bim")

            !------------------------------------------------
            !-----     Apply the resolvent operator     -----
            !------------------------------------------------

            neklab_forcing => dummy_neklab_forcing

            block
               ! GMRES-required variables.
               type(gmres_dp_opts) :: opts
               integer :: info
               class(axpby_linop_rdp), allocatable :: S

               ! Initialize operator S = I - exptA.
               S = axpby_linop_rdp(Id, exptA, 1.0_dp, -1.0_dp)
               ! GMRES default options.
               opts = gmres_dp_opts(kdim=64, verbose=.false., atol=1.0e-8_dp, rtol=1.0e-6_dp, maxiter=10)
               ! GMRES solver.
               call vec_out%zero()
               select type(vec_out)
               type is(nek_zvector)
                  call gmres(S, b%re, vec_out%re, info, options=opts)
               end select
            end block

            ! ----- Imaginary part -----

            ! Set initial condition to zero.
            select type(vec_out)
            type is(nek_zvector)
               call vec2nek(vxp, vyp, vzp, prp, tp, vec_out%re)
            end select
            exptA%tau = exptA%tau / 4.0_dp ; call exptA%init()

            ! Time integration.
            time = 0.0_dp ; lastep = 0 ; fintim = param(10)
            do istep = 1, nsteps
               ! Zero-out the resolvent forcing.
               call resolvent_forcing%zero()
               ! Update the resolvent forcing for the corresponding time.
               call resolvent_forcing%add(vec_in)
               call resolvent_forcing%scal(exp(im_dp*self%omega*time))
               resolvent_forcing%re = resolvent_forcing%im
               ! Nek5000 time-step.
               call nek_advance()
            enddo

            ! Copy the imaginary-part of the response.
            select type(vec_out)
            type is(nek_zvector)
               call nek2vec(vec_out%im, vxp, vyp, vzp, prp, tp)
            end select

            return
         end subroutine resolvent_matvec
 
         subroutine resolvent_rmatvec(self, vec_in, vec_out)
            class(resolvent_linop), intent(in) :: self
            class(abstract_vector_cdp), intent(in) :: vec_in
            class(abstract_vector_cdp), intent(out) :: vec_out
            return
         end subroutine resolvent_rmatvec

         !-----------------------------
         !-----     UTILITIES     -----
         !-----------------------------

         subroutine dummy_neklab_forcing(ffx, ffy, ffz, ix, iy, iz, ieg)
            real(kind=dp), intent(inout) :: ffx, ffy, ffz
            integer      , intent(in)    :: ix, iy, iz, ieg
            return
         end subroutine dummy_neklab_forcing

         subroutine neklab_resolvent_forcing(ffx, ffy, ffz, ix, iy, iz, ieg)
            real(kind=dp), intent(inout) :: ffx, ffy, ffz
            integer      , intent(in)    :: ix, iy, iz, ieg

            ! Internal variables.
            integer :: iel, ip
!            integer, external :: gllel

            ! Local index.
            iel = gllel(ieg)
            ip  = ix + nx1*(iy - 1 + ny1*(iz -1 + nz1*(iel - 1)))
            ! Update forcing.
            associate(force => resolvent_forcing%re)
               ffx = ffx + force%vx(ip)
               ffy = ffy + force%vy(ip)
               if (if3d) ffz = ffz + force%vz(ip)
            end associate
            return
         end subroutine neklab_resolvent_forcing
      end module neklab_linops
