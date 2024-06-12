      module neklab
      
      ! --> Abstract Krylov methods.
         use LightKrylov
      ! --> Definition of the abstract vectors in the Nek framework.
         use neklab_vectors
      ! --> Definitions of the abstract linops in the Nek framework.
         use neklab_linops
      ! -->
         use neklab_analysis
      
         private
      
      !------------------------------------------
      !-----     LIGHTKRYLOV RE-EXPORTS     -----
      !------------------------------------------
      
      ! Global variables.
         public :: dp
      
      ! Krylov factorizations.
         public :: initialize_krylov_subspace
         public :: arnoldi_factorization, krylov_schur_restart
         public :: lanczos_bidiagonalization, lanczos_tridiagonalization
      
      ! Matrix factorizations.
         public :: eigs, eighs, svds
         public :: save_eigenspectrum
      
      !----------------------------------
      !-----     NEKLAB EXPORTS     -----
      !----------------------------------
      
      ! Definition of the abstract vectors in the Nek framework.
         public :: nek_dvector
         public :: nek_zvector
      
      ! Implementation of the standard linear operators.
         public :: exptA_linop
         public :: resolvent_linop
      
      !
         public :: linear_stability_analysis_fixed_point
         public :: transient_growth_analysis_fixed_point
      
      ! Various utilities.
         public :: nek2vec, vec2nek
         public :: neklab_forcing
         public :: outpost_vec
      end module neklab
