
      module icepack_meltpond_cesm_adj

      use icepack_kinds
!      use icepack_parameters, only: c0, c1, c2, p01, puny
!      use icepack_parameters, only: rhofresh, rhoi, rhos, Timelt
!      use icepack_warnings, only: warnstr, icepack_warnings_add
!      use icepack_warnings, only: icepack_warnings_setabort, icepack_warnings_aborted

      implicit none
      real (kind=dbl_kind), parameter, public :: &
         c0   = 0.0_dbl_kind, &
         c1   = 1.0_dbl_kind, &
         c2   = 2.0_dbl_kind, &
         p01  = 0.01_dbl_kind

      real (kind=dbl_kind), public :: &
         puny   = 1.0e-11_dbl_kind

      real (kind=dbl_kind), public :: &
         rhos      = 330.0_dbl_kind   ,&! density of snow (kg/m^3)
         rhoi      = 917.0_dbl_kind   ,&! density of ice (kg/m^3)
         rhow      = 1026.0_dbl_kind  ,&! density of seawater (kg/m^3)
         rhofresh  = 1000.0_dbl_kind    ! density of fresh water (kg/m^3)

      real (kind=dbl_kind), public :: &
         Timelt    = 0.0_dbl_kind     ! melting temperature, ice top surface  (C)

      private
      public :: compute_ponds_cesm_adj

!=======================================================================

      contains

!=======================================================================

      subroutine compute_ponds_cesm_adj(dt,    hi_min,       &
                                    pndaspect9, rfrac,       &
                                    meltt,  melts,           &
                                    frain,                   &
                                    aicen, vicen,            &
                                    Tsfcn,                   &
                                    apnd9,  hpnd9,           &
									pndaspect,               &
									apnd,  hpnd)

      real (kind=dbl_kind):: &
         dt,       & ! time step (s)
         hi_min,   & ! minimum ice thickness allowed for thermo (m)
         pndaspect,   &   ! ratio of pond depth to pond fraction, pertubation
		 pndaspect9   ! ratio of pond depth to pond fraction , basic state

      real (kind=dbl_kind), intent(in) :: &
         rfrac, &    ! water fraction retained for melt ponds
         meltt, &
         melts, &
         frain, &
         aicen, &
         vicen

      real (kind=dbl_kind), intent(in) :: &
         Tsfcn

      real (kind=dbl_kind), intent(inout) :: &
         apnd, &
         hpnd, &
         apnd9, &
         hpnd9

!     local temporary variables

      real (kind=dbl_kind) :: &
         volpn  =c0, &
		 volpn9 =c0

      real (kind=dbl_kind) :: &
         hi  =c0                     , & ! ice thickness (m)
         dTs =c0                     , & ! surface temperature diff for freeze-up (C)
         Tp  =c0                     , & ! pond freezing temperature (C)
         apondn =c0, &
         hpondn =c0, &
!         hi9=c0                     , & ! ice thickness (m)
!         dTs9=c0                    , & ! surface temperature diff for freeze-up (C)
!         Tp9=c0                     , & ! pond freezing temperature (C)
         apondn9 =c0, &
         hpondn9 =c0

!       real (kind=dbl_kind) :: &
!         Td=c0                , & ! temperature difference for freeze-up (C)
!         rexp=c0              , & ! pond contraction scaling
!         dpthhi=c0            ! ratio of pond depth to ice thickness
!      real (kind=dbl_kind), parameter :: &
!
!         Td9       = c2          , & ! temperature difference for freeze-up (C)
!         rexp9     = p01         , & ! pond contraction scaling
!         dpthhi9   = 0.9_dbl_kind  , &  ! ratio of pond depth to ice thickness
!         eps       = 1e-11_dbl_kind
      real (kind=dbl_kind), parameter :: &

         Td       = c2          , & ! temperature difference for freeze-up (C)
         rexp     = p01         , & ! pond contraction scaling
         dpthhi   = 0.9_dbl_kind  , &  ! ratio of pond depth to ice thickness
         eps      = 1e-11_dbl_kind

      character(len=*),parameter :: subname='(compute_ponds_cesm_adj)'


!-----------------------------------------------------------------
! Basic State
!-----------------------------------------------------------------
      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      volpn9 = hpnd9 * apnd9 * aicen
      if (aicen > puny) then

         hi = vicen/aicen

         if (hi < hi_min) then

         !--------------------------------------------------------------
         ! Remove ponds on thin ice
         !--------------------------------------------------------------
            apondn9 = c0
            hpondn9 = c0
            volpn9  = c0

         else
            !-----------------------------------------------------------
            ! Update pond volume
            !-----------------------------------------------------------


            volpn9 = volpn9 &
                  + rfrac/rhofresh*(meltt*rhoi &
                  +                 melts*rhos &
                  +                 frain*  dt)&
                  * aicen

            !-----------------------------------------------------------
            ! Shrink pond volume under freezing conditions
            !-----------------------------------------------------------

			Tp = Timelt - Td

			if (Tp - Tsfcn > c0) then
			  dTs = Tp - Tsfcn
			else
			  dTs = c0
			endif

            volpn9 = volpn9 * exp(rexp*dTs/Tp)

            if (volpn9 > c0) then
			  volpn9 = volpn9
		    else
			  volpn9 = c0
		    endif

            ! fraction of ice covered by ponds
            if (sqrt(volpn9/(pndaspect9*aicen)) < c1) then


              apondn9 = sqrt(volpn9/(pndaspect9*aicen))
		    else

			  apondn9 =c1
			endif


            hpondn9 = pndaspect9 * apondn9

            ! fraction of grid cell covered by ponds

            apondn9 = apondn9 * aicen
            !-----------------------------------------------------------
            ! Limit pond depth
            !-----------------------------------------------------------
			if (hpondn9 > dpthhi*hi) then
			  hpondn9 = dpthhi*hi
			else
			  hpondn9 = hpondn9
			endif


         endif

         !-----------------------------------------------------------
         ! Reload tracer array
         !-----------------------------------------------------------

         apnd9 = apondn9 / aicen
         hpnd9 = hpondn9

      endif
!-----------------------------------------------------------------
! adjoint  state , in reverse mode
!-----------------------------------------------------------------
      !-----------------------------------------------------------------
      ! Identify grid cells where ice can melt
      !-----------------------------------------------------------------
      ! 2
      if (aicen > puny) then
         ! 2-4
         hpondn = hpondn + hpnd
         hpnd = 0
         ! 2-3
         apondn = apondn + c1 / aicen * apnd
		 apnd = 0
         ! 2-1 ( 2-1 and 2-2 are parallel )
         ! nothing to do : hi = vicen/aicen
         ! 2-2
         if (hi < hi_min) then



         !--------------------------------------------------------------
         ! Remove ponds on thin ice
         !--------------------------------------------------------------
         ! 2-2-1 (2-2-1 and 2-2-2~10 are parallel)
            apondn = 0
			hpondn = 0
            volpn  = 0
			apondn9 = c0
            hpondn9 = c0
            volpn9  = c0

         else


        	! 2-2-10
			if (hpondn9 > dpthhi*hi) then
              hpondn = 0
			else
			  hpondn = hpondn
			endif

			! 2-2-9
			apondn = apondn * aicen

            ! 2-2-8
            apondn = apondn &
                   + pndaspect9 * hpondn
            pndaspect = pndaspect &
                      + apondn9 * hpondn
			hpondn = 0

            ! 2-2-7
            if (sqrt(volpn9/(pndaspect9*aicen)) < c1) then
			  volpn = volpn &
			        + c1/c2 &
			        * c1/( sqrt(volpn9/(pndaspect9*aicen))+eps ) &
			        * c1/aicen &
			        * c1/pndaspect9 &
			        * apondn
			  pndaspect = pndaspect &
			        + c1/c2 &
			        * c1/( sqrt(volpn9/(pndaspect9*aicen))+eps ) &
			        * c1/aicen &
			        * -volpn9/(pndaspect9**2) &
			        * apondn
			  apondn = 0
		    else
			  apondn = 0
			endif

            ! 2-2-6
            if (volpn9 > c0) then
			  volpn = volpn
		    else
			  volpn = 0
		    endif

		    ! 2-2-5
!			dTs = dTs
!			Tp = Tp
!			rexp = rexp
!		    volpn = volpn * exp(rexp*dTs/Tp)
!            nothing to do

            ! 2-2-4
!			if (Tp - Tsfcn > c0) then
!			  dTs = Tp - Tsfcn
!			else
!			  dTs = c0
!			endif
!            nothing to do


            ! 2-2-3
!           Td = Td
!			Tp = Tp
!            nothing to do

            !2-2-2
!			volpn = volpn
!            nothing to do

         endif

      endif

      ! 1
      hpnd = hpnd
      apnd = apnd
      hpnd = hpnd + apnd9 * aicen * volpn
      apnd = apnd + hpnd9 * aicen * volpn




      end subroutine compute_ponds_cesm_adj

!=======================================================================

      end module icepack_meltpond_cesm_adj

!=======================================================================
