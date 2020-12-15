!=======================================================================
! Copyright (c) 2018, Los Alamos National Security, LLC
! All rights reserved.
!
! Copyright 2018. Los Alamos National Security, LLC. This software was
! produced under U.S. Government contract DE-AC52-06NA25396 for Los
! Alamos National Laboratory (LANL), which is operated by Los Alamos
! National Security, LLC for the U.S. Department of Energy. The U.S.
! Government has rights to use, reproduce, and distribute this software.
! NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC MAKES ANY
! WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF
! THIS SOFTWARE. If software is modified to produce derivative works,
! such modified software should be clearly marked, so as not to confuse
! it with the version available from LANL.
!
! The full license and distribution policy are available from
! https://github.com/CICE-Consortium
!
!=======================================================================

! Main driver routine for Icepack, the column package for CICE.
! Initializes and steps through the model.
!
! author Elizabeth C. Hunke, LANL
!
      program icedrv

      use icedrv_InitMod
      use icedrv_RunMod
      use icedrv_constants, only: ice_stdout, nu_diag
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
      use icedrv_system, only: icedrv_system_abort
      use icepack_kinds
      use icepack_meltpond_cesm_adj, only: compute_ponds_cesm_adj


      implicit none


      INTEGER,PARAMETER ::  NDIM= 1, MSAVE=7, NWORK=NDIM*(2*MSAVE +1)+2*MSAVE
      DOUBLE PRECISION X(NDIM), G(NDIM), DIAG(NDIM), W(NWORK)
      DOUBLE PRECISION F, EPS, XTOL, GTOL, T1, T2, STPMIN, STPMAX, J
      INTEGER IPRINT(2),IFLAG, ICALL, NN, M, MP, LP
      LOGICAL DIAGCO
      character(len=*), parameter :: subname='(icedrv)'



      real (kind=dbl_kind)  :: &
         pndaspect    ,&! ratio of pond depth to area fraction
         theta_true   ,&
         pndaspect_true ,&
         theta_guess    ,&
         pndaspect_guess

      real (kind=dbl_kind)  :: &
         apond_timeseries(10000) ,&
         hpond_timeseries(10000)

      real (kind=dbl_kind) ,dimension(10000, 5) :: &
         apond_timeseriesn ,&
         hpond_timeseriesn

      real (kind=dbl_kind)  :: &
         ap_o(10000) ,&
         hp_o(10000)

      real (kind=dbl_kind)  :: &
         ap_g(10000) ,&
         hp_g(10000)

      integer (kind=int_kind) :: &
         i            ! time index

      real (kind=dbl_kind), dimension(5) :: &

         aicen       , & ! concentration of ice
         vicen       , & ! volume per unit area of ice          (m)
         Tsfc        , & ! ice/snow surface temperature, Tsfcn
         meltsn      , & ! snow melt                       (m)
         melttn      , & ! top ice melt                    (m)

         apnd        , & ! melt pond area fraction
         hpnd         ! melt pond depth (m)

      real (kind=dbl_kind) :: &
         frain        ! rainfall rate (kg/m^2 s)


      integer (kind=int_kind) :: &
         istep1       , &! counter, number of steps at current timestep
         ncat         , &
         n            , &
         T               ! integration steps

      real (kind=dbl_kind) :: &
         dt              ! thermodynamics timestep (s)

      real (kind=dbl_kind) :: &
         hi_min          ! minimum ice thickness allowed (m)

      real (kind=dbl_kind) :: &
         rfrac

      real (kind=dbl_kind) :: &
         rfracmin  = 0.15_dbl_kind, & ! minimum retained fraction of meltwater
         rfracmax  = 0.85_dbl_kind    ! maximum retained fraction of meltwater

      real (kind=dbl_kind), dimension(5) :: &
         ad_apnd        , & ! melt pond area fraction
         ad_hpnd         ! melt pond depth (m)

      real (kind=dbl_kind), dimension(5) ::  &
         ad_pndaspect

      real (kind=dbl_kind), dimension(10000)  :: &
         ad_pndaspect_timeseries

      integer (kind=int_kind) :: &
         ii



      !-----------------------------------------------------------------
      ! Set lbfgs parameters
      !-----------------------------------------------------------------
!     The driver for LBFGS must always declare LB2 as EXTERNAL
!
      EXTERNAL LB2
      COMMON /LB3/MP,LP,GTOL,STPMIN,STPMAX
!
      NN=1
      M=5
      IPRINT(1)= 1
      IPRINT(2)= 3
!
!     We do not wish to provide the diagonal matrices Hk0, and
!     therefore set DIAGCO to FALSE.
!
      DIAGCO= .FALSE.
      EPS= 1.0D-5
      XTOL= 1.0D-16
      ICALL=0
      IFLAG=0




      !-----------------------------------------------------------------
      ! Run Icepack fwd - obs
      !-----------------------------------------------------------------
      write(12 ,300)
      print *, '!------------------------------------------------------'
      print *, '! Run Icepack fwd - obs '
      print *, '!------------------------------------------------------'
      write(12 ,300)
      !-----------------------------------------------------------------
      ! Initialize Icepack
      !-----------------------------------------------------------------

      call icedrv_initialize

      !-----------------------------------------------------------------
      ! Run Icepack
      !-----------------------------------------------------------------
      ! theta true
      theta_true = 0.8_dbl_kind
      pndaspect_true = theta_true

      call icedrv_run ( pndaspect_true , apond_timeseries, hpond_timeseries, apond_timeseriesn, hpond_timeseriesn)




      ! observation
      ap_o = apond_timeseries
      hp_o = hpond_timeseries
!
!      save observation.mat xobs yobs zobs;
      open (12, file='ap_o.txt')
      write(12 ,200) ap_o
      close(12)
      open (12, file='hp_o.txt')
      write(12 ,200) hp_o
      close(12)
      200 FORMAT(1x, F16.8)



      !-----------------------------------------------------------------
      ! first - guess
      !-----------------------------------------------------------------
      X  = 0.7_dbl_kind

      pndaspect = X(1)



      !-----------------------------------------------------------------
      ! Optimalization LOOP
      !-----------------------------------------------------------------
      write(12 ,300)
      print *, '!------------------------------------------------------'
      print *, '! Optimalization LOOPs '
      print *, '!------------------------------------------------------'
      write(12 ,300)

 20   CONTINUE
      J = 0.D0

      pndaspect = X(1)

      !-----------------------------------------------------------------
      ! Run Icepack fwd - sim
      !-----------------------------------------------------------------
      write(12 ,300)
      print *, '!------------------------------------------------------'
      print *, '! Run Icepack fwd - sim '
      print *, '!------------------------------------------------------'
      write(12 ,300)

      !-----------------------------------------------------------------
      ! Initialize Icepack
      !-----------------------------------------------------------------

      call icedrv_initialize

      !-----------------------------------------------------------------
      ! Run Icepack
      !-----------------------------------------------------------------

      call icedrv_run ( pndaspect , apond_timeseries, hpond_timeseries, apond_timeseriesn, hpond_timeseriesn)

      ! simulation
      ap_g = apond_timeseries
      hp_g = hpond_timeseries

      open (12, file='ap_g.txt')
      write(12 ,200) ap_g
      close(12)
      open (12, file='hp_g.txt')
      write(12 ,200) hp_g
      close(12)


!      do ii = 1,10000
!        ap_g(ii)=apond_timeseriesn(ii,3)
!        hp_g(ii)=hpond_timeseriesn(ii,3)
!      end do
!
!      open (12, file='ap_g.txt')
!      write(12 ,200) ap_g
!      close(12)
!      open (12, file='hp_g.txt')
!      write(12 ,200) hp_g
!      close(12)



      !-----------------------------------------------------------------
      ! Calculate Cost Funtcion J
      !-----------------------------------------------------------------
      300 format( 1x, // )
      write(12 ,300)
      print *, '!------------------------------------------------------'
      print *, '! Calculate Cost Funtcion J'
      print *, '!------------------------------------------------------'
      write(12 ,300)

!      J = 0.0_dbl_kind
!      do i = 1, 10000
!         T1 = ( ap_g(i) -ap_o(i) )**2
!         T2 = ( hp_g(i) -hp_o(i) )**2
!         J = J +T1 +T2
!      end do

      J = 0.0_dbl_kind
      do i = 1, 10000
!         T1 = ( ap_g(i) -ap_o(i) )**2
         T1 = ( hp_g(i) -hp_o(i) )**2
         J = J +T1
      end do

      open (12, file='J.txt')
      write(12 ,200) J
      close(12)


      print *,'J =' ,J



      !-----------------------------------------------------------------
      ! Run cesm pond adj
      !-----------------------------------------------------------------
!      call luyanginput (aicen, vicen, Tsfc,  meltsn,  melttn, frain ,istep1 )
!      call luyanginput (aicen, vicen, Tsfc,  meltsn,  melttn, frain ,istep1 )
!      call icedrv_adj_run ( pndaspect_guess ,   adj_pndaspect )
!
!

      write(12 ,300)
      print *, '!------------------------------------------------------'
      print *, '! Run cesm pond adj '
      print *, '!------------------------------------------------------'
      write(12 ,300)

!      print *,'hi_min =' ,hi_min
!      print *,'n =' ,n

      ! INITIAL CONDITION
      !-----------------------------------------------------------------
      !J transpose
!      ad_apnd (:) = 1.0_dbl_kind
!      ad_apnd (n) = 2*( ap_g(8760) -ap_o(8760) )
!      ad_hpnd (:) = 0.0_dbl_kind
      ad_pndaspect = 0.0_dbl_kind
!                ad_apnd (n) =ad_apnd (n) + 2*( ap_g(T) -ap_o(T) )*adJ
!                ad_hpnd (n) =ad_hpnd (n) + 2*( hp_g(T) -hp_o(T) )*adJ
!                adJ =0
!      print *,'INITIAL CONDITION :'
!      print *,' ad_apnd =' ,ad_apnd
!      print *,' ad_hpnd =' ,ad_hpnd


      ncat = 5
!      print *,'ncat =' ,ncat



      ! time stepping
      !-----------------------------------------------------------------
      do istep1 = 8760, 1, -1

!        print *,'istep1 =' ,istep1

        ! read in intermediate forcing data !ad forcing backward
        call luyanginput (aicen, vicen, Tsfc,  meltsn,  melttn, frain ,istep1 )

!        print *,  '  FORCING :'
!            print *, '  aicen(4) =' ,aicen (4)
!            print *, '  vicen(4) =' ,vicen (4)
!            print *, '  Tsfc(4) ='  ,Tsfc (4)
!            print *, '  meltsn(4) =' ,meltsn (4)
!            print *, '  melttn(4) =' ,melttn (4)
!            print *, '  frain ='  ,frain

        ! additional focing term
!        ad_apnd ( :) = 1.0_dbl_kind
!        ad_hpnd ( :) = 0.0_dbl_kind
!        ad_apnd ( :) = 2*( ap_g(istep1) -ap_o(istep1) )
!        ad_hpnd ( :) = 0.0_dbl_kind
        ad_apnd ( :) = 0.0_dbl_kind
        ad_hpnd ( :) = 2*( hp_g(istep1) -hp_o(istep1) )

!        print *, '  additional focing term:'
!        print *, '  ad_apnd (:)  =' ,ad_apnd

        ! read in pond basic state
        do n = 1, ncat
            apnd (n) = apond_timeseriesn (istep1, n)
            hpnd (n) = hpond_timeseriesn (istep1, n)
        ENDDO

!        print *, '  BASIC STATE :'
!            print *, '  apnd(4) =' ,apnd (4)
!            print *, '  hpnd(4) =' ,hpnd (4)



        do n = 1, ncat

!            print *,'  n =' ,n

                    rfrac = rfracmin + (rfracmax-rfracmin) * aicen(n)

                                   call compute_ponds_cesm_adj  (dt,    hi_min,    &
                                                           pndaspect, rfrac, &
                                                           melttn(n), meltsn(n), &
                                                           frain,                &
                                                           aicen (n), vicen (n), &
                                                           Tsfc  (n), &
                                                           apnd  (n), hpnd  (n),&
                                                           ad_pndaspect (n)         ,&
                                                           ad_apnd  (n), ad_hpnd  (n) &
                                                           )

         enddo



         ! assemble adj time series
         ad_pndaspect_timeseries (istep1) = sum(ad_pndaspect)

!         print *, '  ADJOINT STATE :'
!         print *, '  ad_pndaspect =' , ad_pndaspect
!         print *, '  ad_pndaspect_timeseries (istep1) =' ,  ad_pndaspect_timeseries (istep1)


      enddo
      !-----------------------------------------------------------------
      ! END time stepping


      !ad results

      open (12, file='ad_pndaspect_timeseries.txt')
      write(12 ,200) ad_pndaspect_timeseries
      close(12)



      !-----------------------------------------------------------------
      ! Calculate  Gradient
      !-----------------------------------------------------------------
      ! calculate gradient
      write(12 ,300)
      print *, '!------------------------------------------------------'
      print *, '! Calculate  Gradient G'
      print *, '!------------------------------------------------------'
      write(12 ,300)

      G = ad_pndaspect_timeseries (1)

      open (12, file='G.txt')
      write(12 ,200) G
      close(12)


      !-----------------------------------------------------------------
      ! fminlbfgs
      !-----------------------------------------------------------------
      write(12 ,300)
      print *, '!------------------------------------------------------'
      print *, '! fminlbfgs '
      print *, '!------------------------------------------------------'
      write(12 ,300)

      CALL LBFGS (NN, M, X, J, G, DIAGCO,DIAG,IPRINT,EPS,XTOL,W,IFLAG)

      IF(IFLAG.LE.0) GO TO 50

      ICALL=ICALL + 1

!C     We allow at most 2000 evaluations of F and G
      IF(ICALL.GT.2000) GO TO 50

      GO TO 20

  50  CONTINUE
      !-----------------------------------------------------------------
      ! END Optimalization LOOP
      !-----------------------------------------------------------------
      write(12 ,300)
      print *, '!------------------------------------------------------'
      print *, '! END Optimalization LOOPs '
      print *, '!------------------------------------------------------'
      write(12 ,300)



      !-----------------------------------------------------------------
      !
      !-----------------------------------------------------------------
!      write(12 ,300)
!      print *, '!------------------------------------------------------'
!      print *, '!  '
!      print *, '!------------------------------------------------------'
!      write(12 ,300)

      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__,line= __LINE__)

      write(ice_stdout, *) "ICEPACK COMPLETED SUCCESSFULLY "

      end program icedrv



!=======================================================================

      subroutine luyanginput( aicen, vicen, Tsfc,  meltsn,  melttn, frain ,istep1)

      use icepack_kinds

      implicit none

      character*6::char_count

      character*100::data_dir

      real (kind=dbl_kind), dimension(5) :: &

         aicen       , & ! concentration of ice
         vicen       , & ! volume per unit area of ice          (m)
         Tsfc        , & ! ice/snow surface temperature, Tsfcn
         meltsn      , & ! snow melt                       (m)
         melttn        ! top ice melt                    (m)

      real (kind=dbl_kind) :: &

         frain        ! rainfall rate (kg/m^2 s)

      integer (kind=int_kind) :: &
         istep1     ! counter, number of steps at current timestep

            call int_to_char(6,  istep1   ,char_count)

            data_dir=&
         '/home/luy/icepack/intermediate_forcing_data/dir1/'

            open(12,file=trim(data_dir)//'melttn'//char_count//".txt") !
            read(12,200)melttn
            close(12)

            open(12,file=trim(data_dir)//'meltsn'//char_count//".txt") !
            read(12,200)meltsn
            close(12)

            open(12,file=trim(data_dir)//'frain'//char_count//".txt") !
            read(12,200)frain
            close(12)

            open(12,file=trim(data_dir)//'aicen'//char_count//".txt") !
            read(12,200)aicen
            close(12)

            open(12,file=trim(data_dir)//'vicen'//char_count//".txt") !
            read(12,200)vicen
            close(12)

            open(12,file=trim(data_dir)//'Tsfc'//char_count//".txt") !
            read(12,200)Tsfc
            close(12)

            200 FORMAT(1x, F16.8)


      end subroutine luyanginput

      subroutine int_to_char(char_len,int,chara)
            implicit none
            integer,intent(in)::char_len,int
            character,intent(out)::chara(1:char_len)
            integer::ii,kk

               kk=int
               do ii=char_len,1,-1
                  chara(ii:ii)=char(mod(kk,10)+ichar('0'))
                  kk=kk/10
               end do
      end subroutine int_to_char
