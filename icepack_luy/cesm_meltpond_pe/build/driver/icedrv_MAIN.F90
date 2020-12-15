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

      !-----------------------------------------------------------------
      ! declare
      !-----------------------------------------------------------------
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
      M=3
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
      ! Run Icepack fwd -  "observation"
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

          call icedrv_run ( pndaspect_true , &
                            apond_timeseries, hpond_timeseries, &
                            apond_timeseriesn, hpond_timeseriesn)



          ! observation
          ap_o = apond_timeseries
          hp_o = hpond_timeseries


          ! save observation
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
              ! Run Icepack fwd - simulation
              !-----------------------------------------------------------------
                  write(12 ,300)
                  print *, '!------------------------------------------------------'
                  print *, '! Run Icepack fwd - simulation '
                  print *, '!------------------------------------------------------'
                  write(12 ,300)

                  !-----------------------------------------------------------------
                  ! Initialize Icepack
                  !-----------------------------------------------------------------

                  call icedrv_initialize

                  !-----------------------------------------------------------------
                  ! Run Icepack
                  !-----------------------------------------------------------------

                  call icedrv_run ( pndaspect , &
                                    apond_timeseries, hpond_timeseries, &
                                    apond_timeseriesn, hpond_timeseriesn)


                  ! simulation
                  ap_g = apond_timeseries
                  hp_g = hpond_timeseries


                  ! save simulation
                  open (12, file='ap_g.txt')
                  write(12 ,200) ap_g
                  close(12)

                  open (12, file='hp_g.txt')
                  write(12 ,200) hp_g
                  close(12)



              !-----------------------------------------------------------------
              ! Calculate Cost Function J
              !-----------------------------------------------------------------
                  300 format( 1x, // )
                  write(12 ,300)
                  print *, '!------------------------------------------------------'
                  print *, '! Calculate Cost Funtcion J'
                  print *, '!------------------------------------------------------'
                  write(12 ,300)


                  J = 0.0_dbl_kind
                  do i = 1, 10000
                     T1 = ( ap_g(i) -ap_o(i) )**2
            !         T1 = ( hp_g(i) -hp_o(i) )**2
                     J = J +T1
                  end do


                  ! save Cost Function
                  open (12, file='J.txt')
                  write(12 ,200) J
                  close(12)


                  print *,'J =' ,J



              !-----------------------------------------------------------------
              ! Run cesm pond adj
              !-----------------------------------------------------------------


                  write(12 ,300)
                  print *, '!------------------------------------------------------'
                  print *, '! Run cesm pond adj '
                  print *, '!------------------------------------------------------'
                  write(12 ,300)


                  ! INITIAL CONDITION
                  ad_pndaspect = 0.0_dbl_kind

                  ncat = 5


                  ! time stepping
                  !-----------------------------------------------------------------
                  do istep1 = 8760, 1, -1


                    ! read in intermediate forcing data  {note: adj forcing is given backward}
                    call luyanginput (aicen, vicen, Tsfc,  meltsn,  melttn, frain ,istep1 )


                    ! additional focing term
                    ad_apnd ( :) = 2*( ap_g(istep1) -ap_o(istep1) )
                    ad_hpnd ( :) = 0.0_dbl_kind



                    ! read in pond basic state
                    do n = 1, ncat
                        apnd (n) = apond_timeseriesn (istep1, n)
                        hpnd (n) = hpond_timeseriesn (istep1, n)
                    enddo



                    do n = 1, ncat

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

                  write(12 ,300)
                  print *, '!------------------------------------------------------'
                  print *, '! Calculate  Gradient G'
                  print *, '!------------------------------------------------------'
                  write(12 ,300)

                  ! calculate gradient
                  G = ad_pndaspect_timeseries (1)


                  ! save Cost Funtcion Gradient
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
              IF(ICALL.GT.3000) GO TO 50

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
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__,line= __LINE__)

      write(ice_stdout, *) "ICEPACK COMPLETED SUCCESSFULLY "

      end program icedrv


!=======================================================================
!
!                          Subroutines
!
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
         '/home/luy/icepack/intermediate_forcing_data/dir2/'

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
