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

      implicit none

      real(kind=dbl_kind) :: &
      apond_timeseries(10000), &
      hpond_timeseries(10000), &
      aice_timeseries(10000),  &
      hice_timeseries(10000)

      real(kind=dbl_kind) :: &
      ap_g(10000), &
      hp_g(10000), &
      ai_g(10000), &
      hi_g(10000)

      character(len=*), parameter :: subname='(icedrv)'

      !-----------------------------------------------------------------
      ! Initialize Icepack
      !-----------------------------------------------------------------

      call icedrv_initialize

      !-----------------------------------------------------------------
      ! Run Icepack
      !-----------------------------------------------------------------

      call icedrv_run ( apond_timeseries, hpond_timeseries, aice_timeseries, hice_timeseries )

      !  simulation
      ap_g = apond_timeseries
      hp_g = hpond_timeseries
      ai_g = aice_timeseries
      hi_g = hice_timeseries

      ! save simulation
      200 FORMAT (1x, F16.8)
      open(12, file = 'ap_g.txt')
      write(12,200) ap_g
      close(12)
      open(12, file = 'hp_g.txt')
      write(12,200) hp_g
      close(12)
      open(12, file = 'ai_g.txt')
      write(12,200) ai_g
      close(12)
      open(12, file = 'hi_g.txt')
      write(12,200) hi_g
      close(12)



      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
          file=__FILE__,line= __LINE__)

      write(ice_stdout, *) "ICEPACK COMPLETED SUCCESSFULLY "

      end program icedrv

!=======================================================================
