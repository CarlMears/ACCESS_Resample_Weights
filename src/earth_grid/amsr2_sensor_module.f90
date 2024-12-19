!     9/10/2020 changed 10/31/2020 a slight change was made to alpha_start to make it
!     be exactly zero at mid-scan for a spin rate of 26.5 rpm
!     also sign changed to implement counter clockwise scan

!     5/27/2019 changed 9/10/2020 maxcel and alpha_start changed to implement extended scan

!     jun 1 2018 changed 5/27/2019.  89 ghz now has same eia as other channels,
!     hence the is only one value for beta. jun 1 2018 version saved with date suffix


!     dec 31 2017 change jun 1 2018:
!     integration time changed for 1.350 msec to 1.800 msec
!     nadir angle at 89 ghz changed from 44.3 to 45.37
!     alpha_start changed from -80 to -78.99
!     scan period changed from 1.967213d0 to 2.264151d0
!     maxcel changed from 649 to 553
!     numcel array is removed
!     maxscan_89 is removed, it was not used
!     ncolds, nhots, ntherms, nndoide removed.  they were not used
!     nfreq changed from 6 to 5

    module sensor_module
        implicit none

        public
      
        integer(4), parameter :: maxrev=100000                                 
 
        integer(4), parameter :: maxcel= 243 !regular AMSR2 scan
        integer(4), parameter :: maxscan=4400

        integer(4), parameter :: nfreq=8
        integer(4), parameter :: nch=2*nfreq

        real(8), parameter :: scan_period=60.0d0/40.0d0
    
        real(8), parameter :: sample_time=2.600d-3   
        real(8), dimension(nfreq), parameter :: beta=  &
                    (/47.505d0, 47.505d0, 47.505d0, 47.505d0, 47.505d0, 47.505d0, 47.505d0, 47.115d0/)
        real(8), dimension(nfreq), parameter :: alpha_start = &
                    (/-75.083, -75.083, -75.083, -75.083, -75.083, -75.083, -75.395, -75.725/)


                    
    end module sensor_module
