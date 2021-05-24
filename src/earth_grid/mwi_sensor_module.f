c     9/10/2020 changed 10/31/2020 a slight change was made to alpha_start to make it
c     be exactly zero at mid-scan for a spin rate of 26.5 rpm
c     also sign changed to implement counter clockwise scan

c     5/27/2019 changed 9/10/2020 maxcel and alpha_start changed to implement extended scan

c     jun 1 2018 changed 5/27/2019.  89 ghz now has same eia as other channels,
c     hence the is only one value for beta. jun 1 2018 version saved with date suffix


c     dec 31 2017 change jun 1 2018:
c     integration time changed for 1.350 msec to 1.800 msec
c     nadir angle at 89 ghz changed from 44.3 to 45.37
c     alpha_start changed from -80 to -78.99
c     scan period changed from 1.967213d0 to 2.264151d0
c     maxcel changed from 649 to 553
c     numcel array is removed
c     maxscan_89 is removed, it was not used
c     ncolds, nhots, ntherms, nndoide removed.  they were not used
c     nfreq changed from 6 to 5

      module mwi_sensor_module
      implicit none
      
      integer(4), parameter :: maxrev=100000                                 
 
c     integer(4), parameter :: maxcel= 553 !regular scan
      integer(4), parameter :: maxcel= 571 !extended scan
      integer(4), parameter :: maxscan=4400

      integer(4), parameter :: nfreq=5
      integer(4), parameter :: nch=2*nfreq

	real(8), parameter :: scan_period=2.264151d0 
	
	real(8), parameter :: sample_time=1.800d-3   
	real(8), parameter :: beta= 44.30d0 
c	real(8), parameter :: alpha_start=-78.9912d0  !regular scan 
c	real(8), parameter :: alpha_start=-81.5674d0  !extended scan
	real(8), parameter :: alpha_start= 81.5670d0  !extended scan, fine tune, rpm=26.5, sign changed to implement counter clockwise scan
      end module mwi_sensor_module
