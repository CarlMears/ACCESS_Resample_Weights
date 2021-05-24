      program show_footprintsize
	implicit none

      real(4),    parameter :: rad=0.01745329252
      real(4),    parameter :: re= 6377.933  !this is re for resampling program

	integer(4) itrg
	real(4) eia,eia1,eia2,dis,dis1,dis2

      real(4) thtnadir,beamwidth(2),scalt
	real(4) thtn,thtn1,thtn2,sinnad,sinnad1,sinnad2
      
c      data scalt/825.463/  !this is altitude for resampling program
      data scalt/833./  !this is altitude for resampling program
      data thtnadir/0.0d0/
      
c     data beamwidth/2.637, 1.735/  !determined to give 38 and 25 km footprints
      data beamwidth/2.6131, 1.7195/  !determined to give 38 and 25 km footprints
     
	do itrg=1,2
	
	thtn =thtnadir
	thtn1=thtn - beamwidth(itrg)/2.
	thtn2=thtn + beamwidth(itrg)/2.
	
	sinnad =sind(thtn)
	sinnad1=sind(thtn1)
	sinnad2=sind(thtn2)

      eia =asind((scalt+re)*sinnad/re)
      eia1=asind((scalt+re)*sinnad1/re)
      eia2=asind((scalt+re)*sinnad2/re)
      
      
      dis =re*rad*(eia - thtn)
      dis1=re*rad*(eia1 - thtn1)
      dis2=re*rad*(eia2 - thtn2)
 
      
      write(*,'(i3,10f10.3)')  itrg,beamwidth(itrg),thtn,eia,dis2-dis1
      
      enddo
	stop 'norm end'
      end
 
