	subroutine read_gobs_adjust(itarget,npixel_max,nobs_max,itrg,ifreq,icase,iscan,icel,iobs, npixel,ilat2,ilon2,gobs)
	implicit none
	character(80) filename
	integer(2) ilat2x,ilon2x
	integer(4) itarget,npixel_max,nobs_max,itrg,ifreq,icase,iobs,iscan,icel,i,ipixel

	integer(2) ilat2(npixel_max,0:nobs_max),ilon2(npixel_max,0:nobs_max)
	integer(4) npixel(0:nobs_max)
	real(8) gobsx(4)
	real(8) gobs(npixel_max,0:nobs_max)
	real(8) xsum
	
c     the targets are only for g11.
	if(itarget.eq.1 .and. icase.ne.1) stop 'error in calling read_gobs, pgm stopped'

	if(itarget.eq.0) then
      write(filename,9001) ifreq,iscan+15,icel
 9001 format('O:\mwi\resampling\mwi_gain_arrays\freq',i1,'\s',i2.2,'c',i3.3,'.dat')
      else
      write(filename,9002) itrg,ifreq,iscan,icel
 9002 format('O:\mwi\resampling\mwi_target_arrays_adjust\target',i1,'\freq',i1,'\s',i2.2,'c',i4.4,'.dat')
      endif

 	call openbig(3,filename,'old')
	ipixel=0; xsum=0
	do i=1,1234567890
	
	if(itarget.eq.0) then
	read(3,end=100) ilat2x,ilon2x,gobsx
	else
	read(3,end=100) ilat2x,ilon2x,gobsx(1)
	endif
	
	ipixel=ipixel+1
	if(ipixel.gt.npixel_max) stop 'pgm stopped, ipixel oob'
	ilat2(ipixel,iobs)=ilat2x
	ilon2(ipixel,iobs)=ilon2x

c     if you change the icase definition here, you must change it in fd_norm
	if(icase.eq.1) gobs(ipixel,iobs)=gobsx(1) 
	if(icase.eq.2) gobs(ipixel,iobs)=gobsx(2) 
	if(icase.eq.3) gobs(ipixel,iobs)=gobsx(1) + gobsx(4)
	if(icase.eq.4) gobs(ipixel,iobs)=gobsx(1) - gobsx(4)
	if(icase.eq.5) gobs(ipixel,iobs)=gobsx(2) + gobsx(3)
	if(icase.eq.6) gobs(ipixel,iobs)=gobsx(2) - gobsx(3)
	if(icase.eq.7) gobs(ipixel,iobs)=gobsx(2) + (gobsx(2)-gobsx(1))  ! = 2*gobsx(2)-gobsx(1)
	if(icase.eq.8) gobs(ipixel,iobs)=gobsx(2) - (gobsx(2)-gobsx(1)) ! = gobsx(1)
	
	npixel(      iobs)=ipixel
	xsum=xsum+gobs( ipixel,iobs)
	enddo
  100	continue
      close(3)
      gobs(1:ipixel,iobs)=gobs(1:ipixel,iobs)/xsum
	return
	end
	