c     may 6 2012 version changed on may 25 2012
c     option iwindow set back to be standard run:  iwindow=0
c     it had been set to 4 four the 4 hour runs
c     this change was to make sure the standard option is used for future runs


      program run_check_fit 
      use dflib   
	implicit none

	integer(4), parameter :: numrun=40
	integer(4), parameter :: maxline=1000

	character(120) filename
	character(300) aline
	character(200) bline
	character(120) ascreen(numrun),ascreen_remote(numrun)  
	character(20) acomputer
	integer(4) ilinesv,nline(numrun)
	integer(4) icmd,irun,igo
	integer(4)  iopt_cel,kscan,itrg 
	integer(4) iline
	


      write(*,*) 'enter 1 if you have copied the executible to ops'  
	read(*,*) igo
	if(igo.ne.1) stop 'pgm stopped by user'
	
      acomputer='\\ops2d'
      
      irun=0
      do itrg=1,4
      do kscan=0,8
      if(itrg.eq.2 .or. itrg.eq.4) then
      if(kscan.lt.4 .or. kscan.gt.5) cycle 
      endif
      do iopt_cel=1,2
      if(itrg.eq.2 .or. itrg.eq.4) then
      if(iopt_cel.ne.1) cycle
      endif
      irun=irun+1
      if(irun.gt.numrun) stop 'error with numrun, pgm stopped'

	call sleepqq(10000) !pause 10 seconds to give programs a chance to read in input files without conflicting

      write(aline,9001) iopt_cel,kscan,itrg
 9001 format('mwi\check_fit.exe ',3i6)

      write(ascreen(irun),9002) iopt_cel,kscan,itrg
 9002 format('mwi\check_fit_',3i2.2,'.lis')

	ascreen_remote(irun)=trim(acomputer) // '\ops\batch\lists\' // ascreen(irun) 
c     no need to make changes after this line

      do 100 icmd=1,999

c     do use commands that are currently processing

      write(filename,9003) trim(acomputer),icmd
 9003 format(a,'\ops\batch\command_list\command_queued_',i3.3,'.txt')
      open(4,file=filename,status='new',action='write',share='denyrw',err=100)
	write(4,'(i6.6)') irun
	write(4,'(a300)') aline
	write(4,'(a120)') ascreen(irun)
	close(4)
	goto 110
  100 continue
      stop 'command queue is full, try again later'
  110 continue
	write(*,'(a)') trim(aline)//' queued'
      enddo  !iopt_cel
      enddo  !kscan
      enddo  !itrg
      

      write(*,*) 'all run queued successfully, in 30 seconds screens will begin to display'
	call sleepqq(30000) 

c     ======================================================================================================
c     ===================================== starting writing out screens ===================================
c     ======================================================================================================

      nline=0

  250 continue

      do 300 irun=1,numrun
	call sleepqq(2000) 

      open(2,file=ascreen_remote(irun),action='read',share='denynone')

	do iline=1,123456789
	ilinesv=iline
	read(2,'(a)',end=290) bline
	if(iline.ge.nline(irun)-maxline)  write(*,'(a)') trim(bline)  
	enddo

  290 continue
      close(2)
	nline(irun)=ilinesv
	write(*,1002) irun,nline(irun),maxline
 1002 format('run no. = ',i4,' total lines currently on screen = ',i6,' max lines displayed = ',i6)
	write(*,*) '=============================================================================='

  300 continue

      goto 250
	end

