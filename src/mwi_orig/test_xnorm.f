      integer(4), parameter :: ntrg=4
      integer(4), parameter :: mfreq=6
	real(8) beamwidth11(mfreq),sumgain11(mfreq)
	real(8) beamwidth22(mfreq),sumgain22(mfreq)
      real(8) beamwidth_trg(ntrg),sumgain_trg(ntrg)
	real(8) xnorm_trg,xnorm_obs
	
	
	call openbig(3,'O:\mwi\antenna\data\theoretical_beamwidth_sumgain.dat','old')
	read(3) beamwidth11,sumgain11
	read(3) beamwidth22,sumgain22
	close(3)

	call openbig(3,'O:\mwi\antenna\data\target_beamwidth_sumgain.dat','old')
	read(3) beamwidth_trg,sumgain_trg
	close(3)

      do ifreq=1,mfreq
      do itrg=1,4
      xnorm_trg=sumgain_trg(itrg)
      xnorm_obs=sumgain11(ifreq)
c	weight(jcel,kcel,iscan,ifreq)=(xnorm_trg/xnorm_obs)*a(iobs)
      write(*,*) itrg,ifreq,xnorm_trg/xnorm_obs
      enddo
      enddo
      stop 'norm end'
      end
	include 'x:\syspgm\openbig.for'
