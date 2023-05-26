!    Forked from o:/ssmis/routines/ssmis_geolocation.f

!    2/22/2008 version changed on 12/4/2009.  now geoloc processing is done even when bit 10 (oob therm) is set
!    i dont think this has any impact on past processing

!    the 4/9/2007 version was changed on 2/22/2008. see O:\ssmis\L2_processing\memo1.txt

    subroutine ssmis_geolocation(igrp)
        use ssmis_module
        use ssmis_l2_module
        implicit none
 
        integer(4), parameter:: imode=1 !1 means only due precession; otherwise do both precession and nutation
        real(8), parameter :: alphax=0  !dummy value for alpha
        real(8), parameter :: betax =0  !dummy value for beta
        
        real(8) x_eci(3),y_eci(3),z_eci(3),l_eci(3)
        real(8) scpos_eci(3),rotangle,days
        real(8) sunvec_eci2000(3),sunvec(3),sundis_km,moonvec_eci2000(3),moonvec(3),moondis_km
        real(8) observation_time
        real(8) precession(3,3),nutation(3,3),np(3,3)
        real(8) b(3),refl(3),xlat,xlon,range,rearth,thtinc,azim,pra,suninc8,sunazm8,sunglt8
        real(8) scpos0(3),scvel0(3)
        real(8) time_since_scan_start
        real(8) r,ru(3),v3,sinlat,coslat,deltau,alpha
        real(8) sc_to_sun(3),dot1,dot2,cossunang
        real(4) xlon_rot
 
        integer(4) igrp,iscan,icel,ierror,iflagx

!    =================================================================================================================
!    ================== do igrp=0 that finds scpos, sun, and moon vector at scan start time =========================
!    =================================================================================================================

        if(igrp.eq.0) then

!    =================================================================================================================
!    ================== initialize arrays with zeros to indicate missing data ========================================
!    =================================================================================================================

        scpos_eci0=0
        x_eci0=0
        y_eci0=0
        z_eci0=0
        sunlat=0
        sunlon=0
        sundis=0
        sunVec0=0
        moonVec0=0
        alpha_sun=0
        beta_sun=0
        kflag_sun=0
        scloc=0
 
        do iscan=1,numscan

            !x    if(iscan_flag(iscan).ne.0) cycle

            iflagx=iscan_flag(iscan)
            iflagx=ibclr(iflagx,10)
            if(iflagx.ne.0) cycle


!           =================================================================================================================
!           ================== compute spacecraft (sc) y coordinate vector (y0) which varies very slowly with time ==========
!           =================================================================================================================

            scpos0=scpos(:,iscan)
            scvel0=scvel(:,iscan)
            call cross_norm(scvel0,scpos0, l_eci)
 
 
            observation_time = scan_time(iscan) 
            days=observation_time/86400.d0 -0.5d0
 
            scpos_eci = scpos0 

!           =================================================================================================================
!           ========================================== find spacecraft x, y, and z coordinate vectors =======================
!           =================================================================================================================

            r=sqrt(dot_product(scpos_eci,scpos_eci))
            ru=scpos_eci/r
            sinlat=ru(3)
            coslat=sqrt(ru(1)*ru(1)+ru(2)*ru(2))
            rearth=rp/sqrt(ffac*coslat**2+sinlat**2)  
            deltau=((r-rearth)/re)*sqrt(1 - (1-ffac)*sinlat**2/((ffac*coslat)**2+sinlat**2))
            z_eci(1)=-ru(1)
            z_eci(2)=-ru(2)
            z_eci(3)=-ru(3)*(1+deltau)/(ffac+deltau)
            z_eci=z_eci/sqrt(dot_product(z_eci,z_eci))    

            call cross_norm(l_eci,z_eci, x_eci)
            call cross_norm(z_eci,x_eci, y_eci)
 
            v3=scvel0(3)/sqrt(dot_product(scvel0,scvel0))
            zang(iscan)=atan2d(ru(3),v3) + 90.
            call fixang8(zang(iscan))

!           =================================================================================================================
!           ============ convert vectors for j2000m coordinates to earth-center inertial coordinates (eci) of date ==========
!           ====================== find angle between vernal equinox and greenwich (rotangle) ==============================
!           =================================================================================================================
 
            call get_gm_angle(observation_time,days, rotangle)
 
            if(imode.eq.1) then
                call fd_precession(days, np)
            else
                call fd_precession(days, precession)
                call fd_nutation(  days, nutation)
                np=matmul(nutation,precession)
            endif
 
 
!            =================================================================================================================
!            =============== compute location of spacecraft subpoint (scloc) and sun information  ============================
!            ==========  calling fd_cell_parameters with option 1 computes location of satellite subpoint  ===================
!            =================================================================================================================
 

            scpos_eci0(:,iscan)=scpos_eci
            x_eci0(:,iscan)=x_eci
            y_eci0(:,iscan)=y_eci
            z_eci0(:,iscan)=z_eci

            call  sun_vector(days,  sunvec_eci2000, sundis_km)     !returns vector in j2000 system, dis is km
            sunvec=matmul(np,sunvec_eci2000)   !unit vector pointing from earth to sun

            call moon_vector(days, moonvec_eci2000,moondis_km)     !returns vector in j2000 system, dis is km
            moonvec=matmul(np,moonvec_eci2000) !unit vector pointing from earth to moon

            sunlat(iscan)=asind(sunvec(3))
            sunlon(iscan)=atan2d(sunvec(2),sunvec(1)) - rotangle  !convert from eci to greenwich
            call fixang4( sunlon(iscan))
            sundis(iscan)=sundis_km/149.619d6

            sunVec0(:,iscan)= sunDis_km*sunVec
            moonVec0(:,iscan)=moonDis_km*moonVec

            sc_to_sun=sundis_km*sunvec - 1.d-3*scpos_eci
            sc_to_sun=sc_to_sun/sqrt(dot_product(sc_to_sun,sc_to_sun))

            cossunang=-dot_product(z_eci,sc_to_sun)    !minus sign accounts for dif between windsat x,y,z and ssmis x,y,z
            if(cossunang.gt. 1) cossunang= 1
            if(cossunang.lt.-1) cossunang=-1
            beta_sun(iscan)=acosd(cossunang)

            dot1= dot_product(x_eci,sc_to_sun)
            dot2=-dot_product(y_eci,sc_to_sun)      !minus sign accounts for dif between windsat x,y,z and ssmis x,y,z
            alpha_sun(iscan)=atan2d(dot2,dot1)

            call fd_cell_parameters(3,re,rp,ffac,scpos_eci,x_eci,y_eci,z_eci,betax,alphax,sunvec,sc_to_sun, &
                                    refl,xlat,xlon,range,rearth,thtinc,azim,pra,suninc8,sunazm8,sunglt8,ierror)
            kflag_sun(iscan)=ierror
        
            call fd_cell_parameters(1,re,rp,ffac,scpos_eci,x_eci,y_eci,z_eci,betax,alphax,sunvec,  &
                                    b,refl,xlat,xlon,range,rearth,thtinc,azim,pra,suninc8,sunazm8,sunglt8,ierror)

            if(ierror.ne.0) then  !subpoint not on earth, which should never occur for ssmis
                stop 'satellite subpoint does not intersect earth, pgm stopped'
            else
                xlon_rot=xlon - rotangle  !convert from inertial to earth
                call fixang4( xlon_rot)
                scloc(1,iscan)=xlat      !subpoint geodetic latitude
                scloc(2,iscan)=xlon_rot      !subpoint longitude
                scloc(3,iscan)=range  !altitude (neters)
            endif
        enddo  !iscan
        return 
    endif     !igrp=0



!    =================================================================================================================
!    ================== do igrp.ne.0 that finds obs cell vectors and values ============================================
!    =================================================================================================================


    if(igrp.ne.0) then
 
!    =================================================================================================================
!    ================== initialize arrays with zeros to indicate missing data ========================================
!    =================================================================================================================

        cellat=0
        cellon=0
        celtht=0
        celphi=0
        celsun=0
 
        do 200 iscan=1,numscan

            iflagx=iscan_flag(iscan)
            iflagx=ibclr(iflagx,10)
            if(iflagx.ne.0) cycle


!           =================================================================================================================
!           ================== compute spacecraft (sc) y coordinate vector (y0) which varies very slowly with time ==========
!           =================================================================================================================

            scpos0=scpos(:,iscan)
            scvel0=scvel(:,iscan)
            call cross_norm(scvel0,scpos0, l_eci)

            if(maxval(abs(scpos_eci0(:,iscan)-scpos0)).gt.0.01) stop 'error in ssmis_geolcation, igrp=0 call not done, pgm stopped'

            sunvec=sunVec0(:,iscan)/sqrt(dot_product(sunVec0(:,iscan),sunVec0(:,iscan)))
 
!           =================================================================================================================
!           ======================================== loop thru cells ========================================================
!           =================================================================================================================
 
            do 100 icel=1,maxcel_grp(igrp)
                alpha=alpha_start(igrp) + (icel-1)*(  144.d0/maxcel_grp(igrp))
                time_since_scan_start =   (icel-1)*(0.7596d0/maxcel_grp(igrp))   !0.7596=(144/360)*1.899
                observation_time = scan_time(iscan) + time_since_scan_start
                days=observation_time/86400.d0 -0.5d0
            
                scpos_eci = scpos0 + time_since_scan_start*scvel0
 
!               =================================================================================================================
!               ========================================== find spacecraft x, y, and z coordinate vectors =======================
!               =================================================================================================================

                r=sqrt(dot_product(scpos_eci,scpos_eci))
                ru=scpos_eci/r
                sinlat=ru(3)
                coslat=sqrt(ru(1)*ru(1)+ru(2)*ru(2))
                rearth=rp/sqrt(ffac*coslat**2+sinlat**2)  
                deltau=((r-rearth)/re)*sqrt(1 - (1-ffac)*sinlat**2/((ffac*coslat)**2+sinlat**2))
                z_eci(1)=-ru(1)
                z_eci(2)=-ru(2)
                z_eci(3)=-ru(3)*(1+deltau)/(ffac+deltau)
                z_eci=z_eci/sqrt(dot_product(z_eci,z_eci))    

                call cross_norm(l_eci,z_eci, x_eci)
                call cross_norm(z_eci,x_eci, y_eci)
            
                !    =================================================================================================================
                !    ============ convert vectors for j2000m coordinates to earth-center inertial coordinates (eci) of date ==========
                !    ====================== find angle between vernal equinox and greenwich (rotangle) ==============================
                !    =================================================================================================================
 
                call get_gm_angle(observation_time,days, rotangle)

 
                !    ===========================================================================================================================
                !    == calling fd_cell_parameters with option 2 computes values at the location where the horn boresight intersect the earth ==
                !    =========================================================================================================================== 

                call fd_cell_parameters(2,re,rp,ffac,scpos_eci,x_eci,y_eci,z_eci,beta(igrp),alpha,sunvec, &
                                        b,refl,xlat,xlon,range,rearth,thtinc,azim,pra,suninc8,sunazm8,sunglt8,ierror)
 
 
                if(ierror.ne.0) then   !boresight does not intersect earth, which should never happen with ssmis
                    stop 'boresight does not intersect earth, pgm stopped'
                else 
                    xlon_rot=xlon - rotangle  !convert from inertial to earth
                    call fixang4( xlon_rot)
                    cellat(icel,iscan)=xlat         ! geodetic latitude of cell (deg)
                    cellon(icel,iscan)=xlon_rot  ! longitude of cell (deg)
                    celtht(icel,iscan)=thtin!    ! look incidence angle at cell (deg)
                    celphi(icel,iscan)=azim         ! look azimuth angle relative to north at cell (deg)
                    celsun(icel,iscan)=sunglt8     ! sun glint angle at cell (deg)
                endif
 
            100 continue  !icel
        200 continue  !iscan
 
        return
        endif  !igrp.ne.0

    end subroutine ssmis_geolocation
 
 
 
 
 
 
!    x, y, z  s/c vectors
!    x vector points in general direction of s/c velocity
!    z vector points is s/c nadir and points towards earth for normal attitude
!    y vector is given by z cross x
 
!    icase=1 means compute satellite subpoint location
!    icase=2 means compute cell locations given alpha and beta
!    icase=3 means compute cell locations given boresight
 
    subroutine fd_cell_parameters(icase,re,rp,ffac,scpos,x,y,z,beta,alpha,sunvec, &
                                  b,refl,xlat,xlon,range,rearth,thtinc,azim,pra,suninc,sunazm,sunglt,ierror)
        implicit none
    
        real(8),    intent(in)    :: re,rp,ffac,scpos(3),x(3),y(3),z(3),beta,alpha,sunvec(3)
        integer(4), intent(in)    :: icase
        real(8),    intent(out)   :: refl(3),xlat,xlon,range,rearth,thtinc,azim,pra,suninc,sunazm,sunglt
        integer(4), intent(out)   :: ierror
        real(8),    intent(inout) :: b(3)
    
        real(8) dot_product_unit8
    
        real(8) r,ru(3),sinlat,coslat,costht
        real(8) sinbeta,cosbeta
        real(8) bx,by,bz,cel1,cel2,cel3
        real(8) coslon,sinlon
        real(8) geoid(3),cosazim,cossunglt
        real(8) deltau
        real(8) ha(3),he(3),ve(3),dot1,dot2
 
        r=sqrt(dot_product(scpos,scpos))
        ru=scpos/r
    
        if(icase.eq.1) then     !find pointing vector b(3) to geoid
            sinlat=ru(3)
            coslat=sqrt(ru(1)*ru(1)+ru(2)*ru(2))
            rearth=rp/sqrt(ffac*coslat**2+sinlat**2)  
            deltau=((r-rearth)/re)*sqrt(1 - (1-ffac)*sinlat**2/((ffac*coslat)**2+sinlat**2))
            b(1)=-ru(1)
            b(2)=-ru(2)
            b(3)=-ru(3)*(1+deltau)/(ffac+deltau)
            b=b/sqrt(dot_product(b,b))    
        endif
 
 
        if(icase.eq.2) then     !compute b(3) from alpha and beta
            cosbeta=cosd(beta)
            sinbeta=sind(beta)
            bx= sinbeta*cosd(alpha) !input alpha defined relative to the x axes,  alpha=-90 denotes bx=0, by=1
            by=-sinbeta*sind(alpha) !input alpha defined relative to the x axes , alpha=-90 denotes bx=0, by=1
            bz= cosbeta                
            b=bx*x + by*y + bz*z
        endif
 
 
            !    for icase=3, b(3) is an input
    
        call loccel_range(r,ru(1),ru(2),ru(3),b(1),b(2),b(3),re,rp, cel1,cel2,cel3,rearth,range,ierror)
 
        if(ierror.ne.0) return     !look vector does not intersect earth

        if(icase.eq.3) return  !for ssmis case 3 is to see if s/c can see sun, so there is nothing else to computed
 
        coslat=sqrt(cel1*cel1+cel2*cel2)     !cosine of geocentric lat
        xlat=atand(cel3/(ffac*coslat))       !geodetic cell latitude
        xlon=atan2d(cel2,cel1)
        call fixang8( xlon)
    
        if(icase.eq.1) return  !nothing else is computed for icase=1
 
        !    =================================================================================================================
        !    ======= compute earth inc. ang. (thtinc) and earth azimuth angle wrt. north (azim) ==============================
        !    =================================================================================================================
        
        coslon=cosd(xlon)
        sinlon=sind(xlon)
        coslat=cosd(xlat)
        geoid(1)=coslat*coslon
        geoid(2)=coslat*sinlon
        geoid(3)=sind(xlat)
    
        costht=-dot_product_unit8(b,geoid)
        thtinc=acosd(costht)
 
        if(abs(costht).gt.0.9999999999d0) then    !too close to nadir to get azim
            azim=0
        else
            cosazim=(-sinlon*b(1)+coslon*b(2))/sqrt(1.-costht*costht)
            if(cosazim.lt.-1.) cosazim=-1.
            if(cosazim.gt. 1.) cosazim= 1.
            azim=acosd(cosazim)
            if(b(3)+costht*geoid(3).lt.0.) azim=-azim
            azim=90.-azim     ! convert to relative to clockwise from north
            if(azim.lt.  0.) azim=azim+360.
            if(azim.gt.360.) azim=azim-360.
        endif
 
        !    =================================================================================================================
        !    ======= compute sun glint angle (sunglt), sun inc. ang. (suninc) and sun azimuth angle (sunazm) =================
        !    =================================================================================================================
 
        refl=b + 2*costht*geoid
    
        cossunglt=dot_product_unit8(refl,sunvec)
        sunglt=acosd(cossunglt)
 
        !    compute sun tht and azimuth angles
        costht=dot_product_unit8(sunvec,geoid)
        suninc=acosd(costht)
 
        if(abs(costht).gt.0.9999999999d0) then    !too close to nadir to get sunazm
            sunazm=0
        else
            cosazim=(-sinlon*sunvec(1)+coslon*sunvec(2))/sqrt(1.-costht*costht)
            if(cosazim.lt.-1.) cosazim=-1.
            if(cosazim.gt. 1.) cosazim= 1.
            sunazm=acosd(cosazim)
            if(sunvec(3)-costht*geoid(3).lt.0.) sunazm=-sunazm     !minus sign used here becuase sun is pointing away from surface
            sunazm=90-sunazm     ! convert to relative to clockwise from north
            call fixang8(sunazm)
        endif
 
 
 
        !    =================================================================================================================
        !    ============================= compute polarization rotation angle (pra) =========================================
        !    =================================================================================================================
 
 
        call cross_norm(b,z, ha)      !ha=kx(-z)=(-b)x(-z)=bxz
 
        if(thtinc.gt.0.01) then
            call cross_norm(geoid,b, he) !he=kxg=(-b)xg=gxb
        else
            he=ha
        endif
 
        call cross_norm(b,he, ve)     !ve=hexk=bxhe
    
        dot1=dot_product(ve,ha)
        dot2=dot_product(he,ha)
        pra=atan2d(dot1,dot2)
 
        return
    end subroutine fd_cell_parameters
 
    subroutine loccel_range(r,r1,r2,r3,s1,s2,s3,re,rp, cel1,cel2,cel3,rearth,range,ierror)
        implicit none
    
        real(8), intent(in)  :: r,r1,r2,r3,s1,s2,s3
        real(8), intent(out) :: cel1,cel2,cel3,rearth,range
        real(8) delta,a,cosx,b,c,celmag,arg
        real(8) re,rp
        integer(4) ierror
 
        delta=(re/rp)*(re/rp)-1.
        a=1.+delta*s3*s3
        cosx=-r1*s1-r2*s2-r3*s3
        b=cosx-delta*r3*s3
        c=1.-(re/r)*(re/r)+delta*r3*r3
 
        arg=b*b-a*c
        if(arg.lt.0) then !vector does not intercept earth
            ierror=1
        return
        else
            ierror=0
        endif
 
        range=(b-sqrt(arg))/a
        cel1=r1+range*s1
        cel2=r2+range*s2
        cel3=r3+range*s3
        celmag=sqrt(cel1*cel1+cel2*cel2+cel3*cel3)
        cel1=cel1/celmag
        cel2=cel2/celmag
        cel3=cel3/celmag
        rearth=r*celmag
        range=range*r
        if(range.le.0) ierror=1
        return
    end subroutine loccel_range
 
 
                
 
 
