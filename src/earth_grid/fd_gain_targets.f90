module target_patterns

    private
    public :: fd_gain_targets,fd_gain_targets_km

contains
    
    subroutine fd_gain_targets(beamwidth,tht, gain)
        implicit none
 
        real(8), parameter :: pi=3.141592653589793d0
        real(8), parameter :: rad=pi/180.d0
        real(8), parameter :: ln2=dlog(2.d0)
        real(8), parameter :: c1=4.d-4 
        real(8), parameter :: c2=0.456008940886013d0  !=sqrt(0.3d0*dlog(2.d0))

        real(8) beamwidth,tht, gain
        real(8) ththalf,ainv,x

        ththalf=0.5d0*beamwidth
        ainv=(rad*ththalf)**2/ln2
        x=tht/ththalf
        gain = (exp(-ln2*x*x) +  c1*exp(-c2*x))/(1.0021d0*pi*ainv) 
        !1.0021 is fudge factor to bring norm to a little less than unity
        return
    end

    real(8) function fd_gain_targets_km(beamwidth_km,distance_km)

        implicit none 

        real, intent(in) :: beamwidth_km
        real, intent(in) :: distance_km
       
        real(8), parameter :: pi=3.141592653589793d0
        real(8), parameter :: rad=pi/180.d0
        real(8), parameter :: ln2=dlog(2.d0)
        
        real(8) half_km, ainv, x

        half_km = 0.5d0*beamwidth_km
        ainv=(rad*half_km)**2/ln2
        x=distance_km/half_km
        fd_gain_targets_km = exp(-ln2*x*x)/(1.0021d0*pi*ainv) 
        !1.0021 is fudge factor to bring norm to a little less than unity
        return
    end function fd_gain_targets_km
      
end module target_patterns