!    6/2/2021 converted to f2003, IEEE types
!
!     5/9/2008 version updated 2/4/2014.  routine minang4 and minang8 added

!     11/10/2007 version updated on 5/9/2008.  The fixang routines were made more robust to 
!     floating point inconsistencies

module math_routines
    use, intrinsic :: iso_fortran_env, only: real32, real64, int32, ERROR_UNIT
    implicit none

    private
    public :: cross_norm,dot_product_unit8,invert_3by3,fixang4,fixang8,minang4,minang8
contains

    subroutine cross_norm(x, y, z)
        real(real64) x(3),y(3),z(3),xmag
    
        z(1)=x(2)*y(3)-x(3)*y(2)
        z(2)=x(3)*y(1)-x(1)*y(3)
        z(3)=x(1)*y(2)-x(2)*y(1)
    
        xmag=sqrt(dot_product(z,z))
        z=z/xmag
    end subroutine cross_norm
 

    real(real64) function dot_product_unit8(a,b)

        real(real64) a(3),b(3),dot
        dot=dot_product(a,b)
        if(dot.lt.-1.) dot=-1.
        if(dot.gt. 1.) dot= 1.
        dot_product_unit8=dot
    end functIon dot_product_unit8
 
 
    subroutine invert_3by3(a, ainv,det)

        real(real64),dimension(3,3),intent(in)  ::  a
        real(real64),dimension(3,3),intent(out) ::  ainv
        real(real64),intent(out)                :: det

        real(real64) q1,q2,q3
 
        q1=a(2,2)*a(3,3)-a(3,2)*a(2,3)
        q2=a(3,2)*a(1,3)-a(1,2)*a(3,3)
        q3=a(1,2)*a(2,3)-a(2,2)*a(1,3)
        det=a(1,1)*q1 + a(2,1)*q2 + a(3,1)*q3
        if(abs(det).lt.1.e-30) return
 
        ainv(1,1)=q1
        ainv(1,2)=q2
        ainv(1,3)=q3
    
        ainv(2,1)=a(3,1)*a(2,3)-a(2,1)*a(3,3)
        ainv(2,2)=a(1,1)*a(3,3)-a(3,1)*a(1,3)
        ainv(2,3)=a(2,1)*a(1,3)-a(1,1)*a(2,3)
 
        ainv(3,1)=a(2,1)*a(3,2)-a(3,1)*a(2,2)
        ainv(3,2)=a(3,1)*a(1,2)-a(1,1)*a(3,2)
        ainv(3,3)=a(1,1)*a(2,2)-a(2,1)*a(1,2)
 
        ainv=ainv/det
    end
 
    subroutine fixang4( angle)

        real(real32),intent(inout) :: angle
        angle=angle-360.*floor(angle/360.)  
        if(angle.lt.0 .or. angle.gt.359.9999) angle=0
    end
 
    subroutine fixang8( angle)

        real(real64),intent(inout) :: angle
        angle=angle-360.d0*floor(angle/360.d0)  
        if(angle.lt.0 .or. angle.gt.359.9999999999d0) angle=0
    end    
 
    subroutine minang4( angle)

        real(real32),intent(inout) :: angle
        if(angle.lt.-180.) angle=angle+360.
        if(angle.gt. 180.) angle=angle-360.
    end
 
    subroutine minang8( angle)

        real(real64),intent(inout) :: angle
        if(angle.lt.-180.0D0) angle=angle+360.0D0
        if(angle.gt. 180.0D0) angle=angle-360.0D0
    end

end module math_routines