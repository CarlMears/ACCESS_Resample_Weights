module pythag_module

    use, intrinsic :: iso_fortran_env, only: int16, int32, real32, real64

    public pythag

contains

    real(real64) function pythag(a,b)

        
        real(real64),intent(in) :: a,b
        real(real64) absa,absb
        absa=abs(a)
        absb=abs(b)
        if(absa.gt.absb)then
            pythag=absa*sqrt(1.+(absb/absa)**2)
        else
            if(absb.eq.0.)then
                pythag=0.
            else
                pythag=absb*sqrt(1.+(absa/absb)**2)
            endif
        endif
        return
    end function pythag
end module pythag_module