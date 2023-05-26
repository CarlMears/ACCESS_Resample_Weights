module svd_invert

    use pythag_module, only: pythag

    private
    public svdinv

contains
    subroutine svdinv(a,n, a_inv,ierror)
!	    Finds the inverse of a square matrix using the Singular
!	    Value Decomposition method.  Values of 1/w(j) corresponding
!	    to w(j) smaller than a fixed threshold are set to zero.
        use, intrinsic :: iso_fortran_env, only: int16, int32, real32, real64
        implicit none
        integer(int32),intent(in) ::  n
        real(real64), intent(in)  ::  a(1:n,1:n)
        real(real64), intent(out)  ::  a_inv(1:n,1:n)
        integer(int32),intent(out) ::  ierror

        real(real64) :: u(1:n,1:n),w(1:n),w_inv(1:n)
        real(real64) :: v(1:n,1:n),temp(1:n,1:n),wmin,wmax
        integer (int32) :: i

        u=a
        call svdcmp(u,n,n,n,n,w,v,ierror)

        wmax=maxval(w)
        wmin=wmax*10.0**(-14)
        where(w.ge.wmin)
            w_inv=1.0/w
        elsewhere 
            w_inv=0.0
        endwhere

        do i=1,n
            temp(i,:)=w_inv(i)*u(:,i)
        enddo

        a_inv=matmul(v,temp)

    end subroutine svdinv





!     for this special version iteration now go to 200
!     Changes made to this routine
!     1.  Iterations now go to 100
!     2.  Error flag set if it does not converge
!     3.  Stop statement put in for weird case (I assume it will never happen)
!     4.  rv1 is dimension n rather than NMAX.  I looked at the code, and it appears the index for rv1 never will exceed n
    subroutine svdcmp(a,m,mp,n,np,w,v,ierror)

        use, intrinsic :: iso_fortran_env, only: int16, int32, real32, real64
        use pythag_module, only: pythag

        integer(int32) m,mp,n,np,ierror
        real(real64) a(mp,np),v(np,np),w(np)
        integer(int32) i,its,j,jj,k,l,nm
        real(real64) anorm,c,f,g,h,s,scale,x,y,z,rv1(n)

        ierror=0
        g=0.0
        scale=0.0
        anorm=0.0
        do 25 i=1,n
            l=i+1
            rv1(i)=scale*g
            g=0.0
            s=0.0
            scale=0.0
            if(i.le.m)then
                do 11 k=i,m
                    scale=scale+abs(a(k,i))
11              continue
                if(scale.ne.0.0d0)then
                    do 12 k=i,m
                        a(k,i)=a(k,i)/scale
                        s=s+a(k,i)*a(k,i)
12                  continue
                    f=a(i,i)
                    g=-sign(sqrt(s),f)
                    h=f*g-s
                    a(i,i)=f-g
                    do 15 j=l,n
                        s=0.0
                        do 13 k=i,m
                            s=s+a(k,i)*a(k,j)
13                      continue
                        f=s/h
                        do 14 k=i,m
                            a(k,j)=a(k,j)+f*a(k,i)
14                      continue
15                  continue
                    do 16 k=i,m
                        a(k,i)=scale*a(k,i)
16                  continue
                endif
            endif
            w(i)=scale *g
            g=0.0
            s=0.0
            scale=0.0
            if((i.le.m).and.(i.ne.n))then
                do 17 k=l,n
                    scale=scale+abs(a(i,k))
17              continue
                if(scale.ne.0.0d0)then
                    do 18 k=l,n
                        a(i,k)=a(i,k)/scale
                        s=s+a(i,k)*a(i,k)
18                  continue
                    f=a(i,l)
                    g=-sign(sqrt(s),f)
                    h=f*g-s
                    a(i,l)=f-g
                    do 19 k=l,n
                        rv1(k)=a(i,k)/h
19                  continue
                    do 23 j=l,m
                        s=0.0
                        do 21 k=l,n
                            s=s+a(j,k)*a(i,k)
21                      continue
                        do 22 k=l,n
                            a(j,k)=a(j,k)+s*rv1(k)
22                      continue
23                  continue
                    do 24 k=l,n
                        a(i,k)=scale*a(i,k)
24                  continue
                endif
            endif
            anorm=max(anorm,(abs(w(i))+abs(rv1(i))))
25      continue

        do 32 i=n,1,-1
            if(i.lt.n)then
                if(g.ne.0.0)then
                    do 26 j=l,n
                        v(j,i)=(a(i,j)/a(i,l))/g
26                      continue
                    do 29 j=l,n
                        s=0.0
                        do 27 k=l,n
                            s=s+a(i,k)*v(k,j)
27                          continue
                        do 28 k=l,n
                            v(k,j)=v(k,j)+s*v(k,i)
28                          continue
29                      continue
                endif
                do 31 j=l,n
                    v(i,j)=0.0
                    v(j,i)=0.0
31                  continue
            endif
            v(i,i)=1.0
            g=rv1(i)
            l=i
32      continue

        do 39 i=min(m,n),1,-1
            l=i+1
            g=w(i)
            do 33 j=l,n
                a(i,j)=0.0
33          continue
            if(g.ne.0.0d0)then
                g=1.0/g
                do 36 j=l,n
                    s=0.0
                    do 34 k=l,m
                        s=s+a(k,i)*a(k,j)
34                  continue
                    f=(s/a(i,i))*g
                    do 35 k=i,m
                        a(k,j)=a(k,j)+f*a(k,i)
35                  continue
36              continue
                do 37 j=i,m
                    a(j,i)=a(j,i)*g
37              continue
            else
                do 38 j= i,m
                    a(j,i)=0.0
38              continue
            endif
            a(i,i)=a(i,i)+1.0
39      continue

        do 49 k=n,1,-1
            do 48 its=1,200
                do 41 l=k,1,-1
                    nm=l-1
                    if((abs(rv1(l))+anorm).eq.anorm)  goto 2
                    if((abs(w(nm))+anorm).eq.anorm)  goto 1
41              continue
                write(*,*) l
                stop 'svdcmp stopped, i dont think it should ever ket to this point, else l would be 0'
1               c=0.0
                s=1.0
                do 43 i=l,k
                    f=s*rv1(i)
                    rv1(i)=c*rv1(i)
                    if((abs(f)+anorm).eq.anorm) goto 2
                    g=w(i)
                    h=pythag(f,g)
                    w(i)=h
                    h=1.0/h
                    c= (g*h)
                    s=-(f*h)
                    do 42 j=1,m
                        y=a(j,nm)
                        z=a(j,i)
                        a(j,nm)=(y*c)+(z*s)
                        a(j,i)=-(y*s)+(z*c)
42                  continue
43              continue
2               z=w(k)
                if(l.eq.k)then
                    if(z.lt.0.0)then
                        w(k)=-z
                        do 44 j=1,n
                            v(j,k)=-v(j,k)
44                      continue
                    endif
                    goto 3
                endif

                if(its.eq.300) WRITE(*,*) 'no convergence in svdcmp'
                if(its.eq.300) ierror=1
                x=w(l)
                nm=k-1
                y=w(nm)
                g=rv1(nm)
                h=rv1(k)
                f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y)
                g=pythag(f,1.D0)
                f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
                c=1.0
                s=1.0
                do 47 j=l,nm
                    i=j+1
                    g=rv1(i)
                    y=w(i)
                    h=s*g
                    g=c*g
                    z=pythag(f,h)
                    rv1(j)=z
                    c=f/z
                    s=h/z
                    f= (x*c)+(g*s)
                    g=-(x*s)+(g*c)
                    h=y*s
                    y=y*c
                    do 45 jj=1,n
                        x=v(jj,j)
                        z=v(jj,i)
                        v(jj,j)= (x*c)+(z*s)
                        v(jj,i)=-(x*s)+(z*c)
45                  continue
                    z=pythag(f,h)
                    w(j)=z
                    if(z.ne.0.0)then
                        z=1.0/z
                        c=f*z
                        s=h*z
                    endif
                    f= (c*g)+(s*y)
                    x=-(s*g)+(c*y)
                    do 46 jj=1,m
                        y=a(jj,j)
                        z=a(jj,i)
                        a(jj,j)= (y*c)+(z*s)
                        a(jj,i)=-(y*s)+(z*c)
46                  continue
47              continue
                rv1(l)=0.0
                rv1(k)=f
                w(k)=x
48          continue
3       continue  !goto target
49      continue
        return
    end subroutine svdcmp

    
end module svd_invert

