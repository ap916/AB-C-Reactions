subroutine makeMomentumGrid(Kgrid,npts,length)

    implicit none
    integer*8 npts,i,m
    real*8,dimension(1:npts)::Kgrid
    real*8 length,N,pi

    pi = 3.141592653589
    N = real(npts)
    do i=1,npts
        m = i-1
        if (m.lt.N/2.0) then
            kgrid(i) = (2*pi*m)/length
        else if (m.ge.N/2.0) then
            kgrid(i) = (2*pi*(m-N))/length
        endif
    end do
    
end subroutine makeMomentumGrid