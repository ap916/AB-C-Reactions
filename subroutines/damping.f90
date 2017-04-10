! xd is some critical point 

subroutine damping(Damp,gridpts,isize)
    implicit none
    integer*8 i,isize
    real*8 xd,dx
    real*8,dimension(isize)::Damp
    real*8,dimension(1:isize),intent(in)::gridpts

    dx = 0.01
    xd = gridpts(isize) - gridpts(1) - 3.0

    do i = 1,isize
        if ( gridpts(i) .le. xd ) then
            Damp(i) = 1.0
        else
            Damp(i) = exp(-1.0*dx*(gridpts(i)-xd)**2)
        end if 
    end do

end subroutine  damping
