module smallFunctions
contains

! Function for factorial, Used by PLGNDR Function   
    integer*8 function fact(x)
    integer*8 x,i

    if (x.eq.0) then
        fact=1
        return
    else if (x.eq.1) then
        fact=1
        return
    else
        fact=1
        do i=1,x
        fact=fact*i
        end do
    end if
    
    return
    end function fact

! Function for Double factorial 
    integer*8 function dfact(x)
    integer*8 x,i

    if (x.eq.1) then
        dfact=1
        return
    else
        dfact=1
        do i=1,x,2
        dfact=dfact*i
        end do
    end if

    return
    end function

end module
