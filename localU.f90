module LocalU_mod
    use Multiply_mod
    implicit none
    
    public
    private :: LocalU_metro, phi_new
    
    type(AccCounter) :: Acc_U_local, Acc_U_therm
    real(kind=8) :: phi_new
    
contains
    subroutine LocalU_init()
        call Acc_U_local%init()
        call Acc_U_therm%init()
        return
    end subroutine LocalU_init
    
    subroutine LocalU_clear()
        return
    end subroutine LocalU_clear
    
    subroutine LocalU_reset()
        call Acc_U_local%reset()
        call Acc_U_therm%reset()
        return
    end subroutine LocalU_reset
    
    subroutine LocalU_metro(Op_U, Gr, iseed, nf, ii, ntau)
        use MyMats
! Arguments:
        type(OperatorHubbard), intent(inout) :: Op_U
	    complex(kind=8), dimension(Ndim, Ndim), intent(inout) :: Gr
        integer, intent(inout) :: iseed
        integer, intent(in) :: ii, ntau, nf
!   Local: 
        real(kind=8), external :: ranf
        real(kind=8) :: phi_old, phi_new
        real(kind=8) :: xflip, Xdif, random
        real(kind=8) :: ratio_abs
        complex(kind=8) :: ratio_det, ratio_exp
        integer :: mm, nn

! Local update on space-time (ii, ntau) for auxiliary field flavor (nf)
        phi_old = Conf%phi_list(nf, ii, ntau)
        xflip = ranf(iseed)
        Xdif = dble((xflip - 0.5) * abs(shiftLoc))
        phi_new = phi_old + Xdif
! Calculate Metropolis ratio   
        call Op_U%get_delta(phi_old, phi_new)
        ratio_exp = Op_U%ratio_gaussian
        ratio_det = dcmplx(1.d0,0.d0) + Op_U%Delta * ( dcmplx(1.d0,0.d0) - Gr(ii,ii) )
        ratio_det = dcmplx(1.d0,0.d0) / ratio_det
        ratio_abs = abs(ratio_exp * ratio_det)
! Upgrade Green's function and phi
        random = ranf(iseed)
        if (ratio_abs .gt. random) then
            call Acc_U_local%count(.true.)
            ! Gr(1:Ndim,1:Ndim) = Gr(1:Ndim,1:Ndim) - Gr(1:Ndim,ii) * ratio_det * Op_U%Delta * Gr(ii,1:Ndim)
            call ZGEMM('N', 'N', Ndim, Ndim, 1, -ratio_det * Op_U%Delta, &
                    & Gr(1:Ndim, ii), Ndim, Gr(ii, 1:Ndim), 1, 1.0d0, Gr, Ndim)
            Conf%phi_list(nf, ii, ntau) = phi_new
        else
            call Acc_U_local%count(.false.)
        endif
        return
    end subroutine LocalU_metro
    
    subroutine LocalU_prop_L(Op_U, Prop, iseed, nf, ntau)
        type(OperatorHubbard), intent(inout) :: Op_U
        class(Propagator), intent(inout) :: Prop
        integer, intent(inout) :: iseed
        integer, intent(in) :: ntau, nf
        integer :: ii
        do ii = Ndim, 1, -1
            call LocalU_metro(Op_U, Prop%Gr, iseed, nf, ii, ntau)
            call Op_U%mmult_L(Prop%Gr, Latt, Conf%phi_list(nf, ii, ntau), ii, 1)
            call Op_U%mmult_R(Prop%Gr, Latt, Conf%phi_list(nf, ii, ntau), ii, -1)
        enddo
        ! wrap the left
        do ii = Ndim, 1, -1
            call Op_U%mmult_L(Prop%UUL, Latt, Conf%phi_list(nf, ii, ntau), ii, 1)
        enddo
        return
    end subroutine LocalU_prop_L
    
    subroutine LocalU_prop_R(Op_U, Prop, iseed, nf, ntau)
        type(OperatorHubbard), intent(inout) :: Op_U
        class(Propagator), intent(inout) :: Prop
        integer, intent(inout) :: iseed
        integer, intent(in) :: ntau, nf
        integer :: ii
        do ii = 1, Ndim
            call Op_U%mmult_R(Prop%Gr, Latt, Conf%phi_list(nf, ii, ntau), ii, 1)
            call Op_U%mmult_L(Prop%Gr, Latt, Conf%phi_list(nf, ii, ntau), ii, -1)
            call LocalU_metro(Op_U, Prop%Gr, iseed, nf, ii, ntau)
        enddo
        ! wrap the right
        do ii = 1, Ndim
            call Op_U%mmult_R(Prop%UUR, Latt, Conf%phi_list(nf, ii, ntau), ii, 1)
        enddo
        return
    end subroutine LocalU_prop_R
    
    subroutine LocalU_therm(Op_U, nf, ii, ntau, iseed)
! Arguments: 
        type(OperatorHubbard), intent(inout) :: Op_U
        integer, intent(inout) :: iseed
        integer, intent(in) :: ii, ntau, nf
! Local: 
        real(kind=8), external :: ranf
        real(kind=8) :: phi_old, phi_new
        real(kind=8) :: xflip, Xdif, random
        real(kind=8) :: ratio_abs
        complex(kind=8) :: ratio_exp

! Local update on space-time (ii, ntau) for auxiliary field flavor (nf)
        phi_old = Conf%phi_list(nf, ii, ntau)
        xflip = ranf(iseed)
        Xdif = dble((xflip - 0.5) * abs(shiftLoc))
        phi_new = phi_old + Xdif
! Calculate auxiliary Gaussian ratio   
        call Op_U%get_delta(phi_old, phi_new)
        ratio_exp = Op_U%ratio_gaussian
        ratio_abs = abs(ratio_exp)
! Upgrade phi
        random = ranf(iseed)
        if (ratio_abs .gt. random) then
            call Acc_U_local%count(.true.)
            Conf%phi_list(nf, ii, ntau) = phi_new
        else
            call Acc_U_local%count(.false.)
        endif
        return
    end subroutine LocalU_therm
end module LocalU_mod