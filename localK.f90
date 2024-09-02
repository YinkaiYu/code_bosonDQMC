module LocalK_mod
    use Multiply_mod
    implicit none
    
    public
    private :: LocalK_metro, phi_new
    
    type(AccCounter) :: Acc_U_local, Acc_U_therm
    real(kind=8) :: phi_new
    
contains
    subroutine LocalK_init()
        call Acc_U_local%init()
        call Acc_U_therm%init()
        return
    end subroutine LocalK_init
    
    subroutine LocalK_clear()
        return
    end subroutine LocalK_clear
    
    subroutine LocalK_reset()
        call Acc_U_local%reset()
        call Acc_U_therm%reset()
        return
    end subroutine LocalK_reset
    
    subroutine LocalK_metro(Op_U, Gr, iseed, nf, ii, ntau)
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
    end subroutine LocalK_metro
    
    subroutine LocalK_prop_L(Prop, iseed, nt)
        class(Propagator), intent(inout) :: Prop
        integer, intent(inout) :: iseed
        integer, intent(in) :: nt
        integer :: ii
        do ii = Lq, 1, -1
            call LocalK_metro(Prop%Gr, iseed, ii, nt)
        enddo
        call Op_U%mmult_L(Prop%Gr, Latt, Conf%phi_list, nt, 1)
        call Op_U%mmult_R(Prop%Gr, Latt, Conf%phi_list, nt, -1)
        call Op_U%mmult_L(Prop%UUL, Latt, Conf%phi_list, nt, 1)
        return
    end subroutine LocalK_prop_L
    
    subroutine LocalK_prop_R(Prop, iseed, nt)
        class(Propagator), intent(inout) :: Prop
        integer, intent(inout) :: iseed
        integer, intent(in) :: nt
        integer :: ii
        call Op_U%mmult_R(Prop%Gr, Latt, Conf%phi_list, nt, 1)
        call Op_U%mmult_L(Prop%Gr, Latt, Conf%phi_list, nt, -1)
        do ii = 1, Lq
            call LocalK_metro(Prop%Gr, iseed, ii, nt)
        enddo
        call Op_U%mmult_R(Prop%UUR, Latt, Conf%phi_list, nt, 1)
        return
    end subroutine LocalK_prop_R
    
    subroutine LocalK_therm(ii, ntau, iseed)
! Arguments: 
        integer, intent(inout) :: iseed
        integer, intent(in) :: ii, ntau
! Local: 
        real(kind=8), external :: ranf
        real(kind=8) :: xflip, random, Xdif, ratio_boson
        integer :: ns

        do ns = 1, Naux
            xflip = ranf(iseed)
            Xdif = dble((xflip - 0.5) * abs(shiftWarm(ns)))
            phi_new(ns, ii, ntau) = Conf%phi_list(ns, ii, ntau) + Xdif
        enddo
        ratio_boson = Conf%bosonratio(phi_new, ii, ntau, Latt)
        random = ranf(iseed)
        if (ratio_boson .gt. random) then
            call Acc_U_therm%count(.true.)
            Conf%phi_list(:, ii, ntau) = phi_new(:, ii, ntau)
        else
            call Acc_U_therm%count(.false.)
            phi_new(:, ii, ntau) = Conf%phi_list(:, ii, ntau)
        endif
        return
    end subroutine LocalK_therm
end module LocalK_mod