module OperatorHubbard_mod
    use MyLattice
    implicit none
    public
    
    type :: OperatorHubbard
        complex(kind=8), private :: alpha ! = sqrt(-2UΔτ)
        complex(kind=8), private :: gaussian
        complex(kind=8), private :: expalpha

        complex(kind=8), public :: Delta
        complex(kind=8), public :: ratio_gaussian
    contains
        procedure           :: set          => opK_set
        procedure, private  :: get_exp      => opK_get_exp
        procedure           :: get_delta    => opK_get_delta
        procedure           :: mmult_R      => opK_mmult_R
        procedure           :: mmult_L      => opK_mmult_L
    end type OperatorHubbard
    
contains
    subroutine opK_set(this, RU)
        class(OperatorHubbard), intent(inout) :: this
        real(kind=8), intent(in) :: RU
        this%alpha = dcmplx( 0.d0, 0.d0 )
        if ( RU < -Zero ) this%alpha = dcmplx( sqrt(-2.d0 * RU * Dtau), 0.d0 )
        if ( RU >  Zero ) this%alpha = dcmplx( 0.d0, sqrt( 2.d0 * RU * Dtau) )
        return
    end subroutine opK_set
    
    subroutine opK_get_exp(this, phi, nflag)
        class(OperatorHubbard), intent(inout) :: this
        integer, intent(in) :: nflag ! +1 or -1; propagating direction
        real(kind=8), intent(in) :: phi ! space time local auxiliary field value, for phi_1 or phi_2
        this%gaussian = dcmplx( exp(-0.5d0 * phi * phi), 0.d0 )
        this%expalpha = exp( this%alpha * phi * nflag )
        return
    end subroutine opK_get_exp
    
    subroutine opK_get_delta(this, phi_old, phi_new)
        class(OperatorHubbard), intent(inout) :: this
        real(kind=8), intent(in) :: phi_old, phi_new
        complex(kind=8) :: gaussian_old, gaussian_new
        complex(kind=8) :: expalpha_old, expalpha_new
        call this%get_exp(phi_old, 1)
        gaussian_old = this%gaussian
        expalpha_old = this%expalpha
        call this%get_exp(phi_new, 1)
        gaussian_new = this%gaussian
        expalpha_new = this%expalpha
        this%ratio_gaussian = gaussian_new / gaussian_old
        this%Delta          = expalpha_new / expalpha_old - dcmplx(1.d0,0.d0)
        return
    end subroutine opK_get_delta
    
    subroutine opK_mmult_R(this, Mat, Latt, phi, ntau, nflag)
! Arguments: 
        class(OperatorHubbard), intent(inout) :: this
        complex(kind=8), dimension(Ndim, Ndim), intent(inout) :: Mat
        class(kagomeLattice), intent(in) :: Latt
        real(kind=8), dimension(Naux, Lq, Ltrot), intent(in) :: phi
        integer, intent(in) :: ntau, nflag
! Local: 
        integer :: P(Norb), ii, no, j, sign
        real(kind=8), dimension(Naux) :: vec
        complex(kind=8), dimension(2, Ndim) :: Vhlp

        do ii = 1, Lq
            vec(:) = phi(:, ii, ntau)
            call this%get_exp(vec, nflag) ! output entryC and entryS
            sign = 1
            do no = 1, Norb
                P(no) = Latt%inv_dim_list(ii, no)
            enddo
            Vhlp = dcmplx(0.d0, 0.d0)
            do j = 1, Ndim
                Vhlp(1, j) = this%entryC * Mat(P(1), j) + sign * this%entryS * Mat(P(2), j)
                Vhlp(2, j) = sign * dconjg(this%entryS) * Mat(P(1), j) + this%entryC * Mat(P(2), j)
            enddo
            do j = 1, Ndim
                Mat(P(1), j) = Vhlp(1, j)
                Mat(P(2), j) = Vhlp(2, j)
            enddo
        enddo
        return
    end subroutine opK_mmult_R
    
    subroutine opK_mmult_L(this, Mat, Latt, phi, ntau, nflag) 
! In Mat Out Mat* EXP(D(NF)) * UT(NF)
! Arguments:
        class(OperatorHubbard), intent(inout) :: this
        complex(kind=8), dimension(Ndim, Ndim), intent(inout) :: Mat
        class(kagomeLattice), intent(in) :: Latt
        real(kind=8), dimension(Naux, Lq, Ltrot), intent(in) :: phi
        integer, intent(in) :: ntau, nflag
! Local: 
        integer :: P(Norb), j, ii, no, sign
        real(kind=8), dimension(Naux) :: vec
        complex(kind=8), dimension(Ndim, 2) :: Uhlp

        do ii = 1, Lq
            vec(:) = phi(:, ii, ntau)
            call this%get_exp(vec, nflag) ! output entryC and entryS
            sign = 1
            do no = 1, Norb
                P(no) = Latt%inv_dim_list(ii, no)
            enddo
            Uhlp = dcmplx(0.d0, 0.d0)
            do j = 1, Ndim
                Uhlp(j, 1) = Mat(j, P(1)) * this%entryC + sign * Mat(j, P(2)) * dconjg(this%entryS)
                Uhlp(j, 2) = sign * Mat(j, P(1)) * this%entryS + Mat(j, P(2)) * this%entryC
            enddo
            do j = 1, Ndim
                Mat(j, P(1)) = Uhlp(j, 1)
                Mat(j, P(2)) = Uhlp(j, 2)
            enddo
        enddo
        return
    end subroutine opK_mmult_L
end module OperatorHubbard_mod