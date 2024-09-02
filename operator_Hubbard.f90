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
        procedure           :: set          => opU_set
        procedure, private  :: get_exp      => opU_get_exp
        procedure           :: get_delta    => opU_get_delta
        procedure           :: mmult_R      => opU_mmult_R
        procedure           :: mmult_L      => opU_mmult_L
    end type OperatorHubbard
    
contains
    subroutine opU_set(this, RU)
        class(OperatorHubbard), intent(inout) :: this
        real(kind=8), intent(in) :: RU
        this%alpha = dcmplx( 0.d0, 0.d0 )
        if ( RU < -Zero ) this%alpha = dcmplx( sqrt(-2.d0 * RU * Dtau), 0.d0 )
        if ( RU >  Zero ) this%alpha = dcmplx( 0.d0, sqrt( 2.d0 * RU * Dtau) )
        return
    end subroutine opU_set
    
    subroutine opU_get_exp(this, phi, nflag)
        class(OperatorHubbard), intent(inout) :: this
        integer, intent(in) :: nflag ! +1 or -1; propagating direction
        real(kind=8), intent(in) :: phi ! space time local auxiliary field value, for phi_1 or phi_2
        this%gaussian = dcmplx( exp(-0.5d0 * phi * phi), 0.d0 )
        this%expalpha = exp( this%alpha * phi * nflag )
        return
    end subroutine opU_get_exp
    
    subroutine opU_get_delta(this, phi_old, phi_new)
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
    end subroutine opU_get_delta
    
    subroutine opU_mmult_R(this, Mat, Latt, phi, ii, nflag)
! Arguments: 
        class(OperatorHubbard), intent(inout) :: this
        complex(kind=8), dimension(Ndim, Ndim), intent(inout) :: Mat
        class(kagomeLattice), intent(in) :: Latt
        real(kind=8), intent(in) :: phi
        integer, intent(in) :: ii, nflag
! Local: 
        integer :: jj

        call this%get_exp(phi, nflag)
        do jj = 1, Ndim
            Mat(ii,jj) = this%expalpha * Mat(ii,jj)
        enddo
        return
    end subroutine opU_mmult_R
    
    subroutine opU_mmult_L(this, Mat, Latt, phi, jj, nflag)
! Arguments: 
        class(OperatorHubbard), intent(inout) :: this
        complex(kind=8), dimension(Ndim, Ndim), intent(inout) :: Mat
        class(kagomeLattice), intent(in) :: Latt
        real(kind=8), intent(in) :: phi
        integer, intent(in) :: jj, nflag
! Local: 
        integer :: ii

        call this%get_exp(phi, nflag)
        do ii = 1, Ndim
            Mat(ii,jj) = Mat(ii,jj) * this%expalpha
        enddo
        return
    end subroutine opU_mmult_L
end module OperatorHubbard_mod