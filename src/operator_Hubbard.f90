module OperatorHubbard_mod
    use MyLattice
    implicit none
    public
    
    type :: AccCounter
        real(kind=8), private :: NC_eff_up, ACC_eff_up
        real(kind=8), public :: acc
    contains
        procedure :: init  => Acc_init
        procedure :: reset => Acc_reset
        procedure :: count => Acc_count
        procedure :: ratio => Acc_calc_ratio
    end type AccCounter
    
    type :: OperatorHubbard
        complex(kind=8), private :: alpha
        complex(kind=8), private :: gaussian
        complex(kind=8), private :: expalpha

        complex(kind=8), public :: Delta
        complex(kind=8), public :: ratio_gaussian
        type(AccCounter), public :: Acc_U_local, Acc_U_therm
    contains
        procedure           :: set          => opU_set
        procedure, private  :: get_exp      => opU_get_exp
        procedure           :: get_delta    => opU_get_delta
        procedure           :: mmult_R      => opU_mmult_R
        procedure           :: mmult_L      => opU_mmult_L
    end type OperatorHubbard
    
contains
    real(kind=8) function hs_eta(label) result(eta)
        real(kind=8), intent(in) :: label
        select case (nint(label))
        case (-1)
            eta = -sqrt(2.d0 * (3.d0 - sqrt(6.d0)))
        case (1)
            eta = sqrt(2.d0 * (3.d0 - sqrt(6.d0)))
        case (-2)
            eta = -sqrt(2.d0 * (3.d0 + sqrt(6.d0)))
        case (2)
            eta = sqrt(2.d0 * (3.d0 + sqrt(6.d0)))
        case default
            write(6,*) 'illegal discrete HS label for eta:', label, 'rank', IRANK
            stop
        end select
        return
    end function hs_eta

    real(kind=8) function hs_gamma(label) result(gamma)
        real(kind=8), intent(in) :: label
        select case (abs(nint(label)))
        case (1)
            gamma = 1.d0 + sqrt(6.d0) / 3.d0
        case (2)
            gamma = 1.d0 - sqrt(6.d0) / 3.d0
        case default
            write(6,*) 'illegal discrete HS label for gamma:', label, 'rank', IRANK
            stop
        end select
        return
    end function hs_gamma

    subroutine opU_set(this, RU)
        class(OperatorHubbard), intent(inout) :: this
        real(kind=8), intent(in) :: RU
        this%alpha = dcmplx( 0.d0, 0.d0 )
        if ( RU < -Zero ) this%alpha = dcmplx( sqrt(-RU * Dtau), 0.d0 )
        if ( RU >  Zero ) this%alpha = dcmplx( 0.d0, sqrt( RU * Dtau) )
        return
    end subroutine opU_set
    
    subroutine opU_get_exp(this, phi, nflag)
        class(OperatorHubbard), intent(inout) :: this
        integer, intent(in) :: nflag ! +1 or -1; propagating direction
        real(kind=8), intent(in) :: phi ! space time local auxiliary field value, for phi_1 or phi_2
        this%gaussian = dcmplx(hs_gamma(phi), 0.d0)
        this%expalpha = exp( this%alpha * hs_eta(phi) * nflag )
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
        class(triangularLattice), intent(in) :: Latt
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
        class(triangularLattice), intent(in) :: Latt
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
    
    subroutine Acc_init(this)
        class(AccCounter), intent(inout) :: this
        this%acc = 0.d0
        return
    end subroutine Acc_init
    
    subroutine Acc_reset(this)
        class(AccCounter), intent(inout) :: this
        this%NC_eff_up = 0.d0
        this%ACC_eff_up = 0.d0
        return
    end subroutine Acc_reset
    
    subroutine Acc_count(this, toggle)
        class(AccCounter), intent(inout) :: this
        logical, intent(in) :: toggle
        this%NC_eff_up = this%NC_eff_up + 1
        if (toggle) this%ACC_eff_up = this%ACC_eff_up + 1
        return
    end subroutine Acc_count
    
    subroutine Acc_calc_ratio(this)
        class(AccCounter), intent(inout) :: this
        this%acc = this%acc + this%ACC_eff_up / this%NC_eff_up
        return
    end subroutine Acc_calc_ratio
end module OperatorHubbard_mod
