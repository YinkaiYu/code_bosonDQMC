module ObserEqual_mod
    use ProcessMatrix
    use DQMC_Model_mod
    implicit none
    
    type, public :: ObserEqual
        complex(kind=8), dimension(:,:,:), allocatable  :: den_corr_up, den_corr_do
        real(kind=8)                                    :: density_up,  density_do
    contains
        procedure :: make   => Obs_equal_make
        procedure :: reset  => Obs_equal_reset
        procedure :: ave    => Obs_equal_ave
        procedure :: calc   => Obs_equal_calc
        final :: Obs_equal_clear
    end type ObserEqual
    
contains
    subroutine Obs_equal_make(this)
        class(ObserEqual), intent(inout) :: this
        allocate( this%den_corr_up(Lq, Norb, Norb), this%den_corr_do(Lq, Norb, Norb) )
        return
    end subroutine Obs_equal_make
    
    subroutine Obs_equal_clear(this)
        type(ObserEqual), intent(inout) :: this
        deallocate( this%den_corr_up, this%den_corr_do )
        return
    end subroutine Obs_equal_clear
    
    subroutine Obs_equal_reset(this)
        class(ObserEqual), intent(inout) :: this
        this%den_corr_up = dcmplx(0.d0,0.d0)
        this%den_corr_do = dcmplx(0.d0,0.d0)
        this%density_up  = 0.d0
        this%density_do  = 0.d0
        return
    end subroutine Obs_equal_reset
    
    subroutine Obs_equal_ave(this, Nobs)
        class(ObserEqual), intent(inout) :: this
        integer, intent(in) :: Nobs
        real(kind=8) :: znorm
        znorm = 1.d0 / dble(Nobs)
        this%den_corr_up = this%den_corr_up * znorm
        this%den_corr_do = this%den_corr_do * znorm
        this%density_up  = this%density_up  * znorm
        this%density_do  = this%density_do  * znorm
        return
    end subroutine Obs_equal_ave
    
    subroutine Obs_equal_calc(this, Prop, ntau)
!   Arguments: 
        class(ObserEqual), intent(inout) :: this
        class(Propagator), intent(in) :: Prop
        integer, intent(in) :: ntau
! Local: 
        complex(kind=8), dimension(Ndim, Ndim) :: Grupc, Grup
        complex(kind=8), dimension(Ndim, Ndim) :: Grdoc, Grdo
        integer :: i, j, no1, no2, ii, jj, imj
        
        Grup    = Prop%Gr                           !   Gr(i, j)    = <b_i b^+_j >
        Grupc   = transpose(Grup) - ZKRON           !   Grc(i, j)   = <b^+_i b_j > = <b_i b^+_j > - δ(i,j)
        
        Grdo    = dconjg(Prop%Gr)                   !   Gr(i, j)    = <c_i c^+_j >
        Grdoc   = dconjg(transpose(Grdo)) - ZKRON   !   Grc(i, j)   = <c^+_i c_j > = <c_i c^+_j > - δ(i,j)
        
        do ii = 1, Ndim
            this%density_up = this%density_up + real(Grupc(ii,ii)) / dble(Lq)
            this%density_do = this%density_do + real(Grdoc(ii,ii)) / dble(Lq)
        enddo

        do i = 1, Lq
            do j = 1, Lq
                imj = Latt%imj(i, j)
                do no1 = 1, Norb
                    do no2 = 1, Norb
                        ii = Latt%inv_dim_list(i, no1)
                        jj = Latt%inv_dim_list(j, no2)
                        this%den_corr_up(imj, no1, no2) = this%den_corr_up(imj, no1, no2) + ( Grupc(ii,ii) * Grupc(jj,jj) + Grupc(ii,jj) * Grup(jj,jj) ) / dcmplx(dble(Lq), 0.d0)
                        this%den_corr_do(imj, no1, no2) = this%den_corr_do(imj, no1, no2) + ( Grdoc(ii,ii) * Grdoc(jj,jj) + Grdoc(ii,jj) * Grdo(jj,jj) ) / dcmplx(dble(Lq), 0.d0)
                    enddo
                enddo
            enddo
        enddo

        return
    end subroutine Obs_equal_calc
end module ObserEqual_mod