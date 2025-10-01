!=====================================================================
!                       MODULE HFBCS
!---------------------------------------------------------------------
! Solve Hartree-Fock and BCS equations.
! Skyrme interaction in the mean-field channel.
! Density-dependent delta interaction in the T=1 pairing channel.
! Axial symmetry and time-reversal symmetry imposed.
! Expansion of single-particle wave functions in the cylindrical 
! harmonic-oscillator basis.
!=====================================================================

module hfbcs


!
! Required modules
!

  use angular_momentum
  use database
  use htda
  use io
  use mb_operators
  use tools

!
! Global variables
!

  implicit none
  character(len=5),parameter :: modulename='hfbcs'

!
! Scope control
!

  private
  public :: solve_HFBCS_equations_Broyden, solve_HFBCS_equations, solve_HFHTDA_equations, &
       observables, f, local_one_body_densities, moments, single_particle_sum_rules

  public :: pairing_matrix_elements, BCS, potential, renormalize_densities, Coulomb, energy


contains


  !=====================================================================
  !                SUBROUTINE SOLVE_HFBCS_EQUATIONS_BROYDEN
  !---------------------------------------------------------------------
  ! Solve the Hartree-Fock and BCS equations for a given nucleus, within 
  ! the model given as 2nd argument, restricted to symmetries passed as 
  ! argument 3, by expansion of single-particle states in the given 
  ! basis (argument 4) and with the numerical parameters embedded in 
  ! argument 5.
  ! Iterative resolution of Hartree-Fock and BCS equations using the 
  ! diagonalization method
  ! The solution is contained in the array 'sp_states' (argument 6)
  ! which is of derived type 'sp_wave_function' (embedding all relevant 
  ! properties of the sp states (energy, occupation probability, quantum 
  ! numbers, ...). Pairing properties are embedded in 'nucleus'.
  !=====================================================================

  subroutine solve_HFBCS_equations_Broyden(nucleus, model, symmetries, basis, &
       numerical_parameters, HFfield, sp_states, densities, Coulomb_exchange, &
       Coulomb_exchange_potential)

!
!   Arguments
!
    type(nucleus_properties),intent(inout)              :: nucleus
    type(modelization),intent(inout)                    :: model
    type(symmetry_properties),intent(in)                :: symmetries
    type(cylindrical_basis),intent(inout)               :: basis
    type(numeric) ,intent(inout)                        :: numerical_parameters
    type(field),dimension(2),intent(inout)              :: HFfield
    type(sp_state),dimension(basis%size,2),intent(inout):: sp_states
    type(local_density),dimension(2),intent(inout)      :: densities
    character(len=*),optional,intent(in)                :: Coulomb_exchange
    type(operator),dimension(:),allocatable,optional,intent(inout) :: Coulomb_exchange_potential
!
!   Local variables
!
    character(len=29),parameter :: subroutinename='solve_HFBCS_equations_Broyden'
    integer, dimension(8) :: start, end
    integer :: Nblocks
    integer :: mem_stat
    integer :: iter,isospin
    logical :: converged,finished
    logical :: exact_coul_exch
    real(kind=rk) :: mix

!   Broyden: begin
    integer,parameter :: m=7, nfields=25, itermax=5*(m+3)
    integer :: i,j,k,l
    integer :: NGz,NGr,n,Nitermax
    real(kind=rk),parameter :: alpha=0.5, w0=0.01
    real(kind=rk) :: tol
    real(kind=rk),dimension(m) :: w
    real(kind=rk),dimension(m+1) :: chi
    real(kind=rk),dimension(:,:),allocatable :: V, FPF
    real(kind=rk),dimension(m,m) :: a, beta
    real(kind=rk),dimension(:),allocatable :: tmp, DeltaF, DeltaF_2
    real(kind=rk) :: denominator, denominator1, denominator2, error
    integer :: V_length
!   Broyden: end

!
!   Header message
!
		 
    write(out_unit,'(/80a1)') ('=',i=1,80)
    do i = 1, 3
       write(out_unit,'(a1,78x,a1)') '|','|'
    end do
    write(out_unit,'(a1,5x,a,47x,a1)') '|','HFBCS  I T E R A T I O N S','|'
    do i = 1, 2
       write(out_unit,'(a1,78x,a1)') '|','|'
    end do
    write(out_unit,'(a1,5x,a,46x,a1)') '|','B R O Y D E N   M I X I N G','|'
    do i = 1, 3
       write(out_unit,'(a1,78x,a1)') '|','|'
    end do
    write(out_unit,'(80a1/)') ('=',i=1,80)

    call DATE_AND_TIME(values = start)

!
!   Auxiliary variables
!

    mix=numerical_parameters%convergence%mixing
    Nitermax=numerical_parameters%convergence%Nitermax

    if (.NOT.PRESENT(Coulomb_exchange)) then
       if (debugmode/='low') write(out_unit,'(/a)') &
            '* Treatment of Coulomb exchange: by the Slater approximation'
       exact_coul_exch = .FALSE.
    else
       SELECT CASE (TRIM(Coulomb_exchange))
       CASE('exact')
          write(out_unit,'(/a)') '* Treatment of Coulomb exchange: Exactly'
          exact_coul_exch = .TRUE.
       CASE('slater')
          write(out_unit,'(/a)') '* Treatment of Coulomb exchange: by the Slater approximation'
          exact_coul_exch = .FALSE.
       CASE DEFAULT
          if (debugmode/='low') write(out_unit,'(/a)') &
               '* Treatment of Coulomb exchange: by the Slater approximation'
          exact_coul_exch = .FALSE.
       END SELECT
    end if

    if (nucleus%fix_rank=='n' .and. debugmode=='high') then
       write(out_unit,*)
       write(out_unit,'(a)') '*** Rank of K^pi state not fixed'
       write(out_unit,'(a)') '    Find the K^pi state nearer to Fermi level'
       write(out_unit,*)
    end if

!----------------------
!
!   Memory allocation
!
!----------------------

    NGz = numerical_parameters%integration%NGz
    NGr = numerical_parameters%integration%NGr
    if (symmetries%parity) then
       n = 2 * nfields * NGz/2 * NGr ! total number of field components 
                                     ! (isospin x (4 scalar + 4 vector + 1 tensor fields) x 
                                     ! number of Gauss points)
    else
       n = 2 * nfields * (NGz+1) * NGr
    end if
    if (model%exact_coul_exch) then
       V_length = sum((/(basis%block_size(k)*(basis%block_size(k)+1)/2,k=1,basis%Nblocks)/))
       n = n + V_length
       if (debugmode /= 'low') print*,'V_length=',V_length,' n=',n
    end if
    tol = numerical_parameters%convergence%tolerance
    if (debugmode /= 'low') write(out_unit,'(a,i6,1x/a,e13.6)') &
         "Length of Broyden vectors n=",n,"Tolerance on vector components tol=",tol
    allocate(V(m+2,n), FPF(m+2,n))
    allocate(tmp(n), DeltaF(n), DeltaF_2(n))

    w(:) = 1.0

    V(:,:) = 0.0
    FPF(:,:) = 0.0

!----------------------------------------------------------------------
!
!    "Pre-convergence" using the linear-mixing method
!
!----------------------------------------------------------------------

    write(out_unit,'(/82a1)') ('#',i=1,82)
    write(out_unit,'(//a,i3,a,f4.2/)') '"Pre-convergence" using the linear-mixing method ' // &
         'for',Nitermax,' iterations with mix = ',mix
    write(out_unit,'(/82a1)') ('#',i=1,82)

    iter = 0

    if (model%initial_potential /= 'HF') then
       numerical_parameters%convergence%Nitermax=Nitermax
       numerical_parameters%convergence%mixing=mix
       if (model%exact_coul_exch) then
          call solve_HFBCS_equations(nucleus, model, symmetries, basis, &
               numerical_parameters, HFfield, sp_states, densities, iter, &
               Coulomb_exchange='exact',Coulomb_exchange_potential= &
               Coulomb_exchange_potential,print_header_message=.true.)
       else
          call solve_HFBCS_equations(nucleus, model, symmetries, basis, &
               numerical_parameters, HFfield, sp_states, densities, iter, &
               print_header_message=.true.)
       end if
    end if

!----------------------------------------------------------------------
!
!    Iteration loop (i represents m and j represents n in Baran08, PRC78)
!
!----------------------------------------------------------------------

    numerical_parameters%convergence%mixing=0.0
    numerical_parameters%convergence%Nitermax=1
    model%initial_potential = 'HF'

    write(out_unit,'(/80a1)') ('#',i=1,80)
    write(out_unit,'(//a/)') 'Starting Broyden iteration loop'
    write(out_unit,'(/80a1)') ('#',i=1,80)

!
!   Iteration 1
!

    call HFfield_into_V(HFfield,V(1,:),n,NGz,NGr) ! Initial conditions

    if (model%exact_coul_exch) then
       call solve_HFBCS_equations(nucleus, model, symmetries, basis, &
            numerical_parameters, HFfield, sp_states, densities, iter, &
            Coulomb_exchange='exact',Coulomb_exchange_potential= &
            Coulomb_exchange_potential,print_header_message=.false.)
    else
       call solve_HFBCS_equations(nucleus, model, symmetries, basis, &
            numerical_parameters, HFfield, sp_states, densities, iter, &
            print_header_message=.false.)
    end if
    call HFfield_into_V(HFfield,tmp(:),n,NGz,NGr)
    FPF(1,:) = tmp(:) - V(1,:)
    V(1,:) = V(1,:) + alpha * FPF(1,:)

    write(out_unit,'(/a)') 'Calculation 1 in Broyden iteration 1'

    call V_into_HFfield(V(1,:),n,HFfield,NGz,NGr)
    if (model%exact_coul_exch) then
       call solve_HFBCS_equations(nucleus, model, symmetries, basis, &
            numerical_parameters, HFfield, sp_states, densities, iter, &
            Coulomb_exchange='exact',Coulomb_exchange_potential= &
            Coulomb_exchange_potential,print_header_message=.false.)
    else
       call solve_HFBCS_equations(nucleus, model, symmetries, basis, &
            numerical_parameters, HFfield, sp_states, densities, iter, &
            print_header_message=.false.)
    end if
    call HFfield_into_V(HFfield,tmp(:),n,NGz,NGr)
    FPF(1,:) = tmp(:) - V(1,:)

    error = sqrt(dot_product(FPF(1,:),FPF(1,:)))

    write(out_unit,'(/a,1x,i3,2x,a,f25.15)') 'Calculation 2 in Broyden iteration',1,'error=',error

    if (error <= tol) then

       write(out_unit,'(/a/)') 'Solution found'

    else

!
!      Iteration 2
!
       V(2,:) = V(1,:) + alpha * FPF(1,:)

       call V_into_HFfield(V(2,:),n,HFfield,NGz,NGr)
       if (model%exact_coul_exch) then
          call solve_HFBCS_equations(nucleus, model, symmetries, basis, &
               numerical_parameters, HFfield, sp_states, densities, iter, &
               Coulomb_exchange='exact',Coulomb_exchange_potential= &
               Coulomb_exchange_potential,print_header_message=.false.)
       else
          call solve_HFBCS_equations(nucleus, model, symmetries, basis, &
               numerical_parameters, HFfield, sp_states, densities, iter, &
               print_header_message=.false.)
       end if
       call HFfield_into_V(HFfield,tmp(:),n,NGz,NGr)
       FPF(2,:) = tmp(:) - V(2,:)
       error = sqrt(dot_product(FPF(2,:),FPF(2,:)))

       write(out_unit,'(/a,1x,i3,2x,a,f25.15)') 'Broyden iteration',2,'error=',error

       if (error <= tol) then

          write(out_unit,'(/a/)') 'Solution found'

       else

!
!          Iterations 3 to m + 2
!

          first: do i = 3, m + 2

             chi(:) = 0.0

!            Matrix a (real symmetric)

             a(:,:) = 0.0
             do k = 1, i-2
!               Diagonal matrix element
                a(k,k) = w(k)**2 + w0**2
!               Off-diagonal matrix elements
                DeltaF_2 = FPF(k+1,:) - FPF(k,:)
                 denominator = sqrt(dot_product(DeltaF_2,DeltaF_2))
                 if (denominator <= 0) cycle
                 DeltaF_2 = DeltaF_2/denominator
                 do l = 1, k-1
                    DeltaF = FPF(l+1,:) - FPF(l,:)
                    denominator = sqrt(dot_product(DeltaF,DeltaF))
                    if (denominator <= 0) cycle
                    DeltaF = DeltaF/denominator
                    a(k,l) = w(k) * w(l) * dot_product(DeltaF, DeltaF_2)
                    a(l,k) = a(k,l)
                 end do
              end do

!             Matrix beta (involves matrix inversion)

              beta(:,:) = 0.0
              call inverse(a(1:i-2,1:i-2),beta(1:i-2,1:i-2),i-2)

!
!       Coefficients chi(j) for 1 <= j <= i-1  and 3 <= i <= m+2
!
!               /  gamma(i-1,1)
!              | ---------------                           if j = 1
!              |  ||F(2)-F(1)||
!              |
!              |   gamma(i-1,j)          gamma(i-1,j-1)
!    chi(j) = -| -----------------  -  -----------------   if 2 <= j <= i - 2
!              |  ||F(j+1)-F(j)||       ||F(j)-F(j-1)||
!              |
!              |         gamma(i-1,i-2)
!              | 1  -  -------------------                 if j = i - 1
!               \       ||F(i-1)-F(i-2)||
!
!    with
!                     i-2
!                    ----
!                    \
!    gamma(i-1,j) =  /    beta(k,j) DeltaF(k) . F(i-1)     for 1 <= j <= i-2
!                    ----
!                     k=1

!             j = 1

              tmp(:) = 0.0
              do k = 1, i-2 
                 DeltaF(:) = FPF(k+1,:) - FPF(k,:)
                 denominator = dot_product(DeltaF,DeltaF)
                 if (denominator <= 0) cycle
                 DeltaF(:) = DeltaF/sqrt(denominator)
                 tmp(:) = tmp(:) + beta(k,1) * DeltaF(:)
              end do
              DeltaF(:) = FPF(2,:) - FPF(1,:)
              denominator = dot_product(DeltaF,DeltaF)
              if (denominator <= 0) stop 'denominator <= 0 for j = 1'
              chi(1) = dot_product(tmp,FPF(i-1,:))/sqrt(denominator)

!             2 <= j <= i-2

              do j = 2, i-2
                 tmp(:) = 0.0
                 DeltaF(:) = FPF(j+1,:) - FPF(j,:)
                 denominator1 = dot_product(DeltaF,DeltaF)
                 DeltaF(:) = FPF(j,:) - FPF(j-1,:)
                 denominator2 = dot_product(DeltaF,DeltaF)
                 if (denominator1 <= 0 .or. denominator2 <= 0) cycle
                 do k = 1, i-2
                    DeltaF(:) = FPF(k+1,:) - FPF(k,:)
                    denominator = dot_product(DeltaF,DeltaF)
                    if (denominator <= 0) cycle
                    DeltaF(:) = DeltaF/sqrt(denominator)
                    tmp(:) = tmp(:) + (beta(k,j)/sqrt(denominator1) - &
                         beta(k,j-1)/sqrt(denominator2)) * DeltaF(:)
                 end do
                 chi(j) = dot_product(tmp,FPF(i-1,:))
              end do

!             j = i-1

              tmp(:) = 0.0
              do k = 1, i-2
                 DeltaF(:) = FPF(k+1,:) - FPF(k,:)
                 denominator = dot_product(DeltaF,DeltaF)
                 if (denominator <= 0) cycle
                 DeltaF(:) = DeltaF/sqrt(denominator)
                 tmp(:) = tmp(:) + beta(k,i-2) * DeltaF
              end do
              DeltaF(:) = FPF(i-1,:) - FPF(i-2,:)
              denominator = dot_product(DeltaF,DeltaF)
              if (denominator <= 0) then
                 write(out_unit,'(a,1x,i3,2(2x,a,f25.15))') 'iteration',2, &
                      'denominator =',denominator,'error=',error
                 call fatal(modulename,subroutinename,'denominator <= 0 for k = 0')
              end if
              chi(i-1) = 1.0 - dot_product(tmp,FPF(i-1,:))/sqrt(denominator)

!
!       Next iteration solution:
!
!            i-1
!           ----
!           \
!    v(i) = /    chi(j) * (v(j) + alpha * FPF(j))
!           ----
!            j=1
!

              V(i,:) = 0.0
              do j = 1, i-1
                 V(i,:) = V(i,:) + chi(j)*(V(j,:) + alpha * FPF(j,:))
              end do

              call V_into_HFfield(V(i,:),n,HFfield,NGz,NGr)
              if (model%exact_coul_exch) then
                 call solve_HFBCS_equations(nucleus,model,symmetries,basis, &
                      numerical_parameters,HFfield,sp_states,densities,iter, &
                      Coulomb_exchange='exact',Coulomb_exchange_potential= &
                      Coulomb_exchange_potential,print_header_message=.false.)
              else
                 call solve_HFBCS_equations(nucleus,model,symmetries,basis, &
                      numerical_parameters,HFfield,sp_states,densities,iter, &
                      print_header_message=.false.)
              end if
              call HFfield_into_V(HFfield,tmp(:),n,NGz,NGr)
              FPF(i,:) = tmp(:) - V(i,:)

              error = sqrt(dot_product(FPF(i,:),FPF(i,:)))
              write(out_unit,'(/a,1x,i3,2x,a,f25.15)') 'Broyden iteration',i,'error=',error

              if (error <= tol) then
                 write(out_unit,'(/a/)') 'Solution found'
                 exit first
              end if
     
           end do first

!
!          Iterations i >= m + 3
!

           if (error > tol) then

              error = tol + 1.0
              i = m + 3

              call V_into_HFfield(V(m+2,:),n,HFfield,NGz,NGr)
              if (model%exact_coul_exch) then
                 call solve_HFBCS_equations(nucleus,model,symmetries,basis, &
                      numerical_parameters,HFfield,sp_states,densities,iter, &
                      Coulomb_exchange='exact',Coulomb_exchange_potential= &
                      Coulomb_exchange_potential)
              else
                 call solve_HFBCS_equations(nucleus,model,symmetries,basis, &
                      numerical_parameters,HFfield,sp_states,densities,iter, &
                      print_header_message=.false.)
              end if
              call HFfield_into_V(HFfield,tmp(:),n,NGz,NGr)
              FPF(m+2,:) = tmp(:) - V(m+2,:)

recursion:    do while (error > tol)

                 do j = 1, m+1
                    V(j,:) = V(j+1,:)
                    FPF(j,:) = FPF(j+1,:)
                 end do

!                Matrix a (real symmetric)

                 a(:,:) = 0.0
                 do k = 1, m
!                   Diagonal matrix element
                    a(k,k) = w(k)**2 + w0**2
!                   Off-diagonal matrix elements
                    DeltaF_2 = FPF(k+1,:) - FPF(k,:)
                    denominator = sqrt(dot_product(DeltaF_2,DeltaF_2))
                    if (denominator <= 0) cycle
                    DeltaF_2 = DeltaF_2/denominator
                    do l = 1, k-1
                       DeltaF = FPF(l+1,:) - FPF(l,:)
                       denominator = sqrt(dot_product(DeltaF,DeltaF))
                       if (denominator <= 0) cycle
                       DeltaF = DeltaF/denominator
                       a(k,l) = w(k) * w(l) * dot_product(DeltaF, DeltaF_2)
                       a(l,k) = a(k,l)
                    end do
                 end do

!                Matrix beta (involves matrix inversion)

                 beta(:,:) = 0.0
                 call inverse(a(1:m,1:m),beta(1:m,1:m),m)

!
!       Coefficients chi(j) for i - m - 1 <= j <= i - 1
!
!               /  gamma(i-1,i-m-1)
!              | --------------------                      if j = i - m - 1
!              |  ||F(i-m)-F(i-m-1)||
!              |
!              |   gamma(i-1,j)          gamma(i-1,j-1)
!    chi(j) = -| -----------------  -  -----------------   if i - m <= j <= i - 2
!              |  ||F(j+1)-F(j)||       ||F(j)-F(j-1)||
!              |
!              |         gamma(i-1,i-2)
!              | 1  -  -------------------                 if j = i - 1
!               \       ||F(i-1)-F(i-2)||
!
!    with
!                    i-2
!                    ----
!                    \
!    gamma(i-1,j) =  /    beta(k,j) DeltaF(k) . F(i-1)     for i - m - 1 <= j <= i - 2
!                    ----
!                   k=i-m-1

!                j = i - m - 1

                 tmp(:) = 0.0
                 do k = 1, m
                    DeltaF(:) = FPF(k+1,:) - FPF(k,:)
                    denominator = dot_product(DeltaF,DeltaF)
                    if (denominator <= 0) cycle
                    DeltaF(:) = DeltaF/sqrt(denominator)
                    tmp(:) = tmp(:) + beta(k,1) * DeltaF(:)
                 end do
                 DeltaF(:) = FPF(2,:) - FPF(1,:)
                 denominator = dot_product(DeltaF,DeltaF)
                 if (denominator <= 0) stop 'denominator <= 0 for j = 1'
                 chi(1) = dot_product(tmp,FPF(m+1,:))/sqrt(denominator)

!                i - m <= j <= i - 2

                 do j = 2, m
                    tmp(:) = 0.0
                    DeltaF(:) = FPF(j+1,:) - FPF(j,:)
                    denominator1 = dot_product(DeltaF,DeltaF)
                    DeltaF(:) = FPF(j,:) - FPF(j-1,:)
                    denominator2 = dot_product(DeltaF,DeltaF)
                    if (denominator1 <= 0 .or. denominator2 <= 0) cycle
                    do k = 1, m
                       DeltaF(:) = FPF(k+1,:) - FPF(k,:)
                       denominator = dot_product(DeltaF,DeltaF)
                       if (denominator <= 0) cycle
                       DeltaF(:) = DeltaF/sqrt(denominator)
                       tmp(:) = tmp(:) + (beta(k,j)/sqrt(denominator1) - &
                            beta(k,j-1)/sqrt(denominator2)) * DeltaF(:)
                    end do
                    chi(j) = dot_product(tmp,FPF(m+1,:))
                 end do

!                j = i - 1

                 tmp(:) = 0.0
                 do k = 1, m
                    DeltaF(:) = FPF(k+1,:) - FPF(k,:)
                    denominator = dot_product(DeltaF,DeltaF)
                    if (denominator <= 0) cycle
                    DeltaF(:) = DeltaF/sqrt(denominator)
                    tmp(:) = tmp(:) + beta(k,m) * DeltaF
                 end do
                 DeltaF(:) = FPF(m+1,:) - FPF(m,:)
                 denominator = dot_product(DeltaF,DeltaF)
                 if (denominator <= 0) then
                    write(out_unit,'(a,1x,i3,2(2x,a,f25.15))') 'iteration',2, &
                         'denominator =',denominator,'error=',error
                    call fatal(modulename,subroutinename,'denominator <= 0 for k = 0')
                 end if
                 chi(m+1) = 1.0 - dot_product(tmp,FPF(m+1,:))/sqrt(denominator)

!
!       Next iteration solution:
!
!            i-1
!            ----
!            \
!    v(i) =  /    chi(j) * (v(j) + alpha * FPF(j))
!            ----
!           j=i-m-1
!

                 V(m+2,:) = 0.0
                 do j = 1, m+1
                    V(m+2,:) = V(m+2,:) + chi(j)*(V(j,:) + alpha * FPF(j,:))
                 end do

                 call V_into_HFfield(V(m+2,:),n,HFfield,NGz,NGr)
                 if (model%exact_coul_exch) then
                    call solve_HFBCS_equations(nucleus,model,symmetries,basis, &
                         numerical_parameters,HFfield,sp_states,densities,iter,&
                         Coulomb_exchange='exact',Coulomb_exchange_potential= &
                         Coulomb_exchange_potential,print_header_message=.false.)
                 else
                    call solve_HFBCS_equations(nucleus,model,symmetries,basis, &
                         numerical_parameters, HFfield, sp_states, densities, &
                         iter,print_header_message=.false.)
                 end if
                 call HFfield_into_V(HFfield,tmp(:),n,NGz,NGr)
                 FPF(m+2,:) = tmp(:) - V(m+2,:)

                 tmp(:) = FPF(m+2,:)
                 error = sqrt(dot_product(tmp,tmp))
                 write(out_unit,'(/a,1x,i3,2x,a,f25.15)') 'Broyden iteration',i,'error=',error

                 if (error <= tol) then
                    write(out_unit,'(/a/)') 'Solution found'
                    exit recursion
                 end if

                 if (i >= itermax) then
                    if (model%exact_coul_exch) then
                       call output(nucleus, model, symmetries, basis, numerical_parameters, &
                            HFfield, densities, Coulomb_exchange_potential=Coulomb_exchange_potential)
                    else
                       call output(nucleus, model, symmetries, basis, numerical_parameters, &
                            HFfield, densities)
                    end if
                    call fatal(modulename,subroutinename,'Solution not found so far')
                 end if

                 i = i + 1

              end do recursion

           end if

        end if

     end if

    numerical_parameters%convergence%Nitermax = Nitermax

!---------------------------------------------
!
!   Memory deallocation for local variables
!
!---------------------------------------------

    deallocate(V, FPF, DeltaF, DeltaF_2, tmp)

    call DATE_AND_TIME(values = end)
    call elapsed_time(start,end,"Broyden iterations")

    return

  contains

    subroutine HFfield_into_V(HFfield,V,n,NGz,NGr)

!
!   Arguments
!
      implicit none
      type(field),dimension(2),intent(in) :: HFfield
      integer,intent(in) :: n,NGz,NGr
      real(kind=rk),dimension(n),intent(inout) :: V

!
!   Local variables
!
      integer :: isospin, shift, NGtot, Nz, expected_n
      integer :: V_length, i, j, k, l, block_size, first, last

      if (symmetries%parity) then
         Nz = NGz/2
      else
         Nz = NGz + 1
      end if
      NGtot = Nz*NGr
      expected_n = 2*nfields*NGtot ! nfields = 25 (see last line in the loop over isospin)
      if (model%exact_coul_exch) then
         V_length = sum((/(basis%block_size(k)*(basis%block_size(k)+1)/2,k=1,basis%Nblocks)/))
         expected_n = expected_n + V_length
      end if
      if (n /= expected_n) then
         print*,'n=',n,' 2*nfields*NGtot=',2*nfields*NGtot,' expected_n=',expected_n
         stop "Size problem in HFfield_into_V"
      end if

      do isospin = 1, 2
         shift = nfields*(isospin-1)
         V(shift*NGtot+1:(1+shift)*NGtot) = reshape(HFfield(isospin)%central(iHmin:NGz/2,1:NGr),(/NGtot/))
         V((1+shift)*NGtot+1:(2+shift)*NGtot) = reshape(HFfield(isospin)%spin_orbit(iHmin:NGz/2,1:NGr),(/NGtot/))
         V((2+shift)*NGtot+1:(3+shift)*NGtot) = reshape(HFfield(isospin)%effective_mass(iHmin:NGz/2,1:NGr),(/NGtot/))
         V((3+shift)*NGtot+1:(4+shift)*NGtot) = reshape(HFfield(isospin)%S%z(iHmin:NGz/2,1:NGr),(/NGtot/))
         V((4+shift)*NGtot+1:(5+shift)*NGtot) = reshape(HFfield(isospin)%S%rho(iHmin:NGz/2,1:NGr),(/NGtot/))
         V((5+shift)*NGtot+1:(6+shift)*NGtot) = reshape(HFfield(isospin)%S%phi(iHmin:NGz/2,1:NGr),(/NGtot/))
         V((6+shift)*NGtot+1:(7+shift)*NGtot) = reshape(HFfield(isospin)%div_s(iHmin:NGz/2,1:NGr),(/NGtot/))
         V((7+shift)*NGtot+1:(8+shift)*NGtot) = reshape(HFfield(isospin)%alpha%z(iHmin:NGz/2,1:NGr),(/NGtot/))
         V((8+shift)*NGtot+1:(9+shift)*NGtot) = reshape(HFfield(isospin)%alpha%rho(iHmin:NGz/2,1:NGr),(/NGtot/))
         V((9+shift)*NGtot+1:(10+shift)*NGtot) = reshape(HFfield(isospin)%alpha%phi(iHmin:NGz/2,1:NGr),(/NGtot/))
         V((10+shift)*NGtot+1:(11+shift)*NGtot) = reshape(HFfield(isospin)%C%z(iHmin:NGz/2,1:NGr),(/NGtot/))
         V((11+shift)*NGtot+1:(12+shift)*NGtot) = reshape(HFfield(isospin)%C%rho(iHmin:NGz/2,1:NGr),(/NGtot/))
         V((12+shift)*NGtot+1:(13+shift)*NGtot) = reshape(HFfield(isospin)%C%phi(iHmin:NGz/2,1:NGr),(/NGtot/))
         V((13+shift)*NGtot+1:(14+shift)*NGtot) = reshape(HFfield(isospin)%tensor_so(1)%z(iHmin:NGz/2,1:NGr),(/NGtot/))
         V((14+shift)*NGtot+1:(15+shift)*NGtot) = reshape(HFfield(isospin)%tensor_so(1)%rho(iHmin:NGz/2,1:NGr),(/NGtot/))
         V((15+shift)*NGtot+1:(16+shift)*NGtot) = reshape(HFfield(isospin)%tensor_so(1)%phi(iHmin:NGz/2,1:NGr),(/NGtot/))
         V((16+shift)*NGtot+1:(17+shift)*NGtot) = reshape(HFfield(isospin)%tensor_so(2)%z(iHmin:NGz/2,1:NGr),(/NGtot/))
         V((17+shift)*NGtot+1:(18+shift)*NGtot) = reshape(HFfield(isospin)%tensor_so(2)%rho(iHmin:NGz/2,1:NGr),(/NGtot/))
         V((18+shift)*NGtot+1:(19+shift)*NGtot) = reshape(HFfield(isospin)%tensor_so(2)%phi(iHmin:NGz/2,1:NGr),(/NGtot/))
         V((19+shift)*NGtot+1:(20+shift)*NGtot) = reshape(HFfield(isospin)%tensor_so(3)%z(iHmin:NGz/2,1:NGr),(/NGtot/))
         V((20+shift)*NGtot+1:(21+shift)*NGtot) = reshape(HFfield(isospin)%tensor_so(3)%rho(iHmin:NGz/2,1:NGr),(/NGtot/))
         V((21+shift)*NGtot+1:(22+shift)*NGtot) = reshape(HFfield(isospin)%tensor_so(3)%phi(iHmin:NGz/2,1:NGr),(/NGtot/))
         V((22+shift)*NGtot+1:(23+shift)*NGtot) = reshape(HFfield(isospin)%D%z(iHmin:NGz/2,1:NGr),(/NGtot/))
         V((23+shift)*NGtot+1:(24+shift)*NGtot) = reshape(HFfield(isospin)%D%rho(iHmin:NGz/2,1:NGr),(/NGtot/))
         V((24+shift)*NGtot+1:(25+shift)*NGtot) = reshape(HFfield(isospin)%D%phi(iHmin:NGz/2,1:NGr),(/NGtot/))
      end do

      if (model%exact_coul_exch .and. allocated(Coulomb_exchange_potential)) then
         l = (25+nfields)*NGtot
         do k = 1, basis%Nblocks
            block_size=basis%block_size(k)
            first = MERGE(1,SUM(basis%block_size(1:k-1))+1,(k == 1))
            last = first + block_size - 1
            shift = first - 1
            do i = first, last
               do j = i, last
                  l = l + 1
                  V(l) = Coulomb_exchange_potential(k)%HO(i-shift,j-shift)
               end do
            end do
         end do
      end if

    end subroutine HFfield_into_V
  
    subroutine V_into_HFfield(V,n,HFfield,NGz,NGr)

!
!   Arguments
!
      implicit none
      real(kind=rk),dimension(n) :: V
      integer :: n,NGz,NGr
      type(field),dimension(2) :: HFfield

!
!   Local variables
!
      integer :: isospin, shift, NGtot, Nz, expected_n
      integer :: V_length, i, j, k, l, block_size, first, last

      if (symmetries%parity) then
         Nz = NGz/2
      else
         Nz = NGz + 1
      end if
      NGtot = Nz*NGr
      expected_n = 2*nfields*NGtot ! nfields = 25 (see last line in the loop over isospin)
      if (model%exact_coul_exch) then
         V_length = sum((/(basis%block_size(k)*(basis%block_size(k)+1)/2,k=1,basis%Nblocks)/))
         expected_n = expected_n + V_length
      end if
      if (n /= expected_n) then
         print*,'n=',n,' 2*nfields*NGtot=',2*nfields*NGtot,' expected_n=',expected_n
         stop "Size problem in V_into_HFfield"
      end if

      do isospin = 1, 2
         shift = nfields*(isospin-1)
         HFfield(isospin)%central(iHmin:NGz/2,1:NGr) = reshape(V(shift*NGtot+1:(1+shift)*NGtot),(/Nz,NGr/))
         HFfield(isospin)%spin_orbit(iHmin:NGz/2,1:NGr) = reshape(V((1+shift)*NGtot+1:(2+shift)*NGtot),(/Nz,NGr/))
         HFfield(isospin)%effective_mass(iHmin:NGz/2,1:NGr) = reshape(V((2+shift)*NGtot+1:(3+shift)*NGtot),(/Nz,NGr/))
         HFfield(isospin)%S%z(iHmin:NGz/2,1:NGr) = reshape(V((3+shift)*NGtot+1:(4+shift)*NGtot),(/Nz,NGr/))
         HFfield(isospin)%S%rho(iHmin:NGz/2,1:NGr) = reshape(V((4+shift)*NGtot+1:(5+shift)*NGtot),(/Nz,NGr/))
         HFfield(isospin)%S%phi(iHmin:NGz/2,1:NGr) = reshape(V((5+shift)*NGtot+1:(6+shift)*NGtot),(/Nz,NGr/))
         HFfield(isospin)%div_s(iHmin:NGz/2,1:NGr) = reshape(V((6+shift)*NGtot+1:(7+shift)*NGtot),(/Nz,NGr/))
         HFfield(isospin)%alpha%z(iHmin:NGz/2,1:NGr) = reshape(V((7+shift)*NGtot+1:(8+shift)*NGtot),(/Nz,NGr/))
         HFfield(isospin)%alpha%rho(iHmin:NGz/2,1:NGr) = reshape(V((8+shift)*NGtot+1:(9+shift)*NGtot),(/Nz,NGr/))
         HFfield(isospin)%alpha%phi(iHmin:NGz/2,1:NGr) = reshape(V((9+shift)*NGtot+1:(10+shift)*NGtot),(/Nz,NGr/))
         HFfield(isospin)%C%z(iHmin:NGz/2,1:NGr) = reshape(V((10+shift)*NGtot+1:(11+shift)*NGtot),(/Nz,NGr/))
         HFfield(isospin)%C%rho(iHmin:NGz/2,1:NGr) = reshape(V((11+shift)*NGtot+1:(12+shift)*NGtot),(/Nz,NGr/))
         HFfield(isospin)%C%phi(iHmin:NGz/2,1:NGr) = reshape(V((12+shift)*NGtot+1:(13+shift)*NGtot),(/Nz,NGr/))
         HFfield(isospin)%tensor_so(1)%z(iHmin:NGz/2,1:NGr) = reshape(V((13+shift)*NGtot+1:(14+shift)*NGtot),(/Nz,NGr/))
         HFfield(isospin)%tensor_so(1)%rho(iHmin:NGz/2,1:NGr) = reshape(V((14+shift)*NGtot+1:(15+shift)*NGtot),(/Nz,NGr/))
         HFfield(isospin)%tensor_so(1)%phi(iHmin:NGz/2,1:NGr) = reshape(V((15+shift)*NGtot+1:(16+shift)*NGtot),(/Nz,NGr/))
         HFfield(isospin)%tensor_so(2)%z(iHmin:NGz/2,1:NGr) = reshape(V((16+shift)*NGtot+1:(17+shift)*NGtot),(/Nz,NGr/))
         HFfield(isospin)%tensor_so(2)%rho(iHmin:NGz/2,1:NGr) = reshape(V((17+shift)*NGtot+1:(18+shift)*NGtot),(/Nz,NGr/))
         HFfield(isospin)%tensor_so(2)%phi(iHmin:NGz/2,1:NGr) = reshape(V((18+shift)*NGtot+1:(19+shift)*NGtot),(/Nz,NGr/))
         HFfield(isospin)%tensor_so(3)%z(iHmin:NGz/2,1:NGr) = reshape(V((19+shift)*NGtot+1:(20+shift)*NGtot),(/Nz,NGr/))
         HFfield(isospin)%tensor_so(3)%rho(iHmin:NGz/2,1:NGr) = reshape(V((20+shift)*NGtot+1:(21+shift)*NGtot),(/Nz,NGr/))
         HFfield(isospin)%tensor_so(3)%phi(iHmin:NGz/2,1:NGr) = reshape(V((21+shift)*NGtot+1:(22+shift)*NGtot),(/Nz,NGr/))
         HFfield(isospin)%D%z(iHmin:NGz/2,1:NGr) = reshape(V((22+shift)*NGtot+1:(23+shift)*NGtot),(/Nz,NGr/))
         HFfield(isospin)%D%rho(iHmin:NGz/2,1:NGr) = reshape(V((23+shift)*NGtot+1:(24+shift)*NGtot),(/Nz,NGr/))
         HFfield(isospin)%D%phi(iHmin:NGz/2,1:NGr) = reshape(V((24+shift)*NGtot+1:(25+shift)*NGtot),(/Nz,NGr/))
      end do

      if (model%exact_coul_exch .and. allocated(Coulomb_exchange_potential)) then
         l = (25+nfields)*NGtot
         do k = 1, basis%Nblocks
            block_size=basis%block_size(k)
            first = MERGE(1,SUM(basis%block_size(1:k-1))+1,(k == 1))
            last = first + block_size - 1
            shift = first - 1
            do i = first, last
               do j = i, last
                  l = l + 1
                  Coulomb_exchange_potential(k)%HO(i-shift,j-shift) = V(l)
               end do
            end do
         end do
      end if

    end subroutine V_into_HFfield
  
    subroutine inverse_hfbcs(a,c,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================
      implicit none 
      integer n
      double precision a(n,n), c(n,n)
      double precision L(n,n), U(n,n), b(n), d(n), x(n)
      double precision coeff
      integer i, j, k

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
      L=0.0
      U=0.0
      b=0.0

! step 1: forward elimination
      do k=1, n-1
         do i=k+1,n
            coeff=a(i,k)/a(k,k)
            L(i,k) = coeff
            do j=k+1,n
               a(i,j) = a(i,j)-coeff*a(k,j)
            end do
         end do
      end do

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
      do i=1,n
         L(i,i) = 1.0
      end do
! U matrix is the upper triangular part of A
      do j=1,n
         do i=1,j
            U(i,j) = a(i,j)
         end do
      end do

! Step 3: compute columns of the inverse matrix C
      do k=1,n
         b(k)=1.0
         d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
         do i=2,n
            d(i)=b(i)
            do j=1,i-1
               d(i) = d(i) - L(i,j)*d(j)
            end do
         end do
! Step 3b: Solve Ux=d using the back substitution
         x(n)=d(n)/U(n,n)
         do i = n-1,1,-1
            x(i) = d(i)
            do j=n,i+1,-1
               x(i)=x(i)-U(i,j)*x(j)
            end do
            x(i) = x(i)/u(i,i)
         end do
! Step 3c: fill the solutions x(n) into column k of C
         do i=1,n
            c(i,k) = x(i)
         end do
         b(k)=0.0
      end do

    end subroutine inverse_hfbcs

  end subroutine solve_HFBCS_equations_Broyden


  !=====================================================================
  !                       SUBROUTINE SOLVE_HFBCS_EQUATIONS
  !---------------------------------------------------------------------
  ! Solve the Hartree-Fock and BCS equations for a given nucleus, within 
  ! the model given as 2nd argument, restricted to symmetries passed as 
  ! argument 3, by expansion of single-particle states in the given 
  ! basis (argument 4) and with the numerical parameters embedded in 
  ! argument 5.
  ! Iterative resolution of Hartree-Fock and BCS equations using the 
  ! diagonalization method
  ! The solution is contained in the array 'sp_states' (argument 6)
  ! which is of derived type 'sp_wave_function' (embedding all relevant 
  ! properties of the sp states (energy, occupation probability, quantum 
  ! numbers, ...). Pairing properties are embedded in 'nucleus'.
  !=====================================================================

  subroutine solve_HFBCS_equations(nucleus, model, symmetries, basis, &
       numerical_parameters, HFfield, sp_states, densities, iter, &
       Coulomb_exchange, Coulomb_exchange_potential, print_header_message)

!
!   Arguments
!
    type(nucleus_properties),intent(inout)              :: nucleus
    type(modelization),intent(inout)                    :: model
    type(symmetry_properties),intent(in)                :: symmetries
    type(cylindrical_basis),intent(inout)               :: basis
    type(numeric) ,intent(in)                           :: numerical_parameters
    type(field),dimension(2),intent(inout)              :: HFfield
    type(sp_state),dimension(basis%size,2),intent(inout):: sp_states
    type(local_density),dimension(2),intent(inout)      :: densities
    integer,intent(inout)                               :: iter
    character(len=*),optional,intent(in)                :: Coulomb_exchange
    type(operator),dimension(:),allocatable,optional,intent(inout) :: Coulomb_exchange_potential
    logical,optional,intent(in) :: print_header_message
!
!   Local variables
!
    character(len=21),parameter :: subroutinename='solve_HFBCS_equations'
    integer, dimension(8) :: start, end
    integer, dimension(8) :: start_H, end_H, start_k, end_k, start_isospin, end_isospin
    integer, dimension(8) :: start_potential, end_potential, start_occupations, end_occupations
    integer, dimension(8) :: start_partners, end_partners, start_BCS, end_BCS
    integer, dimension(8) :: start_rank, end_rank
    integer, dimension(8) :: start_Nilsson, end_Nilsson
    character(len=100) :: str_time
    integer :: nb_interval, id
    integer, dimension(2) :: Npart,nb_qp
    integer :: readjustment_counter, readj_count_max
    real(kind=rk) :: error
    logical :: accelerate
    integer :: Nblocks
    integer :: k,n,shift
    integer :: mem_stat
    integer :: isospin
    integer :: i,j,iHF,jHF
    logical :: converged,finished
    logical :: exact_coul_exch
    real(kind=rk) :: mix
    real(kind=rk) :: amplitude

    type(operator),dimension(:),allocatable  :: h_HF        ! Hartree-Fock Hamiltonian
    type(operator),dimension(:),allocatable  :: previous_Coul_exch ! Exchange Coulomb potential
         
    real(kind=rk),dimension(:,:),allocatable :: C           ! Expansion coefficients
    real(kind=rk),dimension(:),allocatable   :: Ctmp
    real(kind=rk),dimension(:,:),allocatable :: sp_energies ! Single-particle energies
    real(kind=rk),dimension(:,:),allocatable :: Vpair       ! Pairing matrix elements

    real(kind=rk) :: Etot1,Etot2,Etot3,Q20_1,Q20_2,Q20_3,Q40_1,Q40_2,Q40_3
    real(kind=rk),dimension(2) :: Epot
    real(kind=rk) :: e_1st_unoccupied,fermi,Vpair_max
    type(expansion) :: psi_bar
    real(kind=rk) :: ovlp,tmp
    integer,dimension(basis%size) :: partner,order
    real(kind=rk),dimension(basis%size/2) :: quasi_pair_energy
    integer,dimension(basis%size/2) :: sp_label,tmp_label
    real(kind=rk),dimension(basis%size,basis%size) :: np_overlaps
    real(kind=rk),dimension(basis%size) :: overlap
    logical :: filled_qp

    integer :: m, max_rank
    real(kind=rk),parameter :: ph_energy_min = 1.0e12_rk
    real(kind=rk),dimension(:),allocatable :: ph_energy ! Particle-hole excitation energy wrt Fermi level

    real(kind=rk) :: r, z, dr, dz, density, t0

    integer :: l
    type(spherical_basis) :: sph_basis
    character(len=1000) :: sph_config, str
    character(len=5),dimension(2*l_max+1) :: str_l_j
    real(kind=rk),dimension(2*l_max+1) :: weight_l_j
    integer,dimension(2*l_max+1) :: permut_sph
    real(kind=rk),parameter :: weight_threshold = 2.0_rk ! in %

!
!   Header message
!

    if (numerical_parameters%convergence%broyden=='n') then

       write(out_unit,'(/80a1)') ('=',i=1,80)
       do i = 1, 3
          write(out_unit,'(a1,78x,a1)') '|','|'
       end do
       write(out_unit,'(a1,5x,a,47x,a1)') '|','HFBCS  I T E R A T I O N S','|'
       do i = 1, 2
          write(out_unit,'(a1,78x,a1)') '|','|'
       end do
       write(out_unit,'(a1,5x,a,48x,a1)') '|','L I N E A R   M I X I N G','|'
       do i = 1, 3
          write(out_unit,'(a1,78x,a1)') '|','|'
       end do
       write(out_unit,'(80a1/)') ('=',i=1,80)

    end if

    if (mod(basis%size,2) == 1) stop 'odd number of basis states => problem for quasi pairs'

    call DATE_AND_TIME(values = start)

!    if (model%h_odd_perturbative) open(unit=27,file='h_odd_perturbative.txt',status='replace')
!    if (model%h_tensor_perturbative) open(unit=37,file='h_tensor_perturbative.txt',status='replace')

!
!   Auxiliary variables
!

    Npart(:)=(/ nucleus%N , nucleus%Z /)
    mix=numerical_parameters%convergence%mixing

    if (.NOT.PRESENT(Coulomb_exchange)) then
       if (debugmode/='low') write(out_unit,'(/a)') &
            '* Treatment of Coulomb exchange: by the Slater approximation'
       exact_coul_exch = .FALSE.
    else
       SELECT CASE (TRIM(Coulomb_exchange))
       CASE('exact')
          write(out_unit,'(/a)') '* Treatment of Coulomb exchange: Exactly'
          exact_coul_exch = .TRUE.
       CASE('slater')
          write(out_unit,'(/a)') '* Treatment of Coulomb exchange: by the Slater approximation'
          exact_coul_exch = .FALSE.
       CASE DEFAULT
          if (debugmode/='low') write(out_unit,'(/a)') &
               '* Treatment of Coulomb exchange: by the Slater approximation'
          exact_coul_exch = .FALSE.
       END SELECT
    end if

    if (nucleus%fix_rank=='n' .and. debugmode=='high') then
       write(out_unit,*)
       write(out_unit,'(a)') '*** Rank of K^pi state not fixed'
       write(out_unit,'(a)') '    Find the K^pi state nearer to Fermi level'
       write(out_unit,*)
    end if

!----------------------
!
!   Memory allocation
!
!----------------------

!
!   Local variables
!

    Nblocks=basis%Nblocks

    allocate(h_HF(Nblocks),stat=mem_stat)
    if (mem_stat/=0) call fatal(modulename, subroutinename, &
         'Problem with allocate(h_HF)')

    do k=1,Nblocks
       n=basis%block_size(k)
       allocate(h_HF(k)%HO(n,n),h_HF(k)%eigenval(n),stat=mem_stat)
       if (mem_stat/=0) call fatal(modulename, subroutinename, &
            'Problem with allocate(h_HF%HO etc...)')
    end do

    allocate(sp_energies(basis%size,2),stat=mem_stat)
    if (mem_stat/=0) call fatal(modulename, subroutinename, &
         'Problem with allocate(sp_energies)')

    allocate(Vpair(basis%size,basis%size),stat=mem_stat)
    if (mem_stat/=0) call fatal(modulename, subroutinename, &
         'Problem with allocate(Vpair)')

    allocate(psi_bar%cylHO(1:basis%size))
    psi_bar%cylHO=0.0

    if (exact_coul_exch) then
       if (.not. allocated(Coulomb_exchange_potential)) then
          allocate(Coulomb_exchange_potential(Nblocks),previous_Coul_exch(Nblocks),stat=mem_stat)
          if (mem_stat/=0) call fatal(modulename, subroutinename, &
               'Problem with allocate(Coulomb)')
          do k=1,Nblocks
             n=basis%block_size(k)
             allocate(Coulomb_exchange_potential(k)%HO(n,n),previous_Coul_exch(k)%HO(n,n),stat=mem_stat)
             if (mem_stat/=0) call fatal(modulename, subroutinename, &
                  'Problem with allocate(Coul%HO)')
             Coulomb_exchange_potential(k)%HO=0.0
             previous_Coul_exch(k)%HO=0.0
          end do
       else
          allocate(previous_Coul_exch(Nblocks),stat=mem_stat)
          if (mem_stat/=0) call fatal(modulename, subroutinename, &
               'Problem with allocate(Coulomb)')
          do k=1,Nblocks
             n=basis%block_size(k)
             allocate(previous_Coul_exch(k)%HO(n,n),stat=mem_stat)
             if (mem_stat/=0) call fatal(modulename, subroutinename, &
                  'Problem with allocate(Coul%HO)')
             previous_Coul_exch(k)%HO=Coulomb_exchange_potential(k)%HO
          end do
       end if
    end if

!---------------------------------
!
!   Beginning of iteration loop
!
!---------------------------------

    do isospin=1,2
       nb_blocked_part(isospin) = 0
       do i=1,nqp
          if (nucleus%K(isospin,i) /= 0) nb_blocked_part(isospin) = nb_blocked_part(isospin) + 1
       end do
    end do

    open(unit=2,file='eigenval',action='write',status='replace')
    open(unit=41,file='cylindrical_expansion',status='replace')

    finished=.false.
    converged=.false.
    Q20_1=1.0_rk ; Q20_2=2.0_rk ; Q20_3=3.0_rk
    Etot1=1.0_rk ; Etot2=2.0_rk ; Etot3=3.0_rk
    Q40_1=1.0_rk ; Q40_2=2.0_rk ; Q40_3=3.0_rk

    readj_count_max = model%readj_count_max
    accelerate = .false.
!    if (nreadj > 0) readj_count_max = 10 ! Number of loops for constraints readjustment

    loop_iteration: do while (.not.finished)

       iter = iter + 1
       
       accelerate = .false.

       write(out_unit,'(/73a1)') ('=',i=1,73)
       write(out_unit,'(a,i4/)') 'Iteration',iter

       loop_readj: do readjustment_counter = 1, readj_count_max

          if (readj_count_max>1) write(out_unit,'(/2(a,i2,2x))') 'iter=',iter,&
               'readjustment_counter=',readjustment_counter

          nucleus%energy%Epot = 0.0_rk

          call DATE_AND_TIME(values = start_isospin)

          loop_isospin: do isospin=1,2

             if (debugmode=='high') call warning(modulename, subroutinename,'Start loop over isospin')

             if (Npart(isospin)==0) cycle

!
!            Initialization of expansion coefficients to 0
!            (array filled in by blocks ==> need to clean everything at every iteration)
!

             do iHF=1,basis%size
                sp_states(iHF,isospin)%coef%cylHO(:)=0.0
             end do
             sp_states(:,isospin)%pi=0

             call DATE_AND_TIME(values = start_k)

             loop_blocks: do k=1,Nblocks

                n=basis%block_size(k)
                shift = 0
                do i = 1, k-1
                   shift = shift + basis%block_size(i)
                end do

!
!               Hamiltonian matrix elements in the Harmonic Oscillator (HO) basis
!
                call DATE_AND_TIME(values = start_H)

                call hamiltonian_cylHO(h_HF(k)%HO(:,:),k,isospin,basis,HFfield(isospin),iter)

                call DATE_AND_TIME(values = end_H)
                if ((24*0+end_H(5))*3600._rk + end_H(6)*60._rk + end_H(7) + &
                     end_H(8)*0.001_rk - (start_H(5)*3600._rk + start_H(6)*60._rk + &
                     start_H(7) + start_H(8)*0.001_rk) > 1.0_rk) then
                   str_time = "isospin=" // trim(int2char(isospin)) // " block #" // &
                        trim(int2char(k)) // " Hamiltonian matrix calculation"
                   call elapsed_time(start_H,end_H,trim(str_time))
                end if

                if (exact_coul_exch .and. isospin == 2) h_HF(k)%HO = h_HF(k)%HO + &
                     mix*previous_Coul_exch(k)%HO + (1.0_rk-mix)*Coulomb_exchange_potential(k)%HO

!
!               Hamiltonian diagonalization
!

                allocate(C(n,n),Ctmp(n))

                call diagon(h_HF(k)%HO(:,:), h_HF(k)%eigenval(:), C, n)
                
                if (debugmode=='high') then
                   call check_orthogonality(C, n)
                   call check_diagonalization(h_HF(k)%HO(:,:), h_HF(k)%eigenval(:), C, n)
                end if

!
!               Fill in arrays with eigenvalues and eigenvectors
!               (largest expansion coeff. chosen to be positive)
!

                sp_energies(shift+1:shift+n,isospin)=h_HF(k)%eigenval(:)

                do i=1,n
                   do j=1,n
                      Ctmp(j)=abs(C(j,i))
                   end do
                   j=maxloc(Ctmp,1)
                   sp_states(i+shift,isospin)%coef%cylHO(1+shift:n+shift)=C(:,i)* &
                        sign(1.0_rk,C(j,i))
!
!                  Phase convention for Omega < 0 states: |-Omega> ~ T |Omega>
!
                   if (basis%Omega(1+shift) < 0) then
                      sp_states(i+shift,isospin)%coef%cylHO(1+shift:n+shift) = &
                           sp_states(i+shift,isospin)%coef%cylHO(1+shift:n+shift) * &
                           (-1)**((basis%Omega(1+shift)+1)/2)
                   end if
                   sp_states(i+shift,isospin)%energy=h_HF(k)%eigenval(i)
                   sp_states(i+shift,isospin)%block_rank=i
                end do

                deallocate(C,Ctmp)

                if (symmetries%parity) sp_states(1+shift:n+shift,isospin)%pi= &
                     basis%parity(1+shift)
                sp_states(1+shift:n+shift,isospin)%Omega=basis%Omega(1+shift)
                sp_states(1+shift:n+shift,isospin)%block=k

             end do loop_blocks

             call date_and_time(values = end_k)
             write(str_time,'(a,i2,1x,a)') 'isospin=',isospin,'loop over blocks'
             call elapsed_time(start_k,end_k,trim(str_time))

!
!            Sort all eigenvalues in increasing order to define hole and particle states
!            and all eigenvectors accordingly
!

             call sort(sp_energies(:,isospin), basis%size, 'increasing', permut = permut(:,isospin))

             call permutation(sp_states(:,isospin), permut(:,isospin), &
                  basis%size, 'normal')

             do iHF = 1, basis%size
                sp_states(iHF,isospin)%HF_state=iHF
             end do

!
!            Single-particle wave functions and their gradients at Gauss--Hermite--Laguerre
!            integration points
!
             call date_and_time(values = start_k)

             call setup_openmp_parameters(1, basis%size, 16, OMP_parameters%set_bounds, &
                  OMP_parameters%nthreads)
             call omp_set_num_threads(OMP_parameters%nthreads)
             !$OMP PARALLEL DO &
             !$OMP PRIVATE(iHF) &
             !$OMP SHARED(isospin,basis)
             do iHF = 1, basis%size
                call sp_wave_functions(sp_states(iHF,isospin), basis)
             end do
             !$OMP END PARALLEL DO
             
             call date_and_time(values = end_k)
             write(str_time,'(a,i2,1x,a)') 'isospin=',isospin,'sp_wave_functions'
             call elapsed_time(start_k,end_k,trim(str_time))

!
!            Occupation numbers: pairing or pure Hartree--Fock
!

             sp_states(1:Npart(isospin),isospin)%occupation=1.0_rk
             sp_states(Npart(isospin)+1:basis%size,isospin)%occupation=0.0
             nucleus%pairing(isospin)%pairing_energy=0.0

!
!            Initializations:
!
!            * chemical potential to Fermi level (half sum of last occ. and 1st unocc. levels)
!            * pairing energy to 0
!            * pairing gaps and qp quantities to 0
!            * occupation probabilities to the Hatree-Fock values
!
             iHF=Npart(isospin)
             e_1st_unoccupied=sp_states(iHF+1,isospin)%energy
             if (iHF/=0) then
                fermi=0.5_rk*(sp_states(iHF,isospin)%energy+e_1st_unoccupied)
             else
                fermi = e_1st_unoccupied
             end if

             nucleus%pairing(isospin)%chemical_potential=fermi
             nucleus%pairing(isospin)%pairing_energy=0.0
             nucleus%pairing(isospin)%average_gap=0.0
             nucleus%pairing(isospin)%gap(:)=0.0
             nucleus%pairing(isospin)%min_qp_index=0
             nucleus%pairing(isospin)%min_qp_energy=0
             nucleus%pairing(isospin)%Vpair(:,:)=0.0

!
!            Determination of quasi-pairs (which would be Kramers degenerate
!            in the equal-filling approximation) and their energy
!

             partner(:)=0
             sp_states(:,isospin)%pair_partner(isospin)=0
             quasi_pair_energy(:)=1e6_rk
             sp_label(:)=0
             n=0

             call DATE_AND_TIME(values = start_partners)

             call setup_openmp_parameters(1, basis%size, 96, OMP_parameters%set_bounds, &
                  OMP_parameters%nthreads)
             call omp_set_num_threads(OMP_parameters%nthreads)
             !$OMP PARALLEL &
             !$OMP SHARED(isospin,sp_states,partner) &
             !$OMP PRIVATE(id,iHF,jHF,psi_bar,ovlp,tmp)
             id = omp_get_thread_num()
             do iHF = OMP_parameters%set_bounds(id)%indices(1), &
                  OMP_parameters%set_bounds(id)%indices(2)
                if (sp_states(iHF,isospin)%Omega<0) cycle
!               definition of time-reversed state
                psi_bar=time_reversal(sp_states(iHF,isospin)%coef,basis)
                ovlp=1e-09_rk
                partner(iHF)=0
!               search for sp state with opposite Omega and largest overlap with psi_bar
                do jHF=1,basis%size
                   if (sp_states(jHF,isospin)%Omega /= -sp_states(iHF,isospin)%Omega .or. &
                        sp_states(jHF,isospin)%pair_partner(isospin) /= 0) cycle ! In case partner already found
                   tmp=dot_product(sp_states(jHF,isospin)%coef%cylHO,psi_bar%cylHO)
                   if (abs(tmp)>abs(ovlp)) then
                      ovlp=tmp
                      partner(iHF)=jHF
                      overlap(iHF)=ovlp
                   end if
                end do
             end do
             !$OMP END PARALLEL

             do iHF = 1, basis%size
                if (sp_states(iHF,isospin)%Omega<0) cycle
                jHF=partner(iHF)
                if (jHF /= 0) then
                   n=n+1
                   if (n>basis%size/2) call fatal(modulename,subroutinename,'Too many unfound partners')
                   quasi_pair_energy(n)=sp_states(iHF,isospin)%energy+sp_states(jHF,isospin)%energy
                   sp_label(n)=iHF
                   sp_states(iHF,isospin)%pair_partner(isospin)=jHF
                   sp_states(jHF,isospin)%pair_partner(isospin)=iHF
                else
                   write(out_unit,'(a,i4)') 'iHF=',iHF
                   call fatal(modulename,subroutinename,'Partner state not found')
                end if
                if (debugmode=='high' .and. sp_states(iHF,isospin)%energy<0) then
                   jHF=partner(iHF)
                   psi_bar=time_reversal(sp_states(iHF,isospin)%coef,basis)
                   ovlp=dot_product(sp_states(jHF,isospin)%coef%cylHO,psi_bar%cylHO)
                   write(out_unit,'(a,i4,1x,a,i3,a2,a1,2x,a,i4,5(a,f12.6,1x))') &
                        'iHF=',iHF,'Omega-pi=',sp_states(iHF,isospin)%Omega,'/2', &
                        signparity(sp_states(iHF,isospin)%pi),' partner jHF=',jHF, &
                        ' <jHF|iHF-bar>=',ovlp,' quasi-pair energy=',quasi_pair_energy(n), &
                        'overlap(iHF)=',overlap(iHF),'<pi>(iHF)=',sp_states(iHF,isospin)%avg_pi, &
                        '<pi>(jHF)=',sp_states(jHF,isospin)%avg_pi
                end if
             end do
						
             call sort(quasi_pair_energy(:), basis%size/2, 'increasing', permut=order(:))

             tmp_label(:)=0
             do i=1,basis%size/2
                tmp_label(i)=sp_label(order(i))
             end do
             sp_label(:)=tmp_label(:)

             call date_and_time(values = end_partners)
             write(str_time,'(a,i2,1x,a)') 'isospin=',isospin,'partners calculation'
             call elapsed_time(start_partners,end_partners,trim(str_time))

!
!            Occupation numbers in the Hartree--Fock approximation
!

             call date_and_time(values = start_rank)
                
!
!            Finding the closest rank of the blocked state to the Fermi level
!

             if (nucleus%fix_rank == 'n') then
                max_rank = maxval(basis%block_size(:),1)
                if (max_rank <= 0) call fatal(modulename,subroutinename,"max_rank <= 0")
                allocate(ph_energy(max_rank))
                if (symmetries%parity) then
                   do j = 1, nqp
                      ph_energy(:) = ph_energy_min
                      if (nucleus%K(isospin,j) == 0) cycle
                      do m = 1, max_rank
                         do jHF = 1, basis%size
                            if (sp_states(jHF,isospin)%Omega /= nucleus%K(isospin,j) .or. &
                                 sp_states(jHF,isospin)%pi /= nucleus%pi(isospin,j)) cycle
                            if (sp_states(jHF,isospin)%block_rank == m) ph_energy(m) = &
                                 abs(sp_states(jHF,isospin)%energy-fermi)
                         end do
                      end do
                      if (minval(ph_energy) == ph_energy_min) call fatal(modulename,subroutinename, &
                           'Rank not found')
                      nucleus%rank(isospin,j) = minloc(ph_energy,dim=1)
                      if (debugmode /= 'low') write(out_unit,'(4x,3(a,i2),2x,a,f9.3)') 'Blocked K^pi state =', &
                           nucleus%K(isospin,j),'/2,  parity = ',nucleus%pi(isospin,j), &
                           ',  Rank = ', nucleus%rank(isospin,j),'Fermi energy =',fermi
                      do m = 1, max_rank
                         if (m == nucleus%rank(isospin,j)) cycle
                         if (abs(ph_energy(m)-ph_energy(nucleus%rank(isospin,j))) <= 0.5_rk) then
                            write(out_unit,'(4x,2(a,i2),a)') 'NOTE: The distances from Fermi level of ranks ', &
                                 nucleus%rank(isospin,j),' and ',m,' differ by less than 500 keV'
                            write(out_unit,'(3x,a,i2,a,f12.6)') 'Distance from Fermi level for rank ', &
                                 nucleus%rank(isospin,j), ' = ', ph_energy(nucleus%rank(isospin,j))
                            write(out_unit,'(3x,a,i2,a,f12.6)') 'Distance from Fermi level for rank ', &
                                 m, ' = ', ph_energy(m)
                         end if
                      end do
                   end do
                else
                   do j = 1, nqp
                      ph_energy(:) = ph_energy_min
                      if (nucleus%K(isospin,j) == 0) cycle
                      do m = 1, max_rank
                         do jHF = 1, basis%size
                            if (sp_states(jHF,isospin)%Omega /= nucleus%K(isospin,j) .or. &
                                 sp_states(jHF,isospin)%avg_pi*nucleus%pi(isospin,j) < 0.0_rk) cycle
                            if (sp_states(jHF,isospin)%block_rank == m) ph_energy(m) = &
                                 abs(sp_states(jHF,isospin)%energy-fermi)
                         end do
                      end do
                      if (minval(ph_energy) == ph_energy_min) call fatal(modulename,subroutinename, &
                           'Rank not found')
                      nucleus%rank(isospin,j) = minloc(ph_energy,dim=1)
                      if (debugmode /= 'low') write(out_unit,'(4x,3(a,i2),2x,a,f9.3)') 'Blocked K^pi state =', &
                           nucleus%K(isospin,j),'/2,  parity = ',nucleus%pi(isospin,j), &
                           ',  Rank = ', nucleus%rank(isospin,j),'Fermi energy =',fermi
                      do m = 1, max_rank
                         if (m == nucleus%rank(isospin,j)) cycle
                         if (abs(ph_energy(m)-ph_energy(nucleus%rank(isospin,j))) <= 0.5_rk) then
                            write(out_unit,'(4x,2(a,i2),a)') 'NOTE: The distances from Fermi level of ranks ', &
                                 nucleus%rank(isospin,j),' and ',m,' differ by less than 500 keV'
                            write(out_unit,'(3x,a,i2,a,f12.6)') 'Distance from Fermi level for rank ', &
                                 nucleus%rank(isospin,j), ' = ', ph_energy(nucleus%rank(isospin,j))
                            write(out_unit,'(3x,a,i2,a,f12.6)') 'Distance from Fermi level for rank ', &
                                 m, ' = ', ph_energy(m)
                         end if
                      end do
                   end do
                end if
                deallocate(ph_energy)
             end if

             call date_and_time(values = end_rank)
             write(str_time,'(a,i2,1x,a)') 'isospin=',isospin,'rank calculation'
             call elapsed_time(start_rank,end_rank,trim(str_time))

!            Determination of the blocked states
             
             if (model%blocking == 'no') then
                i_odd(isospin,:)=0
             else
                do i=1,nqp
                   i_odd(isospin,i)=0
                   if (nucleus%K(isospin,i) == 0) cycle
                   loop_iHF: do iHF=1,basis%size
                      if (sp_states(iHF,isospin)%Omega==nucleus%K(isospin,i) .and. &
                           sp_states(iHF,isospin)%block_rank==nucleus%rank(isospin,i)) then
                         if (symmetries%parity) then
                            if (sp_states(iHF,isospin)%pi==nucleus%pi(isospin,i)) then
                               i_odd(isospin,i)=iHF
                               exit loop_iHF
                            end if
                         else
                            i_odd(isospin,i)=iHF
                            exit loop_iHF
                         end if
                      end if
                   end do loop_iHF
                end do
             end if
             
             if (debugmode /= 'low') then
                do i=1,nqp
                   if (nucleus%K(isospin,i) == 0) cycle
                   write(out_unit,'(4(a,i4,2x))') &
                        'isospin=',isospin,'i=',i,'i_odd=',i_odd(isospin,i), &
                        'partner(i_odd)=',sp_states(i_odd(isospin,i),isospin)%pair_partner(isospin)
                end do
             end if

!            Initialization of occupation numbers and type component

             sp_states(:,isospin)%occupation = 0.0
             sp_states(:,isospin)%type = 0

!            Occupation numbers of the blocked states

             call DATE_AND_TIME(values = start_occupations)
                
             do i=1,nqp
                if (i_odd(isospin,i)/=0) then
                   if (model%blocking == 'SCB') then
                      sp_states(i_odd(isospin,i),isospin)%occupation = 1.0_rk
                   else if (model%blocking == 'EFA') then
                      sp_states(i_odd(isospin,i),isospin)%occupation = 0.5_rk
                      jHF = sp_states(i_odd(isospin,i),isospin)%pair_partner(isospin)
                      sp_states(jHF,isospin)%occupation = 0.5_rk
                   else
                      call fatal(modulename,subroutinename,'Unrecognized blocking scheme ' &
                           //trim(model%blocking))
                   end if
                else if (nucleus%K(isospin,i) /= 0 .and. model%blocking /= 'no') then
                   call fatal(modulename,subroutinename,'Blocked state not found')
                end if
             end do

!            Occupation numbers of remaining, lowest nb_qp(isospin) quasi-pairs 
!            of states set to 1

             nb_qp(isospin) = (Npart(isospin) - nb_blocked_part(isospin))/2
             jHF = 1
             do k = 1, nb_qp(isospin)
                filled_qp = .false.
                loop_fill_qp: do j = jHF, basis%size/2
                   iHF = sp_label(j)
                   do i = 1,nqp
                      if (iHF == i_odd(isospin,i) .or. partner(iHF) == i_odd(isospin,i)) cycle loop_fill_qp
                   end do
                   sp_states(iHF,isospin)%occupation = 1.0_rk
                   sp_states(partner(iHF),isospin)%occupation = 1.0_rk
                   filled_qp = .true.
                   jHF=j+1
                   exit loop_fill_qp
                end do loop_fill_qp
                if (.not.filled_qp) then
                   write(out_unit,'(a,i2,2(1x,a,i4),2x,a,3i5)') &
                        'isospin=',isospin,'k=',k,'jHF=',jHF,'i_odd(isospin,:)=',i_odd(isospin,:)
                   call fatal(modulename,subroutinename,"Quasi-pair not filled")
                end if
                if (debugmode=='high') write(out_unit,'(a1,2x,2(i4,1x),a2,a1,2(i4,1x),f12.6)') &
                     nucleon_type(isospin),k,sp_states(iHF,isospin)%Omega,'/2', &
                     signparity(sp_states(iHF,isospin)%pi),iHF,partner(iHF),quasi_pair_energy(j)
             end do

!            Define HF_level, HF_state and type components of sp_states

             i=0
             do k=1,basis%size
                sp_states(k,isospin)%HF_level = k
                sp_states(k,isospin)%HF_state = k
                if (sp_states(k,isospin)%occupation >= 0.5_rk) then
                   sp_states(k,isospin)%type = 1
                   i=i+1
                else
                   sp_states(k,isospin)%type = 0
                end if
             end do

             if (model%blocking == 'EFA') then
                do k = 1, nqp
                   if (i_odd(isospin,k) == 0) cycle
                   jHF = sp_states(i_odd(isospin,k),isospin)%pair_partner(isospin)
                   sp_states(jHF,isospin)%type = 0
                   i = i - 1
                end do
             end if

             if (i /= Npart(isospin) .and. model%blocking /= 'no') then
                do iHF = 1, basis%size
                   write(out_unit,'(a1,2x,2(i4,1x),a2,a1,3(i4,1x))') &
                        nucleon_type(isospin),iHF,sp_states(iHF,isospin)%Omega,'/2', &
                        signparity(sp_states(iHF,isospin)%pi),iHF,partner(iHF), &
                        sp_states(iHF,isospin)%type
                end do
                write(out_unit,'(/a,i2,1x,a,i4,a,i4)') 'isospin=',isospin,&
                     'number of hole states =', i, ' whereas Npart=', Npart(isospin)
                call fatal(modulename,subroutinename, &
                     'Problem with the particle or hole character of s.p. states')
             end if

             call date_and_time(values = end_occupations)
             write(str_time,'(a,i2,1x,a)') 'isospin=',isospin,'occupations calculation'
             call elapsed_time(start_occupations,end_occupations,trim(str_time))

             call date_and_time(values = start_BCS)
                
!
!            Matrix elements of the effective interaction in the pairing channel
!            (replacing the Skyrme interaction)
!

             call pairing_matrix_elements(sp_states(:,isospin), densities, isospin, &
                  nucleus, model, basis, Vpair, partner, iter)

!
!            Occupation numbers in the BCS approximation
!
!            Check if BCS calculation is needed:
!            * if max{|Vpair(i,j)|,i,j=1..basis%size} is larger than 1 keV, call BCS
!            * if not, do nothing ==> HF solution
!

             Vpair_max = maxval(abs(Vpair))
             if (abs(Vpair_max)>1e-3_rk) then
                call BCS(nucleus, model, basis, sp_states(:,isospin), Vpair, isospin, &
                     i_odd(isospin,:))
             end if

             tmp = 0
             do iHF = 1, basis%size
                tmp = tmp + sp_states(iHF,isospin)%occupation
             end do
             if (isospin == 1) then
                nucleus%avg_particle_number(1) = tmp
             else if (isospin == 2) then
                nucleus%avg_particle_number(2) = tmp
             end if

             call date_and_time(values = end_BCS)
             write(str_time,'(a,i2,1x,a)') 'isospin=',isospin,'BCS calculation'
             call elapsed_time(start_BCS,end_BCS,trim(str_time))

!
!            Local one-body densities in the coordinate space at Gauss--Hermite--Laguerre
!            integration points
!

             call DATE_AND_TIME(values = start_H)

             call local_one_body_densities(sp_states(:,isospin), basis, densities(isospin), &
                  isospin, iter)

             call DATE_AND_TIME(values = end_H)
             write(str_time,'(a,i2,1x,a)') 'isospin=',isospin,'densities calculation'
             call elapsed_time(start_H,end_H,trim(str_time))

!
!            Potential energy
!
             call DATE_AND_TIME(values = start_H)

             call potential_energy(nucleus, sp_states, basis, h_HF, isospin, Epot(isospin))
             nucleus%energy%Epot = nucleus%energy%Epot + Epot(isospin)

             call DATE_AND_TIME(values = end_H)
             write(str_time,'(a,i2,1x,a)') 'isospin=',isospin,'potential_energy calculation'
             call elapsed_time(start_H,end_H,trim(str_time))
 
          end do loop_isospin

          nucleus%avg_particle_number(0) = nucleus%avg_particle_number(1) + &
               nucleus%avg_particle_number(2)

          call DATE_AND_TIME(values = end_isospin)
          call elapsed_time(start_isospin,end_isospin,"loop on isospin")
!
!         Renormalization of local densities
!

          call renormalize_densities(densities, basis, nucleus)

!
!         Multipole moments of the nuclear density (Q20 enough, save time by 
!         calculating the other moments after convergence)
!

          call moments(nucleus, densities, sp_states, basis, symmetries)

!
!         Readjustment of constraints
!

          if (nreadj > 0) then
             call constraints_readjustment(nucleus,basis,sp_states,model,HFfield,error,accelerate)
             write(out_unit,'(a,e15.6)') 'constraint error=',error
             if (error < 1e-3_rk) exit loop_readj
          end if

       end do loop_readj

!
!      Direct Coulomb potential
!

       call date_and_time(values = start_potential)

       call Coulomb(HFfield(2), densities, basis, symmetries)

!
!      Total energy and exact Coulomb exchange potential if requested
!

       if (exact_coul_exch) then
          do k=1,basis%Nblocks
             previous_Coul_exch(k)%HO = Coulomb_exchange_potential(k)%HO
          end do
          call Exact_Coulomb_potential('exc', Coulomb_exchange_potential, basis, sp_states, nucleus%Z)
          call energy(nucleus, model, basis, HFfield(2), densities, sp_states, Coulomb_exchange_potential)
       else
          call energy(nucleus, model, basis, HFfield(2), densities, sp_states)
       end if
       
       call DATE_AND_TIME(values = end_potential)
       write(str_time,'(a)') 'Coulomb potential calculation'
       call elapsed_time(start_potential,end_potential,trim(str_time))

!
!      Hartree-Fock field (central, effective mass and spin-orbite potentials)
!

       call date_and_time(values = start_potential)
 
       call potential(nucleus, model, symmetries, basis, numerical_parameters, &
            HFfield, densities, exact_coul_exch)
       
       call DATE_AND_TIME(values = end_potential)
       write(str_time,'(a)') 'nuclear mean-field potential calculation'
       call elapsed_time(start_potential,end_potential,trim(str_time))

!
!      Dominant Nilsson configurations
!

       call date_and_time(values = start_Nilsson)
 
       call Nilsson_expansion(nucleus, model, symmetries, basis, sp_states)
       
       call DATE_AND_TIME(values = end_Nilsson)
       write(str_time,'(a)') 'dominant Nilsson configurations'
       call elapsed_time(start_Nilsson,end_Nilsson,trim(str_time))

!
!      Displays results on standard output
!

       call display_iteration_results(nucleus, model, symmetries, iter, basis%size, sp_states)

       if (nucleus%spherical_expansion) then

          do isospin = 1, 2
             call cyl2sph(nucleus,symmetries,basis,sp_states(:,isospin), &
                  sph_basis,isospin)
          end do

          write(out_unit,'(/77a1/,/10x,a/,/77a1)') ('=',i=1,77), &
               'S P H E R I C A L   E X P A N S I O N', &
               ('-',i=1,77)

          write(out_unit,'(/a,f7.2,a)') 'Weight display threshold =',weight_threshold,'%'

          do isospin = 1, 2
             if (symmetries%parity) then 
                write(out_unit,'(/77a1/,a/,77a1)') ('-',i=1,77), &
                     'iso. label Omega pi rank   energy  occup.   weights (%)', &
                     ('-',i=1,77)
             else
                write(out_unit,'(/77a1/,a/,77a1)') ('-',i=1,77), &
                     'iso. label Omega <pi> rank   energy  occup.   weights (%)', &
                     ('-',i=1,77)
             end if
             do iHF = 1, basis%size
                if (sp_states(iHF,isospin)%energy > 2*model%truncmax + &
                     sp_states(Npart(isospin),isospin)%energy) exit
                i = 0
                do l = 0, l_max
                   do j = abs(2*l-1), 2*l+1, 2
                      i = i + 1
                      weight_l_j(i) = sp_states(iHF,isospin)%sph_weight(l,j)
                      str_l_j(i)=trim(spectroscopic_notation(l,j))
                   end do
                end do
                if (i /= 2*l_max+1) stop 'Error of counting of (l,j) combinations'
                call sort(weight_l_j,2*l_max+1,'decreasing',permut_sph)
                sph_config=''
                do i = 1, 2*l_max+1
                   if (weight_l_j(i) < weight_threshold) exit ! Keep weights > 1%
                   write(str,'(f5.1,a)') weight_l_j(i),'% ' // trim(str_l_j(permut_sph(i)))
                   sph_config = trim(sph_config) // trim(str)
                end do
                if (symmetries%parity) then
                   write(out_unit,'(a1,2(1x,i3),a2,1x,a1,1x,i4,1x,f10.4,1x,f7.4,2x,a)') &
                        nucleon_type(isospin),iHF,sp_states(iHF,isospin)%Omega,'/2', &
                        parity_sign(sp_states(iHF,isospin)%pi),sp_states(iHF,isospin)%block_rank, &
                        sp_states(iHF,isospin)%energy,sp_states(iHF,isospin)%occupation, &
                        trim(sph_config)
                else
                   write(out_unit,'(a1,2(1x,i3),a2,1x,f7.3,1x,i4,1x,f10.4,1x,f7.4,2x,a)') &
                        nucleon_type(isospin),iHF,sp_states(iHF,isospin)%Omega,'/2', &
                        sp_states(iHF,isospin)%avg_pi,sp_states(iHF,isospin)%block_rank,&
                        sp_states(iHF,isospin)%energy,sp_states(iHF,isospin)%occupation, &
                        trim(sph_config)
                end if
             end do
             write(out_unit,'(77a1)') ('-',i=1,77)
          end do

          write(out_unit,'(/77a1/)') ('=',i=1,77)

       end if

!
!      Calculation of rho local density on an equidistant mesh
!
       ! t0 = (8*model%parameters%B1 + 4*model%parameters%B2)/3.0_rk
       ! isospin = 1
       ! z = 0.00_rk
       ! dr = 0.25_rk
       ! r = 0.0_rk
       ! open(unit=77,file='central_potential_t0.txt',status='replace')
       ! open(unit=78,file='rho.txt',status='replace')
       ! write(77,'(a,f12.3,a)') '# z  r  V_WS  U(t0=',t0,')'
       ! do while(r <= 7.0_rk)
       !    density = rho(r,z,isospin,basis,sp_states(:,isospin))
       !    write(77,'(4(f12.6,1x))') z, r, &
       !         -70.0_rk/(1.0_rk+exp((r-1.35_rk*nucleus%A**(1.0_rk/3.0_rk))/0.7_rk)), &
       !         1.5_rk*density*t0
       !    write(78,'(3(f12.6,1x))') z, r, density
       !   r = r + dr
       ! end do
       ! close(77)
       ! close(78)

!
!      Test for convergence
!

       Etot3=Etot2
       Etot2=Etot1
       Etot1=nucleus%energy%Etot
       Q20_3=Q20_2
       Q20_2=Q20_1
       Q20_1=nucleus%Q20(0)
       Q40_3=Q40_2
       Q40_2=Q40_1
       Q40_1=nucleus%Q40(0)

       converged = ( abs(Etot1-Etot2)+abs(Etot2-Etot3) <= &
            2.0_rk*numerical_parameters%convergence%E_precision ) .and. &
            ( abs(Q20_1-Q20_2)+abs(Q20_2-Q20_3) <= &
            2.0_rk*numerical_parameters%convergence%Q20_precision) .and. &
	    ( abs(Q40_1-Q40_2) + abs(Q40_2-Q40_3) <= &
	    2.0_rk*numerical_parameters%convergence%Q40_precision)
       if (iter >= numerical_parameters%convergence%Nitermax.or.converged) &
            finished=.true.

!
!      Display detailed information about the single-particle states
!

       if (debugmode=='high'.or.finished) then
          do isospin=1,2
             write(41,'(2(a,f8.3,2x))') 'hbar*omega_p =',hb0(isospin)*bp2, &
                  'hbar*omega_z =',hb0(isospin)*bz2
             write(41,'(a)') 'iter nucl. type ' // &
                  'iHF  e  Omega  pi iHO  <iHO|iHF>  nz  np  Lambda  Sigma'  
             do iHF=1,basis%size
                if (sp_states(iHF,isospin)%energy > sp_states(Npart(isospin),isospin)%energy + &
                     2*model%truncmax) exit
                write(2,'(i4,2x,a1,i4,1x,f12.6,1x,i4,a2,1x,a1,1x,2(f12.6,1x),2(1x,i4))') &
                     iter-1,nucleon_type(isospin),iHF, &
                     sp_states(iHF,isospin)%energy, &
                     sp_states(iHF,isospin)%Omega,'/2 ', &
                     signparity(sp_states(iHF,isospin)%pi), &
                     sp_states(iHF,isospin)%occupation, &
                     nucleus%pairing(isospin)%gap(iHF), &
                     sp_states(iHF,isospin)%block_rank, &
                     sp_states(iHF,isospin)%pair_partner(isospin)
                do i=1,basis%size
                   if (100*(sp_states(iHF,isospin)%coef%cylHO(i))**2<1.0_rk) cycle
                   write(41,'(i4,2x,a1,i4,1x,f13.6,1x,i4,a2,1x,a1,1x,i4,2x,e15.8,3(1x,i4),i4,a2)') &
                        iter-1,nucleon_type(isospin),iHF,sp_states(iHF,isospin)%energy, &
                        sp_states(iHF,isospin)%Omega,'/2', &
                        signparity((-1)**(basis%nz(i)+basis%lambda(i))),i, &
                        sp_states(iHF,isospin)%coef%cylHO(i),basis%nz(i), &
                        2*basis%nr(i)+abs(basis%Lambda(i)),basis%Lambda(i),basis%ns(i),'/2'
                end do
             end do
             write(2,*)
          end do
          write(2,*)
       end if

    end do loop_iteration

!---------------------------
!
!   End of iteration loop
!
!---------------------------

!
!   Calculate neutron-proton overlaps
!

    call DATE_AND_TIME(values = start_H)

    np_overlaps(:,:)=0.0
    call calculate_np_overlaps(basis, sp_states, np_overlaps, progress=.FALSE.)
    if (debugmode /= 'low') then
       do iHF=1,basis%size
          if (sp_states(iHF,1)%Omega<0) cycle
          jHF=MAXLOC(abs(np_overlaps(iHF,:)),dim=1)
          i=sp_states(jHF,2)%pair_partner(2)
          j=sp_states(i,2)%pair_partner(2)
          write(out_unit,'(2(a,i4,1x,a,i3,a2,a1,2x),a,f12.6)') &
               'neutron iHF=',iHF,'Omega-pi=',sp_states(iHF,1)%Omega,'/2', &
               signparity(sp_states(iHF,1)%pi),'proton partner jHF=',i, &
               'Omega-pi=',sp_states(i,2)%Omega,'/2',signparity(sp_states(i,2)%pi),&
               '<jHF|iHF>=',np_overlaps(iHF,jHF)
       end do
    end if

    call DATE_AND_TIME(values = end_H)
    write(out_unit,*)
    call elapsed_time(start_H,end_H,"np overlaps calculation")

!
!   Matrix elements of h_HF and sz with final sp_states
!

!    if (model%h_odd_perturbative) call h_odd_perturbative(nucleus, model, symmetries, basis, &
!         HFfield, densities, iter, 27, 'perturbative')

!    if (model%h_odd_perturbative) call h_odd_perturbative(nucleus, model, symmetries, basis, &
!         HFfield, densities, iter, 27, 'exact')

!    if (model%h_tensor_perturbative) call h_tensor_perturbative(nucleus, model, symmetries, basis, &
!         HFfield, densities, iter, 37, 'perturbative')

!    if (model%h_tensor_perturbative) call h_tensor_perturbative(nucleus, model, symmetries, basis, &
!         HFfield, densities, iter, 37, 'exact')

    close(2)
    close(41)
!    if (model%h_odd_perturbative) close(27)
!    if (model%h_tensor_perturbative) close(37)

!---------------------------------------------
!
!   Memory deallocation for local variables
!
!---------------------------------------------

    do k=1,Nblocks
       deallocate(h_HF(k)%HO,h_HF(k)%eigenval)
    end do

    deallocate(h_HF,sp_energies,Vpair)

    deallocate(psi_bar%cylHO)

    call DATE_AND_TIME(values = end)
    call elapsed_time(start,end,"HFBCS calculation")

    return

  contains

    subroutine potential_energy(nucleus, sp_states, basis, h_HF, isospin, Epot)

      implicit none
! Arguments
      type(nucleus_properties),intent(inout) :: nucleus
      type(cylindrical_basis),intent(in) :: basis
      type(sp_state),dimension(basis%size,2),intent(in):: sp_states
      type(operator),dimension(basis%Nblocks),intent(in) :: h_HF
      integer,intent(in) :: isospin
      real(kind=rk),intent(out) :: Epot
! Local variables
      integer :: iblock, iHO, jHO, iHF, n, shift
      integer,dimension(2) :: nz, np, lz, sz
      real(kind=rk) :: v2, vbar

      Epot = 0.0
      do iHF = 1, basis%size
         v2 = sp_states(iHF,isospin)%occupation
         vbar = 0.0_rk
         iblock = sp_states(iHF,isospin)%block
         n = basis%block_size(iblock)
         shift=sum(basis%block_size(1:iblock),1)-n
         do iHO = 1, n
            nz(1) = basis%nz(iHO+shift)
            lz(1) = basis%Lambda(iHO+shift)
            np(1) = 2*basis%nr(iHO+shift) + abs(lz(1))
            sz(1) = basis%ns(iHO+shift)
            do jHO = 1, n
               nz(2) = basis%nz(jHO+shift)
               lz(2) = basis%Lambda(jHO+shift)
               np(2) = 2*basis%nr(jHO+shift) + abs(lz(2))
               sz(2) = basis%ns(jHO+shift)
               vbar = vbar + sp_states(iHF,isospin)%coef%cylHO(iHO+shift) * &
                    sp_states(iHF,isospin)%coef%cylHO(jHO+shift) * &
                    (h_HF(iblock)%HO(iHO,jHO) - HO_psqr(nz,np,lz,sz)*hb0(isospin)*(1.0_rk-c_com/nucleus%A))
            end do
         end do
         Epot = Epot + v2 * vbar
      end do
      Epot = 0.5_rk*Epot

    end subroutine potential_energy

    function rho(r,z,isospin,basis,sp_states)

      implicit none
! Arguments
      real(kind=rk),intent(in) :: r, z
      integer,intent(in) :: isospin
      type(cylindrical_basis),intent(in) :: basis
      type(sp_state),dimension(basis%size),intent(in):: sp_states
! Local variables
      real(kind=rk) :: rho
      real(kind=rk) :: zeta, eta, normcoef1, normcoef2, Hermite1, Hermite2, Laguerre1, Laguerre2
      real(kind=rk) :: s, tmp_plus, tmp_minus
      integer :: iHF, iblock, i, first_HO, last_HO, iHO
      integer :: nzi, nri, li, si

      zeta = bz*z
      eta = bp2*r**2

      s = 0.0_rk
      do iHF = 1, basis%size
         if (sp_states(iHF)%occupation <= 1e-3_rk) cycle
         iblock = sp_states(iHF)%block
         first_HO = sum([(basis%block_size(i),i=1,iblock)]) - basis%block_size(iblock) + 1
         last_HO = first_HO + basis%block_size(iblock) - 1
         tmp_plus = 0.0_rk
         tmp_minus = 0.0_rk
         do iHO = first_HO, last_HO
            nzi = basis%nz(iHO)
            nri = basis%nr(iHO)
            li=basis%Lambda(iHO)
            si=basis%ns(iHO)
            if (si == 1) then
               tmp_plus = tmp_plus + sp_states(iHF)%coef%cylHO(iHO) * &
                    sqrt(fact(nri)/fact(nri+abs(li))/(2**nzi*fact(nzi))) * &
                    sqrt(eta**(abs(li)))*Laguerre_poly(nri,real(abs(li),kind=rk),eta)* &
                    Hermite_poly(nzi,zeta)
            else
               tmp_minus = tmp_minus + sp_states(iHF)%coef%cylHO(iHO) * &
                    sqrt(fact(nri)/fact(nri+abs(li))/(2**nzi*fact(nzi))) * &
                    sqrt(eta**(abs(li)))*Laguerre_poly(nri,real(abs(li),kind=rk),eta)* &
                    Hermite_poly(nzi,zeta)
            end if
         end do
         s = s + sp_states(iHF)%occupation * (tmp_plus**2 + tmp_minus**2)
      end do

      rho = (bz*bp2/pi**1.5_rk) * exp(-(eta+zeta**2)) * s

    end function rho

  end subroutine solve_HFBCS_equations


  !=====================================================================
  !                       SUBROUTINE SOLVE_HFHTDA_EQUATIONS
  !---------------------------------------------------------------------
  ! Solve the Hartree-Fock and HTDA equations for a given nucleus, within 
  ! the model given as 2nd argument, restricted to symmetries passed as 
  ! argument 3, by expansion of single-particle states in the given 
  ! basis (argument 4) and with the numerical parameters embedded in 
  ! argument 5.
  ! Iterative resolution of Hartree-Fock and HTDA equations using the 
  ! diagonalization method
  ! The solution is contained in the array 'sp_states' (argument 6)
  ! which is of derived type 'sp_wave_function' (embedding all relevant 
  ! properties of the sp states (energy, occupation probability, quantum 
  ! numbers, ...). Pairing properties are embedded in 'nucleus'.
  !=====================================================================

  subroutine solve_HFHTDA_equations(nucleus, model, symmetries, basis, &
       numerical_parameters, HTDA_basis, HFfield, sp_states, densities, &
       Coulomb_exchange,Coulomb_exchange_potential)

!
!   Arguments
!
    type(nucleus_properties),intent(inout)             :: nucleus
    type(modelization),intent(in)                      :: model
    type(symmetry_properties),intent(in)               :: symmetries
    type(cylindrical_basis),intent(inout)              :: basis
    type(numeric) ,intent(in)                          :: numerical_parameters
    type(field),dimension(2),intent(inout):: HFfield
    type(sp_state),dimension(basis%size,2)             :: sp_states
    type(local_density),dimension(2),intent(inout)     :: densities
    character(len=*),optional,intent(in)               :: Coulomb_exchange
    type(operator),dimension(:),allocatable,optional,intent(inout) :: Coulomb_exchange_potential
!
!   Local variables
!
    character(len=20),parameter :: subroutinename='solve_HTDA_equations'
    integer, dimension(8) :: start, end
    integer,dimension(2) :: Npart,nb_qp
    integer :: Nblocks
    integer :: k,n,shift
    integer :: mem_stat
    integer :: iter,isospin
    integer :: i,j,iHF,jHF
    logical :: converged,finished
    logical :: exact_coul_exch
    real(kind=rk) :: mix

    type(operator),dimension(:),allocatable  :: h_HF        ! Hartree-Fock Hamiltonian
    type(operator),dimension(:),allocatable  :: previous_Coul_exch ! Exchange Coulomb potential
    real(kind=rk),dimension(:,:),allocatable :: C           ! Expansion coefficients
    real(kind=rk),dimension(:),allocatable   :: Ctmp
    real(kind=rk),dimension(:,:),allocatable :: sp_energies ! Single-particle energies
    real(kind=rk),dimension(:,:),allocatable :: Vpair       ! Pairing matrix elements

    real(kind=rk) :: Etot1,Etot2,Etot3,Q20_1,Q20_2,Q20_3
    real(kind=rk) :: e_1st_unoccupied,fermi
    type(expansion) :: psi_bar
    real(kind=rk) :: ovlp,tmp
    integer,dimension(basis%size) :: partner,order
    real(kind=rk),dimension(basis%size/2) :: quasi_pair_energy
    integer,dimension(basis%size/2) :: sp_label,tmp_label
    real(kind=rk),parameter :: ph_energy_min = 1.0e12_rk
    integer :: m, max_rank
    real(kind=rk),dimension(:),allocatable :: ph_energy ! Particle-hole excitation energy wrt Fermi level
    logical :: filled_qp

!
!   Local variables for HTDA calls
!

  character(len=3)                        :: storage="RAM"! Defines in which way memory will be managed (storage = "RAM" or "HDD").
  type(MB_operator)                       :: Vres         ! Residual interaction for HTDA calculations
  type(MB_operator)                       :: Hiqp         ! Independant quasiparticle operator of the HTDA Hamiltonian construction.
  type(MB_operator)                       :: Hhtda        ! HTDA Hamiltonian operator.
!  type(MB_operator),dimension(-2:2)       :: Q2           ! Q2(mu): Quadrupole operator for nuclear shape studies.
!  type(MB_operator)                       :: r2           ! One-body R square operator
!  type(MB_operator)                       :: J2,Jplus,Jplusbar,Splus,Splusbar
  type(many_body_basis)                   :: HTDA_basis   ! Many Body Basis
  type(useful_window)                     :: RPA_VS, Pairing_VS, Total_window
  real(kind=rk),dimension(:), allocatable :: eigenval_H
  real(kind=rk),dimension(:,:),allocatable:: eigenvec_H
  type(sp_state),dimension(:,:),allocatable  :: Useful_Space
  integer :: nval

!
!   Debugging note
!

    if (debugmode=='medium' .or. debugmode=='high') call debugnote(modulename, &
         subroutinename, 'Beginning of subroutine')

    call DATE_AND_TIME(values = start)

!
!   Auxiliary variables
!

    Npart(:)=(/ nucleus%N , nucleus%Z /)
    mix=numerical_parameters%convergence%mixing

    if (.NOT.PRESENT(Coulomb_exchange)) then
       if (debugmode/='low') write(out_unit,'(/a)') &
            '* Treatment of Coulomb exchange: by the Slater approximation'
       exact_coul_exch = .FALSE.
    else
       SELECT CASE (TRIM(Coulomb_exchange))
       CASE('exact')
          write(out_unit,'(/a)') '* Treatment of Coulomb exchange: Exactly'
          exact_coul_exch = .TRUE.
       CASE('slater')
          write(out_unit,'(/a)') '* Treatment of Coulomb exchange: by the Slater approximation'
          exact_coul_exch = .FALSE.
       CASE DEFAULT
          if (debugmode/='low') write(out_unit,'(/a)') &
               '* Treatment of Coulomb exchange: by the Slater approximation'
          exact_coul_exch = .FALSE.
       END SELECT
    end if

!----------------------
!
!   Memory allocation
!
!----------------------

!
!   Local variables
!

    Nblocks=basis%Nblocks

    allocate(h_HF(Nblocks),stat=mem_stat)
    if (mem_stat/=0) call fatal(modulename, subroutinename, &
         'Problem with allocate(h_HF)')

    do k=1,Nblocks
       n=basis%block_size(k)
       allocate(h_HF(k)%HO(n,n),h_HF(k)%eigenval(n),stat=mem_stat)
       if (mem_stat/=0) call fatal(modulename, subroutinename, &
            'Problem with allocate(h_HF%HO etc...)')
    end do

    allocate(sp_energies(basis%size,2),stat=mem_stat)
    if (mem_stat/=0) call fatal(modulename, subroutinename, &
         'Problem with allocate(sp_energies)')

    allocate(Vpair(basis%size,basis%size),stat=mem_stat)
    if (mem_stat/=0) call fatal(modulename, subroutinename, &
         'Problem with allocate(Vpair)')

    allocate(psi_bar%cylHO(1:basis%size))
    psi_bar%cylHO=0.0

    if (exact_coul_exch) then
       allocate(Coulomb_exchange_potential(Nblocks),previous_Coul_exch(Nblocks),stat=mem_stat)
       if (mem_stat/=0) call fatal(modulename, subroutinename, &
         'Problem with allocate(Coulomb)')
       do k=1,Nblocks
          n=basis%block_size(k)
          allocate(Coulomb_exchange_potential(k)%HO(n,n),previous_Coul_exch(k)%HO(n,n),stat=mem_stat)
          if (mem_stat/=0) call fatal(modulename, subroutinename, &
               'Problem with allocate(Coul%HO etc...)')
          Coulomb_exchange_potential(k)%HO=0.0
          previous_Coul_exch(k)%HO=0.0
       end do
    end if

!---------------------------------
!
!   Beginning of iteration loop
!
!---------------------------------

    do isospin=1,2
       nb_blocked_part(isospin) = 0
       do i=1,nqp
          if (nucleus%K(isospin,i) /= 0) nb_blocked_part(isospin) = nb_blocked_part(isospin) + 1
       end do
    end do

    open(unit=2,file='eigenval',action='write',status='replace')
    open(unit=41,file='cylindrical_expansion',status='replace')

    finished=.false.
    converged=.false.
    Q20_1=1.0_rk ; Q20_2=2.0_rk ; Q20_3=3.0_rk
    Etot1=1.0_rk ; Etot2=2.0_rk ; Etot3=3.0_rk
    iter=1

    loop_iteration: do while (.not.finished)

       loop_isospin: do isospin=1,2

          if (debugmode=='high') call warning(modulename, subroutinename,'Start loop over isospin')

          shift=0

!
!         Initialization of expansion coefficients to 0
!         (array filled in by blocks ==> need to clean everything at every iteration)
!

          do iHF=1,basis%size
             sp_states(iHF,isospin)%coef%cylHO(:)=0.0
          end do
          sp_states(:,isospin)%pi=0

          loop_blocks: do k=1,Nblocks

             n=basis%block_size(k)

!
!            Hamiltonian matrix elements in the Harmonic Oscillator (HO) basis
!

             call hamiltonian_cylHO(h_HF(k)%HO(:,:),k,isospin,basis,HFfield(isospin),iter)

             if (exact_coul_exch .and. isospin == 2) h_HF(k)%HO = h_HF(k)%HO + &
                  mix*previous_Coul_exch(k)%HO + (1.0_rk-mix)*Coulomb_exchange_potential(k)%HO

!
!            Hamiltonian diagonalization
!

             allocate(C(n,n),Ctmp(n))

             call diagon(h_HF(k)%HO(:,:), h_HF(k)%eigenval(:), C, n)

             if (debugmode=='high') then
                call check_orthogonality(C, n)
                call check_diagonalization(h_HF(k)%HO(:,:), h_HF(k)%eigenval(:), C, n)
             end if

!
!            Fill in arrays with eigenvalues and eigenvectors
!            (largest expansion coeff. chosen to be positive)
!

             sp_energies(shift+1:shift+n,isospin)=h_HF(k)%eigenval(:)

             do i=1,n
                do j=1,n
                   Ctmp(j)=abs(C(j,i))
                end do
                j=maxloc(Ctmp,1)
                sp_states(i+shift,isospin)%coef%cylHO(1+shift:n+shift)=C(:,i)* &
                     sign(1.0_rk,C(j,i))
!
!                  Phase convention for Omega < 0 states: |-Omega> ~ T |Omega>
!
                if (basis%Omega(1+shift) < 0) then
                   sp_states(i+shift,isospin)%coef%cylHO(1+shift:n+shift) = &
                        sp_states(i+shift,isospin)%coef%cylHO(1+shift:n+shift) * &
                        (-1)**((basis%Omega(1+shift)+1)/2)
                end if
                sp_states(i+shift,isospin)%energy=h_HF(k)%eigenval(i)
                sp_states(i+shift,isospin)%block_rank=i
             end do

             deallocate(C,Ctmp)

             if (symmetries%parity) sp_states(1+shift:n+shift,isospin)%pi= &
                  basis%parity(1+shift)
             sp_states(1+shift:n+shift,isospin)%Omega=basis%Omega(1+shift)
             sp_states(1+shift:n+shift,isospin)%block=k

             shift=shift+n

          end do loop_blocks

!
!         Sort all eigenvalues in increasing order to define hole and particle states
!         and all eigenvectors accordingly
!

          call sort(sp_energies(:,isospin), basis%size, 'increasing', permut = permut(:,isospin))
         
          call permutation(sp_states(:,isospin), permut(:,isospin), &
               basis%size, 'normal')

!
!         Single-particle wave functions and their gradients at Gauss--Hermite--Laguerre
!         integration points
!

          do iHF = 1, basis%size
             call sp_wave_functions(sp_states(iHF,isospin), basis)
          end do

!
!         Occupation numbers: pairing or pure Hartree--Fock
!

!          if (iter==1 .and. model%initial_potential=='WS' .and. &
!               model%pairing/='seniority') then

             sp_states(1:Npart(isospin),isospin)%occupation=1.0_rk
             sp_states(Npart(isospin)+1:basis%size,isospin)%occupation=0.0
             nucleus%pairing%pairing_energy=0.0

!          else

!
!            Initializations:
!
!            * chemical potential to Fermi level (half sum of last occ. and 1st unocc. levels)
!            * pairing energy to 0
!            * pairing gaps and qp quantities to 0
!            * occupation probabilities to the Hatree-Fock values
!

             iHF=Npart(isospin)
             e_1st_unoccupied=sp_states(iHF+1,isospin)%energy
             if (iHF/=0) then
                fermi=0.5_rk*(sp_states(iHF,isospin)%energy+e_1st_unoccupied)
             else
                fermi = e_1st_unoccupied
             end if
             nucleus%pairing(isospin)%chemical_potential=fermi
             nucleus%pairing(isospin)%pairing_energy=0.0
             nucleus%pairing(isospin)%average_gap=0.0
             nucleus%pairing(isospin)%gap(:)=0.0
             nucleus%pairing(isospin)%min_qp_index=0
             nucleus%pairing(isospin)%min_qp_energy=0
             nucleus%pairing(isospin)%Vpair(:,:)=0.0

!
!            Determination of quasi-pairs (which would be Kramers degenerate
!            in the equal-filling approximation) and their energy
!

             partner(:)=0
             sp_states(:,isospin)%pair_partner(isospin)=0
             quasi_pair_energy(:)=1e6_rk
             sp_label(:)=0
             n=0

             do iHF=1,basis%size
                if (sp_states(iHF,isospin)%Omega<0) cycle
!               definition of time-reversed state
                psi_bar=time_reversal(sp_states(iHF,isospin)%coef,basis)
                ovlp=1e-09_rk
                partner(iHF)=0
!               search for sp state with opposite Omega and largest overlap with psi_bar
                do jHF=1,basis%size
                   if (sp_states(jHF,isospin)%Omega>0) cycle
                   tmp=dot_product(sp_states(jHF,isospin)%coef%cylHO,psi_bar%cylHO)
                   if (abs(tmp)>abs(ovlp)) then
                      ovlp=tmp
                      partner(iHF)=jHF
                   end if
                end do
                if (partner(iHF) == 0) then
                   write(out_unit,*) 'isospn=',isospin,'iHF=',iHF,' partner not found'
                   stop
                end if
                jHF=partner(iHF)
                if (jHF /=0) then
                   n=n+1
                   if (n>basis%size/2) call fatal(modulename,subroutinename,'Too many unfound partners')
                   quasi_pair_energy(n)=sp_states(iHF,isospin)%energy+sp_states(jHF,isospin)%energy
                   sp_label(n)=iHF
                   sp_states(iHF,isospin)%pair_partner(isospin)=jHF
                   sp_states(jHF,isospin)%pair_partner(isospin)=iHF
                else
                   write(out_unit,'(a,i4)') 'iHF=',iHF
                   call fatal(modulename,subroutinename,'Partner state not found')
                end if
                if (debugmode=='high' .and. sp_states(iHF,isospin)%energy<0) &
                     write(out_unit,'(a,i4,1x,a,i3,a2,a1,2x,a,i4,2(a,f12.6,1x))') &
                     'iHF=',iHF,'Omega-pi=',sp_states(iHF,isospin)%Omega,'/2', &
                     signparity(sp_states(iHF,isospin)%pi),' partner jHF=',partner(iHF), &
                     ' <jHF|iHF-bar>=',ovlp,' quasi-pair energy=',quasi_pair_energy(n)
             end do

             call sort(quasi_pair_energy(:), basis%size/2, 'increasing', permut=order(:))

             tmp_label(:)=0
             do i=1,basis%size/2
                tmp_label(i)=sp_label(order(i))
             end do
             sp_label(:)=tmp_label(:)

!
!            Occupation numbers in the Hartree--Fock approximation:
!

!
!               Finding the closest rank of the blocked state to the Fermi level
!
                if (nucleus%fix_rank == 'n') then
                   max_rank = maxval(basis%block_size(:),1)
                   if (max_rank <= 0) call fatal(modulename,subroutinename,"max_rank <= 0")
                   allocate(ph_energy(max_rank))
                   ph_energy(:) = ph_energy_min
                   if (symmetries%parity) then
                      do j = 1, nqp
                         if (nucleus%K(isospin,j) == 0) cycle
                         do m = 1, max_rank
                            do jHF = 1, basis%size
                               if (sp_states(jHF,isospin)%Omega /= nucleus%K(isospin,j) .or. &
                                    sp_states(jHF,isospin)%pi /= nucleus%pi(isospin,j)) cycle
                               if (sp_states(jHF,isospin)%block_rank == m) ph_energy(m) = &
                                    abs(sp_states(jHF,isospin)%energy-fermi)
                            end do
                         end do
                         if (minval(ph_energy) == ph_energy_min) call fatal(modulename,subroutinename, &
                              'Rank not found')
                         nucleus%rank(isospin,j) = minloc(ph_energy,dim=1)
                         if (debugmode /= 'low') write(out_unit,'(4x,3(a,I2))') 'Blocked K^pi state =', &
                              nucleus%K(isospin,j),'/2,  parity = ',nucleus%pi(isospin,j), &
                              ',  Rank = ', nucleus%rank(isospin,j)
                         do m = 1, max_rank
                            if (m == nucleus%rank(isospin,j)) cycle
                            if (abs(ph_energy(m)-ph_energy(nucleus%rank(isospin,j))) <= 0.5_rk) then
                               write(out_unit,'(4x,2(a,i2),a)') 'NOTE: The distances from Fermi level of ranks ', &
                                    nucleus%rank(isospin,j),' and ',m,' differ by less than 500 keV'
                               write(out_unit,'(3x,a,i2,a,f12.6)') 'Distance from Fermi level for rank ', &
                                    nucleus%rank(isospin,j), ' = ', ph_energy(nucleus%rank(isospin,j))
                               write(out_unit,'(3x,a,i2,a,f12.6)') 'Distance from Fermi level for rank ', &
                                    m, ' = ', ph_energy(m)
                            end if
                         end do
                      end do
                   else
                      do j = 1, nqp
                         if (nucleus%K(isospin,j) == 0) cycle
                         do m = 1, max_rank
                            do jHF = 1, basis%size
                               if (sp_states(jHF,isospin)%Omega /= nucleus%K(isospin,j)) cycle
                               if (sp_states(jHF,isospin)%block_rank == m) ph_energy(m) = &
                                    abs(sp_states(jHF,isospin)%energy-fermi)
                            end do
                         end do
                         if (minval(ph_energy) == ph_energy_min) call fatal(modulename,subroutinename, &
                              'Rank not found')
                         nucleus%rank(isospin,j) = minloc(ph_energy,dim=1)
                         if (debugmode /= 'low') write(out_unit,'(4x,3(a,I2))') 'Blocked K^pi state =', &
                              nucleus%K(isospin,j),'/2,  parity = ',nucleus%pi(isospin,j), &
                              ',  Rank = ', nucleus%rank(isospin,j)
                         do m = 1, max_rank
                            if (m == nucleus%rank(isospin,j)) cycle
                            if (abs(ph_energy(m)-ph_energy(nucleus%rank(isospin,j))) <= 0.5_rk) then
                               write(out_unit,'(4x,2(a,i2),a)') 'NOTE: The distances from Fermi level of ranks ', &
                                    nucleus%rank(isospin,j),' and ',m,' differ by less than 500 keV'
                               write(out_unit,'(3x,a,i2,a,f12.6)') 'Distance from Fermi level for rank ', &
                                    nucleus%rank(isospin,j), ' = ', ph_energy(nucleus%rank(isospin,j))
                               write(out_unit,'(3x,a,i2,a,f12.6)') 'Distance from Fermi level for rank ', &
                                    m, ' = ', ph_energy(m)
                            end if
                         end do
                      end do
                   end if
                   deallocate(ph_energy)
                end if
                
!
!            Determination of the blocked states
!

             do i=1,nqp
                i_odd(isospin,i)=0
                if (nucleus%K(isospin,i) == 0) cycle
                loop_iHF: do iHF=1,basis%size
                   if (sp_states(iHF,isospin)%Omega==nucleus%K(isospin,i) .and. &
                        sp_states(iHF,isospin)%block_rank==nucleus%rank(isospin,i)) then
                      i_odd(isospin,i)=iHF
                      if (symmetries%parity) then
                         if (sp_states(iHF,isospin)%pi==nucleus%pi(isospin,i)) then
                            i_odd(isospin,i)=iHF
                            exit loop_iHF
                         end if
                      else
                         exit loop_iHF
                      end if
                   end if
                end do loop_iHF
             end do
             
             if (debugmode /= 'low') then
                do i=1,nqp
                   if (nucleus%K(isospin,i) == 0) cycle
                   write(out_unit,'(4(a,i4,2x))') &
                        'isospin=',isospin,'i=',i,'i_odd=',i_odd(isospin,i), &
                        'partner(i_odd)=',sp_states(i_odd(isospin,i),isospin)%pair_partner(isospin)
                end do
             end if

!               Initialization of occupation numbers and type component

                sp_states(:,isospin)%occupation = 0.0
                sp_states(:,isospin)%type = 0

!               Occupation numbers of the blocked states

                do i=1,nqp
                   if (i_odd(isospin,i)/=0) then
                      if (model%blocking == 'SCB') then
                         sp_states(i_odd(isospin,i),isospin)%occupation = 1.0_rk
                      else if (model%blocking == 'EFA') then
                         sp_states(i_odd(isospin,i),isospin)%occupation = 0.5_rk
                         jHF = sp_states(i_odd(isospin,i),isospin)%pair_partner(isospin)
                         sp_states(jHF,isospin)%occupation = 0.5_rk
                      else
                         call fatal(modulename,subroutinename,'Unrecognized blocking scheme ' &
                              //trim(model%blocking))
                      end if
                   else if (nucleus%K(isospin,i) /= 0) then
                      call fatal(modulename,subroutinename,'Blocked state not found')
                   end if
                end do

!            Occupation numbers of remaining, lowest nb_qp(isospin) quasi-pairs 
!            of states set to 1

             nb_qp(isospin) = (Npart(isospin) - nb_blocked_part(isospin))/2
             jHF = 1
             do k = 1, nb_qp(isospin)
                filled_qp = .false.
                loop_fill_qp: do j = jHF, basis%size/2
                   iHF = sp_label(j)
                   do i = 1,nqp
                      if (iHF == i_odd(isospin,i) .or. partner(iHF) == i_odd(isospin,i)) cycle loop_fill_qp
                   end do
                   sp_states(iHF,isospin)%occupation = 1.0_rk
                   sp_states(partner(iHF),isospin)%occupation = 1.0_rk
                   filled_qp = .true.
                   jHF=j+1
                   exit loop_fill_qp
                end do loop_fill_qp
                if (.not.filled_qp) then
                   write(out_unit,'(a,i2,2(1x,a,i4),2x,a,3i5)') &
                        'isospin=',isospin,'k=',k,'jHF=',jHF,'i_odd(isospin,:)=',i_odd(isospin,:)
                   call fatal(modulename,subroutinename,"Quasi-pair not filled")
                end if
                if (debugmode=='high') write(out_unit,'(a1,2x,2(i4,1x),a2,a1,2(i4,1x),f12.6)') &
                     nucleon_type(isospin),k,sp_states(iHF,isospin)%Omega,'/2', &
                     signparity(sp_states(iHF,isospin)%pi),iHF,partner(iHF),quasi_pair_energy(j)
             end do

!
!            Define HF_level, HF_state and type components of sp_states
!

             i=0
             do k=1,basis%size
                sp_states(k,isospin)%HF_level = k
                sp_states(k,isospin)%HF_state = k
                if (sp_states(k,isospin)%occupation > 0.5_rk) then
                   sp_states(k,isospin)%type = 1
                   i=i+1
                else
                   sp_states(k,isospin)%type = 0
                end if
             end do

             if (model%blocking == 'EFA') then
                do k = 1, nqp
                   if (i_odd(isospin,k) == 0) cycle
                   jHF = sp_states(i_odd(isospin,k),isospin)%pair_partner(isospin)
                   sp_states(jHF,isospin)%type = 0
                   i = i - 1
                end do
             end if

             if (i /= Npart(isospin)) then
                do iHF = 1, basis%size
                   write(out_unit,'(a1,2x,2(i4,1x),a2,a1,3(i4,1x))') &
                        nucleon_type(isospin),iHF,sp_states(iHF,isospin)%Omega,'/2', &
                        signparity(sp_states(iHF,isospin)%pi),iHF,partner(iHF), &
                        sp_states(iHF,isospin)%type
                end do
                write(out_unit,'(/a,i2,1x,a,i4,a,i4)') 'isospin=',isospin,&
                     'number of hole states =', i, ' whereas Npart=', Npart(isospin)
                call fatal(modulename,subroutinename, &
                     'Problem with the particle or hole character of s.p. states')
             end if

!          end if

       end do loop_isospin

!       if (iter>1) then

          write(out_unit,'(/a,i4/)') 'HFHTDA iteration #',iter

!
!         Occupation numbers in the HTDA approach
!

          nval = 1
          write(out_unit,'(a)') 'Starting HTDA calculations'
          call HTDA_PROCESS(storage, nucleus, symmetries, model, basis, sp_states, &
               Vres, Hiqp, Hhtda, HTDA_basis, Useful_Space, RPA_VS, Pairing_VS, &
               Total_window, nval, eigenval_H, eigenvec_H)

          nucleus%pairing(:)%pairing_energy = 0.5_rk*eigenval_H(1) ! By convention so that Epair = Ecorr

!       end if

!
!      Local one-body densities in the coordinate space at Gauss--Hermite--Laguerre
!      integration points
!

       do isospin = 1, 2
          call local_one_body_densities(sp_states(:,isospin), basis, densities(isospin), &
               isospin, iter)
       end do

!
!      Renormalization of local densities
!

       call renormalize_densities(densities, basis, nucleus)

!
!      Direct Coulomb potential
!

       call Coulomb(HFfield(2), densities, basis, symmetries)

!
!      Multipole moments of the nuclear density (Q20 enough, save time by 
!      calculating the other moments after convergence)
!

       call moments(nucleus, densities, sp_states, basis, symmetries)

!
!      beta parameters of an isodensity surface (spherical expansion of the 
!      nuclear radius as a function of theta)
!

!       call nuclear_shape(sp_states, basis, symmetries, nucleus)

!
!      Total energy and exact Coulomb exchange potential if requested
!

       if (exact_coul_exch) then
          do k=1,basis%Nblocks
             previous_Coul_exch(k)%HO = Coulomb_exchange_potential(k)%HO
          end do
          call Exact_Coulomb_potential('exc', Coulomb_exchange_potential, basis, sp_states, nucleus%Z)
          call energy(nucleus, model, basis, HFfield(2), densities, sp_states, Coulomb_exchange_potential)
       else
          call energy(nucleus, model, basis, HFfield(2), densities, sp_states)
       end if

!
!      Hartree-Fock field (central, effective mass and spin-orbite potentials)
!

       call potential(nucleus, model, symmetries, basis, numerical_parameters, &
            HFfield, densities, exact_coul_exch)

!
!      Displays results on standard output
!

       call display_iteration_results(nucleus, model, symmetries, iter, basis%size, sp_states)

!
!      Observables
!

!       call observables(nucleus, model, symmetries, basis, sp_states, densities, HFfield)

!
!      Test for convergence
!

       Etot3=Etot2
       Etot2=Etot1
       Etot1=nucleus%energy%Etot
       Q20_3=Q20_2
       Q20_2=Q20_1
       Q20_1=nucleus%Q20(0)

       converged = ( abs(Etot1-Etot2)+abs(Etot2-Etot3) <= &
            2.0_rk*numerical_parameters%convergence%E_precision ) .and. &
            ( abs(Q20_1-Q20_2)+abs(Q20_2-Q20_3) <= &
            2.0_rk*numerical_parameters%convergence%Q20_precision )
       iter=iter+1
       if (iter > numerical_parameters%convergence%Nitermax .or. converged) &
            finished=.true.

!
!      Display detailed information about the single-particle states
!

       if (debugmode=='high'.or.finished) then
          do isospin=1,2
             do iHF=1,basis%size
                if (sp_states(iHF,isospin)%energy>10.0_rk) exit
                write(2,'(i4,2x,a1,i4,1x,f12.6,1x,i4,a2,1x,a1,1x,2(f12.6,1x),2(1x,i4))') &
                     iter-1,nucleon_type(isospin),iHF, &
                     sp_states(iHF,isospin)%energy, &
                     sp_states(iHF,isospin)%Omega,'/2 ', &
                     signparity(sp_states(iHF,isospin)%pi), &
                     sp_states(iHF,isospin)%occupation, &
                     nucleus%pairing(isospin)%gap(iHF), &
                     sp_states(iHF,isospin)%block, &
                     sp_states(iHF,isospin)%block_rank
                do i=1,basis%size
                   if (abs(sp_states(iHF,isospin)%coef%cylHO(i))<1e-6_rk) cycle
                   write(41,'(i4,2x,a1,i4,1x,f13.6,1x,i4,a2,1x,a1,1x,i4,2x,e15.8,3(1x,i4),i4,a2)') &
                        iter-1,nucleon_type(isospin),iHF,sp_states(iHF,isospin)%energy, &
                        sp_states(iHF,isospin)%Omega,'/2', &
                        signparity((-1)**(basis%nz(i)+basis%lambda(i))),i, &
                        sp_states(iHF,isospin)%coef%cylHO(i),basis%nz(i), &
                        2*basis%nr(i)+abs(basis%Lambda(i)),basis%Lambda(i),basis%ns(i),'/2'
                end do
             end do
             write(2,*)
          end do
          write(2,*)
       end if

    end do loop_iteration

!---------------------------
!
!   End of iteration loop
!
!---------------------------

    close(2)
    close(41)

!---------------------------------------------
!
!   Memory deallocation for local variables
!
!---------------------------------------------

    do k=1,Nblocks
       deallocate(h_HF(k)%HO,h_HF(k)%eigenval)
    end do

    deallocate(h_HF,sp_energies,Vpair)

    deallocate(psi_bar%cylHO)

    call DATE_AND_TIME(values = end)
    call elapsed_time(start,end,"HFBCS-like calculation with HTDA-deduced u's and v's")

!   Observables calculated in the BCS-like wavefunction with HTDA-deduced
!   u's and v's

    call observables(nucleus, model, symmetries, basis, sp_states, &
         densities, HFfield)
       
!   HTDA magnetic moments and spin-quenching factor

    if (nucleus%magnetic_moment) then
       do isospin = 1, 2
          do k = 1, Total_window%max_bound(isospin)
             Useful_Space(isospin,k)%avg_sz = sp_states(k,isospin)%avg_sz
          end do
       end do
       call magnetic_moment(nucleus, symmetries, basis, sp_states, HTDA_basis, Useful_space, &
            Total_window, nval, eigenval_H, eigenvec_H)
    end if

  end subroutine solve_HFHTDA_equations

  !=====================================================================
  !                       SUBROUTINE HAMILTONIAN_CYLHO
  !---------------------------------------------------------------------
  ! Calculates the matrix elements of the Skyrme-Hartree-Fock Hamiltonian 
  ! in the cylindrical harmonic-oscillator basis (in block k). 
  !=====================================================================

  subroutine hamiltonian_cylHO(H, k, isospin, basis, HFfield, iter)

!
!   Arguments
!
    integer,intent(in)                                                  :: k,isospin
    type(cylindrical_basis),intent(in)                                  :: basis
    real(kind=rk), &
         dimension(basis%block_size(k),basis%block_size(k)),intent(out) :: H
    type(field),intent(in)                                 :: HFfield
    integer,optional,intent(in)                                         :: iter

!
!   Local variables
!
    integer :: nthreads
    integer :: n,shift,ih,il
    integer :: i,i1,i2,j,j1,j2,l1,l2,s1,s2,nz1,nz2,nr1,nr2,k1,k2,Delta1,Delta2
    real(kind=rk) :: norm,xi,eta,w,sum_H,sum_U,sum_K,sum_SO,sum_S,S_z,S_rho,S_phi, &
         sum_alpha,alpha_z,alpha_rho,alpha_phi,sum_C,C_z,C_rho,C_phi,sum_tensor_so
    real(kind=rk) :: sum_D,D_z,D_rho,D_phi
    type(spinor) :: chi1,chi2,phi1,phi2,ket_z,ket_rho,ket_phi,ket,bra
    type(vector_spinor) :: grad_phi1,grad_phi2,sigma_grad_phi2
    type(vector_spin_operator) :: sigma
    complex(kind=rk),dimension(2,2) :: Csigma,Jsigma_z,Jsigma_rho,Jsigma_phi

    norm=fz/(2.0_rk*bz*bp2)

!
!   Definition of the cylindrical Pauli matrices
!
    do i=1,2
       sigma%z(i,i)=-cmplx((-1)**i)
       sigma%rho(i,i)=cmplx(0)
       sigma%phi(i,i)=cmplx(0)
    end do

    sigma%z(1,2)=cmplx(0)
    sigma%rho(1,2)=cmplx(1) !*exp(-i * phi) with phi = 0 (axial symmetry)
    sigma%phi(1,2)=-i_cplx

    sigma%z(2,1)=cmplx(0)
    sigma%rho(2,1)=cmplx(1)
    sigma%phi(2,1)=i_cplx

!
!   Number of basis states preceding block # k
!

    n=basis%block_size(k)
    shift=0
    do j=1,k
       shift=shift+basis%block_size(j)
    end do
    shift=shift-n

!
!   Initialization of H to 0
!

    H(:,:)=0.

    nthreads = n
    call setup_openmp_parameters(1, n, nthreads, OMP_parameters%set_bounds, &
            OMP_parameters%nthreads)
    call omp_set_num_threads(OMP_parameters%nthreads)
    !$OMP PARALLEL DO PRIVATE(i1,j1,l1,s1,nz1,nr1,k1,Delta1,chi1) &
    !$OMP PRIVATE(i2,j2,l2,s2,nz2,nr2,k2,Delta2,chi2) &
    !$OMP PRIVATE(sum_H,sum_U,sum_K,sum_SO,sum_S,sum_alpha,sum_C) &
    !$OMP PRIVATE(sum_tensor_so,sum_D,ih,il,xi,eta,w) &
    !$OMP PRIVATE(phi1,phi2,grad_phi1,grad_phi2,sigma_grad_phi2) &
    !$OMP PRIVATE(S_z,S_rho,S_phi,ket_z,ket_rho,ket_phi,ket) &
    !$OMP PRIVATE(alpha_z,alpha_rho,alpha_phi,C_z,C_rho,C_phi,Csigma) &
    !$OMP PRIVATE(Jsigma_z,Jsigma_rho,Jsigma_phi,D_z,D_rho,D_phi,bra)

!
!   Matrix elements in block k
!

    loop_i1: do i1=1,n

       j1=i1+shift
       l1=basis%Lambda(j1)
       s1=basis%ns(j1)
       nz1=basis%nz(j1)
       nr1=basis%nr(j1)
       k1=2*l1+s1
       Delta1=(l1-abs(l1))/2
       chi1=spinor(kronecker(s1,1),kronecker(s1,-1))

       loop_i2: do i2=1,i1

          j2=i2+shift
          l2=basis%Lambda(j2)
          s2=basis%ns(j2)
          nz2=basis%nz(j2)
          nr2=basis%nr(j2)
          k2=2*l2+s2
          Delta2=(l2-abs(l2))/2
          chi2=spinor(kronecker(s2,1),kronecker(s2,-1))

          if (k1/=k2) cycle

          sum_H=0.0
          sum_U=0.0
          sum_K=0.0
          sum_SO=0.0
          sum_S=0.0
          sum_alpha=0.0
          sum_C=0.0
          sum_tensor_so=0.0
          sum_D=0.0

          do ih=iHmin,basis%NGz/2

             if (ih==0) cycle
             xi=basis%xHerm(ih)

             do il=1,basis%NGr

                w=basis%G_Herm(ih)*basis%G_Lag(il)
                eta=basis%xLag(il)

                phi1=sqrt(2.0_rk*bz)*bp*basis%P_Herm(nz1,ih)*basis%P_Lag(nr1,abs(l1),il)* &
                     (-1)**Delta1*sqrt(eta**abs(l1))*exp(-0.5_rk*(xi**2+eta))*chi1
                phi2=sqrt(2.0_rk*bz)*bp*basis%P_Herm(nz2,ih)*basis%P_Lag(nr2,abs(l2),il)* &
                     (-1)**Delta2*sqrt(eta**abs(l2))*exp(-0.5_rk*(xi**2+eta))*chi2
                sum_U=sum_U+w*HFfield%central(ih,il)*real(scalar_product(phi1,phi2))

                grad_phi1=grad(j1,ih,il,basis)
                grad_phi2=grad(j2,ih,il,basis)
                sum_K=sum_K+w*HFfield%effective_mass(ih,il)* &
                     real(scalar_product(grad_phi1,grad_phi2))

                sigma_grad_phi2=vector_product(sigma,grad_phi2)
                sum_SO=sum_SO+w*HFfield%spin_orbit(ih,il)* &
                     aimag(scalar_product(grad_phi1,sigma_grad_phi2))

                S_z=HFfield%S%z(ih,il)
                S_rho=HFfield%S%rho(ih,il)
                S_phi=HFfield%S%phi(ih,il)
                sum_S=sum_S+w*real(scalar_product(phi1,S_z*(sigma%z*phi2))+ &
                     scalar_product(phi1,S_rho*(sigma%rho*phi2))+ &
                     scalar_product(phi1,S_phi*(sigma%phi*phi2)))

                ket_z=spin_operator_action(sigma%z,phi2)
                ket_rho=spin_operator_action(sigma%rho,phi2)
                ket_phi=spin_operator_action(sigma%phi,phi2)
                ket=spin_operator_action(sigma%z,grad_phi2%z) + &
                     spin_operator_action(sigma%rho,grad_phi2%rho) + &
                     spin_operator_action(sigma%phi,grad_phi2%phi)
                sum_S=sum_S+w*HFfield%div_s(ih,il)*(real(scalar_product(grad_phi1%z,ket_z)) + &
                     real(scalar_product(grad_phi1%rho,ket_rho)) + &
                     real(scalar_product(grad_phi1%phi,ket_phi)) + &
                     real(scalar_product(phi1,ket)))
 
                alpha_z=HFfield%alpha%z(ih,il)
                alpha_rho=HFfield%alpha%rho(ih,il)
                alpha_phi=HFfield%alpha%phi(ih,il)
                sum_alpha=sum_alpha+w*aimag(scalar_product(phi1,alpha_z*grad_phi2%z)+ &
                     scalar_product(phi1,alpha_rho*grad_phi2%rho)+ &
                     scalar_product(phi1,alpha_phi*grad_phi2%phi))

                C_z=HFfield%C%z(ih,il)
                C_rho=HFfield%C%rho(ih,il)
                C_phi=HFfield%C%phi(ih,il)
                Csigma=C_z*sigma%z+C_rho*sigma%rho+C_phi*sigma%phi
                sigma_grad_phi2=spin_operator_action(Csigma,grad_phi2)
                sum_C=sum_C+w*real(scalar_product(grad_phi1,sigma_grad_phi2))

                Jsigma_z=HFfield%tensor_so(1)%z(ih,il)*sigma%z + &
                     HFfield%tensor_so(1)%rho(ih,il)*sigma%rho + &
                     HFfield%tensor_so(1)%phi(ih,il)*sigma%phi
                Jsigma_rho=HFfield%tensor_so(2)%z(ih,il)*sigma%z + &
                     HFfield%tensor_so(2)%rho(ih,il)*sigma%rho + &
                     HFfield%tensor_so(2)%phi(ih,il)*sigma%phi
                Jsigma_phi=HFfield%tensor_so(3)%z(ih,il)*sigma%z + &
                     HFfield%tensor_so(3)%rho(ih,il)*sigma%rho + &
                     HFfield%tensor_so(3)%phi(ih,il)*sigma%phi

                ket_z=spin_operator_action(Jsigma_z,grad_phi2%z)
                ket_rho=spin_operator_action(Jsigma_rho,grad_phi2%rho)
                ket_phi=spin_operator_action(Jsigma_phi,grad_phi2%phi)
                sum_tensor_so=sum_tensor_so+w*(aimag(scalar_product(phi1,ket_z)) + &
                     aimag(scalar_product(phi1,ket_rho))+aimag(scalar_product(phi1,ket_phi)))

                ket_z=spin_operator_action(Jsigma_z,phi2)
                ket_rho=spin_operator_action(Jsigma_rho,phi2)
                ket_phi=spin_operator_action(Jsigma_phi,phi2)
                sum_tensor_so=sum_tensor_so-w*(aimag(scalar_product(grad_phi1%z,ket_z)) + &
                     aimag(scalar_product(grad_phi1%rho,ket_rho)) + &
                     aimag(scalar_product(grad_phi1%phi,ket_phi)))

                D_z=HFfield%D%z(ih,il)
                D_rho=HFfield%D%rho(ih,il)
                D_phi=HFfield%D%phi(ih,il)

                ket=D_z*grad_phi2%z+D_rho*grad_phi2%rho+D_phi*grad_phi2%phi
                bra=spin_operator_action(sigma%z,grad_phi1%z)+ &
                     spin_operator_action(sigma%rho,grad_phi1%rho)+ &
                     spin_operator_action(sigma%phi,grad_phi1%phi)
                sum_D=sum_D+0.5_rk*w*real(scalar_product(bra,ket))

                bra=D_z*grad_phi1%z+D_rho*grad_phi1%rho+D_phi*grad_phi1%phi ! D field is real
                ket=spin_operator_action(sigma%z,grad_phi2%z)+ &
                     spin_operator_action(sigma%rho,grad_phi2%rho)+ &
                     spin_operator_action(sigma%phi,grad_phi2%phi)
                sum_D=sum_D+0.5_rk*w*real(scalar_product(bra,ket))

             end do

          end do

          sum_H=norm*(sum_K + sum_U + sum_SO + sum_S + sum_alpha + sum_C + sum_tensor_so + sum_D)
          H(i1,i2)=sum_H ! Lower triangle of matrix H
          if (debugmode=='high'.and.present(iter)) then
             write(40,'(12(i4,2x),f12.6)') iter,isospin,j1,j2,nz1,nr1,l1,s1,nz2,nr2,l2,s2,H(i1,i2)
             write(67,'(3(i4,2x),6(f12.8,1x))') iter,j1,j2,sum_K,sum_U,sum_SO,sum_S,sum_alpha,sum_H
          end if

       end do loop_i2

    end do loop_i1

    !$OMP END PARALLEL DO

    do i1 = 1, n-1
       do i2 = i1+1, n
          H(i1,i2) = H(i2,i1) ! Upper triangle of matrix H
       end do
    end do

    if (debugmode=='high'.and.present(iter)) write(40,*)

    return
    
  end subroutine hamiltonian_cylHO


  !=====================================================================
  !                       SUBROUTINE h_odd_perturbative
  !---------------------------------------------------------------------
  ! Calculates the first-order eigenstates and eigenvalues of 
  !      h_HF = h_even + h_odd
  ! where h_odd is treated as a perturbation with respect to h_even.
  ! Calculates the expectation value of sz to deduce the gs quenching at 
  ! first order in h_odd.
  !=====================================================================

  subroutine h_odd_perturbative(nucleus, model, symmetries, basis, HFfield, &
       densities, iter, iunit, calculation_type)

!
!   Arguments
!
    type(nucleus_properties),intent(inout) :: nucleus
    type(modelization),intent(in) :: model
    type(symmetry_properties),intent(in) :: symmetries
    type(cylindrical_basis),intent(in) :: basis
    type(field),dimension(2),intent(in) :: HFfield
    type(local_density),dimension(2),intent(in) :: densities
    integer,intent(in) :: iter,iunit
    character(len=*) :: calculation_type
!
!   Local variables
!
    character(len=18),parameter :: subroutinename='h_odd_perturbative'
    integer, dimension(2) :: Npart,nb_qp
    integer :: Nsph,ih,il,NGz,NGr
    integer :: Nblocks,iblock
    integer :: k,n,shift
    integer :: mem_stat
    integer :: isospin,i,j,iHF,jHF
    type(field_Bi_contributions) :: h_HF_mtx_elem,quenching_factor
    real(kind=rk) :: sum_sz,expval_sz_unperturbed,ei,ej,avg_sz,avg_lz,hodd_sz,gs_q, &
         expval_sz_phi_i
    type(field_Bi_contributions),dimension(2) :: expval_sz,delta
    type(field_Bi_contributions) :: s_odd
    real(kind=rk),dimension(2) :: expval_sz_core
    real(kind=rk) :: B3, B4, B5, B6, B10, B11, B12, B13, B14, B15, B18, B19, x0, x1, x2, x3
    real(kind=rk) :: delta_S_t0, delta_S_t0x0, delta_S_t1, delta_S_t1x1, delta_S_t2, delta_S_t2x2, &
         delta_S_t3, delta_S_t3x3, delta_S_W
    real(kind=rk) :: delta_A_t1, delta_A_t1x1, delta_A_t2, delta_A_t2x2, delta_A_W
    real(kind=rk) :: delta_C_t1, delta_C_t1x1, delta_C_t2, delta_C_t2x2
    real(kind=rk) :: delta_t0, delta_t0x0, delta_t1, delta_t1x1, delta_t2, delta_t2x2, &
         delta_t3, delta_t3x3, delta_W
    real(kind=rk) :: delta_1S0, delta_3S1, delta_1P1, delta_3P0, delta_3P1, delta_3P2

    type(nucleus_properties) :: local_nucleus
    integer,dimension(2,nqp) :: i_odd_local
    type(field),dimension(2) :: HFfield_even
    integer,dimension(basis%size,2) :: permut_inv
    type(sp_state),dimension(basis%size,2) :: sp_states
    type(operator),dimension(:),allocatable  :: h_HF        ! Hartree-Fock Hamiltonian
    real(kind=rk),dimension(:,:),allocatable :: C           ! Expansion coefficients
    real(kind=rk),dimension(:),allocatable   :: Ctmp
    real(kind=rk),dimension(:,:),allocatable :: sp_energies ! Single-particle energies
    real(kind=rk),dimension(:,:),allocatable :: Vpair       ! Pairing matrix elements

    real(kind=rk) :: e_1st_unoccupied,fermi,Vpair_max
    type(expansion) :: psi_bar
    real(kind=rk) :: ovlp,tmp
    integer,dimension(basis%size) :: partner,order
    real(kind=rk),dimension(basis%size/2) :: quasi_pair_energy
    integer,dimension(basis%size/2) :: sp_label,tmp_label
    real(kind=rk),dimension(basis%size,basis%size) :: np_overlaps
    logical :: filled_qp

    integer :: count, m, max_rank
    real(kind=rk),parameter :: ph_energy_min = 1.0e12_rk
    real(kind=rk),dimension(:),allocatable :: ph_energy ! Particle-hole excitation energy wrt Fermi level


    if (mod(nucleus%N,2) == 0 .and. mod(nucleus%Z,2) == 0) return

    write(iunit,'(a,i3)') 'Iteration # ',iter

    Npart(:)=(/ nucleus%N , nucleus%Z /)
    local_nucleus = nucleus
    i_odd_local = i_odd

!
!   Even part of the HF field
!
    NGz=basis%NGz
    NGr=basis%NGr
    do isospin = 1, 2
       allocate(HFfield_even(isospin)%central(iHmin:NGz/2,NGr), &
            HFfield_even(isospin)%effective_mass(iHmin:NGz/2,NGr), &
            HFfield_even(isospin)%spin_orbit(iHmin:NGz/2,NGr), &
            HFfield_even(isospin)%Coulomb(iHmin:NGz/2,NGr))
       allocate(HFfield_even(isospin)%S%z(iHmin:NGz/2,NGr), &
            HFfield_even(isospin)%S%rho(iHmin:NGz/2,NGr), &
            HFfield_even(isospin)%S%phi(iHmin:NGz/2,NGr))
       allocate(HFfield_even(isospin)%alpha%z(iHmin:NGz/2,NGr), &
            HFfield_even(isospin)%alpha%rho(iHmin:NGz/2,NGr), &
            HFfield_even(isospin)%alpha%phi(iHmin:NGz/2,NGr))
       allocate(HFfield_even(isospin)%C%z(iHmin:NGz/2,NGr), &
            HFfield_even(isospin)%C%rho(iHmin:NGz/2,NGr), &
            HFfield_even(isospin)%C%phi(iHmin:NGz/2,NGr))
       allocate(HFfield_even(isospin)%D%z(iHmin:NGz/2,NGr), &
            HFfield_even(isospin)%D%rho(iHmin:NGz/2,NGr), &
            HFfield_even(isospin)%D%phi(iHmin:NGz/2,NGr))
       do i=1,3
          allocate(HFfield_even(isospin)%tensor_so(i)%z(iHmin:NGz/2,NGr), &
               HFfield_even(isospin)%tensor_so(i)%rho(iHmin:NGz/2,NGr), &
               HFfield_even(isospin)%tensor_so(i)%phi(iHmin:NGz/2,NGr))
       end do
       allocate(HFfield_even(isospin)%div_s(iHmin:NGz/2,NGr))
       HFfield_even(isospin)%central = HFfield(isospin)%central
       HFfield_even(isospin)%Coulomb = HFfield(isospin)%Coulomb
       HFfield_even(isospin)%effective_mass = HFfield(isospin)%effective_mass
       HFfield_even(isospin)%spin_orbit = HFfield(isospin)%spin_orbit
       do i = 1, 3
          HFfield_even(isospin)%tensor_so(i)%z = HFfield(isospin)%tensor_so(i)%z
          HFfield_even(isospin)%tensor_so(i)%rho = HFfield(isospin)%tensor_so(i)%rho
          HFfield_even(isospin)%tensor_so(i)%phi = HFfield(isospin)%tensor_so(i)%phi
       end do
       HFfield_even(isospin)%S%z(:,:) = 0.0
       HFfield_even(isospin)%S%rho(:,:) = 0.0
       HFfield_even(isospin)%S%phi(:,:) = 0.0
       HFfield_even(isospin)%alpha%z(:,:) = 0.0
       HFfield_even(isospin)%alpha%rho(:,:) = 0.0
       HFfield_even(isospin)%alpha%phi(:,:) = 0.0
       HFfield_even(isospin)%C%z(:,:) = 0.0
       HFfield_even(isospin)%C%rho(:,:) = 0.0
       HFfield_even(isospin)%C%phi(:,:) = 0.0
       HFfield_even(isospin)%D%z(:,:) = 0.0
       HFfield_even(isospin)%D%rho(:,:) = 0.0
       HFfield_even(isospin)%D%phi(:,:) = 0.0
       HFfield_even(isospin)%div_s(:,:) = 0.0
    end do

!
!   Memeory allocation
!

    Nsph=(n_max+1)*(l_max+1)**2
    do isospin=1,2
       do i=1,basis%size
          allocate(sp_states(i,isospin)%coef%cylHO(basis%size),stat=mem_stat)
          if (mem_stat/=0) call fatal(modulename, subroutinename, &
               'Problem with allocate(sp_states%coef%cylHO)')
          allocate(sp_states(i,isospin)%coef%sphHO(Nsph),stat=mem_stat)
          if (mem_stat/=0) call fatal(modulename, subroutinename, &
               'Problem with allocate(sp_states%coef%sphHO)')
          allocate(sp_states(i,isospin)%wf(iHmin:basis%NGz/2,basis%NGr),stat=mem_stat)
          if (mem_stat/=0) call fatal(modulename, subroutinename, &
               'Problem with allocate(sp_states%wf)')
          sp_states(i,isospin)%tau=2*MOD(isospin,2) - 1
          sp_states(i,isospin)%Omega=0
          sp_states(i,isospin)%pi=0
          sp_states(i,isospin)%energy=0.0
          sp_states(i,isospin)%occupation=0.0
          sp_states(i,isospin)%avg_j=0.0
          sp_states(i,isospin)%avg_l=0.0
          sp_states(i,isospin)%avg_pi=0.0
          sp_states(i,isospin)%block=0
          sp_states(i,isospin)%block_rank=0
          sp_states(i,isospin)%pair_partner(:)=0
          sp_states(i,isospin)%coef%cylHO(:)=0.0
          sp_states(i,isospin)%coef%sphHO(:)=0.0
          do ih=iHmin,basis%NGz/2
             if (ih==0) cycle
             do il=1,basis%NGr
                sp_states(i,isospin)%wf(ih,il)=spinor(cmplx(0),cmplx(0))
             end do
          end do
       end do
    end do

    Nblocks=basis%Nblocks

    allocate(h_HF(Nblocks),stat=mem_stat)
    if (mem_stat/=0) call fatal(modulename, subroutinename, &
         'Problem with allocate(h_HF)')

    do k=1,Nblocks
       n=basis%block_size(k)
       allocate(h_HF(k)%HO(n,n),h_HF(k)%eigenval(n),stat=mem_stat)
       if (mem_stat/=0) call fatal(modulename, subroutinename, &
            'Problem with allocate(h_HF%HO etc...)')
    end do

    allocate(sp_energies(basis%size,2),stat=mem_stat)
    if (mem_stat/=0) call fatal(modulename, subroutinename, &
         'Problem with allocate(sp_energies)')

    allocate(Vpair(basis%size,basis%size),stat=mem_stat)
    if (mem_stat/=0) call fatal(modulename, subroutinename, &
         'Problem with allocate(Vpair)')

    allocate(psi_bar%cylHO(1:basis%size))
    psi_bar%cylHO=0.0

    do isospin = 1, 2

       shift=0

!
!      Initialization of expansion coefficients to 0
!      (array filled in by blocks ==> need to clean everything at every iteration)
!

       do iHF=1,basis%size
          sp_states(iHF,isospin)%coef%cylHO(:)=0.0
       end do
       sp_states(:,isospin)%pi=0

       loop_blocks: do k=1,Nblocks

          n=basis%block_size(k)

!
!         Hamiltonian matrix elements in the Harmonic Oscillator (HO) basis
!

          if (trim(calculation_type)=='exact') then
             call hamiltonian_cylHO(h_HF(k)%HO(:,:),k,isospin,basis,HFfield(isospin),iter)
          else if (trim(calculation_type)=='perturbative') then
             call hamiltonian_cylHO(h_HF(k)%HO(:,:),k,isospin,basis,HFfield_even(isospin),iter)
          else
             stop 'Wrong calculation type'
          end if

!
!         Hamiltonian diagonalization
!

          allocate(C(n,n),Ctmp(n))

          call diagon(h_HF(k)%HO(:,:), h_HF(k)%eigenval(:), C, n)

!
!         Fill in arrays with eigenvalues and eigenvectors
!         (largest expansion coeff. chosen to be positive)
!

          sp_energies(shift+1:shift+n,isospin)=h_HF(k)%eigenval(:)

          do i=1,n
             do j=1,n
                Ctmp(j)=abs(C(j,i))
             end do
             j=maxloc(Ctmp,1)
             sp_states(i+shift,isospin)%coef%cylHO(1+shift:n+shift)=C(:,i)* &
                  sign(1.0_rk,C(j,i))
             sp_states(i+shift,isospin)%energy=h_HF(k)%eigenval(i)
             sp_states(i+shift,isospin)%block_rank=i
          end do

          deallocate(C,Ctmp)

          if (symmetries%parity) sp_states(1+shift:n+shift,isospin)%pi= &
               basis%parity(1+shift)
          sp_states(1+shift:n+shift,isospin)%Omega=basis%Omega(1+shift)
          sp_states(1+shift:n+shift,isospin)%block=k

          shift=shift+n

       end do loop_blocks

!
!      Sort all eigenvalues in increasing order to define hole and particle states
!      and all eigenvectors accordingly
!

       call sort(sp_energies(:,isospin), basis%size, 'increasing', permut = permut(:,isospin))

       call permutation(sp_states(:,isospin), permut(:,isospin), &
            basis%size, 'normal')

       permut_inv(:,isospin) = 0
       do iHF = 1, basis%size
          sp_states(iHF,isospin)%HF_state=iHF
          permut_inv(permut(iHF,isospin),isospin) = iHF
       end do
!
!      Single-particle wave functions and their gradients at Gauss--Hermite--Laguerre
!      integration points
!

       do iHF = 1, basis%size
          call sp_wave_functions(sp_states(iHF,isospin), basis)
       end do

!
!      Occupation numbers: pairing or pure Hartree--Fock
!

       if (iter==1 .and. model%initial_potential=='WS' .and. model%pairing/='seniority') then

          sp_states(1:Npart(isospin),isospin)%occupation=1.0_rk
          sp_states(Npart(isospin)+1:basis%size,isospin)%occupation=0.0
          local_nucleus%pairing%pairing_energy=0.0

       else

!
!      Initializations:
!
!      * chemical potential to Fermi level (half sum of last occ. and 1st unocc. levels)
!      * pairing energy to 0
!      * pairing gaps and qp quantities to 0
!      * occupation probabilities to the Hatree-Fock values
!
          iHF=Npart(isospin)
          e_1st_unoccupied=sp_states(iHF+1,isospin)%energy
          if (iHF/=0) then
             fermi=0.5_rk*(sp_states(iHF,isospin)%energy+e_1st_unoccupied)
          else
             fermi = e_1st_unoccupied
          end if

          local_nucleus%pairing(isospin)%chemical_potential=fermi
          local_nucleus%pairing(isospin)%pairing_energy=0.0
          local_nucleus%pairing(isospin)%average_gap=0.0
          local_nucleus%pairing(isospin)%gap(:)=0.0
          local_nucleus%pairing(isospin)%min_qp_index=0
          local_nucleus%pairing(isospin)%min_qp_energy=0

!
!         Determination of quasi-pairs (which would be Kramers degenerate
!         in the equal-filling approximation) and their energy
!

          partner(:)=0
          sp_states(:,isospin)%pair_partner(isospin)=0
          quasi_pair_energy(:)=1e6_rk
          sp_label(:)=0
          n=0

          do iHF=1,basis%size
             if (sp_states(iHF,isospin)%Omega<0) cycle
!            definition of time-reversed state
             psi_bar=time_reversal(sp_states(iHF,isospin)%coef,basis)
             ovlp=1e-09_rk
             partner(iHF)=0
!            search for sp state with opposite Omega and largest overlap with psi_bar
             do jHF=1,basis%size
                if (sp_states(jHF,isospin)%Omega /= -sp_states(iHF,isospin)%Omega) cycle
                tmp=dot_product(sp_states(jHF,isospin)%coef%cylHO,psi_bar%cylHO)
                if (abs(tmp)>abs(ovlp)) then
                   ovlp=tmp
                   partner(iHF)=jHF
                end if
             end do
             jHF=partner(iHF)
             if (jHF /= 0) then
                n=n+1
                if (n>basis%size/2) call fatal(modulename,subroutinename,'Too many unfound partners')
                quasi_pair_energy(n)=sp_states(iHF,isospin)%energy+sp_states(jHF,isospin)%energy
                sp_label(n)=iHF
                sp_states(iHF,isospin)%pair_partner(isospin)=jHF
                sp_states(jHF,isospin)%pair_partner(isospin)=iHF
             else
                write(out_unit,'(a,i4)') 'iHF=',iHF
                call fatal(modulename,subroutinename,'Partner state not found')
             end if
          end do

          call sort(quasi_pair_energy(:), basis%size/2, 'increasing', permut=order(:))

          tmp_label(:)=0
          do i=1,basis%size/2
             tmp_label(i)=sp_label(order(i))
          end do
          sp_label(:)=tmp_label(:)

!
!         Occupation numbers in the Hartree--Fock approximation
!

!
!         Finding the closest rank of the blocked state to the Fermi level
!
          if (local_nucleus%fix_rank == 'n') then
             max_rank = maxval(basis%block_size(:),1)
             if (max_rank <= 0) call fatal(modulename,subroutinename,"max_rank <= 0")
             allocate(ph_energy(max_rank))
             ph_energy(:) = ph_energy_min
             if (symmetries%parity) then
                do j = 1, 3
                   if (local_nucleus%K(isospin,j) == 0) cycle
                   do m = 1, max_rank
                      do jHF = 1, basis%size
                         if (sp_states(jHF,isospin)%Omega /= local_nucleus%K(isospin,j) .or. &
                              sp_states(jHF,isospin)%pi /= local_nucleus%pi(isospin,j)) cycle
                         if (sp_states(jHF,isospin)%block_rank == m) ph_energy(m) = &
                              abs(sp_states(jHF,isospin)%energy-fermi)
                      end do
                   end do
                   if (minval(ph_energy) == ph_energy_min) call fatal(modulename,subroutinename, &
                        'Rank not found')
                   local_nucleus%rank(isospin,j) = minloc(ph_energy,dim=1)
                   if (debugmode /= 'low') write(out_unit,'(4x,3(a,I2))') 'Blocked K^pi state =', &
                        local_nucleus%K(isospin,j),'/2,  parity = ',local_nucleus%pi(isospin,j), &
                        ',  Rank = ', local_nucleus%rank(isospin,j)
                   do m = 1, max_rank
                      if (m == local_nucleus%rank(isospin,j)) cycle
                      if (abs(ph_energy(m)-ph_energy(local_nucleus%rank(isospin,j))) <= 0.5_rk) then
                         write(out_unit,'(4x,2(a,i2),a)') 'NOTE: The distances from Fermi level of ranks ', &
                              local_nucleus%rank(isospin,j),' and ',m,' differ by less than 500 keV'
                         write(out_unit,'(3x,a,i2,a,f12.6)') 'Distance from Fermi level for rank ', &
                              local_nucleus%rank(isospin,j), ' = ', ph_energy(local_nucleus%rank(isospin,j))
                         write(out_unit,'(3x,a,i2,a,f12.6)') 'Distance from Fermi level for rank ', &
                              m, ' = ', ph_energy(m)
                      end if
                   end do
                end do
             else
                do j = 1, 3
                   if (nucleus%K(isospin,j) == 0) cycle
                   do m = 1, max_rank
                      do jHF = 1, basis%size
                         if (sp_states(jHF,isospin)%Omega /= nucleus%K(isospin,j)) cycle
                         if (sp_states(jHF,isospin)%block_rank == m) ph_energy(m) = &
                              abs(sp_states(jHF,isospin)%energy-fermi)
                      end do
                   end do
                   if (minval(ph_energy) == ph_energy_min) call fatal(modulename,subroutinename, &
                        'Rank not found')
                   nucleus%rank(isospin,j) = minloc(ph_energy,dim=1)
                   if (debugmode /= 'low') write(out_unit,'(4x,3(a,I2))') 'Blocked K^pi state =', &
                        nucleus%K(isospin,j),'/2,  parity = ',nucleus%pi(isospin,j), &
                        ',  Rank = ', nucleus%rank(isospin,j)
                   do m = 1, max_rank
                      if (m == nucleus%rank(isospin,j)) cycle
                      if (abs(ph_energy(m)-ph_energy(nucleus%rank(isospin,j))) <= 0.5_rk) then
                         write(out_unit,'(4x,2(a,i2),a)') 'NOTE: The distances from Fermi level of ranks ', &
                              nucleus%rank(isospin,j),' and ',m,' differ by less than 500 keV'
                         write(out_unit,'(3x,a,i2,a,f12.6)') 'Distance from Fermi level for rank ', &
                              nucleus%rank(isospin,j), ' = ', ph_energy(nucleus%rank(isospin,j))
                         write(out_unit,'(3x,a,i2,a,f12.6)') 'Distance from Fermi level for rank ', &
                              m, ' = ', ph_energy(m)
                      end if
                   end do
                end do
             end if
             deallocate(ph_energy)
          end if

!         Determination of the blocked states

          do i=1,nqp
             i_odd_local(isospin,i)=0
             if (local_nucleus%K(isospin,i) == 0) cycle
             loop_iHF: do iHF=1,basis%size
                if (sp_states(iHF,isospin)%Omega==local_nucleus%K(isospin,i) .and. &
                     sp_states(iHF,isospin)%block_rank==local_nucleus%rank(isospin,i)) then
                   if (symmetries%parity) then
                      if (sp_states(iHF,isospin)%pi==local_nucleus%pi(isospin,i)) then
                         i_odd_local(isospin,i)=iHF
                         exit loop_iHF
                      end if
                   else
                      i_odd_local(isospin,i)=iHF
                      exit loop_iHF
                   end if
                end if
             end do loop_iHF
          end do
             
!         Initialization of occupation numbers and type component

          sp_states(:,isospin)%occupation = 0.0
          sp_states(:,isospin)%type = 0

!         Occupation numbers of the blocked states

          do i=1,nqp
             if (i_odd_local(isospin,i)/=0) then
                if (model%blocking == 'SCB') then
                   sp_states(i_odd_local(isospin,i),isospin)%occupation = 1.0_rk
                else if (model%blocking == 'EFA') then
                   sp_states(i_odd_local(isospin,i),isospin)%occupation = 0.5_rk
                   jHF = sp_states(i_odd_local(isospin,i),isospin)%pair_partner(isospin)
                   sp_states(jHF,isospin)%occupation = 0.5_rk
                else
                   call fatal(modulename,subroutinename,'Unrecognized blocking scheme ' &
                        //trim(model%blocking))
                end if
             else if (local_nucleus%K(isospin,i) /= 0) then
                call fatal(modulename,subroutinename,'Blocked state not found')
             end if
          end do

!         Occupation numbers of remaining, lowest nb_qp(isospin) quasi-pairs 
!         of states set to 1

          nb_qp(isospin) = (Npart(isospin) - nb_blocked_part(isospin))/2
          jHF = 1
          do k = 1, nb_qp(isospin)
             filled_qp = .false.
             loop_fill_qp: do j = jHF, basis%size/2
                iHF = sp_label(j)
                do i = 1, nqp
                   if (iHF == i_odd_local(isospin,i) .or. partner(iHF) == i_odd_local(isospin,i)) cycle loop_fill_qp
                end do
                sp_states(iHF,isospin)%occupation = 1.0_rk
                sp_states(partner(iHF),isospin)%occupation = 1.0_rk
                filled_qp = .true.
                jHF=j+1
                exit loop_fill_qp
             end do loop_fill_qp
             if (.not.filled_qp) then
                write(out_unit,'(a,i2,2(1x,a,i4),2x,a,3i5)') &
                     'isospin=',isospin,'k=',k,'jHF=',jHF,'i_odd_local(isospin,:)=',i_odd_local(isospin,:)
                call fatal(modulename,subroutinename,"Quasi-pair not filled")
             end if
          end do

!         Define HF_level, HF_state and type components of sp_states

          i=0
          do k=1,basis%size
             sp_states(k,isospin)%HF_level = k
             sp_states(k,isospin)%HF_state = k
             if (sp_states(k,isospin)%occupation >= 0.5_rk) then
                sp_states(k,isospin)%type = 1
                i=i+1
             else
                sp_states(k,isospin)%type = 0
             end if
          end do

          if (model%blocking == 'EFA') then
             do k = 1, nqp
                if (i_odd_local(isospin,k) == 0) cycle
                jHF = sp_states(i_odd_local(isospin,k),isospin)%pair_partner(isospin)
                sp_states(jHF,isospin)%type = 0
                i = i - 1
             end do
          end if

          if (i /= Npart(isospin)) then
             do iHF = 1, basis%size
                write(out_unit,'(a1,2x,2(i4,1x),a2,a1,3(i4,1x))') &
                     nucleon_type(isospin),iHF,sp_states(iHF,isospin)%Omega,'/2', &
                     signparity(sp_states(iHF,isospin)%pi),iHF,partner(iHF), &
                     sp_states(iHF,isospin)%type
             end do
             write(out_unit,'(/a,i2,1x,a,i4,a,i4)') 'isospin=',isospin,&
                  'number of hole states =', i, ' whereas Npart=', Npart(isospin)
             call fatal(modulename,subroutinename, &
                  'Problem with the particle or hole character of s.p. states')
          end if

!
!         Matrix elements of the effective interaction in the pairing channel
!         (replacing the Skyrme interaction)
!

          call pairing_matrix_elements(sp_states(:,isospin), densities, isospin, &
               local_nucleus, model, basis, Vpair, partner, iter)

!
!         Occupation numbers in the BCS approximation
!
!         Check if BCS calculation is needed:
!         * if max{|Vpair(i,j)|,i,j=1..basis%size} is larger than 1 keV, call BCS
!         * if not, do nothing ==> HF solution
!

          Vpair_max = maxval(abs(Vpair))
          if (abs(Vpair_max)>1e-3_rk) then
             call BCS(local_nucleus, model, basis, sp_states(:,isospin), Vpair, isospin, &
                  i_odd_local(isospin,:))
          end if

       end if

    end do

!
!   Checking number of blocked s.p. states
!

    count = 0
    do isospin = 1, 2
       do i = 1, nqp
          if (i_odd_local(isospin,i) /= 0) count = count + 1
       end do
    end do
    if (count > 1) call warning(modulename, subroutinename, 'More than 1 blocked state')

!
!   Calculation of matrix elements and 1st-order eigenstates
!

    s_odd%tot = 0.0
    s_odd%S = 0.0
    s_odd%A = 0.0
    s_odd%C = 0.0
    s_odd%D = 0.0

    expval_sz_phi_i = 0

    do isospin = 1, 2

       write(iunit,'(/a)') 'Particle type: ' // particle_type(isospin)

       expval_sz(isospin)%tot = 0.0
       expval_sz(isospin)%K = 0.0
       expval_sz(isospin)%U = 0.0
       expval_sz(isospin)%so = 0.0
       expval_sz(isospin)%W = 0.0
       expval_sz(isospin)%S = 0.0
       expval_sz(isospin)%A = 0.0
       expval_sz(isospin)%C = 0.0
       expval_sz(isospin)%D = 0.0
  
       do iHF = 1, basis%size

          if (sp_states(iHF,isospin)%occupation < 1e-3_rk) cycle
          ei = sp_states(iHF,isospin)%energy

          call sz_matrix_element(sp_states(iHF,isospin), &
               sp_states(iHF,isospin),basis,expval_sz_unperturbed)
          expval_sz(isospin)%tot = expval_sz(isospin)%tot + sp_states(iHF,isospin)%occupation * &
               expval_sz_unperturbed
          sp_states(iHF,isospin)%avg_sz = expval_sz_unperturbed

          iblock=sp_states(iHF,isospin)%block
          shift=sum(basis%block_size(1:iblock),1)-basis%block_size(iblock)

          do j = 1, basis%block_size(iblock)

             jHF = permut_inv(j + shift,isospin)
             ej = sp_states(jHF,isospin)%energy

             call sz_matrix_element(sp_states(iHF,isospin),sp_states(jHF,isospin),basis,sum_sz)
             if (abs(sum_sz) < 0.01_rk) cycle

             call h_HF_matrix_element(sp_states(iHF,isospin), &
                  sp_states(jHF,isospin),isospin,basis,model,HFfield(isospin), &
                  densities,h_HF_mtx_elem)

             if (jHF /= iHF) then

                hodd_sz = sum_sz*(h_HF_mtx_elem%S(0)+h_HF_mtx_elem%A(0)+h_HF_mtx_elem%C(0)+ &
                     h_HF_mtx_elem%D(0))/(ei-ej)
                expval_sz(isospin)%tot = expval_sz(isospin)%tot + sp_states(iHF,isospin)%occupation * &
                     2*hodd_sz
                do i = 1, 21
                   expval_sz(isospin)%S(i) = expval_sz(isospin)%S(i) + sp_states(iHF,isospin)%occupation * &
                        2*sum_sz*h_HF_mtx_elem%S(i)/(ei-ej)
                   expval_sz(isospin)%A(i) = expval_sz(isospin)%A(i) + sp_states(iHF,isospin)%occupation * &
                        2*sum_sz*h_HF_mtx_elem%A(i)/(ei-ej)
                   expval_sz(isospin)%C(i) = expval_sz(isospin)%C(i) + sp_states(iHF,isospin)%occupation * &
                        2*sum_sz*h_HF_mtx_elem%C(i)/(ei-ej)
                   expval_sz(isospin)%D(i) = expval_sz(isospin)%D(i) + sp_states(iHF,isospin)%occupation * &
                        2*sum_sz*h_HF_mtx_elem%D(i)/(ei-ej)
                end do

                write(iunit,'(/2(a,i4,1x,a,f11.6,2x,i3,a2,a1,2x),3(a,f11.6,1x))') &
                     'iHF=',iHF,' e =', ei, &
                     sp_states(iHF,isospin)%Omega,'/2', &
                     signparity(sp_states(iHF,isospin)%pi), &
                     'jHF=',jHF,' e =', sp_states(jHF,isospin)%energy, &
                     sp_states(jHF,isospin)%Omega,'/2', &
                     signparity(sp_states(jHF,isospin)%pi), &
                     '<iHF|h_HF|jHF> =',h_HF_mtx_elem%tot,'<iHF|s_z|jHF>  =',sum_sz, &
                     '<i|h_odd|j><i|s_z|j>/(e_i-e_j) =',hodd_sz

             else

                write(iunit,'(/2(a,i4,1x,a,f11.6,2x,i3,a2,a1,2x),2(a,f11.6,1x))') &
                     'iHF=',iHF,' e =', ei, &
                     sp_states(iHF,isospin)%Omega,'/2', &
                     signparity(sp_states(iHF,isospin)%pi), &
                     'jHF=',jHF,' e =', sp_states(jHF,isospin)%energy, &
                     sp_states(jHF,isospin)%Omega,'/2', &
                     signparity(sp_states(jHF,isospin)%pi), &
                     '<iHF|h_HF|jHF> =',h_HF_mtx_elem%tot,'<iHF|s_z|jHF>  =',sum_sz

             end if

             write(iunit,'(79x,4(a,f11.6,1x)/,79x,4(a,f11.6,1x))') &
                  'h_K =',h_HF_mtx_elem%K(0),'h_U =',h_HF_mtx_elem%U(0), &
                  'h_so =', h_HF_mtx_elem%so(0),'h_W =',h_HF_mtx_elem%W(0), &
                  'h_S =',h_HF_mtx_elem%S(0),'h_A =',h_HF_mtx_elem%A(0), &
                  'h_C  =',h_HF_mtx_elem%C(0),'h_D =',h_HF_mtx_elem%D(0)
             write(iunit,'(79x,3(a,f11.6,1x)/,79x,4(a,f11.6,1x))') &
                  'h_S(B9)=',h_HF_mtx_elem%S(9),'h_S(B10)=',h_HF_mtx_elem%S(10), &
                  'h_S(B11)=',h_HF_mtx_elem%S(11),'h_S(B12)=',h_HF_mtx_elem%S(12), &
                  'h_S(B13)=',h_HF_mtx_elem%S(13),'h_S(B14)=',h_HF_mtx_elem%S(14), &
                  'h_S(B15)=',h_HF_mtx_elem%S(15)
             write(iunit,'(79x,3(a,f11.6,1x))') &
                  'h_A(B3)=',h_HF_mtx_elem%A(3),'h_A(B4)=',h_HF_mtx_elem%A(4), &
                  'h_A(B9)=',h_HF_mtx_elem%A(9)

          end do

       end do

       expval_sz(isospin)%S(0) = sum(expval_sz(isospin)%S(1:21),1)
       expval_sz(isospin)%A(0) = sum(expval_sz(isospin)%A(1:21),1)
       expval_sz(isospin)%C(0) = sum(expval_sz(isospin)%C(1:21),1)
       expval_sz(isospin)%D(0) = sum(expval_sz(isospin)%D(1:21),1)
       expval_sz_core(isospin) = expval_sz(isospin)%tot

       if (mod(Npart(isospin),2) == 1) then

          gs_q = gs(isospin)

          iHF = i_odd_local(isospin,1)
          iblock=sp_states(iHF,isospin)%block
          shift=sum(basis%block_size(1:iblock))-basis%block_size(iblock)
          expval_sz_phi_i = sp_states(iHF,isospin)%avg_sz
          s_odd%tot = sp_states(iHF,isospin)%avg_sz
          do j = 1, basis%block_size(iblock)
             jHF = permut_inv(j + shift,isospin)
             if (jHF == iHF) cycle
             call h_HF_matrix_element(sp_states(iHF,isospin), &
                  sp_states(jHF,isospin),isospin,basis,model,HFfield(isospin), &
                  densities,h_HF_mtx_elem)
             call sz_matrix_element(sp_states(iHF,isospin),sp_states(jHF,isospin),basis,sum_sz)
             ej = sp_states(jHF,isospin)%energy
             hodd_sz = sum_sz*(h_HF_mtx_elem%S(0)+h_HF_mtx_elem%A(0)+h_HF_mtx_elem%C(0)+ &
                  h_HF_mtx_elem%D(0))/(ei-ej)
             s_odd%tot = s_odd%tot + sp_states(iHF,isospin)%occupation * 2*hodd_sz  ! s_odd = <psi_i|sz|psi_i>
             s_odd%S(:) = s_odd%S(:) + sp_states(iHF,isospin)%occupation * &
                  2 * sum_sz*h_HF_mtx_elem%S(:)/(ei-ej)
             s_odd%A(:) = s_odd%A(:) + sp_states(iHF,isospin)%occupation * &
                  2 * sum_sz*h_HF_mtx_elem%A(:)/(ei-ej)
             s_odd%C(:) = s_odd%C(:) + sp_states(iHF,isospin)%occupation * &
                  2 * sum_sz*h_HF_mtx_elem%C(:)/(ei-ej)
             s_odd%D(:) = s_odd%D(:) + sp_states(iHF,isospin)%occupation * &
                  2 * sum_sz*h_HF_mtx_elem%D(:)/(ei-ej)
          end do

! s_odd = <psi|sz|psi_i> = <phi|sz|phi_i>  + order-1 term in h_odd

          expval_sz_core(isospin) = expval_sz_core(isospin) - s_odd%tot
          expval_sz(isospin)%S(:) = expval_sz(isospin)%S(:) - s_odd%S(:)
          expval_sz(isospin)%A(:) = expval_sz(isospin)%A(:) - s_odd%A(:)
          expval_sz(isospin)%C(:) = expval_sz(isospin)%C(:) - s_odd%C(:)
          expval_sz(isospin)%D(:) = expval_sz(isospin)%D(:) - s_odd%D(:)

          write(iunit,'(6(a,f10.6,1x))') &
            '<sz>' // nucleon_type(isospin) // '=', expval_sz(isospin)%tot, &
            's_odd =',s_odd%tot, '<phi_i|s_z|phi_i> =',sp_states(iHF,isospin)%avg_sz, &
            '<sz>' // nucleon_type(isospin) // '(core,S)=', expval_sz(isospin)%S(0), &
            '<sz>' // nucleon_type(isospin) // '(core,A)=', expval_sz(isospin)%A(0), &
            '<sz>' // nucleon_type(isospin) // '(core,C)=', expval_sz(isospin)%C(0)

       else

          write(iunit,'(4(a,f10.6,1x))') &
            '<sz>' // nucleon_type(isospin) // '=', expval_sz(isospin)%tot, &
            '<sz>' // nucleon_type(isospin) // '(core,S)=', expval_sz(isospin)%S(0), &
            '<sz>' // nucleon_type(isospin) // '(core,A)=', expval_sz(isospin)%A(0), &
            '<sz>' // nucleon_type(isospin) // '(core,C)=', expval_sz(isospin)%C(0)

       end if

    end do

!
!   Perturbative gs quenching (with respect to h_even eigenstates): delta at order 1 in h_odd
!

    if (abs(expval_sz_phi_i) < 1e-6_rk) call fatal(modulename,subroutinename, &
         's_odd too small ==> division by 0')
    do isospin = 1,2
       delta(isospin)%tot = (gl(isospin) - gs(isospin))/gs_q * expval_sz_core(isospin)/expval_sz_phi_i
       delta(isospin)%S(:) = (gl(isospin) - gs(isospin))/gs_q * expval_sz(isospin)%S(:)/expval_sz_phi_i
       delta(isospin)%A(:) = (gl(isospin) - gs(isospin))/gs_q * expval_sz(isospin)%A(:)/expval_sz_phi_i
       delta(isospin)%C(:) = (gl(isospin) - gs(isospin))/gs_q * expval_sz(isospin)%C(:)/expval_sz_phi_i
       delta(isospin)%D(:) = (gl(isospin) - gs(isospin))/gs_q * expval_sz(isospin)%D(:)/expval_sz_phi_i
    end do

    quenching_factor%tot = 1.0_rk - (delta(1)%tot+delta(2)%tot)

    write(iunit,'(/2(a,f10.6,1x)/)') &
         'Perturbative gs quenching factor: eta =',quenching_factor%tot, &
         '==> 1 - eta = delta =',delta(1)%tot+delta(2)%tot

    write(iunit,'(4(a,f10.6,1x)//,20x,5(a,f10.6,1x)/,20x,4(a,f10.6,1x)/)') &
         'delta(n)=',delta(1)%tot,'delta(n,S)=',delta(1)%S(0),'delta(n,A)=',delta(1)%A(0), &
         'delta(n,C)=',delta(1)%C(0),'delta(n,S,B9) =',delta(1)%S(9), &
         'delta(n,S,B10)=',delta(1)%S(10),'delta(n,S,B11)=',delta(1)%S(11), &
         'delta(n,S,B12)=',delta(1)%S(12),'delta(n,S,B13)=',delta(1)%S(13), &
         'delta(n,S,B14)=',delta(1)%S(14),'delta(n,S,B15)=',delta(1)%S(15), &
         'delta(n,S,B18)=',delta(1)%S(18),'delta(n,S,B19)=',delta(1)%S(19)
    write(iunit,'(20x,3(a,f10.6,1x)/)') &
         'delta(n,A,B3) =',delta(1)%A(3),'delta(n,A,B4) =',delta(1)%A(4), &
         'delta(n,A,B9) =',delta(1)%A(9)
    write(iunit,'(20x,2(a,f10.6,1x)/)') &
         'delta(n,C,B14) =',delta(1)%C(14),'delta(n,C,B15) =',delta(1)%C(15)

    write(iunit,'(4(a,f10.6,1x)//,20x,5(a,f10.6,1x)/,20x,4(a,f10.6,1x)/)') &
         'delta(p)=',delta(2)%tot,'delta(p,S)=',delta(2)%S(0),'delta(p,A)=',delta(2)%A(0), &
         'delta(p,C)=',delta(2)%C(0),'delta(p,S,B9) =',delta(2)%S(9), &
         'delta(p,S,B10)=',delta(2)%S(10),'delta(p,S,B11)=',delta(2)%S(11), &
         'delta(p,S,B12)=',delta(2)%S(12),'delta(p,S,B13)=',delta(2)%S(13), &
         'delta(p,S,B14)=',delta(2)%S(14),'delta(p,S,B15)=',delta(2)%S(15), &
         'delta(p,S,B18)=',delta(2)%S(18),'delta(p,S,B19)=',delta(2)%S(19)
    write(iunit,'(20x,3(a,f10.6,1x)/)') &
         'delta(p,A,B3) =',delta(2)%A(3),'delta(p,A,B4) =',delta(2)%A(4), &
         'delta(p,A,B9) =',delta(2)%A(9)
    write(iunit,'(20x,2(a,f10.6,1x)/)') &
         'delta(p,C,B14) =',delta(2)%C(14),'delta(p,C,B15) =',delta(2)%C(15)
    write(iunit,'(4(a,f10.6,1x))') &
         'delta=',delta(1)%tot+delta(2)%tot,'delta(S)=',delta(1)%S(0)+delta(2)%S(0), &
         'delta(A)=',delta(1)%A(0)+delta(2)%A(0), &
         'delta(C)=',delta(1)%C(0)+delta(2)%C(0)

    B3 = model%parameters%B3
    B4 = model%parameters%B4
    B5 = model%parameters%B5
    B6 = model%parameters%B6
    B10 = model%parameters%B10
    B11 = model%parameters%B11
    B12 = model%parameters%B12
    B13 = model%parameters%B13
    B14 = model%parameters%B14
    B15 = model%parameters%B15
    B18 = model%parameters%B18
    B19 = model%parameters%B19
    write(iunit,'(/5(2(a,f13.6,1x)/))') &
         'B3 =',B3,'B4 =',B4,&
         'B10=',B10,'B11=',B11, &
         'B12=',B12,'B13=',B13, &
         'B14=',B14,'B15=',B15, &
         'B18=',B18,'B19=',B19
    write(iunit,'(/a)') &
         'Decomposition of S, A, and C contributions in terms of coupling constants:'
    write(iunit,'(/5(a,f10.6,1x)/,4(a,f10.6,1x))') &
         'delta(S,B9)=',delta(1)%S(9)+delta(2)%S(9), &
         'delta(S,B10)=',delta(1)%S(10)+delta(2)%S(10),'delta(S,B11)=',delta(1)%S(11)+delta(2)%S(11), &
         'delta(S,B12)=',delta(1)%S(12)+delta(2)%S(12),'delta(S,B13)=',delta(1)%S(13)+delta(2)%S(13), &
         'delta(S,B14)=',delta(1)%S(14)+delta(2)%S(14),'delta(S,B15)=',delta(1)%S(15)+delta(2)%S(15), &
         'delta(S,B18)=',delta(1)%S(18)+delta(2)%S(18),'delta(S,B19)=',delta(1)%S(19)+delta(2)%S(19)
    write(iunit,'(a,f10.6,1x)') 'delta(S,B9)+...+delta(S,B19)=',delta(1)%S(9)+delta(2)%S(9)+ &
         delta(1)%S(10)+delta(2)%S(10) + delta(1)%S(11)+delta(2)%S(11) + &
         delta(1)%S(12)+delta(2)%S(12) + delta(1)%S(13)+delta(2)%S(13) + &
         delta(1)%S(14)+delta(2)%S(14) + delta(1)%S(15)+delta(2)%S(15) + &
         delta(1)%S(18)+delta(2)%S(18) + delta(1)%S(19)+delta(2)%S(19)
    write(iunit,'(/3(a,f10.6,1x))') &
         'delta(A,B3)=',delta(1)%A(3)+delta(2)%A(3),'delta(A,B4)=',delta(1)%A(4)+delta(2)%A(4), &
         'delta(A,B9)=',delta(1)%A(9)+delta(2)%A(9)
    write(iunit,'(a,f10.6,1x)') 'delta(A,B3)+delta(A,B4)+delta(A,B9)=',delta(1)%A(3)+delta(2)%A(3)+ &
         delta(1)%A(4)+delta(2)%A(4)+ delta(1)%A(9)+delta(2)%A(9)
    write(iunit,'(/2(a,f10.6,1x))') &
         'delta(C,B14)=',delta(1)%C(14)+delta(2)%C(14),'delta(C,B15)=',delta(1)%C(15)+delta(2)%C(15)
    write(iunit,'(a,f10.6,1x)') 'delta(C,B14)+delta(C,B15)=',delta(1)%C(14)+delta(2)%C(14)+ &
         delta(1)%C(15)+delta(2)%C(15)

    write(iunit,'(/a)') &
         'Decomposition of S, A and C contributions in terms of Skyrme parameters:'
    delta_S_t1 = 0.0_rk
    delta_S_t1x1 = 0.0_rk
    delta_S_t2 = 0.0_rk
    delta_S_t2x2 = 0.0_rk
    if (abs(B14) > 1e-6_rk) then
       delta_S_t1x1 = delta_S_t1x1 + (0.25_rk + B18/B14) * (delta(1)%S(14)+delta(2)%S(14))
       delta_S_t2x2 = delta_S_t2x2 + (0.75_rk - B18/B14) * (delta(1)%S(14)+delta(2)%S(14))
    end if
    if (abs(B15) > 1e-6_rk) then
       delta_S_t1 = delta_S_t1 + (0.25_rk + B19/B15) * (delta(1)%S(15)+delta(2)%S(15))
       delta_S_t2 = delta_S_t2 + (0.75_rk - B19/B15) * (delta(1)%S(15)+delta(2)%S(15))
    end if
    if (abs(B18) > 1e-6_rk) then
       delta_S_t1x1 = delta_S_t1x1 + (0.75_rk + 3.0_rk*B14/(16.0_rk*B18)) * &
            (delta(1)%S(18)+delta(2)%S(18))
       delta_S_t2x2 = delta_S_t2x2 + (0.25_rk - 3.0_rk*B14/(16.0_rk*B18)) * &
            (delta(1)%S(18)+delta(2)%S(18))
    end if
    if (abs(B19) > 1e-6_rk) then
       delta_S_t1 = delta_S_t1 + (0.75_rk + 3.0_rk*B15/(16.0_rk*B19)) * &
            (delta(1)%S(19)+delta(2)%S(19))
       delta_S_t2 = delta_S_t2 + (0.25_rk - 3.0_rk*B15/(16.0_rk*B19)) * &
            (delta(1)%S(19)+delta(2)%S(19))
    end if
    delta_S_t0=delta(1)%S(11)+delta(2)%S(11)
    delta_S_t0x0=delta(1)%S(10)+delta(2)%S(10)
    delta_S_t3=delta(1)%S(13)+delta(2)%S(13)
    delta_S_t3x3=delta(1)%S(12)+delta(2)%S(12)
    delta_S_W=delta(1)%S(9)+delta(2)%S(9)
    write(iunit,'(/a,f10.6/,4(2(a,f10.6,1x)/),2(a,f10.6,1x))') &
         'delta(S)=',delta(1)%S(0)+delta(2)%S(0), &
         'delta(S,t0)=',delta_S_t0,'delta(S,t0x0)=',delta_S_t0x0, &
         'delta(S,t1)=',delta_S_t1,'delta(S,t1x1)=',delta_S_t1x1, &
         'delta(S,t2)=',delta_S_t2,'delta(S,t2x2)=',delta_S_t2x2, &
         'delta(S,t3)=',delta_S_t3,'delta(S,t3x3)=',delta_S_t3x3, &
         'delta(S,W) =',delta_S_W
    write(iunit,'(2(a,f10.6,1x))') 'delta(S,t1)+delta(S,t2)=',delta_S_t1+delta_S_t2, &
         'delta(S,B15)+delta(S,B19)=',delta(1)%S(15)+delta(2)%S(15)+delta(1)%S(19)+delta(2)%S(19)
    write(iunit,'(2(a,f10.6,1x))') 'delta(S,t1x1)+delta(S,t2x2)=',delta_S_t1x1+delta_S_t2x2, &
         'delta(S,B14)+delta(S,B18)=',delta(1)%S(14)+delta(2)%S(14)+delta(1)%S(18)+delta(2)%S(18)
    write(iunit,'(2(a,f10.6,1x))') 'delta(S,t0)+...+delta(S,W)=', &
         delta_S_t0+delta_S_t0x0+delta_S_t1+delta_S_t1x1+delta_S_t2+delta_S_t2x2+delta_S_t3+delta_S_t3x3+delta_S_W,&
         'delta(S)=',delta(1)%S(0)+delta(2)%S(0)

    delta_A_t1 = 0.0_rk
    delta_A_t1x1 = 0.0_rk
    delta_A_t2 = 0.0_rk
    delta_A_t2x2 = 0.0_rk
    if (abs(B3) > 1e-6_rk) then
       delta_A_t1 = delta_A_t1 + 0.25_rk * (23*B3+4*B4-76*B5-48*B6)/15.0_rk * &
            (delta(1)%A(3)+delta(2)%A(3))/B3
       delta_A_t1x1 = delta_A_t1x1 + 0.125_rk * (-16*B3-8*B4+32*B5+96*B6)/15.0_rk * &
            (delta(1)%A(3)+delta(2)%A(3))/B3
       delta_A_t2 = delta_A_t2 + 0.25_rk * (21*B3-12*B4+28*B5-16*B6)/5.0_rk * &
            (delta(1)%A(3)+delta(2)%A(3))/B3
       delta_A_t2x2 = delta_A_t2x2 + 0.125_rk * (-12*B3+24*B4-16*B5+32*B6)/5.0_rk * &
            (delta(1)%A(3)+delta(2)%A(3))/B3
    end if
    if (abs(B4) > 1e-6_rk) then
       delta_A_t1 = delta_A_t1 - 0.125_rk * (23*B3+4*B4-76*B5-48*B6)/15.0_rk * &
            (delta(1)%A(4)+delta(2)%A(4))/B4
       delta_A_t1x1 = delta_A_t1x1 - 0.25_rk * (-16*B3-8*B4+32*B5+96*B6)/15.0_rk * &
            (delta(1)%A(4)+delta(2)%A(4))/B4
       delta_A_t2 = delta_A_t2 + 0.125_rk * (21*B3-12*B4+28*B5-16*B6)/5.0_rk * &
            (delta(1)%A(4)+delta(2)%A(4))/B4
       delta_A_t2x2 = delta_A_t2x2 + 0.25_rk * (-12*B3+24*B4-16*B5+32*B6)/5.0_rk * &
            (delta(1)%A(4)+delta(2)%A(4))/B4
    end if
    delta_A_W=delta(1)%A(9)+delta(2)%A(9)
    write(iunit,'(/a,f10.6/,2(2(a,f10.6,1x)/),2(a,f10.6,1x))') &
         'delta(A)=',delta(1)%A(0)+delta(2)%A(0), &
         'delta(A,t1)=',delta_A_t1,'delta(A,t1x1)=',delta_A_t1x1, &
         'delta(A,t2)=',delta_A_t2,'delta(A,t2x2)=',delta_A_t2x2, &
         'delta(A,W) =',delta_A_W
    write(iunit,'(2(a,f10.6,1x))') &
         'delta(A,t1)+delta(A,t1x1)+delta(A,t2)+delta(A,t2x2)=',delta_A_t1+delta_A_t1x1+delta_A_t2+delta_A_t2x2, &
         'delta(A,B3)+delta(A,B4)=',delta(1)%A(3)+delta(2)%A(3)+delta(1)%A(4)+delta(2)%A(4)

    delta_C_t1 = 0.0_rk
    delta_C_t1x1 = 0.0_rk
    delta_C_t2 = 0.0_rk
    delta_C_t2x2 = 0.0_rk
    if (abs(B14) > 1e-6_rk) then
       delta_C_t1x1 = delta_C_t1x1 + (0.25_rk + B18/B14) * (delta(1)%C(14)+delta(2)%C(14))
       delta_C_t2x2 = delta_C_t2x2 + (0.75_rk - B18/B14) * (delta(1)%C(14)+delta(2)%C(14))
    end if
    if (abs(B15) > 1e-6_rk) then
       delta_C_t1 = delta_C_t1 + (0.25_rk + B19/B15) * (delta(1)%C(15)+delta(2)%C(15))
       delta_C_t2 = delta_C_t2 + (0.75_rk - B19/B15) * (delta(1)%C(15)+delta(2)%C(15))
    end if
    write(iunit,'(/a,f10.6/,2(2(a,f10.6,1x)/))') 'delta(C)=',delta(1)%C(0)+delta(2)%C(0), &
         'delta(C,t1)=',delta_C_t1,'delta(C,t1x1)=',delta_C_t1x1, &
         'delta(C,t2)=',delta_C_t2,'delta(C,t2x2)=',delta_C_t2x2

    write(iunit,'(a)') &
         'Decomposition of delta in terms of Skyrme parameters:'
    delta_t0 = delta_S_t0
    delta_t0x0 = delta_S_t0x0
    delta_t1 = delta_S_t1+delta_A_t1+delta_C_t1
    delta_t1x1 = delta_S_t1x1+delta_A_t1x1+delta_C_t1x1
    delta_t2 = delta_S_t2+delta_A_t2+delta_C_t2
    delta_t2x2 = delta_S_t2x2+delta_A_t2x2+delta_C_t2x2
    delta_t3 = delta_S_t3
    delta_t3x3 = delta_S_t3x3
    delta_W = delta_S_W+delta_A_W
    write(iunit,'(/a,f10.6/, 4(2(a,f10.6,1x)/),a,f10.6)') &
         'delta=',delta(1)%tot+delta(2)%tot, &
         'delta(t0)=',delta_t0,'delta(t0x0)=',delta_t0x0, &
         'delta(t1)=',delta_t1,'delta(t1x1)=',delta_t1x1, &
         'delta(t2)=',delta_t2,'delta(t2x2)=',delta_t2x2, &
         'delta(t3)=',delta_t3,'delta(t3x3)=',delta_t3x3, &
         'delta(W) =',delta_W
    write(iunit,'(a,f10.6)') 'delta(t0)+...+delta(W)=', &
         delta_t0+delta_t0x0+delta_t1+delta_t1x1+ &
         delta_t2+delta_t2x2+delta_t3+delta_t3x3+delta_W

    write(iunit,'(/a)') &
         'Decomposition of delta in terms of partial waves:'
    if (abs(B11) > 1e-6_rk) then
       x0 = -B10/B11
       delta_1S0 = 0.5_rk * (1.0_rk-x0) * (delta_t0 - delta_t0x0/x0)
       delta_3S1 = 0.5_rk * (1.0_rk+x0) * (delta_t0 + delta_t0x0/x0)
    else
       delta_1S0 = 0.5_rk * delta_t0
       delta_3S1 = 0.5_rk * delta_t0
    end if
    if (abs(B13) > 1e-6_rk .and. abs(B12) > 1e-6_rk) then
       x3 = -B12/B13
       delta_1S0 = delta_1S0 + 0.5_rk * (1.0_rk-x3) * (delta_t3 - delta_t3x3/x3)
       delta_3S1 = delta_3S1 + 0.5_rk * (1.0_rk+x3) * (delta_t3 + delta_t3x3/x3)
    else
       delta_1S0 = delta_1S0 + 0.5_rk * delta_t3
       delta_3S1 = delta_3S1 + 0.5_rk * delta_t3
    end if
    if (abs(B15+4*B19) > 1e-6_rk .and. abs(B14+4*B18) > 1e-6_rk) then
       x1 = -(B14+4*B18)/(B15+4*B19)
       delta_1S0 = delta_1S0 + 0.5_rk * (1.0_rk-x1) * (delta_t1 - delta_t1x1/x1)
       delta_3S1 = delta_3S1 + 0.5_rk * (1.0_rk+x1) * (delta_t1 + delta_t1x1/x1)
    else
       delta_1S0 = delta_1S0 + 0.5_rk * delta_t1
       delta_3S1 = delta_3S1 + 0.5_rk * delta_t1
    end if
    if (abs(-3*B15+4*B19) > 1e-6_rk .and. abs(-3*B14+4*B18) > 1e-6_rk) then
       x2 = (-3*B14+4*B18)/(-3*B15+4*B19)
       delta_1P1 = 0.5_rk * (1.0_rk-x2) * (delta_t2 - delta_t2x2/x2)
       delta_3P0 = 0.5_rk * (1.0_rk+x2) * (delta_t2 + delta_t2x2/x2)
    else
       delta_1P1 = 0.5_rk * delta_t2
       delta_3P0 = 0.5_rk * delta_t2
    end if
    delta_3P0 = delta_3P0/3.0_rk
    delta_3P1 = delta_3P0
    delta_3P2 = delta_3P0
    delta_3P0 = delta_3P0 + delta_W
    delta_3P1 = delta_3P1 - 0.5_rk*delta_W
    delta_3P2 = delta_3P2 + 0.5_rk*delta_W
    write(iunit,'(/a,f10.6/, 2(a,f10.6,1x)/, a,f10.6/, 3(a,f10.6,1x))') &
         'delta=',delta(1)%tot+delta(2)%tot, &
         'delta(1S0)=',delta_1S0,'delta(3S1)=',delta_3S1, &
         'delta(1P1)=',delta_1P1, &
         'delta(3P0)=',delta_3P0,'delta(3P1)=',delta_3P1,'delta(3P2)=',delta_3P2
    write(iunit,'(2(a,f10.6,1x))') 'delta(1S0)+delta(3S1)=',delta_1S0+delta_3S1, &
         'delta(t0)+delta(t0x0)+delta(t1)+delta(t1x1)+delta(t3)+delta(t3x3)=', &
         delta_t0+delta_t0x0+delta_t1+delta_t1x1+delta_t3+delta_t3x3
    write(iunit,'(2(a,f10.6,1x))') 'delta(1P1)+delta(3P0)+...+delta(3P2)=', &
         delta_1P1+delta_3P0+delta_3P1+delta_3P2, &
         'delta(t2)+delta(t2x2)+delta(W)=',delta_t2+delta_t2x2+delta_W

!
!   Memory deallocation
!

    do k=1,Nblocks
       deallocate(h_HF(k)%HO,h_HF(k)%eigenval)
    end do

    deallocate(h_HF,sp_energies,Vpair)

    deallocate(psi_bar%cylHO)

    do isospin = 1, 2
       deallocate(HFfield_even(isospin)%central)
       deallocate(HFfield_even(isospin)%effective_mass)
       deallocate(HFfield_even(isospin)%spin_orbit)
       deallocate(HFfield_even(isospin)%Coulomb)
       do i = 1, 3
          deallocate(HFfield_even(isospin)%tensor_so(i)%z)
          deallocate(HFfield_even(isospin)%tensor_so(i)%rho)
          deallocate(HFfield_even(isospin)%tensor_so(i)%phi)
       end do
       deallocate(HFfield_even(isospin)%S%z)
       deallocate(HFfield_even(isospin)%S%rho)
       deallocate(HFfield_even(isospin)%S%phi)
       deallocate(HFfield_even(isospin)%alpha%z)
       deallocate(HFfield_even(isospin)%alpha%rho)
       deallocate(HFfield_even(isospin)%alpha%phi)
       deallocate(HFfield_even(isospin)%C%z)
       deallocate(HFfield_even(isospin)%C%rho)
       deallocate(HFfield_even(isospin)%C%phi)
       deallocate(HFfield_even(isospin)%D%z)
       deallocate(HFfield_even(isospin)%D%rho)
       deallocate(HFfield_even(isospin)%D%phi)
       deallocate(HFfield_even(isospin)%div_s)
    end do

  end subroutine h_odd_perturbative


  !=====================================================================
  !                       SUBROUTINE h_tensor_perturbative
  !---------------------------------------------------------------------
  ! Calculates the first-order eigenstates and eigenvalues of 
  !      h_HF = h_HF(te=to=0) + h_HF(te,to)
  ! where h_HF(te,to) is treated as a perturbation with respect to 
  ! h_HF(te=to=0).
  !=====================================================================

  subroutine h_tensor_perturbative(nucleus, model, symmetries, basis, HFfield, &
       densities, iter, iunit, calculation_type)

!
!   Arguments
!
    type(nucleus_properties),intent(inout) :: nucleus
    type(modelization),intent(in) :: model
    type(symmetry_properties),intent(in) :: symmetries
    type(cylindrical_basis),intent(in) :: basis
    type(field),dimension(2),intent(in) :: HFfield
    type(local_density),dimension(2),intent(in) :: densities
    integer,intent(in) :: iter,iunit
    character(len=*) :: calculation_type
!
!   Local variables
!
    character(len=18),parameter :: subroutinename='h_odd_perturbative'
    integer, dimension(2) :: Npart,nb_qp
    integer :: Nsph,ih,il,NGz,NGr
    integer :: Nblocks,iblock
    integer :: k,n,shift
    integer :: mem_stat
    integer :: isospin,i,j,iHF,jHF
    type(field_Bi_contributions) :: h_HF_mtx_elem,quenching_factor
    real(kind=rk) :: sum_sz,expval_sz_unperturbed,ei,ej,avg_sz,avg_lz,hodd_sz,gs_q, &
         expval_sz_phi_i
    type(field_Bi_contributions),dimension(2) :: expval_sz,delta
    type(field_Bi_contributions) :: s_odd
    real(kind=rk),dimension(2) :: expval_sz_core
    real(kind=rk) :: B14, B15, B16, B17, B18, B19, B20, B21
    real(kind=rk) :: alpha_C, alpha_T, beta_C, beta_T, &
         h_W_alpha_C, h_W_alpha_T, h_W_beta_C, h_W_beta_T

    type(nucleus_properties) :: local_nucleus
    integer,dimension(2,nqp) :: i_odd_local
    type(field),dimension(2) :: HFfield_no_tensor
    integer,dimension(basis%size,2) :: permut_inv
    type(sp_state),dimension(basis%size,2) :: sp_states
    type(operator),dimension(:),allocatable  :: h_HF        ! Hartree-Fock Hamiltonian
    real(kind=rk),dimension(:,:),allocatable :: C           ! Expansion coefficients
    real(kind=rk),dimension(:),allocatable   :: Ctmp
    real(kind=rk),dimension(:,:),allocatable :: sp_energies ! Single-particle energies
    real(kind=rk),dimension(:,:),allocatable :: Vpair       ! Pairing matrix elements

    real(kind=rk) :: e_1st_unoccupied,fermi,Vpair_max
    type(expansion) :: psi_bar
    real(kind=rk) :: ovlp,tmp
    integer,dimension(basis%size) :: partner,order
    real(kind=rk),dimension(basis%size/2) :: quasi_pair_energy
    integer,dimension(basis%size/2) :: sp_label,tmp_label
    real(kind=rk),dimension(basis%size,basis%size) :: np_overlaps
    logical :: filled_qp

    integer :: count, m, max_rank
    real(kind=rk),parameter :: ph_energy_min = 1.0e12_rk
    real(kind=rk),dimension(:),allocatable :: ph_energy ! Particle-hole excitation energy wrt Fermi level


!
!   Auxiliary variables
!

    B14=model%parameters%B14
    B15=model%parameters%B15
    B16=model%parameters%B16
    B17=model%parameters%B17
    B18=model%parameters%B18
    B19=model%parameters%B19
    B20=model%parameters%B20
    B21=model%parameters%B21

    write(iunit,'(a,1x,a/)') 'Skyrme parametrization:',model%mean_field

    write(iunit,'(a)') 'Skyrme coupling constants involving tensor parameters:'
    write(iunit,'(a)') 'B14 = -(t1*x1+t2*x2)/8 + (te+to)/4        B15 = (t1-t2)/8 - (te-to)/4'
    write(iunit,'(a)') 'B16 = -3*(te+to)/8                        B17 = 3*(te-to)/8'
    write(iunit,'(a)') 'B18 = -(3*t1*x1-t2*x2)/32 + (3*te-to)/16  B19 = (3*t1+t2)/32 - (3*te+to)/16'
    write(iunit,'(a)') 'B20 = 3*(3*te-to)/16                      B21 = -3*(3*te+to)/16'
    write(iunit,'(a)') 'where 3*te = T and 3*to = U from original Skyrme paper Nucl. Phys. 9, 615 (1958)'
    write(iunit,'(a,f12.6,2x)',advance='no') 'B14 =',B14
    write(iunit,'(a,f12.6)') 'B15 =',B15
    write(iunit,'(a,f12.6,2x)',advance='no') 'B16 =',B16
    write(iunit,'(a,f12.6)') 'B17 =',B17
    write(iunit,'(a,f12.6,2x)',advance='no') 'B18 =',B18
    write(iunit,'(a,f12.6)') 'B19 =',B19
    write(iunit,'(a,f12.6,2x)',advance='no') 'B20 =',B20
    write(iunit,'(a,f12.6)') 'B21 =',B21
    write(iunit,'(/a)') 'Coupling constants of J^2 terms in Hamiltonian density:'
    write(iunit,'(a)') '* alpha = alpha_C + alpha_T = B14+B15-(B16+B17)'
    write(iunit,'(a)') '   alpha_C = (t1-t2)/8 - (t1*x1+t2*x2)/8 = B14+B15+2*(B16+B17)/3'
    write(iunit,'(a)') '   alpha_T = 5*to/4 = -5*(B16+B17)/3  (like-particle, T=1)'
    write(iunit,'(a)') '* beta = beta_C + beta_T = B14-B16'
    write(iunit,'(a)') '   beta_C = -(t1*x1+t2*x2)/8 = B14+2*B16/3'
    write(iunit,'(a)') '   beta_T = 5*(te+to)/8 = -5*B16/3  (np, T=0)'
    alpha_C = B14+B15+2*(B16+B17)/3
    alpha_T = -5*(B16+B17)/3
    beta_C = B14+2*B16/3
    beta_T = -5*B16/3
    write(iunit,'(3(a,f12.6,2x))') 'alpha =',alpha_C+alpha_T,'alpha_C =',alpha_C, &
         'alpha_T =',alpha_T
    write(iunit,'(3(a,f12.6,2x))') 'beta  =',beta_C+beta_T, &
         'beta_C  =',beta_C, 'beta_T  =',beta_T
    
    write(iunit,'(/a,i3)') 'Iteration # ',iter

    Npart(:)=(/ nucleus%N , nucleus%Z /)
    local_nucleus = nucleus
    i_odd_local = i_odd

!
!   No-tensor part of the HF field
!
    NGz=basis%NGz
    NGr=basis%NGr
    do isospin = 1, 2
       allocate(HFfield_no_tensor(isospin)%central(iHmin:NGz/2,NGr), &
            HFfield_no_tensor(isospin)%effective_mass(iHmin:NGz/2,NGr), &
            HFfield_no_tensor(isospin)%spin_orbit(iHmin:NGz/2,NGr), &
            HFfield_no_tensor(isospin)%Coulomb(iHmin:NGz/2,NGr))
       allocate(HFfield_no_tensor(isospin)%S%z(iHmin:NGz/2,NGr), &
            HFfield_no_tensor(isospin)%S%rho(iHmin:NGz/2,NGr), &
            HFfield_no_tensor(isospin)%S%phi(iHmin:NGz/2,NGr))
       allocate(HFfield_no_tensor(isospin)%alpha%z(iHmin:NGz/2,NGr), &
            HFfield_no_tensor(isospin)%alpha%rho(iHmin:NGz/2,NGr), &
            HFfield_no_tensor(isospin)%alpha%phi(iHmin:NGz/2,NGr))
       allocate(HFfield_no_tensor(isospin)%C%z(iHmin:NGz/2,NGr), &
            HFfield_no_tensor(isospin)%C%rho(iHmin:NGz/2,NGr), &
            HFfield_no_tensor(isospin)%C%phi(iHmin:NGz/2,NGr))
       allocate(HFfield_no_tensor(isospin)%D%z(iHmin:NGz/2,NGr), &
            HFfield_no_tensor(isospin)%D%rho(iHmin:NGz/2,NGr), &
            HFfield_no_tensor(isospin)%D%phi(iHmin:NGz/2,NGr))
       do i=1,3
          allocate(HFfield_no_tensor(isospin)%tensor_so(i)%z(iHmin:NGz/2,NGr), &
               HFfield_no_tensor(isospin)%tensor_so(i)%rho(iHmin:NGz/2,NGr), &
               HFfield_no_tensor(isospin)%tensor_so(i)%phi(iHmin:NGz/2,NGr))
       end do
       allocate(HFfield_no_tensor(isospin)%div_s(iHmin:NGz/2,NGr))
       HFfield_no_tensor(isospin)%central = HFfield(isospin)%central
       HFfield_no_tensor(isospin)%Coulomb = HFfield(isospin)%Coulomb
       HFfield_no_tensor(isospin)%effective_mass = HFfield(isospin)%effective_mass
       HFfield_no_tensor(isospin)%spin_orbit = HFfield(isospin)%spin_orbit
       do i = 1, 3
          HFfield_no_tensor(isospin)%tensor_so(i)%z = 0
          HFfield_no_tensor(isospin)%tensor_so(i)%rho = 0
          HFfield_no_tensor(isospin)%tensor_so(i)%phi = 0
       end do
       HFfield_no_tensor(isospin)%S%z(:,:) = 0.0
       HFfield_no_tensor(isospin)%S%rho(:,:) = 0.0
       HFfield_no_tensor(isospin)%S%phi(:,:) = 0.0
       HFfield_no_tensor(isospin)%alpha%z(:,:) = 0.0
       HFfield_no_tensor(isospin)%alpha%rho(:,:) = 0.0
       HFfield_no_tensor(isospin)%alpha%phi(:,:) = 0.0
       HFfield_no_tensor(isospin)%C%z(:,:) = 0.0
       HFfield_no_tensor(isospin)%C%rho(:,:) = 0.0
       HFfield_no_tensor(isospin)%C%phi(:,:) = 0.0
       HFfield_no_tensor(isospin)%D%z(:,:) = 0.0
       HFfield_no_tensor(isospin)%D%rho(:,:) = 0.0
       HFfield_no_tensor(isospin)%D%phi(:,:) = 0.0
       HFfield_no_tensor(isospin)%div_s(:,:) = 0.0
    end do

!
!   Memeory allocation
!

    Nsph=(n_max+1)*(l_max+1)**2
    do isospin=1,2
       do i=1,basis%size
          allocate(sp_states(i,isospin)%coef%cylHO(basis%size),stat=mem_stat)
          if (mem_stat/=0) call fatal(modulename, subroutinename, &
               'Problem with allocate(sp_states%coef%cylHO)')
          allocate(sp_states(i,isospin)%coef%sphHO(Nsph),stat=mem_stat)
          if (mem_stat/=0) call fatal(modulename, subroutinename, &
               'Problem with allocate(sp_states%coef%sphHO)')
          allocate(sp_states(i,isospin)%wf(iHmin:basis%NGz/2,basis%NGr),stat=mem_stat)
          if (mem_stat/=0) call fatal(modulename, subroutinename, &
               'Problem with allocate(sp_states%wf)')
          sp_states(i,isospin)%tau=2*MOD(isospin,2) - 1
          sp_states(i,isospin)%Omega=0
          sp_states(i,isospin)%pi=0
          sp_states(i,isospin)%energy=0.0
          sp_states(i,isospin)%occupation=0.0
          sp_states(i,isospin)%avg_j=0.0
          sp_states(i,isospin)%avg_l=0.0
          sp_states(i,isospin)%avg_pi=0.0
          sp_states(i,isospin)%block=0
          sp_states(i,isospin)%block_rank=0
          sp_states(i,isospin)%pair_partner(:)=0
          sp_states(i,isospin)%coef%cylHO(:)=0.0
          sp_states(i,isospin)%coef%sphHO(:)=0.0
          do ih=iHmin,basis%NGz/2
             if (ih==0) cycle
             do il=1,basis%NGr
                sp_states(i,isospin)%wf(ih,il)=spinor(cmplx(0),cmplx(0))
             end do
          end do
       end do
    end do

    Nblocks=basis%Nblocks

    allocate(h_HF(Nblocks),stat=mem_stat)
    if (mem_stat/=0) call fatal(modulename, subroutinename, &
         'Problem with allocate(h_HF)')

    do k=1,Nblocks
       n=basis%block_size(k)
       allocate(h_HF(k)%HO(n,n),h_HF(k)%eigenval(n),stat=mem_stat)
       if (mem_stat/=0) call fatal(modulename, subroutinename, &
            'Problem with allocate(h_HF%HO etc...)')
    end do

    allocate(sp_energies(basis%size,2),stat=mem_stat)
    if (mem_stat/=0) call fatal(modulename, subroutinename, &
         'Problem with allocate(sp_energies)')

    allocate(Vpair(basis%size,basis%size),stat=mem_stat)
    if (mem_stat/=0) call fatal(modulename, subroutinename, &
         'Problem with allocate(Vpair)')

    allocate(psi_bar%cylHO(1:basis%size))
    psi_bar%cylHO=0.0

    do isospin = 1, 2

       shift=0

!
!      Initialization of expansion coefficients to 0
!      (array filled in by blocks ==> need to clean everything at every iteration)
!

       do iHF=1,basis%size
          sp_states(iHF,isospin)%coef%cylHO(:)=0.0
       end do
       sp_states(:,isospin)%pi=0

       loop_blocks: do k=1,Nblocks

          n=basis%block_size(k)

!
!         Hamiltonian matrix elements in the Harmonic Oscillator (HO) basis
!

          if (trim(calculation_type)=='exact') then
             call hamiltonian_cylHO(h_HF(k)%HO(:,:),k,isospin,basis,HFfield(isospin),iter)
          else if (trim(calculation_type)=='perturbative') then
             call hamiltonian_cylHO(h_HF(k)%HO(:,:),k,isospin,basis,HFfield_no_tensor(isospin),iter)
          else
             stop 'Wrong calculation type'
          end if

!
!         Hamiltonian diagonalization
!

          allocate(C(n,n),Ctmp(n))

          call diagon(h_HF(k)%HO(:,:), h_HF(k)%eigenval(:), C, n)

!
!         Fill in arrays with eigenvalues and eigenvectors
!         (largest expansion coeff. chosen to be positive)
!

          sp_energies(shift+1:shift+n,isospin)=h_HF(k)%eigenval(:)

          do i=1,n
             do j=1,n
                Ctmp(j)=abs(C(j,i))
             end do
             j=maxloc(Ctmp,1)
             sp_states(i+shift,isospin)%coef%cylHO(1+shift:n+shift)=C(:,i)* &
                  sign(1.0_rk,C(j,i))
             sp_states(i+shift,isospin)%energy=h_HF(k)%eigenval(i)
             sp_states(i+shift,isospin)%block_rank=i
          end do

          deallocate(C,Ctmp)

          if (symmetries%parity) sp_states(1+shift:n+shift,isospin)%pi= &
               basis%parity(1+shift)
          sp_states(1+shift:n+shift,isospin)%Omega=basis%Omega(1+shift)
          sp_states(1+shift:n+shift,isospin)%block=k

          shift=shift+n

       end do loop_blocks

!
!      Sort all eigenvalues in increasing order to define hole and particle states
!      and all eigenvectors accordingly
!

       call sort(sp_energies(:,isospin), basis%size, 'increasing', permut = permut(:,isospin))

       call permutation(sp_states(:,isospin), permut(:,isospin), &
            basis%size, 'normal')

       permut_inv(:,isospin) = 0
       do iHF = 1, basis%size
          sp_states(iHF,isospin)%HF_state=iHF
          permut_inv(permut(iHF,isospin),isospin) = iHF
       end do
!
!      Single-particle wave functions and their gradients at Gauss--Hermite--Laguerre
!      integration points
!

       do iHF = 1, basis%size
          call sp_wave_functions(sp_states(iHF,isospin), basis)
       end do

!
!      Occupation numbers: pairing or pure Hartree--Fock
!

       if (iter==1 .and. model%initial_potential=='WS' .and. model%pairing/='seniority') then

          sp_states(1:Npart(isospin),isospin)%occupation=1.0_rk
          sp_states(Npart(isospin)+1:basis%size,isospin)%occupation=0.0
          local_nucleus%pairing%pairing_energy=0.0

       else

!
!      Initializations:
!
!      * chemical potential to Fermi level (half sum of last occ. and 1st unocc. levels)
!      * pairing energy to 0
!      * pairing gaps and qp quantities to 0
!      * occupation probabilities to the Hatree-Fock values
!
          iHF=Npart(isospin)
          e_1st_unoccupied=sp_states(iHF+1,isospin)%energy
          if (iHF/=0) then
             fermi=0.5_rk*(sp_states(iHF,isospin)%energy+e_1st_unoccupied)
          else
             fermi = e_1st_unoccupied
          end if

          local_nucleus%pairing(isospin)%chemical_potential=fermi
          local_nucleus%pairing(isospin)%pairing_energy=0.0
          local_nucleus%pairing(isospin)%average_gap=0.0
          local_nucleus%pairing(isospin)%gap(:)=0.0
          local_nucleus%pairing(isospin)%min_qp_index=0
          local_nucleus%pairing(isospin)%min_qp_energy=0

!
!         Determination of quasi-pairs (which would be Kramers degenerate
!         in the equal-filling approximation) and their energy
!

          partner(:)=0
          sp_states(:,isospin)%pair_partner(isospin)=0
          quasi_pair_energy(:)=1e6_rk
          sp_label(:)=0
          n=0

          do iHF=1,basis%size
             if (sp_states(iHF,isospin)%Omega<0) cycle
!            definition of time-reversed state
             psi_bar=time_reversal(sp_states(iHF,isospin)%coef,basis)
             ovlp=1e-09_rk
             partner(iHF)=0
!            search for sp state with opposite Omega and largest overlap with psi_bar
             do jHF=1,basis%size
                if (sp_states(jHF,isospin)%Omega /= -sp_states(iHF,isospin)%Omega) cycle
                tmp=dot_product(sp_states(jHF,isospin)%coef%cylHO,psi_bar%cylHO)
                if (abs(tmp)>abs(ovlp)) then
                   ovlp=tmp
                   partner(iHF)=jHF
                end if
             end do
             jHF=partner(iHF)
             if (jHF /= 0) then
                n=n+1
                if (n>basis%size/2) call fatal(modulename,subroutinename,'Too many unfound partners')
                quasi_pair_energy(n)=sp_states(iHF,isospin)%energy+sp_states(jHF,isospin)%energy
                sp_label(n)=iHF
                sp_states(iHF,isospin)%pair_partner(isospin)=jHF
                sp_states(jHF,isospin)%pair_partner(isospin)=iHF
             else
                write(out_unit,'(a,i4)') 'iHF=',iHF
                call fatal(modulename,subroutinename,'Partner state not found')
             end if
          end do

          call sort(quasi_pair_energy(:), basis%size/2, 'increasing', permut=order(:))

          tmp_label(:)=0
          do i=1,basis%size/2
             tmp_label(i)=sp_label(order(i))
          end do
          sp_label(:)=tmp_label(:)

!
!         Occupation numbers in the Hartree--Fock approximation
!

!
!         Finding the closest rank of the blocked state to the Fermi level
!
          if (local_nucleus%fix_rank == 'n') then
             max_rank = maxval(basis%block_size(:),1)
             if (max_rank <= 0) call fatal(modulename,subroutinename,"max_rank <= 0")
             allocate(ph_energy(max_rank))
             ph_energy(:) = ph_energy_min
             if (symmetries%parity) then
                do j = 1, 3
                   if (local_nucleus%K(isospin,j) == 0) cycle
                   do m = 1, max_rank
                      do jHF = 1, basis%size
                         if (sp_states(jHF,isospin)%Omega /= local_nucleus%K(isospin,j) .or. &
                              sp_states(jHF,isospin)%pi /= local_nucleus%pi(isospin,j)) cycle
                         if (sp_states(jHF,isospin)%block_rank == m) ph_energy(m) = &
                              abs(sp_states(jHF,isospin)%energy-fermi)
                      end do
                   end do
                   if (minval(ph_energy) == ph_energy_min) call fatal(modulename,subroutinename, &
                        'Rank not found')
                   local_nucleus%rank(isospin,j) = minloc(ph_energy,dim=1)
                   if (debugmode /= 'low') write(out_unit,'(4x,3(a,I2))') 'Blocked K^pi state =', &
                        local_nucleus%K(isospin,j),'/2,  parity = ',local_nucleus%pi(isospin,j), &
                        ',  Rank = ', local_nucleus%rank(isospin,j)
                   do m = 1, max_rank
                      if (m == local_nucleus%rank(isospin,j)) cycle
                      if (abs(ph_energy(m)-ph_energy(local_nucleus%rank(isospin,j))) <= 0.5_rk) then
                         write(out_unit,'(4x,2(a,i2),a)') 'NOTE: The distances from Fermi level of ranks ', &
                              local_nucleus%rank(isospin,j),' and ',m,' differ by less than 500 keV'
                         write(out_unit,'(3x,a,i2,a,f12.6)') 'Distance from Fermi level for rank ', &
                              local_nucleus%rank(isospin,j), ' = ', ph_energy(local_nucleus%rank(isospin,j))
                         write(out_unit,'(3x,a,i2,a,f12.6)') 'Distance from Fermi level for rank ', &
                              m, ' = ', ph_energy(m)
                      end if
                   end do
                end do
             else
                do j = 1, 3
                   if (nucleus%K(isospin,j) == 0) cycle
                   do m = 1, max_rank
                      do jHF = 1, basis%size
                         if (sp_states(jHF,isospin)%Omega /= nucleus%K(isospin,j)) cycle
                         if (sp_states(jHF,isospin)%block_rank == m) ph_energy(m) = &
                              abs(sp_states(jHF,isospin)%energy-fermi)
                      end do
                   end do
                   if (minval(ph_energy) == ph_energy_min) call fatal(modulename,subroutinename, &
                        'Rank not found')
                   nucleus%rank(isospin,j) = minloc(ph_energy,dim=1)
                   if (debugmode /= 'low') write(out_unit,'(4x,3(a,I2))') 'Blocked K^pi state =', &
                        nucleus%K(isospin,j),'/2,  parity = ',nucleus%pi(isospin,j), &
                        ',  Rank = ', nucleus%rank(isospin,j)
                   do m = 1, max_rank
                      if (m == nucleus%rank(isospin,j)) cycle
                      if (abs(ph_energy(m)-ph_energy(nucleus%rank(isospin,j))) <= 0.5_rk) then
                         write(out_unit,'(4x,2(a,i2),a)') 'NOTE: The distances from Fermi level of ranks ', &
                              nucleus%rank(isospin,j),' and ',m,' differ by less than 500 keV'
                         write(out_unit,'(3x,a,i2,a,f12.6)') 'Distance from Fermi level for rank ', &
                              nucleus%rank(isospin,j), ' = ', ph_energy(nucleus%rank(isospin,j))
                         write(out_unit,'(3x,a,i2,a,f12.6)') 'Distance from Fermi level for rank ', &
                              m, ' = ', ph_energy(m)
                      end if
                   end do
                end do
             end if
             deallocate(ph_energy)
          end if

!         Determination of the blocked states

          do i=1,nqp
             i_odd_local(isospin,i)=0
             if (local_nucleus%K(isospin,i) == 0) cycle
             loop_iHF: do iHF=1,basis%size
                if (sp_states(iHF,isospin)%Omega==local_nucleus%K(isospin,i) .and. &
                     sp_states(iHF,isospin)%block_rank==local_nucleus%rank(isospin,i)) then
                   if (symmetries%parity) then
                      if (sp_states(iHF,isospin)%pi==local_nucleus%pi(isospin,i)) then
                         i_odd_local(isospin,i)=iHF
                         exit loop_iHF
                      end if
                   else
                      i_odd_local(isospin,i)=iHF
                      exit loop_iHF
                   end if
                end if
             end do loop_iHF
          end do
             
!         Initialization of occupation numbers and type component

          sp_states(:,isospin)%occupation = 0.0
          sp_states(:,isospin)%type = 0

!         Occupation numbers of the blocked states

          do i=1,nqp
             if (i_odd_local(isospin,i)/=0) then
                if (model%blocking == 'SCB') then
                   sp_states(i_odd_local(isospin,i),isospin)%occupation = 1.0_rk
                else if (model%blocking == 'EFA') then
                   sp_states(i_odd_local(isospin,i),isospin)%occupation = 0.5_rk
                   jHF = sp_states(i_odd_local(isospin,i),isospin)%pair_partner(isospin)
                   sp_states(jHF,isospin)%occupation = 0.5_rk
                else
                   call fatal(modulename,subroutinename,'Unrecognized blocking scheme ' &
                        //trim(model%blocking))
                end if
             else if (local_nucleus%K(isospin,i) /= 0) then
                call fatal(modulename,subroutinename,'Blocked state not found')
             end if
          end do

!         Occupation numbers of remaining, lowest nb_qp(isospin) quasi-pairs 
!         of states set to 1

          nb_qp(isospin) = (Npart(isospin) - nb_blocked_part(isospin))/2
          jHF = 1
          do k = 1, nb_qp(isospin)
             filled_qp = .false.
             loop_fill_qp: do j = jHF, basis%size/2
                iHF = sp_label(j)
                do i = 1, nqp
                   if (iHF == i_odd_local(isospin,i) .or. partner(iHF) == i_odd_local(isospin,i)) cycle loop_fill_qp
                end do
                sp_states(iHF,isospin)%occupation = 1.0_rk
                sp_states(partner(iHF),isospin)%occupation = 1.0_rk
                filled_qp = .true.
                jHF=j+1
                exit loop_fill_qp
             end do loop_fill_qp
             if (.not.filled_qp) then
                write(out_unit,'(a,i2,2(1x,a,i4),2x,a,3i5)') &
                     'isospin=',isospin,'k=',k,'jHF=',jHF,'i_odd_local(isospin,:)=',i_odd_local(isospin,:)
                call fatal(modulename,subroutinename,"Quasi-pair not filled")
             end if
          end do

!         Define HF_level, HF_state and type components of sp_states

          i=0
          do k=1,basis%size
             sp_states(k,isospin)%HF_level = k
             sp_states(k,isospin)%HF_state = k
             if (sp_states(k,isospin)%occupation >= 0.5_rk) then
                sp_states(k,isospin)%type = 1
                i=i+1
             else
                sp_states(k,isospin)%type = 0
             end if
          end do

          if (model%blocking == 'EFA') then
             do k = 1, nqp
                if (i_odd_local(isospin,k) == 0) cycle
                jHF = sp_states(i_odd_local(isospin,k),isospin)%pair_partner(isospin)
                sp_states(jHF,isospin)%type = 0
                i = i - 1
             end do
          end if

          if (i /= Npart(isospin)) then
             do iHF = 1, basis%size
                write(out_unit,'(a1,2x,2(i4,1x),a2,a1,3(i4,1x))') &
                     nucleon_type(isospin),iHF,sp_states(iHF,isospin)%Omega,'/2', &
                     signparity(sp_states(iHF,isospin)%pi),iHF,partner(iHF), &
                     sp_states(iHF,isospin)%type
             end do
             write(out_unit,'(/a,i2,1x,a,i4,a,i4)') 'isospin=',isospin,&
                  'number of hole states =', i, ' whereas Npart=', Npart(isospin)
             call fatal(modulename,subroutinename, &
                  'Problem with the particle or hole character of s.p. states')
          end if

!
!         Matrix elements of the effective interaction in the pairing channel
!         (replacing the Skyrme interaction)
!

          call pairing_matrix_elements(sp_states(:,isospin), densities, isospin, &
               local_nucleus, model, basis, Vpair, partner, iter)

!
!         Occupation numbers in the BCS approximation
!
!         Check if BCS calculation is needed:
!         * if max{|Vpair(i,j)|,i,j=1..basis%size} is larger than 1 keV, call BCS
!         * if not, do nothing ==> HF solution
!

          Vpair_max = maxval(abs(Vpair))
          if (abs(Vpair_max)>1e-3_rk) then
             call BCS(local_nucleus, model, basis, sp_states(:,isospin), Vpair, isospin, &
                  i_odd_local(isospin,:))
          end if

       end if

    end do

!
!   Checking number of blocked s.p. states
!

    count = 0
    do isospin = 1, 2
       do i = 1, nqp
          if (i_odd_local(isospin,i) /= 0) count = count + 1
       end do
    end do
    if (count > 1) call warning(modulename, subroutinename, 'More than 1 blocked state')

!
!   Calculation of matrix elements and 1st-order eigenstates
!

    s_odd%tot = 0.0
    s_odd%S = 0.0
    s_odd%A = 0.0
    s_odd%C = 0.0
    s_odd%D = 0.0

    expval_sz_phi_i = 0

    do isospin = 1, 2

       write(iunit,'(/a)') 'Particle type: ' // particle_type(isospin)

       expval_sz(isospin)%tot = 0.0
       expval_sz(isospin)%K = 0.0
       expval_sz(isospin)%U = 0.0
       expval_sz(isospin)%so = 0.0
       expval_sz(isospin)%W = 0.0
       expval_sz(isospin)%S = 0.0
       expval_sz(isospin)%A = 0.0
       expval_sz(isospin)%C = 0.0
       expval_sz(isospin)%D = 0.0
  
       do iHF = 1, basis%size

!          if (sp_states(iHF,isospin)%occupation < 1e-3_rk) cycle
          ei = sp_states(iHF,isospin)%energy
          if (ei >= 0.0) cycle

          call sz_matrix_element(sp_states(iHF,isospin), &
               sp_states(iHF,isospin),basis,expval_sz_unperturbed)
          expval_sz(isospin)%tot = expval_sz(isospin)%tot + sp_states(iHF,isospin)%occupation * &
               expval_sz_unperturbed
          sp_states(iHF,isospin)%avg_sz = expval_sz_unperturbed

          iblock=sp_states(iHF,isospin)%block
          shift=sum(basis%block_size(1:iblock),1)-basis%block_size(iblock)

          do j = 1, basis%block_size(iblock)

             jHF = permut_inv(j + shift,isospin)
             if (jHF /= iHF) cycle ! Diagonal h_HF matrix elements in HF basis only
             ej = sp_states(jHF,isospin)%energy

             call sz_matrix_element(sp_states(iHF,isospin),sp_states(jHF,isospin),basis,sum_sz)
             if (abs(sum_sz) < 0.01_rk) cycle

             call h_HF_matrix_element(sp_states(iHF,isospin), &
                  sp_states(jHF,isospin),isospin,basis,model,HFfield(isospin), &
                  densities,h_HF_mtx_elem)

             if (jHF /= iHF) then

                hodd_sz = sum_sz*(h_HF_mtx_elem%S(0)+h_HF_mtx_elem%A(0)+h_HF_mtx_elem%C(0)+ &
                     h_HF_mtx_elem%D(0))/(ei-ej)
                expval_sz(isospin)%tot = expval_sz(isospin)%tot + sp_states(iHF,isospin)%occupation * &
                     2*hodd_sz
                do i = 1, 21
                   expval_sz(isospin)%S(i) = expval_sz(isospin)%S(i) + sp_states(iHF,isospin)%occupation * &
                        2*sum_sz*h_HF_mtx_elem%S(i)/(ei-ej)
                   expval_sz(isospin)%A(i) = expval_sz(isospin)%A(i) + sp_states(iHF,isospin)%occupation * &
                        2*sum_sz*h_HF_mtx_elem%A(i)/(ei-ej)
                   expval_sz(isospin)%C(i) = expval_sz(isospin)%C(i) + sp_states(iHF,isospin)%occupation * &
                        2*sum_sz*h_HF_mtx_elem%C(i)/(ei-ej)
                   expval_sz(isospin)%D(i) = expval_sz(isospin)%D(i) + sp_states(iHF,isospin)%occupation * &
                        2*sum_sz*h_HF_mtx_elem%D(i)/(ei-ej)
                end do

                write(iunit,'(/2(a,i4,1x,a,f11.6,2x,i3,a2,a1,2x),3(a,f11.6,1x))') &
                     'iHF=',iHF,' e =', ei, &
                     sp_states(iHF,isospin)%Omega,'/2', &
                     signparity(sp_states(iHF,isospin)%pi), &
                     'jHF=',jHF,' e =', sp_states(jHF,isospin)%energy, &
                     sp_states(jHF,isospin)%Omega,'/2', &
                     signparity(sp_states(jHF,isospin)%pi), &
                     '<iHF|h_HF|jHF> =',h_HF_mtx_elem%tot,'<iHF|s_z|jHF>  =',sum_sz, &
                     '<i|h_odd|j><i|s_z|j>/(e_i-e_j) =',hodd_sz

             else

                write(iunit,'(/1(a,i4,1x,a,f11.6,2x,i3,a2,a1,2x),2(a,f11.6,1x))') &
                     'iHF=',iHF,' e =', ei, &
                     sp_states(iHF,isospin)%Omega,'/2', &
                     signparity(sp_states(iHF,isospin)%pi), &
                     '<iHF|h_HF|jHF> =',h_HF_mtx_elem%tot,'<iHF|s_z|jHF>  =',sum_sz

             end if

             h_W_alpha_C = alpha_C * h_HF_mtx_elem%W(15)/B15
             h_W_alpha_T = alpha_T/5 * (2*h_HF_mtx_elem%W(15)/B15 &
                  - 3 * h_HF_mtx_elem%W(17)/B17)
             h_W_beta_C = beta_C * (h_HF_mtx_elem%W(14)/B14 - h_HF_mtx_elem%W(15)/B15)
             h_W_beta_T = beta_T/5 * (2*(h_HF_mtx_elem%W(14)/B14 - h_HF_mtx_elem%W(15)/B15) &
                  - 3 * (h_HF_mtx_elem%W(16)/B16 - h_HF_mtx_elem%W(17)/B17))

             write(iunit,'(8x,4(a,f11.6,1x)/,8x,4(a,f11.6,1x))') &
                  'h_K =',h_HF_mtx_elem%K(0),'h_U =',h_HF_mtx_elem%U(0), &
                  'h_so =', h_HF_mtx_elem%so(0),'h_W =',h_HF_mtx_elem%W(0), &
                  'h_S =',h_HF_mtx_elem%S(0),'h_A =',h_HF_mtx_elem%A(0), &
                  'h_C  =',h_HF_mtx_elem%C(0),'h_D =',h_HF_mtx_elem%D(0)
             write(iunit,'(8x,2(a,f11.6,1x)/,8x,2(a,f11.6,1x))') &
                  'h_W(B14)=',h_HF_mtx_elem%W(14),'h_W(B15)=',h_HF_mtx_elem%W(15), &
                  'h_W(B16)=',h_HF_mtx_elem%W(16),'h_W(B17)=',h_HF_mtx_elem%W(17)
             write(iunit,'(8x,2(a,f11.6,1x)/,8x,2(a,f11.6,1x))') &
                  'h_W(alpha_C)=',h_W_alpha_C,'h_W(alpha_T)=',h_W_alpha_T, &
                  'h_W(beta_C) =',h_W_beta_C, 'h_W(beta_T) =',h_W_beta_T
             write(iunit,'(8x,3(a,f11.6,1x)/, 2(8x,4(a,f11.6,1x)/),8x,2(a,f11.6,1x))') &
                  'h_S(B9) =',h_HF_mtx_elem%S(9),'h_S(B10)=',h_HF_mtx_elem%S(10), &
                  'h_S(B11)=',h_HF_mtx_elem%S(11),'h_S(B12)=',h_HF_mtx_elem%S(12), &
                  'h_S(B13)=',h_HF_mtx_elem%S(13),'h_S(B14)=',h_HF_mtx_elem%S(14), &
                  'h_S(B15)=',h_HF_mtx_elem%S(15), &
                  'h_S(B16)=',h_HF_mtx_elem%S(16),'h_S(B17)=',h_HF_mtx_elem%S(17), &
                  'h_S(B18)=',h_HF_mtx_elem%S(18),'h_S(B19)=',h_HF_mtx_elem%S(19), &
                  'h_S(B20)=',h_HF_mtx_elem%S(20),'h_S(B21)=',h_HF_mtx_elem%S(21)
             write(iunit,'(8x,3(a,f11.6,1x))') &
                  'h_A(B3) =',h_HF_mtx_elem%A(3),'h_A(B4) =',h_HF_mtx_elem%A(4), &
                  'h_A(B9) =',h_HF_mtx_elem%A(9)
             write(iunit,'(8x,2(a,f11.6,1x))') &
                  'h_C(B14)=',h_HF_mtx_elem%C(14),'h_C(B15)=',h_HF_mtx_elem%C(15)
             write(iunit,'(8x,2(a,f11.6,1x))') &
                  'h_D(B16)=',h_HF_mtx_elem%D(16),'h_D(B17)=',h_HF_mtx_elem%D(17)
          end do

       end do

       expval_sz(isospin)%S(0) = sum(expval_sz(isospin)%S(1:21),1)
       expval_sz(isospin)%A(0) = sum(expval_sz(isospin)%A(1:21),1)
       expval_sz(isospin)%C(0) = sum(expval_sz(isospin)%C(1:21),1)
       expval_sz(isospin)%D(0) = sum(expval_sz(isospin)%D(1:21),1)
       expval_sz_core(isospin) = expval_sz(isospin)%tot

       if (mod(Npart(isospin),2) == 1) then

          gs_q = gs(isospin)

          iHF = i_odd_local(isospin,1)
          iblock=sp_states(iHF,isospin)%block
          shift=sum(basis%block_size(1:iblock))-basis%block_size(iblock)
          expval_sz_phi_i = sp_states(iHF,isospin)%avg_sz
          s_odd%tot = sp_states(iHF,isospin)%avg_sz
          do j = 1, basis%block_size(iblock)
             jHF = permut_inv(j + shift,isospin)
             if (jHF == iHF) cycle
             call h_HF_matrix_element(sp_states(iHF,isospin), &
                  sp_states(jHF,isospin),isospin,basis,model,HFfield(isospin), &
                  densities,h_HF_mtx_elem)
             call sz_matrix_element(sp_states(iHF,isospin),sp_states(jHF,isospin),basis,sum_sz)
             ej = sp_states(jHF,isospin)%energy
             hodd_sz = sum_sz*(h_HF_mtx_elem%S(0)+h_HF_mtx_elem%A(0)+h_HF_mtx_elem%C(0)+ &
                  h_HF_mtx_elem%D(0))/(ei-ej)
             s_odd%tot = s_odd%tot + sp_states(iHF,isospin)%occupation * 2*hodd_sz  ! s_odd = <psi_i|sz|psi_i>
             s_odd%S(:) = s_odd%S(:) + sp_states(iHF,isospin)%occupation * &
                  2 * sum_sz*h_HF_mtx_elem%S(:)/(ei-ej)
             s_odd%A(:) = s_odd%A(:) + sp_states(iHF,isospin)%occupation * &
                  2 * sum_sz*h_HF_mtx_elem%A(:)/(ei-ej)
             s_odd%C(:) = s_odd%C(:) + sp_states(iHF,isospin)%occupation * &
                  2 * sum_sz*h_HF_mtx_elem%C(:)/(ei-ej)
             s_odd%D(:) = s_odd%D(:) + sp_states(iHF,isospin)%occupation * &
                  2 * sum_sz*h_HF_mtx_elem%D(:)/(ei-ej)
          end do

! s_odd = <psi|sz|psi_i> = <phi|sz|phi_i>  + order-1 term in h_odd

          expval_sz_core(isospin) = expval_sz_core(isospin) - s_odd%tot
          expval_sz(isospin)%S(:) = expval_sz(isospin)%S(:) - s_odd%S(:)
          expval_sz(isospin)%A(:) = expval_sz(isospin)%A(:) - s_odd%A(:)
          expval_sz(isospin)%C(:) = expval_sz(isospin)%C(:) - s_odd%C(:)
          expval_sz(isospin)%D(:) = expval_sz(isospin)%D(:) - s_odd%D(:)

          write(iunit,'(6(a,f10.6,1x))') &
            '<sz>' // nucleon_type(isospin) // '=', expval_sz(isospin)%tot, &
            's_odd =',s_odd%tot, '<phi_i|s_z|phi_i> =',sp_states(iHF,isospin)%avg_sz, &
            '<sz>' // nucleon_type(isospin) // '(core,S)=', expval_sz(isospin)%S(0), &
            '<sz>' // nucleon_type(isospin) // '(core,A)=', expval_sz(isospin)%A(0), &
            '<sz>' // nucleon_type(isospin) // '(core,C)=', expval_sz(isospin)%C(0)

       else

          write(iunit,'(4(a,f10.6,1x))') &
            '<sz>' // nucleon_type(isospin) // '=', expval_sz(isospin)%tot, &
            '<sz>' // nucleon_type(isospin) // '(core,S)=', expval_sz(isospin)%S(0), &
            '<sz>' // nucleon_type(isospin) // '(core,A)=', expval_sz(isospin)%A(0), &
            '<sz>' // nucleon_type(isospin) // '(core,C)=', expval_sz(isospin)%C(0)

       end if

    end do

!
!   Memory deallocation
!

    do k=1,Nblocks
       deallocate(h_HF(k)%HO,h_HF(k)%eigenval)
    end do

    deallocate(h_HF,sp_energies,Vpair)

    deallocate(psi_bar%cylHO)

    do isospin = 1, 2
       deallocate(HFfield_no_tensor(isospin)%central)
       deallocate(HFfield_no_tensor(isospin)%effective_mass)
       deallocate(HFfield_no_tensor(isospin)%spin_orbit)
       deallocate(HFfield_no_tensor(isospin)%Coulomb)
       do i = 1, 3
          deallocate(HFfield_no_tensor(isospin)%tensor_so(i)%z)
          deallocate(HFfield_no_tensor(isospin)%tensor_so(i)%rho)
          deallocate(HFfield_no_tensor(isospin)%tensor_so(i)%phi)
       end do
       deallocate(HFfield_no_tensor(isospin)%S%z)
       deallocate(HFfield_no_tensor(isospin)%S%rho)
       deallocate(HFfield_no_tensor(isospin)%S%phi)
       deallocate(HFfield_no_tensor(isospin)%alpha%z)
       deallocate(HFfield_no_tensor(isospin)%alpha%rho)
       deallocate(HFfield_no_tensor(isospin)%alpha%phi)
       deallocate(HFfield_no_tensor(isospin)%C%z)
       deallocate(HFfield_no_tensor(isospin)%C%rho)
       deallocate(HFfield_no_tensor(isospin)%C%phi)
       deallocate(HFfield_no_tensor(isospin)%D%z)
       deallocate(HFfield_no_tensor(isospin)%D%rho)
       deallocate(HFfield_no_tensor(isospin)%D%phi)
       deallocate(HFfield_no_tensor(isospin)%div_s)
    end do

  end subroutine h_tensor_perturbative


  !=====================================================================
  !                       SUBROUTINE h_HF_matrix_element
  !---------------------------------------------------------------------
  ! Calculates the matrix elements of the Skyrme-Hartree-Fock Hamiltonian 
  ! between two HF states (decomposed in the HO basis). 
  !=====================================================================

  subroutine h_HF_matrix_element(psi1, psi2, isospin, basis, model, HFfield, &
       densities,h_HF)

!
!   Arguments
!
    type(sp_state),intent(in) :: psi1, psi2
    integer,intent(in) :: isospin
    type(cylindrical_basis),intent(in) :: basis
    type(modelization),intent(in) :: model
    type(field),intent(in) :: HFfield
    type(local_density),dimension(2),intent(in) :: densities
    type(field_Bi_contributions),intent(out) :: h_HF

!
!   Local variables
!
    integer :: ih,il
    integer :: i,j,i1,i2,k,j1,j2,shift,n
    integer :: l1,l2,s1,s2,nz1,nz2,nr1,nr2,k1,k2,Delta1,Delta2
    real(kind=rk) :: C1,C2
    real(kind=rk) :: norm,xi,eta,w
    real(kind=rk) :: B1,B2,B3,B4,B5,B6,B7,B8,B9,B10,B11,B12,B13,B14,B15,B16,B17,&
         B18,B19,B20,B21,alpha
    real(kind=rk) :: rho_q,rho_n,rho_p,rho_tot,rho_a, &
         tau_q,tau_tot,Delta_rho_q,Delta_rho_tot,div_J_q,div_J_tot, &
         s_z_tot,s_rho_tot,s_phi_tot
    real(kind=rk) :: j_tot_z,j_tot_rho,j_tot_phi,rot_j_tot_z,rot_j_tot_rho,rot_j_tot_phi, &
         s_tot_z,s_tot_rho,s_tot_phi,rot_s_tot_z,rot_s_tot_rho,rot_s_tot_phi
    real(kind=rk) :: t_tot_z,t_tot_rho,t_tot_phi,Delta_s_tot_z,Delta_s_tot_rho,Delta_s_tot_phi, &
         F_tot_z,F_tot_rho,F_tot_phi,trace_J_q,trace_J_tot
    real(kind=rk) :: div_s, div_s_q
    real(kind=rk) :: rho_qbar,s2_qbar
    type(vector) :: s_qbar
    real(kind=rk),dimension(1:21) :: sum_S,sum_A,sum_C,sum_W,sum_D
    real(kind=rk) :: sum_S_tot,sum_A_tot,sum_C_tot,sum_U,sum_K,sum_SO,S_z,S_rho,S_phi, &
         alpha_z,alpha_rho,alpha_phi,C_z,C_rho,C_phi,sum_tensor_so
    real(kind=rk) :: sum_D_tot,D_z,D_rho,D_phi
    type(spinor) :: chi1,chi2,phi1,phi2,ket_z,ket_rho,ket_phi,ket,bra
    type(vector_spinor) :: grad_phi1,grad_phi2,sigma_grad_phi2
    type(vector_spin_operator) :: sigma
    complex(kind=rk),dimension(2,2) :: Csigma,Jsigma_z,Jsigma_rho,Jsigma_phi

    norm=fz/(2.0_rk*bz*bp2)

!
!   Definition of the cylindrical Pauli matrices
!
    do i=1,2
       sigma%z(i,i)=-cmplx((-1)**i)
       sigma%rho(i,i)=cmplx(0)
       sigma%phi(i,i)=cmplx(0)
    end do

    sigma%z(1,2)=cmplx(0)
    sigma%rho(1,2)=cmplx(1) !*exp(-i * phi) with phi = 0 (axial symmetry)
    sigma%phi(1,2)=-i_cplx

    sigma%z(2,1)=cmplx(0)
    sigma%rho(2,1)=cmplx(1)
    sigma%phi(2,1)=i_cplx

!
!   Number of basis states preceding block # k
!

    k=psi1%block
    n=basis%block_size(k)
    shift=0
    do j=1,k
       shift=shift+basis%block_size(j)
    end do
    shift=shift-n

!
!   Auxiliary variables
!

    B1=model%parameters%B1
    B2=model%parameters%B2
    B3=model%parameters%B3
    B4=model%parameters%B4
    B5=model%parameters%B5
    B6=model%parameters%B6
    B7=model%parameters%B7
    B8=model%parameters%B8
    B9=model%parameters%B9
    B10=model%parameters%B10
    B11=model%parameters%B11
    B12=model%parameters%B12
    B13=model%parameters%B13
    B14=model%parameters%B14
    B15=model%parameters%B15
    B16=model%parameters%B16
    B17=model%parameters%B17
    B18=model%parameters%B18
    B19=model%parameters%B19
    B20=model%parameters%B20
    B21=model%parameters%B21
    alpha=model%parameters%alpha

!
!   Matrix element
!

    h_HF%tot = 0.0
    h_HF%K = 0.0
    h_HF%U = 0.0
    h_HF%SO = 0.0
    h_HF%S = 0.0
    h_HF%A = 0.0
    h_HF%C = 0.0
    h_HF%W = 0.0
    h_HF%D = 0.0

    if (psi1%block /= psi2%block) return ! Axial symmetry

    loop_i1: do i1=1,n

       j1=i1+shift
       l1=basis%Lambda(j1)
       s1=basis%ns(j1)
       nz1=basis%nz(j1)
       nr1=basis%nr(j1)
       k1=2*l1+s1
       Delta1=(l1-abs(l1))/2
       chi1=spinor(kronecker(s1,1),kronecker(s1,-1))
       C1 = psi1%coef%cylHO(j1)
       if (abs(C1) < TBME_tol) cycle

       loop_i2: do i2=1,n

          j2=i2+shift
          l2=basis%Lambda(j2)
          s2=basis%ns(j2)
          nz2=basis%nz(j2)
          nr2=basis%nr(j2)
          k2=2*l2+s2
          Delta2=(l2-abs(l2))/2
          chi2=spinor(kronecker(s2,1),kronecker(s2,-1))
          C2 = psi2%coef%cylHO(j2)
          if (abs(C2) < TBME_tol .or. abs(C1*C2) < TBME_tol**2) cycle

          sum_U=0.0
          sum_K=0.0
          sum_SO=0.0
          sum_S=0.0
          sum_S_tot=0.0
          sum_A=0.0
          sum_A_tot=0.0
          sum_C_tot=0.0
          sum_C=0.0
          sum_tensor_so=0.0
          sum_W=0.0
          sum_D=0.0
          sum_D_tot=0.0

          do ih=iHmin,basis%NGz/2

             if (ih==0) cycle
             xi=basis%xHerm(ih)

             do il=1,basis%NGr

                w=basis%G_Herm(ih)*basis%G_Lag(il)
                eta=basis%xLag(il)

!
!               Local densities
!
                rho_tot=densities(1)%rho(ih,il)+densities(2)%rho(ih,il)
                s_tot_z=densities(1)%s%z(ih,il)+densities(2)%s%z(ih,il)
                s_tot_rho=densities(1)%s%rho(ih,il)+densities(2)%s%rho(ih,il)
                s_tot_phi=densities(1)%s%phi(ih,il)+densities(2)%s%phi(ih,il)
                rot_j_tot_z=densities(1)%rot_j%z(ih,il)+densities(2)%rot_j%z(ih,il)
                rot_j_tot_rho=densities(1)%rot_j%rho(ih,il)+densities(2)%rot_j%rho(ih,il)
                rot_j_tot_phi=densities(1)%rot_j%phi(ih,il)+densities(2)%rot_j%phi(ih,il)
                t_tot_z=densities(1)%t%z(ih,il)+densities(2)%t%z(ih,il)
                t_tot_rho=densities(1)%t%rho(ih,il)+densities(2)%t%rho(ih,il)
                t_tot_phi=densities(1)%t%phi(ih,il)+densities(2)%t%phi(ih,il)
                j_tot_z=densities(1)%j%z(ih,il)+densities(2)%j%z(ih,il)
                j_tot_rho=densities(1)%j%rho(ih,il)+densities(2)%j%rho(ih,il)
                j_tot_phi=densities(1)%j%phi(ih,il)+densities(2)%j%phi(ih,il)
                rot_s_tot_z=densities(1)%rot_s%z(ih,il)+densities(2)%rot_s%z(ih,il)
                rot_s_tot_rho=densities(1)%rot_s%rho(ih,il)+densities(2)%rot_s%rho(ih,il)
                rot_s_tot_phi=densities(1)%rot_s%phi(ih,il)+densities(2)%rot_s%phi(ih,il)
                Delta_s_tot_z=densities(1)%Delta_s%z(ih,il)+densities(2)%Delta_s%z(ih,il)
                Delta_s_tot_rho=densities(1)%Delta_s%rho(ih,il)+densities(2)%Delta_s%rho(ih,il)
                Delta_s_tot_phi=densities(1)%Delta_s%phi(ih,il)+densities(2)%Delta_s%phi(ih,il)
                F_tot_z=densities(1)%F%z(ih,il)+densities(2)%F%z(ih,il)
                F_tot_rho=densities(1)%F%rho(ih,il)+densities(2)%F%rho(ih,il)
                F_tot_phi=densities(1)%F%phi(ih,il)+densities(2)%F%phi(ih,il)

!
!               HO wave functions and their derivatives
!
                phi1=sqrt(2.0_rk*bz)*bp*basis%P_Herm(nz1,ih)*basis%P_Lag(nr1,abs(l1),il)* &
                     (-1)**Delta1*sqrt(eta**abs(l1))*exp(-0.5_rk*(xi**2+eta))*chi1
                phi2=sqrt(2.0_rk*bz)*bp*basis%P_Herm(nz2,ih)*basis%P_Lag(nr2,abs(l2),il)* &
                     (-1)**Delta2*sqrt(eta**abs(l2))*exp(-0.5_rk*(xi**2+eta))*chi2
                sum_U=sum_U+w*HFfield%central(ih,il)*real(scalar_product(phi1,phi2))

                grad_phi1=grad(j1,ih,il,basis)
                grad_phi2=grad(j2,ih,il,basis)
                sum_K=sum_K+w*HFfield%effective_mass(ih,il)* &
                     real(scalar_product(grad_phi1,grad_phi2))

                sigma_grad_phi2=vector_product(sigma,grad_phi2)
                sum_SO=sum_SO+w*HFfield%spin_orbit(ih,il)* &
                     aimag(scalar_product(grad_phi1,sigma_grad_phi2))

!
!               Spin-field contributions
!
                S_z=HFfield%S%z(ih,il)
                S_rho=HFfield%S%rho(ih,il)
                S_phi=HFfield%S%phi(ih,il)
                sum_S_tot=sum_S_tot+w*real(scalar_product(phi1,S_z*(sigma%z*phi2))+ &
                     scalar_product(phi1,S_rho*(sigma%rho*phi2))+ &
                     scalar_product(phi1,S_phi*(sigma%phi*phi2)))

                S_z = B9*(rot_j_tot_z+densities(isospin)%rot_j%z(ih,il))
                S_rho = B9*(rot_j_tot_rho+densities(isospin)%rot_j%rho(ih,il))
                sum_S(9) = sum_S(9) + w*real(scalar_product(phi1,S_z*(sigma%z*phi2))+ &
                     scalar_product(phi1,S_rho*(sigma%rho*phi2)))

                S_z = 2.0_rk * B10 * s_tot_z
                S_rho = 2.0_rk * B10 * s_tot_rho
                sum_S(10) = sum_S(10) + w*real(scalar_product(phi1,S_z*(sigma%z*phi2))+ &
                     scalar_product(phi1,S_rho*(sigma%rho*phi2)))

                S_z = 2.0_rk * B11 * densities(isospin)%s%z(ih,il)
                S_rho = 2.0_rk * B11 * densities(isospin)%s%rho(ih,il)
                sum_S(11) = sum_S(11) + w*real(scalar_product(phi1,S_z*(sigma%z*phi2))+ &
                     scalar_product(phi1,S_rho*(sigma%rho*phi2)))

                S_z = 2.0_rk * B12 * rho_tot**alpha * s_tot_z
                S_rho = 2.0_rk * B12 * rho_tot**alpha * s_tot_rho
                sum_S(12) = sum_S(12) + w*real(scalar_product(phi1,S_z*(sigma%z*phi2))+ &
                     scalar_product(phi1,S_rho*(sigma%rho*phi2)))

                S_z = 2.0_rk * B13 * rho_tot**alpha * densities(isospin)%s%z(ih,il)
                S_rho = 2.0_rk * B13 * rho_tot**alpha * densities(isospin)%s%rho(ih,il)
                sum_S(13) = sum_S(13) + w*real(scalar_product(phi1,S_z*(sigma%z*phi2))+ &
                     scalar_product(phi1,S_rho*(sigma%rho*phi2)))

                S_z = -B14*t_tot_z
                S_rho = -B14*t_tot_rho
                sum_S(14) = sum_S(14) + w*real(scalar_product(phi1,S_z*(sigma%z*phi2))+ &
                     scalar_product(phi1,S_rho*(sigma%rho*phi2)))

                S_z = -B15*densities(isospin)%t%z(ih,il)
                S_rho = -B15*densities(isospin)%t%rho(ih,il)
                sum_S(15) = sum_S(15) + w*real(scalar_product(phi1,S_z*(sigma%z*phi2))+ &
                     scalar_product(phi1,S_rho*(sigma%rho*phi2)))

                S_z = -2*B16*densities(isospin)%F%z(ih,il)
                S_rho = -2*B16*densities(isospin)%F%rho(ih,il)
                sum_S(16) = sum_S(16) + w*real(scalar_product(phi1,S_z*(sigma%z*phi2))+ &
                     scalar_product(phi1,S_rho*(sigma%rho*phi2)))

                S_z = -2*B17*densities(isospin)%F%z(ih,il)
                S_rho = -2*B17*densities(isospin)%F%rho(ih,il)
                sum_S(17) = sum_S(17) + w*real(scalar_product(phi1,S_z*(sigma%z*phi2))+ &
                     scalar_product(phi1,S_rho*(sigma%rho*phi2)))

                S_z = 2*B18*Delta_s_tot_z
                S_rho = 2*B18*Delta_s_tot_rho
                sum_S(18) = sum_S(18) + w*real(scalar_product(phi1,S_z*(sigma%z*phi2))+ &
                     scalar_product(phi1,S_rho*(sigma%rho*phi2)))

                S_z = 2*B19*densities(isospin)%Delta_s%z(ih,il)
                S_rho = 2*B19*densities(isospin)%Delta_s%rho(ih,il)
                sum_S(19) = sum_S(19) + w*real(scalar_product(phi1,S_z*(sigma%z*phi2))+ &
                     scalar_product(phi1,S_rho*(sigma%rho*phi2)))

                ket_z=spin_operator_action(sigma%z,phi2)
                ket_rho=spin_operator_action(sigma%rho,phi2)
                ket_phi=spin_operator_action(sigma%phi,phi2)
                ket=spin_operator_action(sigma%z,grad_phi2%z) + &
                     spin_operator_action(sigma%rho,grad_phi2%rho) + &
                     spin_operator_action(sigma%phi,grad_phi2%phi)

                div_s = densities(1)%div_s(ih,il)+densities(2)%div_s(ih,il)
                sum_S(20) = sum_S(20) + w*2*B20*div_s*(real(scalar_product(grad_phi1%z,ket_z)) + &
                     real(scalar_product(grad_phi1%rho,ket_rho)) + &
                     real(scalar_product(grad_phi1%phi,ket_phi)) + &
                     real(scalar_product(phi1,ket)))

                div_s_q = densities(isospin)%div_s(ih,il)
                sum_S(21) = sum_S(21) + w*2*B21*div_s_q*(real(scalar_product(grad_phi1%z,ket_z)) + &
                     real(scalar_product(grad_phi1%rho,ket_rho)) + &
                     real(scalar_product(grad_phi1%phi,ket_phi)) + &
                     real(scalar_product(phi1,ket)))

                sum_S_tot=sum_S_tot+w*HFfield%div_s(ih,il)*(real(scalar_product(grad_phi1%z,ket_z)) + &
                     real(scalar_product(grad_phi1%rho,ket_rho)) + &
                     real(scalar_product(grad_phi1%phi,ket_phi)) + &
                     real(scalar_product(phi1,ket)))

                alpha_z=HFfield%alpha%z(ih,il)
                alpha_rho=HFfield%alpha%rho(ih,il)
                alpha_phi=HFfield%alpha%phi(ih,il)
                sum_A_tot=sum_A_tot+w*aimag(scalar_product(phi1,alpha_z*grad_phi2%z)+ &
                     scalar_product(phi1,alpha_rho*grad_phi2%rho)+ &
                     scalar_product(phi1,alpha_phi*grad_phi2%phi))

                alpha_z=-2.0_rk*B3*j_tot_z
                alpha_rho=-2.0_rk*B3*j_tot_rho
                alpha_phi=-2.0_rk*B3*j_tot_phi
                sum_A(3)=sum_A(3)+w*aimag(scalar_product(phi1,alpha_z*grad_phi2%z)+ &
                     scalar_product(phi1,alpha_rho*grad_phi2%rho) + &
                     scalar_product(phi1,alpha_phi*grad_phi2%phi))

                alpha_z=-2.0_rk*B4*densities(isospin)%j%z(ih,il)
                alpha_rho=-2.0_rk*B4*densities(isospin)%j%rho(ih,il)
                alpha_phi=-2.0_rk*B4*densities(isospin)%j%phi(ih,il)
                sum_A(4)=sum_A(4)+w*aimag(scalar_product(phi1,alpha_z*grad_phi2%z)+ &
                     scalar_product(phi1,alpha_rho*grad_phi2%rho) + &
                     scalar_product(phi1,alpha_phi*grad_phi2%phi))

                alpha_z=B9*(rot_s_tot_z+densities(isospin)%rot_s%z(ih,il))
                alpha_rho=B9*(rot_s_tot_rho+densities(isospin)%rot_s%rho(ih,il))
                alpha_phi=B9*(rot_s_tot_phi+densities(isospin)%rot_s%phi(ih,il))
                sum_A(9)=sum_A(9)+w*aimag(scalar_product(phi1,alpha_z*grad_phi2%z)+ &
                     scalar_product(phi1,alpha_rho*grad_phi2%rho) + &
                     scalar_product(phi1,alpha_phi*grad_phi2%phi))

                C_z=HFfield%C%z(ih,il)
                C_rho=HFfield%C%rho(ih,il)
                C_phi=HFfield%C%phi(ih,il)
                Csigma=C_z*sigma%z+C_rho*sigma%rho+C_phi*sigma%phi
                sigma_grad_phi2=spin_operator_action(Csigma,grad_phi2)
                sum_C_tot=sum_C_tot+w*real(scalar_product(grad_phi1,sigma_grad_phi2))

                C_z=-B14*s_tot_z
                C_rho=-B14*s_tot_rho
                C_phi=-B14*s_tot_phi
                Csigma=C_z*sigma%z+C_rho*sigma%rho+C_phi*sigma%phi
                sigma_grad_phi2=spin_operator_action(Csigma,grad_phi2)
                sum_C(14)=sum_C(14)+w*real(scalar_product(grad_phi1,sigma_grad_phi2))

                C_z=-B15*densities(isospin)%s%z(ih,il)
                C_rho=-B15*densities(isospin)%s%rho(ih,il)
                C_phi=-B15*densities(isospin)%s%phi(ih,il)
                Csigma=C_z*sigma%z+C_rho*sigma%rho+C_phi*sigma%phi
                sigma_grad_phi2=spin_operator_action(Csigma,grad_phi2)
                sum_C(15)=sum_C(15)+w*real(scalar_product(grad_phi1,sigma_grad_phi2))

                Jsigma_z=HFfield%tensor_so(1)%z(ih,il)*sigma%z + &
                     HFfield%tensor_so(1)%rho(ih,il)*sigma%rho + &
                     HFfield%tensor_so(1)%phi(ih,il)*sigma%phi
                Jsigma_rho=HFfield%tensor_so(2)%z(ih,il)*sigma%z + &
                     HFfield%tensor_so(2)%rho(ih,il)*sigma%rho + &
                     HFfield%tensor_so(2)%phi(ih,il)*sigma%phi
                Jsigma_phi=HFfield%tensor_so(3)%z(ih,il)*sigma%z + &
                     HFfield%tensor_so(3)%rho(ih,il)*sigma%rho + &
                     HFfield%tensor_so(3)%phi(ih,il)*sigma%phi
                ket_z=spin_operator_action(Jsigma_z,grad_phi2%z)
                ket_rho=spin_operator_action(Jsigma_rho,grad_phi2%rho)
                ket_phi=spin_operator_action(Jsigma_phi,grad_phi2%phi)
                sum_tensor_so=sum_tensor_so+w*(aimag(scalar_product(phi1,ket_z)) + &
                     aimag(scalar_product(phi1,ket_rho))+aimag(scalar_product(phi1,ket_phi)))
                ket_z=spin_operator_action(Jsigma_z,phi2)
                ket_rho=spin_operator_action(Jsigma_rho,phi2)
                ket_phi=spin_operator_action(Jsigma_phi,phi2)
                sum_tensor_so=sum_tensor_so-w*(aimag(scalar_product(grad_phi1%z,ket_z)) + &
                     aimag(scalar_product(grad_phi1%rho,ket_rho)) + &
                     aimag(scalar_product(grad_phi1%phi,ket_phi)))

                Jsigma_z = B14 * ( &
                     (densities(1)%tensor_J(1)%z(ih,il)+densities(2)%tensor_J(1)%z(ih,il))*sigma%z + &
                     (densities(1)%tensor_J(1)%rho(ih,il)+densities(2)%tensor_J(1)%rho(ih,il))*sigma%rho + &
                     (densities(1)%tensor_J(1)%phi(ih,il)+densities(2)%tensor_J(1)%phi(ih,il))*sigma%phi)
                Jsigma_rho = B14 * ( &
                     (densities(1)%tensor_J(2)%z(ih,il)+densities(2)%tensor_J(2)%z(ih,il))*sigma%z + &
                     (densities(1)%tensor_J(2)%rho(ih,il)+densities(2)%tensor_J(2)%rho(ih,il))*sigma%rho + &
                     (densities(1)%tensor_J(2)%phi(ih,il)+densities(2)%tensor_J(2)%phi(ih,il))*sigma%phi)
                Jsigma_phi = B14 * ( &
                     (densities(1)%tensor_J(3)%z(ih,il)+densities(2)%tensor_J(3)%z(ih,il))*sigma%z + &
                     (densities(1)%tensor_J(3)%rho(ih,il)+densities(2)%tensor_J(3)%rho(ih,il))*sigma%rho + &
                     (densities(1)%tensor_J(3)%phi(ih,il)+densities(2)%tensor_J(3)%phi(ih,il))*sigma%phi)
                ket_z=spin_operator_action(Jsigma_z,grad_phi2%z)
                ket_rho=spin_operator_action(Jsigma_rho,grad_phi2%rho)
                ket_phi=spin_operator_action(Jsigma_phi,grad_phi2%phi)
                sum_W(14) = sum_W(14) + w*(aimag(scalar_product(phi1,ket_z)) + &
                     aimag(scalar_product(phi1,ket_rho))+aimag(scalar_product(phi1,ket_phi)))
                ket_z=spin_operator_action(Jsigma_z,phi2)
                ket_rho=spin_operator_action(Jsigma_rho,phi2)
                ket_phi=spin_operator_action(Jsigma_phi,phi2)
                sum_W(14)=sum_W(14)-w*(aimag(scalar_product(grad_phi1%z,ket_z)) + &
                     aimag(scalar_product(grad_phi1%rho,ket_rho)) + &
                     aimag(scalar_product(grad_phi1%phi,ket_phi)))

                Jsigma_z = B15 * ( &
                     densities(isospin)%tensor_J(1)%z(ih,il)*sigma%z + &
                     densities(isospin)%tensor_J(1)%rho(ih,il)*sigma%rho + &
                     densities(isospin)%tensor_J(1)%phi(ih,il)*sigma%phi)
                Jsigma_rho = B15 * ( &
                     densities(isospin)%tensor_J(2)%z(ih,il)*sigma%z + &
                     densities(isospin)%tensor_J(2)%rho(ih,il)*sigma%rho + &
                     densities(isospin)%tensor_J(2)%phi(ih,il)*sigma%phi)
                Jsigma_phi = B15 * ( &
                     densities(isospin)%tensor_J(3)%z(ih,il)*sigma%z + &
                     densities(isospin)%tensor_J(3)%rho(ih,il)*sigma%rho + &
                     densities(isospin)%tensor_J(3)%phi(ih,il)*sigma%phi)
                ket_z=spin_operator_action(Jsigma_z,grad_phi2%z)
                ket_rho=spin_operator_action(Jsigma_rho,grad_phi2%rho)
                ket_phi=spin_operator_action(Jsigma_phi,grad_phi2%phi)
                sum_W(15) = sum_W(15) + w*(aimag(scalar_product(phi1,ket_z)) + &
                     aimag(scalar_product(phi1,ket_rho))+aimag(scalar_product(phi1,ket_phi)))
                ket_z=spin_operator_action(Jsigma_z,phi2)
                ket_rho=spin_operator_action(Jsigma_rho,phi2)
                ket_phi=spin_operator_action(Jsigma_phi,phi2)
                sum_W(15)=sum_W(15)-w*(aimag(scalar_product(grad_phi1%z,ket_z)) + &
                     aimag(scalar_product(grad_phi1%rho,ket_rho)) + &
                     aimag(scalar_product(grad_phi1%phi,ket_phi)))

                trace_J_tot = densities(1)%tensor_J(1)%z(ih,il)+densities(2)%tensor_J(1)%z(ih,il) + &
                     densities(1)%tensor_J(2)%rho(ih,il)+densities(2)%tensor_J(2)%rho(ih,il) + &
                     densities(1)%tensor_J(3)%phi(ih,il)+densities(2)%tensor_J(3)%phi(ih,il)
                Jsigma_z = B16 * ( &
                     (densities(1)%tensor_J(1)%z(ih,il)+densities(2)%tensor_J(1)%z(ih,il)+trace_J_tot)*sigma%z + &
                     (densities(1)%tensor_J(2)%z(ih,il)+densities(2)%tensor_J(2)%z(ih,il))*sigma%rho + &
                     (densities(1)%tensor_J(3)%z(ih,il)+densities(2)%tensor_J(3)%z(ih,il))*sigma%phi)
                Jsigma_rho = B16 * ( &
                     (densities(1)%tensor_J(1)%rho(ih,il)+densities(2)%tensor_J(1)%rho(ih,il))*sigma%z + &
                     (densities(1)%tensor_J(2)%rho(ih,il)+densities(2)%tensor_J(2)%rho(ih,il)+trace_J_tot)*sigma%rho + &
                     (densities(1)%tensor_J(3)%rho(ih,il)+densities(2)%tensor_J(3)%rho(ih,il))*sigma%phi)
                Jsigma_phi = B16 * ( &
                     (densities(1)%tensor_J(1)%phi(ih,il)+densities(2)%tensor_J(1)%phi(ih,il))*sigma%z + &
                     (densities(1)%tensor_J(2)%phi(ih,il)+densities(2)%tensor_J(2)%phi(ih,il))*sigma%rho + &
                     (densities(1)%tensor_J(3)%phi(ih,il)+densities(2)%tensor_J(3)%phi(ih,il)+trace_J_tot)*sigma%phi)
                ket_z=spin_operator_action(Jsigma_z,grad_phi2%z)
                ket_rho=spin_operator_action(Jsigma_rho,grad_phi2%rho)
                ket_phi=spin_operator_action(Jsigma_phi,grad_phi2%phi)
                sum_W(16) = sum_W(16) + w*(aimag(scalar_product(phi1,ket_z)) + &
                     aimag(scalar_product(phi1,ket_rho))+aimag(scalar_product(phi1,ket_phi)))
                ket_z=spin_operator_action(Jsigma_z,phi2)
                ket_rho=spin_operator_action(Jsigma_rho,phi2)
                ket_phi=spin_operator_action(Jsigma_phi,phi2)
                sum_W(16)=sum_W(16)-w*(aimag(scalar_product(grad_phi1%z,ket_z)) + &
                     aimag(scalar_product(grad_phi1%rho,ket_rho)) + &
                     aimag(scalar_product(grad_phi1%phi,ket_phi)))

                trace_J_q = densities(isospin)%tensor_J(1)%z(ih,il) + &
                     densities(isospin)%tensor_J(2)%rho(ih,il) + &
                     densities(isospin)%tensor_J(3)%phi(ih,il)
                Jsigma_z = B17 * ( &
                     (densities(isospin)%tensor_J(1)%z(ih,il)+trace_J_q)*sigma%z + &
                     densities(isospin)%tensor_J(2)%z(ih,il)*sigma%rho + &
                     densities(isospin)%tensor_J(3)%z(ih,il)*sigma%phi)
                Jsigma_rho = B17 * ( &
                     densities(isospin)%tensor_J(1)%rho(ih,il)*sigma%z + &
                     (densities(isospin)%tensor_J(2)%rho(ih,il)+trace_J_q)*sigma%rho + &
                     densities(isospin)%tensor_J(3)%rho(ih,il)*sigma%phi)
                Jsigma_phi = B17 * ( &
                     densities(isospin)%tensor_J(1)%phi(ih,il)*sigma%z + &
                     densities(isospin)%tensor_J(2)%phi(ih,il)*sigma%rho + &
                     (densities(isospin)%tensor_J(3)%phi(ih,il)+trace_J_q)*sigma%phi)
                ket_z=spin_operator_action(Jsigma_z,grad_phi2%z)
                ket_rho=spin_operator_action(Jsigma_rho,grad_phi2%rho)
                ket_phi=spin_operator_action(Jsigma_phi,grad_phi2%phi)
                sum_W(17) = sum_W(17) + w*(aimag(scalar_product(phi1,ket_z)) + &
                     aimag(scalar_product(phi1,ket_rho))+aimag(scalar_product(phi1,ket_phi)))
                ket_z=spin_operator_action(Jsigma_z,phi2)
                ket_rho=spin_operator_action(Jsigma_rho,phi2)
                ket_phi=spin_operator_action(Jsigma_phi,phi2)
                sum_W(17)=sum_W(17)-w*(aimag(scalar_product(grad_phi1%z,ket_z)) + &
                     aimag(scalar_product(grad_phi1%rho,ket_rho)) + &
                     aimag(scalar_product(grad_phi1%phi,ket_phi)))

                D_z=HFfield%D%z(ih,il)
                D_rho=HFfield%D%rho(ih,il)
                D_phi=HFfield%D%phi(ih,il)
                ket=D_z*grad_phi2%z+D_rho*grad_phi2%rho+D_phi*grad_phi2%phi
                bra=spin_operator_action(sigma%z,grad_phi1%z)+ &
                     spin_operator_action(sigma%rho,grad_phi1%rho)+ &
                     spin_operator_action(sigma%phi,grad_phi1%phi)
                sum_D_tot=sum_D_tot+0.5_rk*w*real(scalar_product(bra,ket))
                bra=D_z*grad_phi1%z+D_rho*grad_phi1%rho+D_phi*grad_phi1%phi ! D field is real
                ket=spin_operator_action(sigma%z,grad_phi2%z)+ &
                     spin_operator_action(sigma%rho,grad_phi2%rho)+ &
                     spin_operator_action(sigma%phi,grad_phi2%phi)
                sum_D_tot=sum_D_tot+0.5_rk*w*real(scalar_product(bra,ket))

                D_z = -2*B16*s_tot_z
                D_rho = -2*B16*s_tot_rho
                D_phi = -2*B16*s_tot_phi
                ket=D_z*grad_phi2%z+D_rho*grad_phi2%rho+D_phi*grad_phi2%phi
                bra=spin_operator_action(sigma%z,grad_phi1%z)+ &
                     spin_operator_action(sigma%rho,grad_phi1%rho)+ &
                     spin_operator_action(sigma%phi,grad_phi1%phi)
                sum_D(16)=sum_D(16)+0.5_rk*w*real(scalar_product(bra,ket))
                bra=D_z*grad_phi1%z+D_rho*grad_phi1%rho+D_phi*grad_phi1%phi ! D field is real
                ket=spin_operator_action(sigma%z,grad_phi2%z)+ &
                     spin_operator_action(sigma%rho,grad_phi2%rho)+ &
                     spin_operator_action(sigma%phi,grad_phi2%phi)
                sum_D(16)=sum_D(16)+0.5_rk*w*real(scalar_product(bra,ket))

                D_z = -2*B17*densities(isospin)%s%z(ih,il)
                D_rho = -2*B17*densities(isospin)%s%rho(ih,il)
                D_phi = -2*B17*densities(isospin)%s%phi(ih,il)
                ket=D_z*grad_phi2%z+D_rho*grad_phi2%rho+D_phi*grad_phi2%phi
                bra=spin_operator_action(sigma%z,grad_phi1%z)+ &
                     spin_operator_action(sigma%rho,grad_phi1%rho)+ &
                     spin_operator_action(sigma%phi,grad_phi1%phi)
                sum_D(17)=sum_D(17)+0.5_rk*w*real(scalar_product(bra,ket))
                bra=D_z*grad_phi1%z+D_rho*grad_phi1%rho+D_phi*grad_phi1%phi ! D field is real
                ket=spin_operator_action(sigma%z,grad_phi2%z)+ &
                     spin_operator_action(sigma%rho,grad_phi2%rho)+ &
                     spin_operator_action(sigma%phi,grad_phi2%phi)
                sum_D(17)=sum_D(17)+0.5_rk*w*real(scalar_product(bra,ket))

             end do

          end do

          h_HF%K(0) = h_HF%K(0) + C1 * C2 * sum_K
          h_HF%U(0) = h_HF%U(0) + C1 * C2 * sum_U
          h_HF%SO(0) = h_HF%SO(0) + C1 * C2 * sum_SO
          h_HF%S(0) = h_HF%S(0) + C1 * C2 * sum_S_tot
          h_HF%A(0) = h_HF%A(0) + C1 * C2 * sum_A_tot
          h_HF%C(0) = h_HF%C(0) + C1 * C2 * sum_C_tot
          do i = 1, 21
             h_HF%S(i) = h_HF%S(i) + C1 * C2 * sum_S(i)
             h_HF%A(i) = h_HF%A(i) + C1 * C2 * sum_A(i)
             h_HF%C(i) = h_HF%C(i) + C1 * C2 * sum_C(i)
             h_HF%W(i) = h_HF%W(i) + C1 * C2 * sum_W(i)
             h_HF%D(i) = h_HF%D(i) + C1 * C2 * sum_D(i)
          end do
          if (abs(h_HF%S(0) - sum(h_HF%S(1:21),1)) > 1e-1_rk) stop 'Inconsistency for S-field contribution'
          if (abs(h_HF%A(0) - sum(h_HF%A(1:21),1)) > 1e-1_rk) stop 'Inconsistency for A-field contribution'
          if (abs(h_HF%C(0) - sum(h_HF%C(1:21),1)) > 1e-1_rk) stop 'Inconsistency for C-field contribution'
          h_HF%W(0) = h_HF%W(0) + C1 * C2 * sum_tensor_so
          h_HF%D(0) = h_HF%D(0) + C1 * C2 * sum_D_tot
          if (abs(h_HF%W(0) - sum(h_HF%W(1:21),1)) > 1e-1_rk) stop 'Inconsistency for W-field contribution'
          if (abs(h_HF%D(0) - sum(h_HF%D(1:21),1)) > 1e-1_rk) stop 'Inconsistency for W-field contribution'

       end do loop_i2

    end do loop_i1

    do i = 0, 21
       h_HF%K(i) = norm * h_HF%K(i)
       h_HF%U(i) = norm * h_HF%U(i)
       h_HF%SO(i) = norm * h_HF%SO(i)
       h_HF%S(i) = norm * h_HF%S(i)
       h_HF%A(i) = norm * h_HF%A(i)
       h_HF%C(i) = norm * h_HF%C(i)
       h_HF%W(i) = norm * h_HF%W(i)
       h_HF%D(i) = norm * h_HF%D(i)
    end do

    h_HF%tot = h_HF%K(0) + h_HF%U(0) + h_HF%SO(0) + h_HF%S(0) + h_HF%A(0) + h_HF%C(0) + h_HF%D(0) + h_HF%W(0)

    return

  end subroutine h_HF_matrix_element


  !=====================================================================
  !                       SUBROUTINE sz_matrix_element
  !---------------------------------------------------------------------
  ! Calculates the matrix elements of the z-component of the spin 
  ! operator on the symmetry z-axis between two HF states (decomposed 
  ! in the HO basis).
  !=====================================================================

  subroutine sz_matrix_element(psi1, psi2, basis, sum_sz)

!
!   Arguments
!
    type(sp_state),intent(in) :: psi1, psi2
    type(cylindrical_basis),intent(in) :: basis
    real(kind=rk),intent(out) :: sum_sz

!
!   Local variables
!
    integer :: i,j,i1,i2,k,j1,j2,shift,n
    integer :: s1,s2,nz1,nz2,nr1,nr2
    real(kind=rk) :: C1,C2

!
!   Number of basis states preceding block # k
!

    k=psi1%block
    n=basis%block_size(k)
    shift=0
    do j=1,k
       shift=shift+basis%block_size(j)
    end do
    shift=shift-n

!
!   Matrix element
!

    sum_sz=0.0

    if (psi1%block /= psi2%block) return ! Axial symmetry
    if (psi1%tau /= psi2%tau) return ! Tz symmetry

    loop_i1: do i1=1,n

       j1=i1+shift
       s1=basis%ns(j1)
       nz1=basis%nz(j1)
       nr1=basis%nr(j1)
       C1 = psi1%coef%cylHO(j1)

       loop_i2: do i2=1,n

          j2=i2+shift
          s2=basis%ns(j2)
          nz2=basis%nz(j2)
          nr2=basis%nr(j2)
          C2 = psi2%coef%cylHO(j2)

          sum_sz = sum_sz + C1 * C2 * kronecker(nz1,nz2) * kronecker(nr1,nr2) * &
               kronecker(s1,s2) * s1

       end do loop_i2

    end do loop_i1

    sum_sz = 0.5_rk * sum_sz

    return

  end subroutine sz_matrix_element


  !=====================================================================
  !                       SUBROUTINE PAIRING_MATRIX_ELEMENTS
  !---------------------------------------------------------------------
  ! Calculates the matrix elements of the density-dependent delta 
  ! interaction (DDDI) or seniority force.
  !=====================================================================

  subroutine pairing_matrix_elements(sp_states, densities, isospin, &
       nucleus, model, basis, Vpair, partner, iter)

!
!   Arguments
!
    type(cylindrical_basis),intent(in)                         :: basis
    type(sp_state),dimension(basis%size),intent(in)            :: sp_states
    type(local_density),dimension(2),intent(in)                :: densities
    integer,intent(in)                                         :: isospin
    type(nucleus_properties),intent(in)                        :: nucleus
!    type(nucleus_properties),intent(inout)                     :: nucleus
    type(modelization),intent(in)                              :: model
    real(kind=rk),dimension(basis%size,basis%size),intent(out) :: Vpair
    integer,dimension(basis%size),intent(in)                   :: partner
    integer,optional,intent(in)                                :: iter

!
!   Local variables
!
    character(len=23),parameter :: subroutinename='pairing_matrix_elements'
    integer,dimension(2) :: Npart
    real(kind=rk) :: mu,center,width,fermi
    real(kind=rk),dimension(:,:,:),allocatable :: rho_sp
    integer :: iHF,jHF,icyl
    real(kind=rk) :: ei,fi,ej,fj
    integer :: iBlk,shift
    integer :: iHerm,iLag
    real(kind=rk) :: xi,eta
    real(kind=rk) :: sum_up,sum_down,C,power_eta
    integer :: nz,nr,Lambda
    real(kind=rk) :: V,Gz,rho_isoscalar
    type(sp_state),dimension(4) :: HF_states
    integer :: iHF_bar, jHF_bar

!
!   Debugging note
!
    if (debugmode=='high') call debugnote(modulename, subroutinename)

!
!   Auxiliary variables
!

    Npart(:)=(/ nucleus%N , nucleus%Z /)

!---------------------
!
!   Seniority force
!
!---------------------

    if (model%pairing=='seniority') then

       Vpair(:,:)=-abs(model%strength(isospin))/(11.0_rk+Npart(isospin))

!------------------------------------------
!
!   Density-dependent delta interaction
!
!------------------------------------------

    else if (model%pairing=='DDDI') then

!
!      Auxiliary variables
!

       mu=model%diffuseness
       center=model%truncmax-model%truncmin
       width=0.5_rk*(model%truncmax+model%truncmin)
       fermi=0.5_rk*(sp_states(Npart(isospin))%energy+ &
            sp_states(Npart(isospin)+1)%energy)

!
!      Single-particle densities: sum_s |phi_i(z,r,s)|^2
!     

       allocate(rho_sp(basis%size,iHmin:basis%NGz/2,basis%NGr))
       rho_sp = 0.0

       do iHF=1,basis%size

          if (sp_states(iHF)%Omega <=0) cycle

!
!         Speed up calculation by skipping small contributions
!
          ei=sp_states(iHF)%energy
          fi=f((abs(ei-fermi-0.5_rk*center)-width)/mu)
!!          if ((model%fit_pairing .ne. 'no') .and. (model%pairing .ne. 'seniority')) fi = 1.0_rk
          if ((model%fit_pairing .eq. 'no') .and. (abs(fi)<1e-6_rk)) then	!hock
!          if (abs(fi)<1e-6_rk) then		!original
             rho_sp(iHF,:,:)=0.0
             cycle
          end if

          iBlk=sp_states(iHF)%block
          shift=merge(0,sum(basis%block_size(1:iBlk-1)),iBlk==1)
!          shift=0
!          do k=1,iBlk
!             shift=shift+basis%block_size(k)
!          end do
!          shift=shift-basis%block_size(iBlk)

          do iHerm=iHmin,basis%NGz/2

             if (iHerm==0) cycle
             xi=basis%xHerm(iHerm)

             do iLag=1,basis%NGr

                eta=basis%xLag(iLag)

                sum_up=0.
                sum_down=0.
                do icyl=1+shift,basis%block_size(iBlk)+shift
                   C=sp_states(iHF)%coef%cylHO(icyl)
                   if (abs(C)<1e-06_rk) cycle
                   nz=basis%nz(icyl)
                   nr=basis%nr(icyl)
                   Lambda=basis%Lambda(icyl)
                   if (Lambda==0) then
                      power_eta=1.0_rk
                   else if (eta>=0.0) then
                      power_eta=sqrt(eta**Lambda)
                   else
                      call fatal(modulename, subroutinename, 'eta < 0')
                   end if
                   if (basis%ns(icyl)==1) then
                      sum_up=sum_up+C*power_eta* &
                           basis%P_Herm(nz,iHerm)*basis%P_Lag(nr,Lambda,iLag)
                   else
                      sum_down=sum_down+C*power_eta* &
                           basis%P_Herm(nz,iHerm)*basis%P_Lag(nr,Lambda,iLag)
                   end if
                end do

                rho_sp(iHF,iHerm,iLag)=(sum_up**2+sum_down**2)*&
                     bz*bp2*exp(-(eta+xi**2))/pi

                if (abs(rho_sp(iHF,iHerm,iLag))<1e-14_rk) then
                   rho_sp(iHF,iHerm,iLag)=0.0
                else if (debugmode=='high'.and.present(iter)) then
                   write(88,'(5(i4,1x),e15.8)') iter,isospin,iHF,iHerm,iLag, &
                        rho_sp(iHF,iHerm,iLag)
                end if

             end do

          end do

       end do

       if (debugmode=='high') write(out_unit,'(2x,a/)') 'sp densities: [OK]'

!---------------------------
!
!      Matrix elements
!
!---------------------------

       Vpair(:,:)=0.0

       do iHF=1,basis%size

          if (sp_states(iHF)%Omega < 0) cycle

!
!         Speed up calculation by skipping small contributions
!
          ei=sp_states(iHF)%energy
          fi=f((abs(ei-fermi-0.5_rk*center)-width)/mu)
!!          if ((model%fit_pairing .ne. 'no') .and. (model%pairing .ne. 'seniority')) fi = 1.0_rk
          if ((model%fit_pairing .eq. 'no') .and. (abs(fi)<1e-6_rk)) then	!hock
!          if (abs(fi)<1e-6_rk) then		!original
             Vpair(iHF,:)=0.0
             cycle
          end if

          do jHF=iHF,basis%size

             if (sp_states(jHF)%Omega < 0) cycle

             ej=sp_states(jHF)%energy
             fj=f((abs(ej-fermi-0.5_rk*center)-width)/mu)
!!             if ((model%fit_pairing .ne. 'no') .and. (model%pairing .ne. 'seniority')) fj = 1.0_rk
!             if ((model%fit_pairing .eq. 'no') .and. (abs(fj)<1e-6_rk)) cycle
             if (abs(fj)<1e-6_rk) cycle	!hock

             V=0.0
             do iHerm=iHmin,basis%NGz/2
                if (iHerm==0) cycle
                Gz=basis%G_Herm(iHerm)
                do iLag=1,basis%NGr
                   rho_isoscalar= & ! zero at 1st iter.
                        densities(1)%rho(iHerm,iLag)+densities(2)%rho(iHerm,iLag)
                   V=V+Gz*basis%G_Lag(iLag)* &
                        rho_sp(iHF,iHerm,iLag)*rho_sp(jHF,iHerm,iLag)* &
                        (1.0_rk-model%eta*(rho_isoscalar/model%rho_c)**model%beta)

!                   psi1 = sp_states(iHF)%wf(iHerm,iLag)
!                   psi2 = sp_states(partner(iHF))%wf(iHerm,iLag)
!                   psi3 = sp_states(jHF)%wf(iHerm,iLag)
!                   psi4 = sp_states(partner(jHF))%wf(iHerm,iLag)

!                   V=V+Gz*basis%G_Lag(iLag)*(1.0_rk-model%eta*(rho_isoscalar/model%rho_c)**model%beta)* &
!                        real(scalar_product(psi1,psi2)*scalar_product(psi3,psi4))

                end do
             end do
             V=fz*V


             Vpair(iHF,jHF)=model%strength(isospin)*V*pi/(bz*bp2)
             Vpair(jHF,iHF)=Vpair(iHF,jHF)                        ! Vpair is a symmetric matrix
!             Vpair(partner(iHF),partner(jHF))=Vpair(iHF,jHF)      ! Vdelta is time-reversal invariant
!             Vpair(partner(jHF),partner(iHF))=Vpair(iHF,jHF)      ! Vpair is a symmetric matrix

             if (debugmode=='high' .and. abs(Vpair(iHF,jHF))>1e-6_rk) &
                  write(99,'(3(i4,1x),1(e20.12,1x))') isospin,iHF,jHF,&
                  Vpair(iHF,jHF)

          end do

       end do

       deallocate(rho_sp)    

       if (debugmode=='high') write(out_unit,'(2x,a)') 'matrix elements: [OK]'

    else if (model%pairing=='PSG') then

       Vpair(:,:)=0.0
!       write(out_unit,*) 'a_psg=',a_psg,' psg_strength=',psg_strength
!       write(out_unit,*) 'bz=',bz,' bp=',bp

       do iHF=1,basis%size

          if (sp_states(iHF)%Omega < 0) cycle

          HF_states(1) = sp_states(iHF)
          iHF_bar = partner(iHF)
          if (iHF_bar < 0 .or. iHF_bar > basis%size) call fatal(modulename, subroutinename, &
               'Error with partner(iHF)')
          HF_states(2) = sp_states(iHF_bar)

          do jHF=iHF,basis%size

             if (sp_states(jHF)%Omega < 0) cycle

             HF_states(3) = sp_states(jHF)
             jHF_bar = partner(jHF)
             if (jHF_bar < 0 .or. jHF_bar > basis%size) call fatal(modulename, subroutinename, &
                  'Error with partner(jHF)')
             HF_states(4) = sp_states(jHF_bar)

             Vpair(iHF,jHF) = HF_PSG_antisym(HF_states,basis) ! =HF_PSG_antisym(HF_states,basis)
!             Vpair(iHF,jHF) = HF_DELTA(HF_states,-abs(model%strength(:)),basis)
             Vpair(jHF,iHF) = Vpair(iHF,jHF)
!             v = HF_delta(HF_states,(/0.0_rk,-psg_strength/),basis)
!             xi= -HF_gauss(HF_states,basis,0,1)*psg_strength
!             if (abs(V)>1e-2_rk) write(out_unit,'(2(i4,1x),1(f10.6,1x))') &
!                  iHF,jHF,Vpair(iHF,jHF)/v!,xi/v

             if (debugmode=='high' .and. abs(Vpair(iHF,jHF))>1e-6_rk) &
                  write(99,'(3(i4,1x),2(e20.12,1x))') isospin,iHF,jHF,&
                  0.5_rk*Vpair(iHF,jHF),HF_PSG_antisym(HF_states,basis)/HF_DELTA(HF_states, -abs(model%strength(:)), basis)

          end do

       end do

    else

       call fatal(modulename, subroutinename, 'Wrong name for pairing model')

    end if

    return

  end subroutine pairing_matrix_elements


  !=====================================================================
  !                       SUBROUTINE BCS
  !---------------------------------------------------------------------
  ! Solves the BCS equations for the single-particle spectrum given in 
  ! input (through argument 1). Returns the state-dependent pairing gaps, 
  ! the chemical potential, the single-particle states occupation 
  ! probabilities and the pairing energy.
  !=====================================================================

  subroutine BCS(nucleus, model, basis, sp_states, Vpair, isospin, i_odd)

!
!   Arguments
!
    type(nucleus_properties),intent(inout)                    :: nucleus
    type(modelization),intent(in)                             :: model
    type(cylindrical_basis),intent(in)                        :: basis
    type(sp_state),dimension(basis%size),intent(inout)        :: sp_states
    real(kind=rk),dimension(basis%size,basis%size),intent(in) :: Vpair
    integer,intent(in)                                        :: isospin
    integer,dimension(nqp),intent(in)                           :: i_odd

!
!   Local variables
!
    character(len=3),parameter :: subroutinename='BCS'
    integer,dimension(2) :: Npart
    real(kind=rk) :: mu,center,width,fermi,lambda,e_1st_unoccupied,Vpair_max
    integer :: i,j,n,iHF,jHF,kHF,Nblocked
    integer :: iterBCS
    integer,parameter :: iterBCSmax=500
    real(kind=rk),parameter :: mix=0.75_rk , tol=1e-6_rk
    real(kind=rk) :: sum_gap,sum_Delta,sum_e,denom,sum_1
    real(kind=rk) :: ei,fi,ej,fj,Delta_i
    real(kind=rk),dimension(basis%size) :: Delta,Delta_previous
    real(kind=rk) :: lambda_previous
    real(kind=rk) :: Epair,avg_gap
    real(kind=rk) :: Eqp,Eqpmin,v2,uv,tmp
    integer :: iqpmin
    integer,dimension(nqp) :: j_odd

!
!   Debugging note
!

    if (debugmode=='high') call debugnote(modulename, subroutinename)

!
!   Auxiliary variables
!

    Npart(:)=(/ nucleus%N , nucleus%Z /)
    mu=model%diffuseness
    if (mu<1e-5_rk) then
       call warning(modulename,subroutinename,'mu < 1e-5_rk ==> mu reset to 0.001')
       mu=1e-3_rk
    end if
    center=model%truncmax-model%truncmin
    width=0.5_rk*(model%truncmax+model%truncmin)

!
!   Initializations
!   * chemical potential to Fermi level
!   * pairing gaps to arbitrary values
!   * occupancies to the Hatree-Fock values
!   * pairing energy to 0
!

    n=Npart(isospin)-mod(Npart(isospin),2)
    e_1st_unoccupied=sp_states(n+1)%energy
    if (n/=Npart(isospin)) then
       fermi=0.5_rk*(sp_states(n)%energy+e_1st_unoccupied)
    else
       fermi = e_1st_unoccupied
    end if
!    fermi = sp_states(Npart(isospin))%energy
    lambda=fermi
    Delta(:)=1.0_rk
    Delta_previous(:)=0.0
    nucleus%pairing(isospin)%pairing_energy=0.0
    nucleus%pairing(isospin)%chemical_potential=lambda
    nucleus%pairing(isospin)%average_gap=0.0
    nucleus%pairing(isospin)%gap(:)=0.0
    nucleus%pairing(isospin)%min_qp_index=0
    nucleus%pairing(isospin)%min_qp_energy=0

    j_odd(:)=0
    do i=1,nqp
       if (i_odd(i)/=0) j_odd(i)=sp_states(i_odd(i))%pair_partner(isospin)
    end do

!
!   Check if BCS calculation is needed:
!   yes if max{|Vpair(i,j)|,i,j=1..basis%size} is larger than 1 keV
!   no otherwise
!

!     Vpair_max = maxval(abs(Vpair))
!     do iHF=1,basis%size
!        if (sp_states(iHF)%Omega < 0) cycle
!        do jHF=iHF,basis%size
!           if (sp_states(jHF)%Omega < 0) cycle
!           if (abs(Vpair(iHF,jHF) > 1e-3_rk) &
!                write(54,'(3(a,i4,1x),a,f12.6)') &
!                'isospin =',isospin,'iHF =',iHF,'jHF =',jHF, &
!                '<iHF T(iHF)|Vpair(1-P12)|jHF T(jHF)> =',Vpair(iHF,jHF)
! !          if (Vpair_max<abs(Vpair(iHF,jHF))) Vpair_max=abs(Vpair(iHF,jHF))
!        end do
!     end do
!     if (abs(Vpair_max)<=1e-3_rk) return

!---------------------------
!
!   Solving BCS equations
!
!---------------------------

    do iterBCS=1,iterBCSmax

!
!      Pairing gaps
!

       loop_iHF_1: do iHF=1,basis%size

          do i=1,nqp
             if (iHF == i_odd(i) .or. iHF == j_odd(i)) cycle loop_iHF_1
          end do
          if (sp_states(iHF)%Omega<0) cycle
          jHF=sp_states(iHF)%pair_partner(isospin)
          ei=0.5_rk*(sp_states(iHF)%energy+sp_states(jHF)%energy)
          fi=f((abs(ei-lambda-0.5_rk*center)-width)/mu)
!!          if ((model%fit_pairing .ne. 'no') .and. (model%pairing .ne. 'seniority')) fi = 1.0_rk	!hock
          if (ei>lambda+model%truncmax.and.abs(fi)<1e-5_rk) exit	!original
!          if ((model%fit_pairing .eq. 'no') .and. (ei>lambda+model%truncmax.and.abs(fi)<1e-5_rk)) exit

          sum_gap=0.
          loop_jHF: do jHF=1,basis%size
             do j=1,nqp
                if (jHF == i_odd(j) .or. jHF == j_odd(j)) cycle loop_jHF
             end do
             if (sp_states(jHF)%Omega<0) cycle
             kHF=sp_states(jHF)%pair_partner(isospin)
             ej=0.5_rk*(sp_states(jHF)%energy+sp_states(kHF)%energy)
             fj=f((abs(ej-lambda-0.5_rk*center)-width)/mu)
!!             if ((model%fit_pairing .ne. 'no') .and. (model%pairing .ne. 'seniority')) fj = 1.0_rk
             if (ej>lambda+model%truncmax.and.abs(fj)<1e-5_rk) exit	!original
!             if ((model%fit_pairing .eq. 'no') .and. (ej>lambda+model%truncmax.and.abs(fj)<1e-5_rk)) exit

             Delta_i=Delta(jHF)
             Eqp=sqrt((ej-lambda)**2+fj*Delta_i**2)
             sum_gap=sum_gap+Vpair(iHF,jHF)*fj*Delta_i/Eqp
          end do loop_jHF

          Delta_previous(iHF)=Delta(iHF)
!          if (sum_gap>0.0) call fatal(modulename,subroutinename,'sum_gap should be <0')
          Delta(iHF)=-mix*0.5_rk*sum_gap + (1.0_rk-mix)*Delta_previous(iHF)

       end do loop_iHF_1

!
!      Chemical potential
!

       sum_Delta=0.
       sum_e=0.
       denom=0.
       sum_1=0.
       loop_iHF_2: do iHF=1,basis%size
          do i=1,nqp
             if (iHF == i_odd(i) .or. iHF == j_odd(i)) cycle loop_iHF_2
          end do
          if (sp_states(iHF)%Omega<0) cycle
          jHF=sp_states(iHF)%pair_partner(isospin)
          ei=0.5_rk*(sp_states(iHF)%energy+sp_states(jHF)%energy)
          fi=f((abs(ei-lambda-0.5_rk*center)-width)/mu)
!!          if ((model%fit_pairing .ne. 'no') .and. (model%pairing .ne. 'seniority')) fi = 1.0_rk
          if (ei>lambda+model%truncmax.and.abs(fi)<1e-5_rk) exit	!original
!          if ((model%fit_pairing .eq. 'no') .and. (ei>lambda+model%truncmax.and.abs(fi)<1e-5_rk)) exit
          Delta_i=Delta(iHF)
          Eqp=sqrt((ei-lambda)**2+fi*Delta_i**2)   
          denom=denom+1.0_rk/Eqp
          sum_e=sum_e+(1.0_rk-ei/Eqp)
          sum_Delta=sum_Delta+abs(Delta_i-Delta_previous(iHF))
          sum_1=sum_1+1.0_rk
       end do loop_iHF_2

       if (denom == 0.0) call fatal(modulename, subroutinename, 'denom=0')
       lambda_previous=lambda 
!       lambda=mix*(Npart(isospin)-3+kronecker(i_odd(1),0)+kronecker(i_odd(2),0)+kronecker(i_odd(3),0)-sum_e)/denom &
!            + (1.0_rk-mix)*lambda_previous 
       Nblocked = 0
       do i=1,nqp
          Nblocked = Nblocked + kronecker(i_odd(i),0)
       end do
       Nblocked = nqp - Nblocked
       lambda = mix * (Npart(isospin) - Nblocked - sum_e)/denom + (1.0_rk-mix)*lambda_previous 

!
!      Convergence control
!

       if (abs(lambda-lambda_previous)+sum_Delta < (1.0_rk+sum_1)*tol) exit

    end do

    if (iterBCS==iterBCSmax) call warning(modulename, subroutinename, &
         'BCS solution not converged')

!-----------------------------------------------------------------
!
!   Ocupation numbers, pairing energy and lowest quasi-particle
!
!-----------------------------------------------------------------

    Epair=0.
    avg_gap=0.0
    denom=0.0
    Eqpmin=1e+3_rk
    sum_gap=0.0
    tmp=0.0
    loop_iHF_3: do iHF=1,basis%size
       do i=1,nqp
          if (iHF == i_odd(i) .or. iHF == j_odd(i)) cycle loop_iHF_3
       end do
       if (sp_states(iHF)%Omega<0) cycle
       jHF=sp_states(iHF)%pair_partner(isospin)
       ei=0.5_rk*(sp_states(iHF)%energy+sp_states(jHF)%energy)
       fi=f((abs(ei-lambda-0.5_rk*center)-width)/mu)
!!       if ((model%fit_pairing .ne. 'no') .and. (model%pairing .ne. 'seniority')) fi = 1.0_rk
!       if (ei>lambda+model%truncmax.and.abs(fi)<1e-5_rk) Delta(iHF)=0.0	!original
       if ((model%fit_pairing .eq. 'no') .and. (ei>lambda+model%truncmax.and.abs(fi)<1e-5_rk)) Delta(iHF)=0.0	!hock

       Delta_i=Delta(iHF)
       Eqp=sqrt((ei-lambda)**2+fi*Delta_i**2)
       if (abs(fi*Delta_i**2)<=1e-6) then
          v2=0.5_rk*(1.0_rk-sign(1.0_rk,ei-lambda))
       else
          v2=0.5_rk*(1.0_rk-(ei-lambda)/Eqp)
       end if
       if (v2>1.0_rk .or. v2<0.0) call fatal(modulename, subroutinename, &
            'v2>1 or v2<0')
       sp_states(iHF)%occupation=v2
       sp_states(jHF)%occupation=v2
       Delta(jHF)=Delta_i
       Epair=Epair-0.5_rk*fi*Delta_i**2/Eqp
       uv=sqrt((1.0_rk-v2)*v2)
       avg_gap=avg_gap+Delta_i*uv
       denom=denom+uv
       if (Eqp<Eqpmin) then
          iqpmin=iHF
          Eqpmin=Eqp
       end if
       sum_gap=sum_gap+v2
       tmp=tmp+v2*2
    end do loop_iHF_3

    if (abs(denom)<=1e-6_rk) then
       avg_gap=0.0
    else
       avg_gap=avg_gap/denom
    end if

!
!   Updating components of argument nucleus
!

    nucleus%pairing(isospin)%chemical_potential=lambda
    nucleus%pairing(isospin)%average_gap=avg_gap
    nucleus%pairing(isospin)%gap(:)=Delta(:)
    nucleus%pairing(isospin)%Vpair(:,:) = Vpair(:,:)

    do i=1,nqp
       if (i_odd(i)/=0) then
          nucleus%pairing(isospin)%gap(i_odd(i))=0.0
          if (j_odd(i)/=0) nucleus%pairing(isospin)%gap(j_odd(i))=0.0
       end if
    end do
    nucleus%pairing(isospin)%pairing_energy=Epair
    nucleus%pairing(isospin)%min_qp_index=iqpmin
    nucleus%pairing(isospin)%min_qp_energy=Eqpmin

    do i =1, basis%size
       write(264,*) isospin, sp_states(i)%energy, nucleus%pairing(isospin)%gap(i)
    end do
    write(264,*)
    return

  end subroutine BCS


  !=====================================================================
  !                       FUNCTION F
  !---------------------------------------------------------------------
  ! Smooth cut-off function to limit pairing action to a window around
  ! Fermi level.
  !=====================================================================

  function f(x)

    implicit none
    real(kind=rk) :: x,f

    if (x > 50.0_rk) then
       f = 0.0
    else if (x < -50.0_rk) then
       f = 1.0_rk
    else
       f = 1.0_rk/sqrt(1.0_rk+exp(x))
    end if

    return

  end function f


  !=====================================================================!
  !                     SUBROUTINE SP_WAVE_FUNCTIONS                    !
  !---------------------------------------------------------------------!
  ! Evaluates a single-particle wave function and its gradient at       !
  ! Gauss--Hermite and Gauss--Laguerre integration points.              !
  !                                                                     !
  !   sp_states(iHF)%wf(ih,il)%up = f_iHF^+(z,rho)                      !
  !   sp_states(iHF)%wf(ih,il)%down = f_iHF^-(z,rho)                    !
  !                                                                     !
  !=====================================================================!

  subroutine sp_wave_functions(psi, basis) 

    implicit none
!
!   Arguments
!
    type(sp_state),intent(inout) :: psi
!    type(vector_spinor_function),intent(out) :: grad_psi
    type(cylindrical_basis),intent(in) :: basis

!
!   Local variables
!
    integer :: iblock,shift,ih,il,i,nzi,nri,li,nsi,Delta
    real(kind=rk) :: xi,eta,sum_fp,sum_fm,tmp,norm

!
!   Calculation of the single-particle wave function at the Gauss--Hermite
!   and Gauss--Laguerre integration points
!
    psi%wf(:,:)=spinor(cmplx(0),cmplx(0))
    norm=sqrt(2.0_rk*bz)*bp

    iblock=psi%block
    shift=0
    if (iblock>1) shift=sum(basis%block_size(1:iblock-1))

    do ih=iHmin,basis%NGz/2
       if (ih==0) cycle
       xi=basis%xHerm(ih)
       do il=1,basis%NGr
          eta=basis%xLag(il)
          sum_fp=0.0_rk
          sum_fm=0.0_rk
          do i=1+shift,basis%block_size(iblock)+shift
             nzi=basis%nz(i)
             nri=basis%nr(i)
             li=basis%Lambda(i)
             nsi=basis%ns(i)
             Delta=(li-abs(li))/2
             tmp=psi%coef%cylHO(i)*basis%P_Herm(nzi,ih)*basis%P_Lag(nri,abs(li),il)* &
                  exp(-0.5_rk*(xi**2+eta))*(-1)**Delta
             if (li/=0) tmp=sqrt(eta**abs(li))*tmp
             if (nsi==1) then
                sum_fp=sum_fp+tmp
             else
                sum_fm=sum_fm+tmp
             end if
          end do
          psi%wf(ih,il)=spinor(cmplx(sum_fp*norm),cmplx(sum_fm*norm))
       end do
    end do

    psi%avg_pi = 0
    do i = 1, basis%size
       if (2*basis%Lambda(i) + basis%ns(i) /= psi%Omega) cycle
       psi%avg_pi = psi%avg_pi + psi%coef%cylHO(i)**2 * (-1)**(basis%nz(i)+basis%Lambda(i))
    end do

!
!   Calculation of the gradient of the single-particle wave function
!   at the Gauss--Hermite and Gauss--Laguerre integration points
!
!    if (associated(grad_psi%z)) nullify(grad_psi%z)
!    if (associated(grad_psi%rho)) nullify(grad_psi%rho)
!    if (associated(grad_psi%phi)) nullify(grad_psi%phi)

!    allocate(grad_psi%z(iHmin:basis%NGz/2,basis%NGr),grad_psi%rho(iHmin:basis%NGz/2,basis%NGr), &
!         grad_psi%phi(iHmin:basis%NGz/2,basis%NGr))

!    grad_psi%z(:,:)=spinor(cmplx(0),cmplx(0))
!    grad_psi%rho(:,:)=spinor(cmplx(0),cmplx(0))
!    grad_psi%phi(:,:)=spinor(cmplx(0),cmplx(0))

!    do ih=iHmin,basis%NGz/2
!       if (ih==0) cycle
!       do il=1,basis%NGr
!          g=grad(psi%coef,ih,il,basis)
!          grad_psi%z(ih,il)=g%z
!          grad_psi%rho(ih,il)=g%rho
!          grad_psi%phi(ih,il)=g%phi
!       end do
!    end do

    return

  end subroutine sp_wave_functions


  !=====================================================================!
  !                 SUBROUTINE LOCAL_ONE_BODY_DENSITIES                 !
  !---------------------------------------------------------------------!
  ! Calculates local one-body densities rho (nucleon density), Delta_rho! 
  ! (Laplacian of rho), tau (kinetic density), div_J (divergence of the !
  ! spin-orbit current J), j (current density), div_j (divergence of j),!
  ! rot_j (curl of j), s (spin density), div_s (divergence of  s) and   !
  ! rot_s (curl of s) at Gauss--Hermite--Laguerre integration points.   !
  !=====================================================================!

  subroutine local_one_body_densities(sp_states, basis, densities, isospin, iter)

    implicit none
!
!   Arguments
!

    type(cylindrical_basis),intent(in)     :: basis
    type(sp_state),dimension(basis%size),intent(in) :: sp_states
    type(local_density),intent(inout)      :: densities
    integer,optional,intent(in)            :: isospin,iter
!
!   Local variables
!
    character(len=24),parameter :: subroutinename='local_one_body_densities'
    integer :: nthreads
    integer :: i,ih,il,iHF
    real(kind=rk) :: rho
    complex(kind=rk) :: sum_rho,sum_tau,sum_Delta_rho,sum_div_J,sum_div_s
    type(vector) :: sum_j,sum_rot_j,sum_s,sum_Delta_s,sum_rot_s,sum_T,sum_F, &
         psi_grad_psi,psi_sigma_psi,psi_sigma_grad_psi,Delta_psi_sigma_psi, &
         grad_psi_sigma_grad_psi,psi_sigma_Delta_psi
    type(spinor) :: psi,Delta_psi,sigma_z_Delta_psi,sigma_rho_Delta_psi, &
         sigma_phi_Delta_psi,grad_phi_sigma_rho_psi,grad_phi_sigma_phi_psi, &
         sigma_dot_grad_psi
    type(vector_spinor) :: grad_psi,sigma_grad_psi,sigma_psi, &
         sigma_z_grad_psi,sigma_rho_grad_psi,sigma_phi_grad_psi
    type(vector_spin_operator) :: sigma
    type(vector),dimension(3) :: sum_tensor_J
    real(kind=rk) :: w, tmp_rho, tmp_tau, tmp_sz, tmp_jphi
    integer :: k, shift
    real(kind=rk) :: c
    type(vector_spinor) :: grad_phiHO

!
!   Debugging note
!
    if (debugmode=='high') call debugnote(modulename, subroutinename)

!
!   Definition of the cylindrical Pauli matrices
!
    do i=1,2
       sigma%z(i,i)=-cmplx((-1)**i)
       sigma%rho(i,i)=cmplx(0)
       sigma%phi(i,i)=cmplx(0)
    end do

    sigma%z(1,2)=cmplx(0)
    sigma%rho(1,2)=cmplx(1) !*exp(-i * phi) with phi = 0 (axial symmetry)
    sigma%phi(1,2)=-i_cplx

    sigma%z(2,1)=cmplx(0)
    sigma%rho(2,1)=cmplx(1)
    sigma%phi(2,1)=i_cplx

!
!   Filling-in local densities at mesh points (ih,il)
!

!    int_J2 = 0.0

    nthreads = basis%NGz/2
    call setup_openmp_parameters(iHmin,basis%NGz/2, nthreads, OMP_parameters%set_bounds, &
            OMP_parameters%nthreads)
    call omp_set_num_threads(OMP_parameters%nthreads)

    !$OMP PARALLEL DO &
    !$OMP PRIVATE(ih,il,sum_rho,sum_tau,sum_Delta_rho,sum_div_J,sum_j) &
    !$OMP PRIVATE(sum_rot_j,sum_s,sum_rot_s,sum_T,sum_Delta_s) &
    !$OMP PRIVATE(sum_tensor_J,sum_F,sum_div_s,rho) &
    !$OMP PRIVATE(iHF,psi,Delta_psi,grad_psi,sigma_grad_psi,psi_grad_psi) &
    !$OMP PRIVATE(sigma_psi,psi_sigma_psi,psi_sigma_grad_psi) &
    !$OMP PRIVATE(sigma_z_grad_psi,sigma_rho_grad_psi,sigma_phi_grad_psi) &
    !$OMP PRIVATE(Delta_psi_sigma_psi,grad_psi_sigma_grad_psi) &
    !$OMP PRIVATE(sigma_z_Delta_psi,sigma_rho_Delta_psi,sigma_phi_Delta_psi) &
    !$OMP PRIVATE(psi_sigma_Delta_psi,grad_phi_sigma_rho_psi,grad_phi_sigma_phi_psi) &
    !$OMP PRIVATE(sigma_dot_grad_psi)

    loop_ih: do ih=iHmin,basis%NGz/2

       if (ih==0) cycle

       loop_il: do il=1,basis%NGr

          sum_rho=cmplx(0)
          sum_tau=cmplx(0)
          sum_Delta_rho=cmplx(0)
          sum_div_J=cmplx(0)
          sum_j=vector(cmplx(0),cmplx(0),cmplx(0))
          sum_rot_j=vector(cmplx(0),cmplx(0),cmplx(0))
          sum_s=vector(cmplx(0),cmplx(0),cmplx(0))
          sum_rot_s=vector(cmplx(0),cmplx(0),cmplx(0))
          sum_T=vector(cmplx(0),cmplx(0),cmplx(0))
          sum_Delta_s=vector(cmplx(0),cmplx(0),cmplx(0))
          do i=1,3
             sum_tensor_J(i)=vector(cmplx(0),cmplx(0),cmplx(0))
          end do
          sum_F=vector(cmplx(0),cmplx(0),cmplx(0))
          sum_div_s=cmplx(0)
          rho=sqrt(basis%xLag(il))/bp

          loop_iHF: do iHF=1,basis%size

             psi=sp_states(iHF)%wf(ih,il)
             sum_rho=sum_rho+sp_states(iHF)%occupation * scalar_product(psi,psi)

             Delta_psi=Laplacian(sp_states(iHF)%coef,ih,il,basis)
             sum_Delta_rho=sum_Delta_rho+sp_states(iHF)%occupation * scalar_product(psi,Delta_psi)

             grad_psi%z=spinor(cmplx(0),cmplx(0))
             grad_psi%rho=spinor(cmplx(0),cmplx(0))
             grad_psi%phi=spinor(cmplx(0),cmplx(0))
             grad_phiHO%z=spinor(cmplx(0),cmplx(0))
             grad_phiHO%rho=spinor(cmplx(0),cmplx(0))
             grad_phiHO%phi=spinor(cmplx(0),cmplx(0))
             do i = 1, basis%size
                c = sp_states(iHF)%coef%cylHO(i)
                if (abs(c) < 1e-8_rk) cycle
                grad_phiHO = grad_HO(i,ih,il,basis)
                grad_psi%z   = grad_psi%z   + c * grad_phiHO%z
                grad_psi%rho = grad_psi%rho + c * grad_phiHO%rho
                grad_psi%phi = grad_psi%phi + c * grad_phiHO%phi
             end do
!             grad_psi=grad(sp_states(iHF)%coef,ih,il,basis)
             grad_psi = grad_HF_in_situ(sp_states(iHF)%coef%cylHO,basis%size,ih,il,basis)
             sum_tau=sum_tau+sp_states(iHF)%occupation * scalar_product(grad_psi,grad_psi)

             sigma_grad_psi=vector_product(sigma,grad_psi)
             sum_div_J=sum_div_J+sp_states(iHF)%occupation * scalar_product(grad_psi,sigma_grad_psi)

             psi_grad_psi=vector(scalar_product(psi,grad_psi%z), &
                  scalar_product(psi,grad_psi%rho),scalar_product(psi,grad_psi%phi))
             sum_j=sum_j+sp_states(iHF)%occupation * psi_grad_psi

             sum_rot_j=sum_rot_j+sp_states(iHF)%occupation * vector_product(grad_psi,grad_psi)

             sigma_psi=vector_spinor(sigma%z*psi,sigma%rho*psi,sigma%phi*psi)
             psi_sigma_psi=vector(scalar_product(psi,sigma_psi%z),scalar_product(psi,sigma_psi%rho), &
                  scalar_product(psi,sigma_psi%phi))
             sum_s=sum_s+sp_states(iHF)%occupation * psi_sigma_psi

             psi_sigma_grad_psi=vector(scalar_product(psi,sigma_grad_psi%z), &
                  scalar_product(psi,sigma_grad_psi%rho),scalar_product(psi,sigma_grad_psi%phi))
             sum_rot_s=sum_rot_s+sp_states(iHF)%occupation * (vector_product(grad_psi,sigma_psi)- &
                  psi_sigma_grad_psi)

             sigma_z_grad_psi=spin_operator_action(sigma%z,grad_psi)
             sigma_rho_grad_psi=spin_operator_action(sigma%rho,grad_psi)
             sigma_phi_grad_psi=spin_operator_action(sigma%phi,grad_psi)
             sum_T=sum_T+sp_states(iHF)%occupation * vector( &
                  scalar_product(grad_psi,sigma_z_grad_psi), &
                  scalar_product(grad_psi,sigma_rho_grad_psi), &
                  scalar_product(grad_psi,sigma_phi_grad_psi))

             Delta_psi_sigma_psi=vector(scalar_product(Delta_psi,sigma_psi%z), &
                  scalar_product(Delta_psi,sigma_psi%rho), &
                  scalar_product(Delta_psi,sigma_psi%phi))
             grad_psi_sigma_grad_psi=vector(scalar_product(grad_psi,sigma_z_grad_psi), &
                  scalar_product(grad_psi,sigma_rho_grad_psi), &
                  scalar_product(grad_psi,sigma_phi_grad_psi))
             sigma_z_Delta_psi=spin_operator_action(sigma%z,Delta_psi)
             sigma_rho_Delta_psi=spin_operator_action(sigma%rho,Delta_psi)
             sigma_phi_Delta_psi=spin_operator_action(sigma%phi,Delta_psi)
             psi_sigma_Delta_psi=vector(scalar_product(psi,sigma_z_Delta_psi), &
                  scalar_product(psi,sigma_rho_Delta_psi), &
                  scalar_product(psi,sigma_phi_Delta_psi))
             sum_Delta_s=sum_Delta_s+sp_states(iHF)%occupation*(Delta_psi_sigma_psi+ &
                  2.0_rk*grad_psi_sigma_grad_psi+psi_sigma_Delta_psi)
             grad_phi_sigma_rho_psi=(1.0_rk/rho)*spin_operator_action(sigma%phi,psi)
             grad_phi_sigma_phi_psi=(-1.0_rk/rho)*spin_operator_action(sigma%rho,psi)
             sum_Delta_s%rho=sum_Delta_s%rho+2.0_rk*sp_states(iHF)%occupation* &
                  scalar_product(grad_psi%phi,grad_phi_sigma_rho_psi)
             sum_Delta_s%phi=sum_Delta_s%phi+2.0_rk*sp_states(iHF)%occupation* &
                  scalar_product(grad_psi%phi,grad_phi_sigma_phi_psi)

             sum_tensor_J(1)=sum_tensor_J(1)+sp_states(iHF)%occupation* &
                  vector(scalar_product(psi,sigma_z_grad_psi%z), &
                  scalar_product(psi,sigma_rho_grad_psi%z), &
                  scalar_product(psi,sigma_phi_grad_psi%z))
             sum_tensor_J(2)=sum_tensor_J(2)+sp_states(iHF)%occupation* &
                  vector(scalar_product(psi,sigma_z_grad_psi%rho), &
                  scalar_product(psi,sigma_rho_grad_psi%rho), &
                  scalar_product(psi,sigma_phi_grad_psi%rho))
             sum_tensor_J(3)=sum_tensor_J(3)+sp_states(iHF)%occupation* &
                  vector(scalar_product(psi,sigma_z_grad_psi%phi), &
                  scalar_product(psi,sigma_rho_grad_psi%phi), &
                  scalar_product(psi,sigma_phi_grad_psi%phi))

             sigma_dot_grad_psi=spin_operator_action(sigma%z,grad_psi%z)+ &
                  spin_operator_action(sigma%rho,grad_psi%rho)+ &
                  spin_operator_action(sigma%phi,grad_psi%phi)
             sum_F=sum_F+sp_states(iHF)%occupation * vector( &
                  scalar_product(grad_psi%z,sigma_dot_grad_psi), &
                  scalar_product(grad_psi%rho,sigma_dot_grad_psi), &
                  scalar_product(grad_psi%phi,sigma_dot_grad_psi))

             sum_div_s=sum_div_s+sp_states(iHF)%occupation * &
                  (scalar_product(grad_psi,sigma_psi)+ &
                  scalar_product(psi,sigma_dot_grad_psi))
                  
          end do loop_iHF

          densities%rho(ih,il) = real(sum_rho)/(2.0_rk*pi)
          densities%tau(ih,il) = real(sum_tau)/(2.0_rk*pi)
          densities%Delta_rho(ih,il) = real(sum_Delta_rho)/pi+ &
               2.0_rk*densities%tau(ih,il)
          densities%div_J(ih,il) = -aimag(sum_div_J)/(2.0_rk*pi)

          densities%grad_rho%z(ih,il) = 2*real(sum_j%z)/(2.0_rk*pi)
          densities%grad_rho%rho(ih,il) = 2*real(sum_j%rho)/(2.0_rk*pi)
          densities%grad_rho%phi(ih,il) = 2*real(sum_j%phi)/(2.0_rk*pi)

          densities%j%z(ih,il) = aimag(sum_j%z)/(2.0_rk*pi)
          densities%j%rho(ih,il) = aimag(sum_j%rho)/(2.0_rk*pi)
          densities%j%phi(ih,il) = aimag(sum_j%phi)/(2.0_rk*pi)

          densities%rot_j%z(ih,il) = aimag(sum_rot_j%z)/(2.0_rk*pi)
          densities%rot_j%rho(ih,il) = aimag(sum_rot_j%rho)/(2.0_rk*pi)
          densities%rot_j%phi(ih,il) = aimag(sum_rot_j%phi)/(2.0_rk*pi)

          densities%s%z(ih,il) = real(sum_s%z)/(2.0_rk*pi)
          densities%s%rho(ih,il) = real(sum_s%rho)/(2.0_rk*pi)
          densities%s%phi(ih,il) = real(sum_s%phi)/(2.0_rk*pi)

          densities%rot_s%z(ih,il) = real(sum_rot_s%z)/(2.0_rk*pi)
          densities%rot_s%rho(ih,il) = real(sum_rot_s%rho)/(2.0_rk*pi)
          densities%rot_s%phi(ih,il) = real(sum_rot_s%phi)/(2.0_rk*pi)

          densities%t%z(ih,il) = real(sum_T%z)/(2.0_rk*pi)
          densities%t%rho(ih,il) = real(sum_T%rho)/(2.0_rk*pi)
          densities%t%phi(ih,il) = real(sum_T%phi)/(2.0_rk*pi)

          densities%Delta_s%z(ih,il) = real(sum_Delta_s%z)/(2.0_rk*pi)
          densities%Delta_s%rho(ih,il) = real(sum_Delta_s%rho)/(2.0_rk*pi)- &
               densities%s%rho(ih,il)/rho**2 ! No phi-derivative term (axial symmetry)
          densities%Delta_s%phi(ih,il) = real(sum_Delta_s%phi)/(2.0_rk*pi)- &
               densities%s%phi(ih,il)/rho**2 ! No phi-derivative term (axial symmetry)

          do i=1,3
             densities%tensor_J(i)%z(ih,il) = aimag(sum_tensor_J(i)%z)/(2.0_rk*pi)
             densities%tensor_J(i)%rho(ih,il) = aimag(sum_tensor_J(i)%rho)/(2.0_rk*pi)
             densities%tensor_J(i)%phi(ih,il) = aimag(sum_tensor_J(i)%phi)/(2.0_rk*pi)
          end do

          densities%F%z(ih,il) = real(sum_F%z)/(2.0_rk*pi)
          densities%F%rho(ih,il) = real(sum_F%rho)/(2.0_rk*pi)
          densities%F%phi(ih,il) = real(sum_F%phi)/(2.0_rk*pi)

          densities%div_s(ih,il) = real(sum_div_s)/(2.0_rk*pi) + &
               densities%s%rho(ih,il)/rho

          if (debugmode=='high') write(66,'(3(i3,2x),12(f12.6,2x))') isospin,ih,il, &
               densities%j%z(ih,il),densities%j%rho(ih,il),densities%j%phi(ih,il), &
               densities%rot_s%z(ih,il),densities%rot_s%rho(ih,il),densities%rot_s%phi(ih,il), &
               densities%s%z(ih,il),densities%s%rho(ih,il),densities%s%phi(ih,il), &
               densities%rot_j%z(ih,il),densities%rot_j%rho(ih,il),densities%rot_j%phi(ih,il)

       end do loop_il

    end do loop_ih

    !$OMP END PARALLEL DO

    ! tmp_sz=0
    ! tmp_jphi=0
    ! tmp_tau=0
    ! tmp_rho=0
    ! do ih=iHmin,basis%NGz/2
    !    if (ih==0) cycle
    !    do il=1,basis%NGr
    !       w=basis%G_Herm(ih)*basis%G_Lag(il)
    !       tmp_rho=tmp_rho+w*densities%rho(ih,il)
    !       tmp_tau=tmp_tau+w*densities%tau(ih,il)
    !       tmp_jphi = tmp_jphi+w*densities%j%phi(ih,il)
    !       tmp_sz = tmp_sz+w*densities%s%z(ih,il)
    !    end do
    ! end do
    ! tmp_rho = tmp_rho*fz*pi/(bz*bp2)
    ! tmp_tau = tmp_tau*fz*pi/(bz*bp2)
    ! tmp_jphi = tmp_jphi*fz*pi/(bz*bp2)
    ! tmp_sz = tmp_sz*fz*pi/(bz*bp2)
    ! write(6,'(a,i3,5(2x,a,f9.4))') 'Isospin=',isospin, &
    !      'int(rho) =',tmp_rho, &
    !      'hb^2/(2*m)*int(tau) =',hb0(isospin)*tmp_tau, &
    !      'int(j_phi) =',tmp_jphi, &
    !      'int(s_z) =',tmp_sz

  contains

    function grad_HF_in_situ(coef,n,ih,il,basis) result(grad_psi)

      implicit none
      integer,intent(in)  :: n
      real(kind=rk),dimension(n),intent(in) :: coef
      type(cylindrical_basis),intent(in) :: basis
      type(vector_spinor) :: grad_psi
      integer :: ih,il,i
      real(kind=rk) :: c
      type(vector_spinor) :: grad_phiHO

      grad_psi%z=spinor(cmplx(0),cmplx(0))
      grad_psi%rho=spinor(cmplx(0),cmplx(0))
      grad_psi%phi=spinor(cmplx(0),cmplx(0))
      grad_phiHO%z=spinor(cmplx(0),cmplx(0))
      grad_phiHO%rho=spinor(cmplx(0),cmplx(0))
      grad_phiHO%phi=spinor(cmplx(0),cmplx(0))
      do i = 1, basis%size
         c = coef(i)
         if (abs(c) < 1e-8_rk) cycle
         grad_phiHO = grad_HO(i,ih,il,basis)
         grad_psi%z   = grad_psi%z   + c * grad_phiHO%z
         grad_psi%rho = grad_psi%rho + c * grad_phiHO%rho
         grad_psi%phi = grad_psi%phi + c * grad_phiHO%phi
      end do

    end function grad_HF_in_situ


  end subroutine local_one_body_densities


  !=====================================================================
  !                       SUBROUTINE RENORMALIZE_DENSITIES
  !---------------------------------------------------------------------
  ! Calculates the renormalization factor and apply it to all local 
  ! densities.
  !=====================================================================

  subroutine renormalize_densities(densities, basis, nucleus)

!
!   Arguments
!

    type(local_density),dimension(2),intent(inout) :: densities
    type(cylindrical_basis),intent(in)             :: basis
    type(nucleus_properties),intent(inout)         :: nucleus
!
!   Local variables
!
    character(len=21),parameter :: subroutinename='renormalize_densities'
    integer,dimension(2) :: Npart
    integer :: isospin
    integer :: iHerm,iLag
    real(kind=rk) :: Gz,Gr
    real(kind=rk) :: norm,renorm,weight
    real(kind=rk) :: sum_sz

!
!   Debugging note
!

    if (debugmode=='high') call debugnote(modulename, subroutinename)

!
!   Auxiliary variables
!

    Npart(:)=(/ nucleus%N, nucleus%Z /)
    norm=fz*pi/(bz*bp2)

!
!   Calculating renormalization factor and apply it to all local densities
!

    densities(:)%renorm=0.0

    do isospin=1,2

       do iHerm=iHmin,basis%NGz/2
          if (iHerm==0) cycle
          Gz=basis%G_Herm(iHerm)
          do iLag=1,basis%NGr
             Gr=basis%G_Lag(iLag)
             weight=Gz*Gr*densities(isospin)%rho(iHerm,iLag)
             densities(isospin)%renorm=densities(isospin)%renorm+weight
          end do
       end do

       if (Npart(isospin)/=0) then
          densities(isospin)%renorm=Npart(isospin)/(norm*densities(isospin)%renorm)
       else
          densities(isospin)%renorm=1.0
       end if
       renorm=densities(isospin)%renorm

       sum_sz=0.0
       do iHerm=iHmin,basis%NGz/2
          if (iHerm==0) cycle
          do iLag=1,basis%NGr
             densities(isospin)%rho(iHerm,iLag)=renorm*densities(isospin)%rho(iHerm,iLag)
             densities(isospin)%tau(iHerm,iLag)=renorm*densities(isospin)%tau(iHerm,iLag)
             densities(isospin)%Delta_rho(iHerm,iLag)=renorm*densities(isospin)%Delta_rho(iHerm,iLag)
             densities(isospin)%div_J(iHerm,iLag)=renorm*densities(isospin)%div_J(iHerm,iLag)

             densities(isospin)%grad_rho%z(iHerm,iLag)=renorm*densities(isospin)%grad_rho%z(iHerm,iLag)
             densities(isospin)%grad_rho%rho(iHerm,iLag)=renorm*densities(isospin)%grad_rho%rho(iHerm,iLag)
             densities(isospin)%grad_rho%phi(iHerm,iLag)=renorm*densities(isospin)%grad_rho%phi(iHerm,iLag)

             densities(isospin)%j%z(iHerm,iLag)=renorm*densities(isospin)%j%z(iHerm,iLag)
             densities(isospin)%j%rho(iHerm,iLag)=renorm*densities(isospin)%j%rho(iHerm,iLag)
             densities(isospin)%j%phi(iHerm,iLag)=renorm*densities(isospin)%j%phi(iHerm,iLag)

             densities(isospin)%s%z(iHerm,iLag)=renorm*densities(isospin)%s%z(iHerm,iLag)
             densities(isospin)%s%rho(iHerm,iLag)=renorm*densities(isospin)%s%rho(iHerm,iLag)
             densities(isospin)%s%phi(iHerm,iLag)=renorm*densities(isospin)%s%phi(iHerm,iLag)

             densities(isospin)%rot_j%z(iHerm,iLag)=renorm*densities(isospin)%rot_j%z(iHerm,iLag)
             densities(isospin)%rot_j%rho(iHerm,iLag)=renorm*densities(isospin)%rot_j%rho(iHerm,iLag)
             densities(isospin)%rot_j%phi(iHerm,iLag)=renorm*densities(isospin)%rot_j%phi(iHerm,iLag)

             densities(isospin)%rot_s%z(iHerm,iLag)=renorm*densities(isospin)%rot_s%z(iHerm,iLag)
             densities(isospin)%rot_s%rho(iHerm,iLag)=renorm*densities(isospin)%rot_s%rho(iHerm,iLag)
             densities(isospin)%rot_s%phi(iHerm,iLag)=renorm*densities(isospin)%rot_s%phi(iHerm,iLag)

             densities(isospin)%Delta_s%z(iHerm,iLag)=renorm*densities(isospin)%Delta_s%z(iHerm,iLag)
             densities(isospin)%Delta_s%rho(iHerm,iLag)=renorm*densities(isospin)%Delta_s%rho(iHerm,iLag)
             densities(isospin)%Delta_s%phi(iHerm,iLag)=renorm*densities(isospin)%Delta_s%phi(iHerm,iLag)

             densities(isospin)%T%z(iHerm,iLag)=renorm*densities(isospin)%T%z(iHerm,iLag)
             densities(isospin)%T%rho(iHerm,iLag)=renorm*densities(isospin)%T%rho(iHerm,iLag)
             densities(isospin)%T%phi(iHerm,iLag)=renorm*densities(isospin)%T%phi(iHerm,iLag)

             densities(isospin)%F%z(iHerm,iLag)=renorm*densities(isospin)%F%z(iHerm,iLag)
             densities(isospin)%F%rho(iHerm,iLag)=renorm*densities(isospin)%F%rho(iHerm,iLag)
             densities(isospin)%F%phi(iHerm,iLag)=renorm*densities(isospin)%F%phi(iHerm,iLag)

             densities(isospin)%div_s(iHerm,iLag)=renorm*densities(isospin)%div_s(iHerm,iLag)

             densities(isospin)%tensor_J(3)%z(iHerm,iLag) = renorm * densities(isospin)%tensor_J(3)%z(iHerm,iLag)
             densities(isospin)%tensor_J(3)%rho(iHerm,iLag) = renorm * densities(isospin)%tensor_J(3)%rho(iHerm,iLag)

             sum_sz=sum_sz+basis%G_Herm(iHerm)*basis%G_Lag(iLag)*densities(isospin)%s%z(iHerm,iLag)
          end do
       end do
       nucleus%Sz(isospin)=sum_sz*norm*0.5_rk

    end do

    nucleus%Sz(0)=nucleus%Sz(1)+nucleus%Sz(2)

    if (debugmode=='high') write(out_unit,'(2(a,1x,f12.6/))') &
         'Neutron renormalization factor =',densities(1)%renorm, &
         'Proton renormalization factor  =',densities(2)%renorm

    return

  end subroutine renormalize_densities


  !=====================================================================
  !                       SUBROUTINE CONSTRAINTS_READJUSTMENT
  !---------------------------------------------------------------------
  ! Automatic readjustment of linear constraints using perturbation 
  ! theory and Inglis-Belyaev mass parameters.
  !=====================================================================

  subroutine constraints_readjustment(nucleus,basis,sp_states,model,HFfield,error,accelerate)

    implicit none
!
!   Arguments
!
    type(nucleus_properties),intent(in) :: nucleus
    type(cylindrical_basis),intent(in) :: basis
    type(sp_state),dimension(basis%size,2),intent(in) :: sp_states
    type(modelization),intent(inout) :: model
    type(field),dimension(2),intent(inout) :: HFfield
    real(kind=rk),intent(out) :: error
    logical,intent(inout) :: accelerate

!
!   Local variables
!
    character(len=24),parameter :: subroutinename = 'constraints_readjustment'
    integer :: i, j, isospin, iHerm, iLag
    real(kind=rk),dimension(:),allocatable :: Qdiff,dlambda
    real(kind=rk),dimension(Nconstr) :: Qdiff_max
    real(kind=rk),dimension(:,:),allocatable :: Qmass,Qmass_inv
    real(kind=rk) :: current_constraint, z, r

    if (nreadj <= 0 .or. .not.allocated(readjusted_constraint_index)) return

    if (debugmode /= 'low') call warning(modulename,subroutinename,'Readjusting selected constraints')

!----------------------------------------------------------------------------------------------
!
!   First order variation of Lagrange multipliers lambda_i for readjusted, linear constraints
!
!   Note: lambda is represented by model%constraints(i)%stiffness
!
!----------------------------------------------------------------------------------------------

    allocate(Qdiff(nreadj),Qmass(nreadj,nreadj),Qmass_inv(nreadj,nreadj),dlambda(nreadj))

!
!   Difference between current and targeted values of readjusted constraints
!
    Qdiff_max(:) = 1.0_rk
    Qdiff_max(1) = 1.0_rk
    Qdiff_max(2) = 100.0_rk
    Qdiff_max(3) = 100.0_rk
    Qdiff_max(4) = 100.0_rk
    error = 0.0
    do i = 1, nreadj
       j = readjusted_constraint_index(i)
       select case(j)
       case (1)
          current_constraint = nucleus%zcm(0)
       case (2)
          current_constraint = nucleus%Q20(0)
       case (3)
          current_constraint = nucleus%Q30(0)
       case (4)
          current_constraint = nucleus%Q40(0)
       case default
          call fatal(modulename,subroutinename,'Unimplemented readjusted constraint')
       end select
       Qdiff(i) = model%constraints(j)%value - current_constraint
!       if (abs(Qdiff(i)) > Qdiff_max(j)) Qdiff(i) = sign(Qdiff_max(j),Qdiff(i))
       error = error + abs(Qdiff(i)/max(1.0_rk,abs(model%constraints(j)%value)))
    end do

    call mass_parameters
    call inverse(Qmass,Qmass_inv,nreadj)

!
!   Variations dlambda_i of Lagrange multipliers
!
    do i = 1, nreadj
       dlambda(i) = 0
       do j = 1, nreadj
          dlambda(i) = dlambda(i) + Qmass_inv(i,j)*Qdiff(j)
          write(out_unit,'(2(i2,1x),4(a,e15.6,2x))') &
               i,j,'Qmass=',Qmass(i,j),'Qmass_inv=',Qmass_inv(i,j),&
               'Qdiff=',Qdiff(j),'dlambda=',dlambda(i)
       end do
    end do

!
!   Update of Lagrange multipliers
!
    do i = 1, nreadj
       j = readjusted_constraint_index(i)
       model%constraints(j)%stiffness = model%constraints(j)%stiffness + dlambda(i)
       write(out_unit,'(i2,3(2x,a,e15.6))') j,'old:',model%constraints(j)%stiffness-dlambda(i), &
            'new:',model%constraints(j)%stiffness,'dlambda=',dlambda(i)
    end do

!
!   Update central potential of the HF field
!

    do isospin = 1, 2
       do iHerm = iHmin, basis%NGz/2
          if (iHerm == 0) cycle
          z = basis%xHerm(iHerm)/bz
          do iLag = 1, basis%NGr
             r = sqrt(basis%xlag(iLag))/bp
             do i = 1, nreadj
                j = readjusted_constraint_index(i)
                HFfield(isospin)%central(iHerm,iLag) = HFfield(isospin)%central(iHerm,iLag) - &
                     dlambda(i) * fconstr(j,z,r)
             end do
          end do
       end do
    end do

    deallocate(dlambda,Qdiff,Qmass,Qmass_inv)

    contains

      subroutine mass_parameters

        implicit none
        integer :: i, k1, j, k2, isospin, iHF, jHF
        integer :: iHerm, iLag
        integer :: block_index, block_size, shift
        integer :: iHO, jHO
        real(kind=rk) :: Ci, ui, vi, Eqpi, Cj, uj, vj, Eqpj
        real(kind=rk) :: norm, z, wz, r, weight, phi
        real(kind=rk),dimension(Nconstr,basis%size,basis%size,2) :: Qmx
        real(kind=rk),dimension(2) :: particle_number

        norm = 0.5_rk*fz/(bz*bp**2)

        Qmx = 0.0

        do isospin = 1, 2

           particle_number(isospin)=0.0

           do iHF = 1, basis%size

              block_index = sp_states(iHF,isospin)%block
              block_size = basis%block_size(block_index)
              shift=sum(basis%block_size(1:block_index))-basis%block_size(block_index)

              do iHerm = iHmin, basis%NGz/2
                 if (iHerm == 0) cycle
                 wz = basis%G_Herm(iHerm)
                 do iLag = 1, basis%NGr
                    weight = wz * basis%G_Lag(iLag)
                    particle_number(isospin) = particle_number(isospin) + norm*weight * &
                         sp_states(iHF,isospin)%occupation * &
                         (sp_states(iHF,isospin)%wf(iHerm,iLag)%up**2 + &
                         sp_states(iHF,isospin)%wf(iHerm,iLag)%down**2)
                 end do
              end do

              do jHF = iHF, basis%size

                 Qmx(:,iHF,jHF,isospin) = 0.0
                 if (sp_states(jHF,isospin)%block /= block_index) cycle ! constraint operators commute with Jz (and parity if not broken)

!                Calculation in coordinate representation
!                Note: constraint operators do not depend on spin

                 do iHerm = iHmin, basis%NGz/2

                    if (iHerm == 0) cycle

                    z = basis%xHerm(iHerm)/bz
                    wz = basis%G_Herm(iHerm)

                    do iLag = 1, basis%NGr

                       r = sqrt(basis%xLag(iLag))/bp
                       weight = wz * basis%G_Lag(iLag)
                       phi = sp_states(iHF,isospin)%wf(iHerm,iLag)%up * &
                            sp_states(jHF,isospin)%wf(iHerm,iLag)%up + &
                            sp_states(iHF,isospin)%wf(iHerm,iLag)%down * &
                            sp_states(jHF,isospin)%wf(iHerm,iLag)%down
                       do i = 1, nreadj
                          j = readjusted_constraint_index(i)
                          Qmx(j,iHF,jHF,isospin) = Qmx(j,iHF,jHF,isospin) + weight * phi * fconstr(j,z,r)
                       end do

                    end do

                 end do

                 Qmx(:,iHF,jHF,isospin) = norm * Qmx(:,iHF,jHF,isospin)
!                 Qmx(1,iHF,jHF,isospin) = Qmx(1,iHF,jHF,isospin) / nucleus%A ! op. Q(1) = z, not zcm
                 if (jHF /= iHF) then
                    do i = 1, Nconstr
                       Qmx(i,jHF,iHF,isospin) = Qmx(i,iHF,jHF,isospin)
                    end do
                 end if

              end do

           end do

        end do

!       Trick to accelerate convergence of constraint on zcm

        Qmx(1,:,:,:) = Qmx(1,:,:,:) / model%readj_acceleration_factor

!        if (.not.accelerate .and. trim(model%constraints(1)%readj) == 'y') then
!           if (abs(Qdiff(1)) < 0.1_rk .and. &
!                abs(error-abs(Qdiff(1)/max(1.0_rk,abs(model%constraints(1)%value)))) < 0.01_rk) then
!              accelerate = .true.
!              write(out_unit,'(a)') 'Accelerating convergence of constraint on zcm'
!              Qmx(1,:,:,:) = Qmx(1,:,:,:) / model%readj_acceleration_factor
!           end if
!        end if

!        print*,nucleus%N,particle_number(1)
!        print*,nucleus%Z,particle_number(2)

!        z = 0.0
!        do isospin = 1, 2
!           do iHF = 1, basis%size
!              z = z + sp_states(iHF,isospin)%occupation * Qmx(1,iHF,iHF,isospin)
!           end do
!        end do
!        print*,'zcm=',z,' nucleus%zcm=',nucleus%zcm(0)

        do i = 1, nreadj

           k1 = readjusted_constraint_index(i)

           do j = 1, i

              k2 = readjusted_constraint_index(j)

! Belyaev expression

              Qmass(i,j) = 0.0

              do isospin = 1, 2
                 do iHF = 1, basis%size
                    block_index = sp_states(iHF,isospin)%block
                    vi = sqrt(sp_states(iHF,isospin)%occupation)
                    ui = sqrt(1.0_rk-sp_states(iHF,isospin)%occupation)
                    Eqpi = sqrt((sp_states(iHF,isospin)%energy-nucleus%pairing(isospin)%chemical_potential)**2 + &
                         nucleus%pairing(isospin)%gap(iHF)**2)
                    do jHF = 1, basis%size
                       if (sp_states(jHF,isospin)%block /= block_index) cycle
                       vj = sqrt(sp_states(jHF,isospin)%occupation)
                       uj = sqrt(1.0_rk-sp_states(jHF,isospin)%occupation)
                       Eqpj = sqrt((sp_states(jHF,isospin)%energy-nucleus%pairing(isospin)%chemical_potential)**2 + &
                            nucleus%pairing(isospin)%gap(jHF)**2)
                       Qmass(i,j) = Qmass(i,j) + 2.0_rk*(ui*vj-uj*vi)**2/(Eqpi+Eqpj) * &
                            Qmx(k1,iHF,jHF,isospin)*Qmx(k2,iHF,jHF,isospin)
!                       if (abs(Qmx(k1,iHF,jHF,isospin)*Qmx(k2,iHF,jHF,isospin))*(ui*vj-uj*vi)**2>1e-6_rk) &
!                            write(out_unit,'(5(i5),2(e15.6,1x))') &
!                            k1,k1,isospin,iHF,jHF,Qmx(k1,iHF,jHF,isospin)*Qmx(k2,iHF,jHF,isospin),&
!                            (ui*vj-uj*vi)**2
                    end do
                 end do
              end do

              if (j /= i) Qmass(j,i) = Qmass(i,j)

           end do

        end do
        
      end subroutine mass_parameters

  end subroutine constraints_readjustment



  !=====================================================================
  !                       SUBROUTINE COULOMB
  !---------------------------------------------------------------------
  ! Calculates the direct part of the Coulomb potential.
  !=====================================================================

  subroutine Coulomb(HFfield, densities, basis, symmetries)

!
!   Arguments
!
    type(field),intent(inout)      :: HFfield
    type(local_density),dimension(2),intent(in) :: densities
    type(cylindrical_basis),intent(in)          :: basis
    type(symmetry_properties),intent(in)        :: symmetries
!
!   Local variables
!
    real(kind=rk),dimension(:),allocatable :: rLag
    integer :: iHerm,iLag,iz,ir
    real(kind=rk) :: z,z_prime,r,r_prime,Gz,Gr
    real(kind=rk) :: arg,elliptic_integral,dist
    real(kind=rk) :: V

!
!   Auxiliary variables
!

    allocate(rLag(basis%NGr))
    do iLag=1,basis%NGr
       rLag(iLag)=sqrt(basis%xLag(iLag))/bp
    end do

!
!   Calculating direct Coulomb potential
!

    do iHerm=iHmin,basis%NGz/2

       if (iHerm==0) cycle

       z=basis%xHerm(iHerm)/bz

       do iLag=1,basis%NGr

          r=rLag(iLag)

          V=0.
          do iz=iHmin,basis%NGz/2

             if (iz==0) cycle

             z_prime=basis%xHerm(iz)/bz
             Gz=basis%G_Herm(iz)

             do ir=1,basis%NGr

                r_prime=rLag(ir)
                Gr=basis%G_Lag(ir)
                dist=sqrt((z-z_prime)**2+(r+r_prime)**2)
                arg=1.0_rk-4.0_rk*r*r_prime/dist**2
                elliptic_integral=1.0_rk+arg*(0.4630151_rk+0.1077812_rk*arg)
                if (arg > 1e-6_rk) elliptic_integral=elliptic_integral-&
                     arg*(0.2452727_rk+0.0412496_rk*arg)*log(arg)
                V=V+Gz*Gr*densities(2)%Delta_rho(iz,ir)*elliptic_integral*dist

                if (symmetries%parity) then
                   dist=sqrt((z+z_prime)**2+(r+r_prime)**2)
                   arg=1.0_rk-4.0_rk*r*r_prime/dist**2
                   elliptic_integral=1.0_rk+arg*(0.4630151_rk+0.1077812_rk*arg)
                   if (arg > 1e-6_rk) elliptic_integral=elliptic_integral-&
                     arg*(0.2452727_rk+0.0412496_rk*arg)*log(arg)
                   V=V+Gz*Gr*densities(2)%Delta_rho(iz,ir)*elliptic_integral*dist
                end if

             end do

          end do

          HFfield%Coulomb(iHerm,iLag)=e2*V/(bz*bp2)*c_Coulomb

       end do

    end do

    deallocate(rLag)

    return

  end subroutine Coulomb


  !=====================================================================
  !                       SUBROUTINE MOMENTS
  !---------------------------------------------------------------------
  ! Calculates moments of the nuclear density.
  !=====================================================================

  subroutine moments(nucleus, densities, sp_states, basis, symmetries)

!
!   Arguments
!

    type(nucleus_properties),intent(inout)            :: nucleus
    type(local_density),dimension(2),intent(inout)    :: densities
    type(cylindrical_basis),intent(in)                :: basis
    type(sp_state),dimension(basis%size,2),intent(in) :: sp_states
    type(symmetry_properties),intent(in)              :: symmetries
!
!   Local variables
!
    character(len=7),parameter :: subroutinename='moments'
    integer :: isospin,iHF
    integer :: iHerm,iLag
    real(kind=rk) :: z,Gz,r,Gr
    real(kind=rk) :: norm,weight,tmp

!
!   Auxiliary variables
!

    norm=fz*pi/(bz*bp2)

!
!   Calculating moments of the nuclear density for neutrons, protons and total
!
!   <zcm> and <Q30> are identically zero in the case of parity symmetry
!

    nucleus%zcm(:)=0.0
    nucleus%r2(:)=0.0
    nucleus%Q20(:)=0.0
    nucleus%Q30(:)=0.0
    nucleus%Q40(:)=0.0
    nucleus%QN(:)=0.0

    do isospin=1,2

       do iHerm=iHmin,basis%NGz/2
          if (iHerm==0) cycle
          z=basis%xHerm(iHerm)/bz
          Gz=basis%G_Herm(iHerm)
          do iLag=1,basis%NGr
             r=sqrt(basis%xLag(iLag))/bp
             Gr=basis%G_Lag(iLag)
             weight=Gz*Gr*densities(isospin)%rho(iHerm,iLag)
             nucleus%r2(isospin)=nucleus%r2(isospin)+weight*(z**2+r**2)
             nucleus%Q20(isospin)=nucleus%Q20(isospin)+weight*(2.0_rk*z**2-r**2)
             nucleus%Q40(isospin)=nucleus%Q40(isospin)+weight*(z**4-3.0_rk*z**2*r**2+0.375_rk*r**4)
             nucleus%QN(isospin)=nucleus%QN(isospin)+weight*exp(-(z/a_QN)**2) ! z_min=0
          end do
       end do

       nucleus%r2(isospin)=norm*nucleus%r2(isospin)
       nucleus%Q20(isospin)=norm*nucleus%Q20(isospin)
       nucleus%Q40(isospin)=sqrt(9.0_rk/(4.0_rk*pi))*norm*nucleus%Q40(isospin)
       nucleus%QN(isospin)=norm*nucleus%QN(isospin)

       if (.not.symmetries%parity) then

          tmp = 0.0
          do iHerm=iHmin,basis%NGz/2
             if (iHerm==0) cycle
             z=basis%xHerm(iHerm)/bz
             Gz=basis%G_Herm(iHerm)
             do iLag=1,basis%NGr
                r=sqrt(basis%xLag(iLag))/bp
                Gr=basis%G_Lag(iLag)
                weight=Gz*Gr*densities(isospin)%rho(iHerm,iLag)
                densities(isospin)%renorm=densities(isospin)%renorm+weight
                tmp=tmp+weight
                nucleus%zcm(isospin)=nucleus%zcm(isospin)+weight*z
                nucleus%Q30(isospin)=nucleus%Q30(isospin)+weight*z*(z**2-1.5_rk*r**2)
             end do
          end do

!          nucleus%zcm(isospin)=norm*nucleus%zcm(isospin)/densities(isospin)%renorm
          nucleus%zcm(isospin)=nucleus%zcm(isospin)/tmp
          nucleus%Q30(isospin)=sqrt(7.0_rk/(4.0_rk*pi))*norm*nucleus%Q30(isospin)

          nucleus%avg_pi(isospin) = 0.0
          do iHF = 1, basis%size
             nucleus%avg_pi(isospin) = nucleus%avg_pi(isospin) * &
                  sp_states(iHF,isospin)%occupation*sp_states(iHF,isospin)%avg_pi
          end do

       end if

!       nucleus%zcm(0)=nucleus%zcm(0)+nucleus%zcm(isospin)
       nucleus%r2(0)=nucleus%r2(0)+nucleus%r2(isospin)
       nucleus%Q20(0)=nucleus%Q20(0)+nucleus%Q20(isospin)
       nucleus%Q30(0)=nucleus%Q30(0)+nucleus%Q30(isospin)
       nucleus%Q40(0)=nucleus%Q40(0)+nucleus%Q40(isospin)
       nucleus%QN(0)=nucleus%QN(0)+nucleus%QN(isospin)

       tmp=0.0
       do iHF=1,basis%size
          tmp=tmp+sp_states(iHF,isospin)%occupation*sp_states(iHF,isospin)%Omega/2.0_rk
       end do
       nucleus%Jz(isospin)=tmp

    end do

    nucleus%zcm(0)=(nucleus%N*nucleus%zcm(1)+nucleus%Z*nucleus%zcm(2))/nucleus%A
    nucleus%Jz(0)=nucleus%Jz(1)+nucleus%Jz(2)
    nucleus%avg_pi(0) = nucleus%avg_pi(1) * nucleus%avg_pi(2)

    nucleus%beta2 = sqrt(pi/5.0_rk)*nucleus%Q20(0)/nucleus%r2(0)
    nucleus%beta3 = 4*pi*nucleus%Q30(0)/(3*nucleus%A**2*1.2_rk**3)

    if (debugmode=='high' .or. nreadj > 0) write(out_unit,'(6(a,1x,f12.6,1x,a/))') &
         'zcm=',nucleus%zcm(0),'fm','r2=',nucleus%r2(0),'fm^2', &
         'Q20=',nucleus%Q20(0),'fm^2','Q30=',nucleus%Q30(0),'fm^3', &
         'Q40=',nucleus%Q40(0),'fm^4','QN=',nucleus%QN(0),' '

    return

  end subroutine moments

  subroutine nuclear_shape(sp_states, basis, symmetries, nucleus)

!
!   Arguments
!
    type(cylindrical_basis),intent(in)                :: basis
    type(sp_state),dimension(basis%size,2),intent(in) :: sp_states
    type(symmetry_properties),intent(in)              :: symmetries
    type(nucleus_properties),intent(in)               :: nucleus

!
!   Local variables
!
    real(kind=rk),parameter :: rho_sat = 0.16_rk ! saturation density
    real(kind=rk),parameter :: dr = 0.01_rk, dz = dr ! steps in r and z
    real(kind=rk) :: rho_contour, rho1, rho2
    real(kind=rk) :: r, rmax, z, zmin, zmax
    real(kind=rk) :: theta_max, dtheta, theta, rt, cos_theta, sin_theta
    real(kind=rk) :: vol, R0, c, zcm, newc, newzcm
    integer :: i, l, iH
    integer,parameter :: ntheta = 50, nbeta = 4
    real(kind=rk),dimension(ntheta) :: newR, alpha
    real(kind=rk),dimension(nbeta) :: beta
    integer,dimension(8) :: start, end


    call DATE_AND_TIME(values = start)

!--------------------------------------------------------------------------
!
!               DETERMINATION OF ISODENSITY CONTOUR
!               AND ITS DEFORMATION CHARACTERISTICS
!
!--------------------------------------------------------------------------

    rho_contour = rho_sat/2
    theta_max = pi
    if (symmetries%parity) theta_max = pi/2
    dtheta = theta_max/(ntheta-1)

!
!   Volume, R0 parameter, scale factor c and center-of-mass position
!   (by trapezoidal rule)
!
    theta = 0.0_rk
    vol = 0.0_rk
    c = 0.0_rk
    zcm = 0.0_rk
    do i = 2, ntheta-1
       rt = rtheta(theta)
       vol = vol + sin(theta)*rt**3*dtheta
       c = c + sin(theta)*rt*dtheta
       zcm = zcm + sin(theta)*cos(theta)*rt**4*dtheta
       theta = theta+dtheta
    end do
    vol = vol * 2*pi/3
    zcm = 0.5_rk*pi*zcm/vol
    if (symmetries%parity) then
       vol = vol * 2
       c = c * 2
       zcm = 0.0_rk
    end if
    R0 = (3*vol/(4*pi))**(1.0_rk/3.0_rk)
    c = c/(2*R0)
    write(out_unit,'(/3(a,f9.3,1x,a/))') &
         'Volume         =',vol,'fm^3', &
         'R0             =',R0,'fm', &
         'Scale factor c =',c
    write(out_unit,'(/a,f8.3,1x,a)') 'Center of mass position: z_cm =',zcm,'fm'
    if (.not. symmetries%parity) then
       iH = -basis%NGz/2
       do while(zcm > basis%xHerm(iH)/bz)
          iH = iH + 1
          if (iH == basis%NGz/2) exit
       end do
       write(out_unit,'(3x,2(a,i3,1x,a,f8.3,2x))') &
            'iH=',iH-1,'z_iH=',basis%xHerm(iH-1)/bz, &
            'iH+1=',iH,'z_{iH+1}=',basis%xHerm(iH)/bz
    end if

!
!   beta_l parameters by integration of R(theta)*Leg(l,cos(theta))*sin(theta)
!
    write(out_unit,'(/a,f5.2)') 'Bohr parameters from ' // &
         'isodensity contour rho/rho_sat = ',rho_contour/rho_sat
    do l = 1, nbeta
       beta(l) = 0
       theta = 0
       do i = 2, ntheta-1
          rt = rtheta(theta)
          beta(l) = beta(l) + sin(theta)*rt*sqrt((2*l+1)/(4*pi))*Leg_poly(l,cos(theta))*dtheta
          theta = theta+dtheta
       end do
       beta(l) = beta(l) * 2*pi/(c*R0)
       write(out_unit,'(a,i1,a,f8.3)') 'beta',l,'=',beta(l)
    end do
 
!
!   beta_l parameters from multipole moments
!
    beta(2) = nucleus%beta2
    write(out_unit,'(/a,/a,f8.3)') 'Bohr parameters from multipole moments:', &
         'beta2 =',beta(2)

    if (abs(zcm) < 0.005_rk) then
       call DATE_AND_TIME(values = end)
       call elapsed_time(start,end,"nuclear shape")
       return
    end if

!--------------------------------------------------------------------------
!
!               SHIFTING CENTER OF MASS AND 
!           RECALCULATION OF ISODENSITY CONTOUR
!
!           new z = old z - zcm
!
!--------------------------------------------------------------------------

    open(unit=2,file='density_contour',status='replace')

    theta = 0.0_rk
    alpha(1) = theta
    rt = rtheta(theta)
    newr(1) = rt - zcm
    write(out_unit,'(4(a,f8.3,1x))') 'theta=',theta,'R=',rt, &
         'alpha=',alpha(1)*180/pi,'newR=',newR(1)
    write(2,'(4(f8.3,1x))') rt*cos(theta),rt*sin(theta), &
         newR(1)*cos(alpha(1)),newR(1)*sin(alpha(1))

    theta = dtheta
    newc = 0.0_rk
    vol = 0.0_rk
    newzcm = 0.0_rk
    do i = 2, ntheta-1
       rt = rtheta(theta)
       newR(i) = sqrt(rt**2+zcm**2-2*zcm*rt*cos(theta))
       alpha(i) = acos((rt*cos(theta)-zcm)/newR(i))
       newc = newc + 0.5_rk*(sin(alpha(i-1))*newR(i-1) &
            + sin(alpha(i))*newR(i)) * (alpha(i)-alpha(i-1))
       vol = vol + 0.5_rk*(sin(alpha(i-1))*newR(i-1)**3 &
            + sin(alpha(i))*newR(i)**3) * (alpha(i)-alpha(i-1))
       newzcm = newzcm + 0.5_rk*(sin(alpha(i-1))*cos(alpha(i-1))*newR(i-1)**4 &
            + sin(alpha(i))*cos(alpha(i))*newR(i)**4) * (alpha(i)-alpha(i-1))
       write(out_unit,'(4(a,f8.3,1x))') 'theta=',theta*180/pi,'R=',rt, &
            'alpha=',alpha(i)*180/pi,'newR=',newR(i)
       write(2,'(4(f8.3,1x))') rt*cos(theta),rt*sin(theta), &
            newR(i)*cos(alpha(i)),newR(i)*sin(alpha(i))
       theta = theta + dtheta
    end do
    newc = newc/(2*R0)
    vol = vol * 2*pi/3.0_rk
    newzcm = 0.5_rk*pi*newzcm/vol

    theta = pi
    rt = rtheta(theta)
    alpha(ntheta) = theta
    newr(ntheta) = rt + zcm
    write(out_unit,'(4(a,f8.3,1x))') 'theta=',theta,'R=',rt, &
         'alpha=',alpha(ntheta)*180/pi,'newR=',newR(ntheta)
    write(2,'(4(f8.3,1x))') rt*cos(theta),rt*sin(theta), &
         newR(ntheta)*cos(alpha(ntheta)),newR(ntheta)*sin(alpha(ntheta))

    close(2)

    write(out_unit,'(/a,f8.3,1x,a)') 'Volume of shifted nuclear shape =',vol,'fm^3'
    write(out_unit,'(a,f8.3)') 'Scale factor of shifted nuclear shape =',newc
    write(out_unit,'(a,f8.3,1x,a)') 'Center of mass of shifted nuclear shape =',newzcm,'fm'

    write(out_unit,'(a)') 'beta_l values of shifted nuclear shape:'
    do l = 1, nbeta
       beta(l) = 0.0_rk
       do i = 2, ntheta
          beta(l) = beta(l) + 0.5_rk*(sin(alpha(i-1))*newR(i-1)* &
               sqrt((2*l+1)/(4*pi))*Leg_poly(l,cos(alpha(i-1))) &
               + sin(alpha(i))*newR(i)* &
               sqrt((2*l+1)/(4*pi))*Leg_poly(l,cos(alpha(i)))) * &
               (alpha(i)-alpha(i-1))
       end do
       beta(l) = beta(l) * 2*pi/(newc*R0)
       write(out_unit,'(a,i1,a,f8.3)') 'beta',l,'=',beta(l)
    end do
    
    write(out_unit,*)
    
    call DATE_AND_TIME(values = end)
    call elapsed_time(start,end,"nuclear shape")

  contains

    function rho(r,z)

      implicit none
      real(kind=rk) :: r,z, rho
      real(kind=rk) :: zeta, eta, s, tmp_plus, tmp_minus
      integer :: isospin, iHF, iblock, i, first_HO, last_HO, iHO
      integer :: nzi, nri, li, si

      zeta = z*bz
      eta = r**2*bp2
      s = 0.0_rk
      do isospin = 1, 2
         do iHF = 1, basis%size
            if (sp_states(iHF,isospin)%occupation <= 1e-4_rk) cycle
            iblock = sp_states(iHF,isospin)%block
            first_HO = sum([(basis%block_size(i),i=1,iblock)]) - basis%block_size(iblock) + 1
            last_HO = first_HO + basis%block_size(iblock) - 1
            tmp_plus = 0.0_rk
            tmp_minus = 0.0_rk
            do iHO = first_HO, last_HO
               nzi = basis%nz(iHO)
               nri = basis%nr(iHO)
               li=basis%Lambda(iHO)
               si=basis%ns(iHO)
               if (si == 1) then
                  tmp_plus = tmp_plus + sp_states(iHF,isospin)%coef%cylHO(iHO) * &
                       sqrt(fact(nri)/fact(nri+abs(li))/(2**nzi*fact(nzi))) * &
                       sqrt(eta**(abs(li)))*Laguerre_poly(nri,real(abs(li),kind=rk),eta)* &
                       Hermite_poly(nzi,zeta)
               else
                  tmp_minus = tmp_minus + sp_states(iHF,isospin)%coef%cylHO(iHO) * &
                       sqrt(fact(nri)/fact(nri+abs(li))/(2**nzi*fact(nzi))) * &
                       sqrt(eta**(abs(li)))*Laguerre_poly(nri,real(abs(li),kind=rk),eta)* &
                       Hermite_poly(nzi,zeta)
               end if
            end do
            s = s + sp_states(iHF,isospin)%occupation * (tmp_plus**2 + tmp_minus**2)
         end do
      end do

      rho = (bz*bp2/pi**1.5_rk) * exp(-(eta+zeta**2)) * s

    end function rho

    function rtheta(theta)

      implicit none
      real(kind=rk) :: theta, rtheta
      real(kind=rk) :: rho1, rho2, u, du
      integer :: i

      du = dr
      u = 0.0_rk
      rho2 = rho(u*sin(theta),u*cos(theta))
      rho1 = rho(du*sin(theta),du*cos(theta))
      i = 2
      do while((rho_contour-rho1)*(rho_contour-rho2) > 0.0_rk)
         u = u+du
         rho2 = rho1
         rho1 = rho(u*sin(theta),u*cos(theta))
         i = i + 1
         if (i == 1000) stop 'Isodensity contour not found'
      end do
      rtheta = u-0.5_rk*du

    end function rtheta

  end subroutine nuclear_shape


  !=====================================================================
  !                       SUBROUTINE ENERGY
  !---------------------------------------------------------------------
  ! Calculates the components of the total Hartree-Fock binding energy.
  !=====================================================================

  subroutine energy(nucleus, model, basis, HFfield, densities, sp_states, Coul_exch)

!
!   Arguments
!
    type(nucleus_properties),intent(inout)     :: nucleus
    type(modelization),intent(in)              :: model
    type(cylindrical_basis),intent(in)         :: basis
    type(field),intent(in)        :: HFfield
    type(local_density),dimension(2),intent(in):: densities
    type(sp_state),dimension(basis%size,2),intent(in) :: sp_states
    type(operator),dimension(basis%Nblocks),optional,intent(inout) :: Coul_exch

!
!   Local variables
!
    character(len=6),parameter :: subroutinename='energy'
    integer :: iHerm,iLag
    real(kind=rk) :: B1,B2,B3,B4,B5,B6,B7,B8,B9,B10,B11,B12,B13,B14,B15,B16,B17,&
         B18,B19,B20,B21,alpha,u0,v0
    real(kind=rk) :: Etot,Ekin,Evol,Esurf,Eso,Eodd,Ecoul_dir,Ecoul_exch,Ecoul,Epair,Epot
    real(kind=rk) :: Gz,Gr,weight,norm,r
    real(kind=rk) :: rho_n,rho_p,rho_tot
    real(kind=rk) :: tau_n,tau_p,tau_tot
    real(kind=rk) :: Delta_rho_n,Delta_rho_p,Delta_rho_tot
    real(kind=rk) :: div_J_n,div_J_p,div_J_tot
    real(kind=rk) :: rot_s_j,j2_tot,j2_n,j2_p,s2_tot,s2_n,s2_p
    real(kind=rk) :: sT_n,sT_p,sT_tot,sDelta_s_n,sDelta_s_p,sDelta_s_tot,sF_n,sF_p,sF_tot, &
         div_s_n,div_s_p,div_s_tot
    real(kind=rk) :: tensor_J2_n,tensor_J2_p,tensor_J2_tot,trace_J_n,trace_J_p,trace_J_tot, &
         crossed_J_n,crossed_J_p,crossed_J_tot
    type(vector) :: j_n,j_p,j_tot
    type(vector) :: s_n,s_p,s_tot,rot_s_n,rot_s_p,rot_s_tot, &
         Delta_s_n,Delta_s_p,Delta_s_tot
    type(vector) :: T_n,T_p,T_tot
    type(vector) :: F_n,F_p,F_tot
    integer :: iHF,iblock,n,shift,i,j
    real(kind=rk) :: Ci,Cj
    real(kind=rk) :: sum_e,E_B1,E_B2,E_B3,E_B4,E_B5,E_B6,E_B7,E_B8,E_B9,E_B10,E_B11, &
         E_B12,E_B13,E_B14,E_B15,E_B16,E_B17,E_B18,E_B19,E_B20,E_B21,E_u0,E_v0
    real(kind=rk) :: sum_jnjp, sum_snsp, sum_rho_snsp

!
!   Debugging note
!

    if (debugmode=='high') call debugnote(modulename, subroutinename)

!
!   Auxiliary variables
!

    B1=model%parameters%B1
    B2=model%parameters%B2
    B3=model%parameters%B3
    B4=model%parameters%B4
    B5=model%parameters%B5
    B6=model%parameters%B6
    B7=model%parameters%B7
    B8=model%parameters%B8
    B9=model%parameters%B9
    B10=model%parameters%B10
    B11=model%parameters%B11
    B12=model%parameters%B12
    B13=model%parameters%B13
    B14=model%parameters%B14
    B15=model%parameters%B15
    B16=model%parameters%B16
    B17=model%parameters%B17
    B18=model%parameters%B18
    B19=model%parameters%B19
    B20=model%parameters%B20
    B21=model%parameters%B21
    alpha=model%parameters%alpha
    u0=model%parameters%u0
    v0=model%parameters%v0

    norm=fz*pi/(bz*bp2)

!
!   Calculating components of Skyrme energy functional
!

    Ekin=0.0
    Evol=0.0
    Esurf=0.0
    Eso=0.0
    Eodd=0.0
    Ecoul_dir=0.0
    E_B1=0.0
    E_B2=0.0
    E_B3=0.0
    E_B4=0.0
    E_B5=0.0
    E_B6=0.0
    E_B7=0.0
    E_B8=0.0
    E_B9=0.0
    E_B10=0.0
    E_B11=0.0
    E_B12=0.0
    E_B13=0.0
    E_B14=0.0
    E_B15=0.0
    E_B16=0.0
    E_B17=0.0
    E_B18=0.0
    E_B19=0.0
    E_B20=0.0
    E_B21=0.0
    E_u0=0.0
    E_v0=0.0

    sum_jnjp=0.0
    sum_snsp=0.0

    do iHerm=iHmin,basis%NGz/2

       if (iHerm==0) cycle
       Gz=basis%G_Herm(iHerm)

       do iLag=1,basis%NGr

          Gr=basis%G_Lag(iLag)
          r=sqrt(basis%xLag(iLag))/bp
          weight=Gz*Gr

          rho_n=densities(1)%rho(iHerm,iLag)
          rho_p=densities(2)%rho(iHerm,iLag)
          rho_tot=rho_n+rho_p
          if (rho_tot<0) then
             rho_tot=0.
             call warning (modulename, subroutinename, &
                  'rho_tot found < 0 ==> set to 0')
          end if

          tau_n=densities(1)%tau(iHerm,iLag)
          tau_p=densities(2)%tau(iHerm,iLag)
          tau_tot=tau_n+tau_p

          Delta_rho_n=densities(1)%Delta_rho(iHerm,iLag)
          Delta_rho_p=densities(2)%Delta_rho(iHerm,iLag)
          Delta_rho_tot=Delta_rho_n+Delta_rho_p

          div_J_n=densities(1)%div_J(iHerm,iLag)
          div_J_p=densities(2)%div_J(iHerm,iLag)
          div_J_tot=div_J_n+div_J_p

          j_n=vector(densities(1)%j%z(iHerm,iLag),densities(1)%j%rho(iHerm,iLag), &
               densities(1)%j%phi(iHerm,iLag))
          j_p=vector(densities(2)%j%z(iHerm,iLag),densities(2)%j%rho(iHerm,iLag), &
               densities(2)%j%phi(iHerm,iLag))
          j_tot=j_n+j_p

          s_n=vector(densities(1)%s%z(iHerm,iLag),densities(1)%s%rho(iHerm,iLag), &
               densities(1)%s%phi(iHerm,iLag))
          s_p=vector(densities(2)%s%z(iHerm,iLag),densities(2)%s%rho(iHerm,iLag), &
               densities(2)%s%phi(iHerm,iLag))
          s_tot=s_n+s_p

          rot_s_n=vector(densities(1)%rot_s%z(iHerm,iLag),densities(1)%rot_s%rho(iHerm,iLag), &
               densities(1)%rot_s%phi(iHerm,iLag))
          rot_s_p=vector(densities(2)%rot_s%z(iHerm,iLag),densities(2)%rot_s%rho(iHerm,iLag), &
               densities(2)%rot_s%phi(iHerm,iLag))
          rot_s_tot=rot_s_n+rot_s_p

          Ekin=Ekin+weight*(hb0(1)*tau_n+hb0(2)*tau_p)

          Esurf=Esurf+weight*(B5*rho_tot*Delta_rho_tot+B6*(rho_n*Delta_rho_n+rho_p*Delta_rho_p))

          Eso=Eso+weight*B9*(rho_tot*div_J_tot+rho_n*div_J_n+rho_p*div_J_p)

          rot_s_j=real(scalar_product(rot_s_tot,j_tot)) + real(scalar_product(rot_s_n,j_n)) + &
               real(scalar_product(rot_s_p,j_p))

          j2_tot=real(scalar_product(j_tot,j_tot))
          j2_n=real(scalar_product(j_n,j_n))
          j2_p=real(scalar_product(j_p,j_p))

          s2_tot=real(scalar_product(s_tot,s_tot))
          s2_n=real(scalar_product(s_n,s_n))
          s2_p=real(scalar_product(s_p,s_p))

          T_n=vector(densities(1)%t%z(iHerm,iLag),densities(1)%t%rho(iHerm,iLag), &
               densities(1)%t%phi(iHerm,iLag))
          T_p=vector(densities(2)%t%z(iHerm,iLag),densities(2)%t%rho(iHerm,iLag), &
               densities(2)%t%phi(iHerm,iLag))
          T_tot=T_n+T_p
          sT_n=real(scalar_product(s_n,T_n))
          sT_p=real(scalar_product(s_p,T_p))
          sT_tot=real(scalar_product(s_tot,T_tot))

          Delta_s_n=vector(densities(1)%Delta_s%z(iHerm,iLag), &
               densities(1)%Delta_s%rho(iHerm,iLag), &
               densities(1)%Delta_s%phi(iHerm,iLag))
          Delta_s_p=vector(densities(2)%Delta_s%z(iHerm,iLag), &
               densities(2)%Delta_s%rho(iHerm,iLag), &
               densities(2)%Delta_s%phi(iHerm,iLag))
          Delta_s_tot=Delta_s_n+Delta_s_p
          sDelta_s_n=real(scalar_product(s_n,Delta_s_n))
          sDelta_s_p=real(scalar_product(s_p,Delta_s_p))
          sDelta_s_tot=real(scalar_product(s_tot,Delta_s_tot))

          tensor_J2_n=0.0
          do i=1,3
             tensor_J2_n = tensor_J2_n + densities(1)%tensor_J(i)%z(iHerm,iLag)**2 + &
                  densities(1)%tensor_J(i)%rho(iHerm,iLag)**2 + &
                  densities(1)%tensor_J(i)%phi(iHerm,iLag)**2
          end do
          
          tensor_J2_p=0.0
          do i=1,3
             tensor_J2_p = tensor_J2_p + densities(2)%tensor_J(i)%z(iHerm,iLag)**2 + &
                  densities(2)%tensor_J(i)%rho(iHerm,iLag)**2 + &
                  densities(2)%tensor_J(i)%phi(iHerm,iLag)**2
          end do

          tensor_J2_tot=0.0
          do i=1,3
             tensor_J2_tot = tensor_J2_tot + (densities(1)%tensor_J(i)%z(iHerm,iLag)+ &
                  densities(2)%tensor_J(i)%z(iHerm,iLag))**2 + &
                  (densities(1)%tensor_J(i)%rho(iHerm,iLag)+densities(2)%tensor_J(i)%rho(iHerm,iLag))**2 + &
                  (densities(1)%tensor_J(i)%phi(iHerm,iLag)+densities(2)%tensor_J(i)%phi(iHerm,iLag))**2
          end do

          trace_J_n = densities(1)%tensor_J(1)%z(iHerm,iLag) + densities(1)%tensor_J(2)%rho(iHerm,iLag) + &
               densities(1)%tensor_J(3)%phi(iHerm,iLag)
          trace_J_p = densities(2)%tensor_J(1)%z(iHerm,iLag) + densities(2)%tensor_J(2)%rho(iHerm,iLag) + &
               densities(2)%tensor_J(3)%phi(iHerm,iLag)
          trace_J_tot = trace_J_n + trace_J_p

          crossed_J_n = densities(1)%tensor_J(1)%z(iHerm,iLag)**2 + &
               densities(1)%tensor_J(1)%rho(iHerm,iLag)*densities(1)%tensor_J(2)%z(iHerm,iLag) + &
               densities(1)%tensor_J(1)%phi(iHerm,iLag)*densities(1)%tensor_J(3)%z(iHerm,iLag)
          crossed_J_n = crossed_J_n + &
               densities(1)%tensor_J(2)%z(iHerm,iLag)*densities(1)%tensor_J(1)%rho(iHerm,iLag) + &
               densities(1)%tensor_J(2)%rho(iHerm,iLag)**2 + &
               densities(1)%tensor_J(2)%phi(iHerm,iLag)*densities(1)%tensor_J(3)%rho(iHerm,iLag)
          crossed_J_n = crossed_J_n + &
               densities(1)%tensor_J(3)%z(iHerm,iLag)*densities(1)%tensor_J(1)%phi(iHerm,iLag) + &
               densities(1)%tensor_J(3)%rho(iHerm,iLag)*densities(1)%tensor_J(2)%phi(iHerm,iLag) + &
               densities(1)%tensor_J(3)%phi(iHerm,iLag)**2

          crossed_J_p = densities(2)%tensor_J(1)%z(iHerm,iLag)**2 + &
               densities(2)%tensor_J(1)%rho(iHerm,iLag)*densities(2)%tensor_J(2)%z(iHerm,iLag) + &
               densities(2)%tensor_J(1)%phi(iHerm,iLag)*densities(2)%tensor_J(3)%z(iHerm,iLag)
          crossed_J_p = crossed_J_p + &
               densities(2)%tensor_J(2)%z(iHerm,iLag)*densities(2)%tensor_J(1)%rho(iHerm,iLag) + &
               densities(2)%tensor_J(2)%rho(iHerm,iLag)**2 + &
               densities(2)%tensor_J(2)%phi(iHerm,iLag)*densities(2)%tensor_J(3)%rho(iHerm,iLag)
          crossed_J_p = crossed_J_p + &
               densities(2)%tensor_J(3)%z(iHerm,iLag)*densities(2)%tensor_J(1)%phi(iHerm,iLag) + &
               densities(2)%tensor_J(3)%rho(iHerm,iLag)*densities(2)%tensor_J(2)%phi(iHerm,iLag) + &
               densities(2)%tensor_J(3)%phi(iHerm,iLag)**2

          crossed_J_tot = &
               (densities(1)%tensor_J(1)%z(iHerm,iLag)  +densities(2)%tensor_J(1)%z(iHerm,iLag))**2 + &
               (densities(1)%tensor_J(1)%rho(iHerm,iLag)+densities(2)%tensor_J(1)%rho(iHerm,iLag)) * &
               (densities(1)%tensor_J(2)%z(iHerm,iLag)  +densities(2)%tensor_J(2)%z(iHerm,iLag)) + &
               (densities(1)%tensor_J(1)%phi(iHerm,iLag)+densities(2)%tensor_J(1)%phi(iHerm,iLag)) * &
               (densities(1)%tensor_J(3)%z(iHerm,iLag)  +densities(2)%tensor_J(3)%z(iHerm,iLag))
          crossed_J_tot = crossed_J_tot + &
               (densities(1)%tensor_J(2)%z(iHerm,iLag)  +densities(2)%tensor_J(2)%z(iHerm,iLag)) * &
               (densities(1)%tensor_J(1)%rho(iHerm,iLag)+densities(2)%tensor_J(1)%rho(iHerm,iLag)) + &
               (densities(1)%tensor_J(2)%rho(iHerm,iLag)+densities(2)%tensor_J(2)%rho(iHerm,iLag))**2 + &
               (densities(1)%tensor_J(2)%phi(iHerm,iLag)+densities(2)%tensor_J(2)%phi(iHerm,iLag)) * &
               (densities(1)%tensor_J(3)%rho(iHerm,iLag)+densities(2)%tensor_J(3)%rho(iHerm,iLag))
          crossed_J_tot = crossed_J_tot + &
               (densities(1)%tensor_J(3)%z(iHerm,iLag)  +densities(2)%tensor_J(3)%z(iHerm,iLag)) * &
               (densities(1)%tensor_J(1)%phi(iHerm,iLag)+densities(2)%tensor_J(1)%phi(iHerm,iLag)) + &
               (densities(1)%tensor_J(3)%rho(iHerm,iLag)+densities(2)%tensor_J(3)%rho(iHerm,iLag)) * &
               (densities(1)%tensor_J(2)%phi(iHerm,iLag)+densities(2)%tensor_J(2)%phi(iHerm,iLag)) + &
               (densities(1)%tensor_J(3)%phi(iHerm,iLag)+densities(2)%tensor_J(3)%phi(iHerm,iLag))**2

          F_n=vector(densities(1)%F%z(iHerm,iLag),densities(1)%F%rho(iHerm,iLag), &
               densities(1)%F%phi(iHerm,iLag))
          F_p=vector(densities(2)%F%z(iHerm,iLag),densities(2)%F%rho(iHerm,iLag), &
               densities(2)%F%phi(iHerm,iLag))
          F_tot=F_n+F_p

          sF_n=real(scalar_product(s_n,F_n))
          sF_p=real(scalar_product(s_p,F_p))
          sF_tot=real(scalar_product(s_tot,F_tot))

          div_s_n=densities(1)%div_s(iHerm,iLag)
          div_s_p=densities(2)%div_s(iHerm,iLag)
          div_s_tot=div_s_n+div_s_p

          Evol=Evol+weight*(B1*rho_tot**2 + B2*(rho_n**2+rho_p**2) + B3*rho_tot*tau_tot + &
               B4*(rho_n*tau_n+rho_p*tau_p) + rho_tot**alpha*(B7*rho_tot**2 + &
               B8*(rho_n**2+rho_p**2)) + B14*(tensor_J2_tot-st_tot) + &
               B15*(tensor_J2_n+tensor_J2_p-st_n-st_p) + B16*(trace_J_tot**2 + crossed_J_tot) + &
               B17*(trace_J_n**2 + crossed_J_n + trace_J_p**2 + crossed_J_p))

          Eodd=Eodd+weight*(B9*rot_s_j - B3*j2_tot - B4*(j2_n+j2_p) + rho_tot**alpha * &
               (B12*s2_tot+B13*(s2_n+s2_p)) + B10*s2_tot + B11*(s2_n+s2_p) + &
               B18*sDelta_s_tot + B19*(sDelta_s_n+sDelta_s_p) - 2.0_rk * &
               (B16*sF_tot + B17*(sF_n +sF_p)) + &
               B20*div_s_tot**2 + B21*(div_s_n**2+div_s_p**2))

          E_B1=E_B1+weight*B1*rho_tot**2
          E_B2=E_B2+weight*B2*(rho_n**2+rho_p**2)
          E_B3=E_B3+weight*B3*(rho_tot*tau_tot - j2_tot)
          E_B4=E_B4+weight*B4*(rho_n*tau_n - j2_n + rho_p*tau_p - j2_p)
          E_B5=E_B5+weight*B5*rho_tot*Delta_rho_tot
          E_B6=E_B6+weight*B6*(rho_n*Delta_rho_n+rho_p*Delta_rho_p)
          E_B7=E_B7+weight*B7*rho_tot**alpha*rho_tot**2
          E_B8=E_B8+weight*B8*rho_tot**alpha*(rho_n**2+rho_p**2)
          E_B9=E_B9+weight*B9*(rho_tot*div_J_tot+rho_n*div_J_n+rho_p*div_J_p+rot_s_j)
          E_B10=E_B10+weight*B10*s2_tot
          E_B11=E_B11+weight*B11*(s2_n+s2_p)
          E_B12=E_B12+weight*B12*rho_tot**alpha*s2_tot
          E_B13=E_B13+weight*B13*rho_tot**alpha*(s2_n+s2_p)
          E_B14=E_B14+weight*B14*(tensor_J2_tot-st_tot)
          E_B15=E_B15+weight*B15*(tensor_J2_n+tensor_J2_p-st_n-st_p)
          E_B16=E_B16+weight*B16*(trace_J_tot**2 + crossed_J_tot - 2.0_rk*sF_tot)
          E_B17=E_B17+weight*B17*(trace_J_n**2 + crossed_J_n - 2.0_rk*sF_n + &
               trace_J_p**2 + crossed_J_p - 2.0_rk*sF_p)
          E_B18=E_B18+weight*B18*sDelta_s_tot
          E_B19=E_B19+weight*B19*(sDelta_s_n+sDelta_s_p)
          E_B20=E_B20+weight*B20*div_s_tot**2
          E_B21=E_B21+weight*B21*(div_s_n**2+div_s_p**2)
          E_u0=E_u0+weight*0.75_rk*u0*((rho_n**2-s2_n)*rho_p+(rho_p**2-s2_p)*rho_n)
          E_v0=E_v0+weight*0.75_rk*v0*(rho_n**2-s2_n)*(rho_p**2-s2_p)

          Ecoul_dir=Ecoul_dir+weight*0.5_rk*rho_p*HFfield%Coulomb(iHerm,iLag)

          sum_jnjp=sum_jnjp+weight*real(scalar_product(j_n,j_p))
          sum_snsp=sum_snsp+weight*real(scalar_product(s_n,s_p))
          sum_rho_snsp=sum_rho_snsp+weight*rho_tot**alpha*real(scalar_product(s_n,s_p))
          
       end do

    end do

    Ekin=norm*Ekin*(1.0_rk-c_com*1.0_rk/nucleus%A) ! One-body center of mass correction
    Evol=norm*Evol
    Esurf=norm*Esurf
    Eso=norm*Eso
    Eodd=norm*Eodd
    Ecoul_dir=norm*Ecoul_dir*c_Coulomb
    Epair=sum(nucleus%pairing(:)%pairing_energy)

    E_B1=norm*E_B1
    E_B2=norm*E_B2
    E_B3=norm*E_B3
    E_B4=norm*E_B4
    E_B5=norm*E_B5
    E_B6=norm*E_B6
    E_B7=norm*E_B7
    E_B8=norm*E_B8
    E_B9=norm*E_B9
    E_B10=norm*E_B10
    E_B11=norm*E_B11
    E_B12=norm*E_B12
    E_B13=norm*E_B13
    E_B14=norm*E_B14
    E_B15=norm*E_B15
    E_B16=norm*E_B16
    E_B17=norm*E_B17
    E_B18=norm*E_B18
    E_B19=norm*E_B19
    E_B20=norm*E_B20
    E_B21=norm*E_B21
    E_u0=norm*E_u0
    E_v0=norm*E_v0

    sum_jnjp=sum_jnjp*norm
    sum_snsp=sum_snsp*norm
    sum_rho_snsp=sum_rho_snsp*norm

    Epot = nucleus%energy%Epot

!
!   Exchange Coulomb energy
!

    Ecoul_exch=0.0

    if (PRESENT(Coul_exch) .and. abs(c_Coulomb) > 0.0_rk) then

!      Exact calculation

       do iHF=1,basis%size
          if (abs(sp_states(iHF,2)%occupation) < 1e-6_rk) cycle
          iblock=sp_states(iHF,2)%block
          n=basis%block_size(iblock)
          shift=0
          if (iblock>1) shift=sum(basis%block_size(1:iblock-1))
          do i=1,n
             Ci=sp_states(iHF,2)%coef%cylHO(i+shift)
             do j=1,n
                Cj=sp_states(iHF,2)%coef%cylHO(j+shift)
                Ecoul_exch=Ecoul_exch+0.5_rk*sp_states(iHF,2)%occupation*Coul_exch(iblock)%HO(i,j)*Ci*Cj
             end do
          end do
       end do

    else if (abs(c_Slater) > 0.0_rk .and. abs(c_Coulomb) > 0.0_rk) then

!      Slater approximation

       do iHerm=iHmin,basis%NGz/2
          if (iHerm==0) cycle
          Gz=basis%G_Herm(iHerm)
          do iLag=1,basis%NGr
             Gr=basis%G_Lag(iLag)
             weight=Gz*Gr
             rho_p=densities(2)%rho(iHerm,iLag)
             if (rho_p<0) then
                rho_p=0.
                call warning (modulename, subroutinename, 'rho_p found < 0 ==> set to 0')
             end if
             Ecoul_exch = Ecoul_exch + weight*rho_p**(4.0_rk/3.0_rk)
          end do
       end do

       Ecoul_exch = -norm*Ecoul_exch*0.75_rk*e2*(3.0_rk/pi)**(1.0_rk/3.0_rk)

    else

       Ecoul_exch = 0.0_rk

    end if

!
!   Update nucleus energy components
!

    Ecoul=Ecoul_dir+Ecoul_exch
    Etot=Ekin+Evol+Esurf+Eso+Eodd+Ecoul+Epair+E_u0+E_v0
    nucleus%energy=energy_properties(Etot,Ekin,Evol,Esurf,Eso,Eodd,Ecoul,Epair,&
         (/E_B1,E_B2,E_B3,E_B4,E_B5,E_B6,E_B7,E_B8,E_B9,E_B10,E_B11,E_B12,E_B13, &
         E_B14,E_B15,E_B16,E_B17,E_B18,E_B19,E_B20,E_B21,E_u0,E_v0/),Epot)

!   Test cancellations of separate contributions to Etot when A=1 (no self-interaction condition)
!    write(out_unit,'(/4(a,f13.6,2x))') 'E_B1 =',E_B1,'E_B2 =',E_B2,'E_B3 =',E_B3,'E_B4 =',E_B4
!    write(out_unit,'(3(a,f13.6,2x))') 'E_B5 =',E_B5,'E_B6 =',E_B6,'E_B9 =',E_B9
!    write(out_unit,'(3(4(a,f13.6,2x)/),2(a,f13.6,2x)/)') &
!         'E_B10=',E_B10,'E_B11=',E_B11,'E_B12=',E_B12,'E_B13=',E_B13, &
!         'E_B14=',E_B14,'E_B15=',E_B15,'E_B16=',E_B16,'E_B17=',E_B17, &
!         'E_B18=',E_B18,'E_B19=',E_B19,'E_B20=',E_B20,'E_B21=',E_B21, &
!         'E_u0=',E_u0,'E_v0=',E_v0
!    write(out_unit,'(1(a,f13.6,2x)/)') 'E_B1 + ... + E_B21 + E_u0 + E_v0 =',&
!         E_B1+E_B2+E_B3+E_B4+E_B5+E_B6+E_B9+E_B10+E_B11+E_B12+E_B13+E_B14+&
!         E_B15+E_B16+E_B17+E_B18+E_B19+E_B20+E_B21+E_u0+E_v0

    if (nucleus%A == 1) then
       if (B1/=0 .and. B10/=0) write(out_unit,'(1(a,f10.6,2x))') 'Int(rho^2-s^2)=',E_B1/B1-E_B10/B10
       if (B3/=0 .and. B5/=0) write(out_unit,'(2(a,f10.6,2x))') &
            'Int(rho.tau-j^2)=',E_B3/B3,'-1/4*Int(rho.Delta rho)=',-0.25_rk*E_B5/B5
       if (B15/=0 .and. B19/=0) write(out_unit,'(2(a,f10.6,2x))') &
            'Int(s.T-J^2)=',-E_B15/B15,'-1/4*Int(s.Delta s)=',-0.25_rk*E_B19/B19
       if (B9/=0) write(out_unit,'(1(a,f10.6,2x))') 'Int(rho.div(J)+s.rot(j))=',E_B9/B9
       if (B16/=0 .and. B20/=0) write(out_unit,'(2(a,f10.6,2x))') &
            'Int(s.F-1/2*[(Tr J)^2+Tr(J.J)])=',-0.5_rk*E_B16/B16, &
            '1/4*Int((div s)^2)=',0.25_rk*E_B20/B20
       write(out_unit,*)
    end if

    if (alpha==0.0 .or. B7 == 0.0) then
       sum_e=0.0
       do i=1,2
          do iHF=1,basis%size
             sum_e = sum_e + sp_states(iHF,i)%occupation*sp_states(iHF,i)%energy
          end do
       end do
       write(out_unit,'(a/,2(1x,a,f13.6)/)') 'Set t3=0 and no Coulomb or exact Coulomb:',&
            'sum_i e_i = ',sum_e,' 2*Etot - Ekin =',2*Etot - Ekin
    end if

    if (mod(nucleus%N,2) == 1 .and. mod(nucleus%Z,2) == 1) then
       write(out_unit,'(/a/)') '--------------------------------------------------------'
       write(out_unit,'(a/)') ' GM splitting analysis (valid for fixed s.p. spectrum):'
       write(out_unit,'(a/)') '--------------------------------------------------------'
       write(out_unit,'(a,f10.2,1x,a)') '-4*B3*int(jn.jp) =',-4*B3*sum_jnjp*1e3_rk,'keV'
       write(out_unit,'(a,f10.2,1x,a)') '4*B10*int(sn.sp) =',4*B10*sum_snsp*1e3_rk,'keV'
       write(out_unit,'(a,f10.2,1x,a)') '4*B12*int(rho^alpha*sn.sp) =',4*B12*sum_rho_snsp*1e3_rk,'keV'
    end if

    return

  end subroutine energy


  !=====================================================================
  !                       SUBROUTINE POTENTIAL
  !---------------------------------------------------------------------
  ! Calculates the components of the Hartree-Fock potential: effective
  ! mass, central, spin-orbit, spin, current, tensor and odd.
  !=====================================================================

  subroutine potential(nucleus, model, symmetries, basis, numerical_parameters, &
       HFfield, densities, exact_coul_exch)

!
!   Arguments
!
    type(nucleus_properties),intent(in)            :: nucleus
    type(modelization),intent(in)                  :: model
    type(symmetry_properties),intent(in)           :: symmetries
    type(cylindrical_basis),intent(in)             :: basis
    type(numeric),intent(in)                       :: numerical_parameters
    type(field),dimension(2),intent(inout)         :: HFfield
    type(local_density),dimension(2),intent(inout) :: densities
    logical,intent(in)                             :: exact_coul_exch

!
!   Local variables
!
    character(len=9),parameter :: subroutinename='potential'
    real(kind=rk) :: B1,B2,B3,B4,B5,B6,B7,B8,B9,B10,B11,B12,B13,B14,B15,B16,B17,&
         B18,B19,B20,B21,alpha,u0,v0
    integer :: isospin,iHerm,iLag,k
    real(kind=rk) :: rho_q,rho_n,rho_p,rho_tot,rho_a, &
         tau_q,tau_tot,Delta_rho_q,Delta_rho_tot,div_J_q,div_J_tot, &
         s_z_tot,s_rho_tot,s_phi_tot
    real(kind=rk) :: j_tot_z,j_tot_rho,j_tot_phi,rot_j_tot_z,rot_j_tot_rho,rot_j_tot_phi, &
         s_tot_z,s_tot_rho,s_tot_phi,rot_s_tot_z,rot_s_tot_rho,rot_s_tot_phi
    real(kind=rk) :: t_tot_z,t_tot_rho,t_tot_phi,Delta_s_tot_z,Delta_s_tot_rho,Delta_s_tot_phi, &
         F_tot_z,F_tot_rho,F_tot_phi,trace_J_q,trace_J_tot
    real(kind=rk) :: rho_qbar,s2_qbar
    type(vector) :: s_qbar
    real(kind=rk) :: Vcentr,Vconstr,z,r
    real(kind=rk),dimension(Nconstr) :: stiffness,Qtarget,Qvalue
    real(kind=rk) :: mix
    integer,parameter :: field_unit = 15
    real(kind=rk) :: norm, Gz, weight
    real(kind=rk),dimension(2) :: integral_spin_field, integral_s_contribution, &
         integral_T_contribution, integral_so_contribution, integral_Delta_s_contribution

!
!   Debugging note
!

    if (debugmode=='medium'.or.debugmode=='high') call debugnote(modulename, subroutinename)

!
!   Auxiliary variables
!

    B1=model%parameters%B1
    B2=model%parameters%B2
    B3=model%parameters%B3
    B4=model%parameters%B4
    B5=model%parameters%B5
    B6=model%parameters%B6
    B7=model%parameters%B7
    B8=model%parameters%B8
    B9=model%parameters%B9
    B10=model%parameters%B10
    B11=model%parameters%B11
    B12=model%parameters%B12
    B13=model%parameters%B13
    B14=model%parameters%B14
    B15=model%parameters%B15
    B16=model%parameters%B16
    B17=model%parameters%B17
    B18=model%parameters%B18
    B19=model%parameters%B19
    B20=model%parameters%B20
    B21=model%parameters%B21
    alpha=model%parameters%alpha
    u0=model%parameters%u0
    v0=model%parameters%v0

    Qvalue(1)=nucleus%zcm(0)
    Qvalue(2)=nucleus%Q20(0)
    Qvalue(3)=nucleus%Q30(0)
    Qvalue(4)=nucleus%Q40(0)
    Qvalue(5)=nucleus%QN(0)
    Qvalue(6)=nucleus%r2(0)

    Qtarget(:)=model%constraints(:)%value

    stiffness(:) = model%constraints(:)%stiffness
    do k = 1, Nconstr
       if (trim(model%constraints(k)%readj) == 'n') then
          stiffness(k) = stiffness(k) * (Qvalue(k) - Qtarget(k)) ! quadratic constraint
       else
          stiffness(k) = -stiffness(k) ! linear constraint
       end if
    end do
    if (symmetries%parity) then
       stiffness(1)=0.0 ; Qtarget(1)=0.0
       stiffness(3)=0.0 ; Qtarget(3)=0.0
    else
       stiffness(1) = stiffness(1)/nucleus%A ! constraining function should be z/A, but argument A missing in fconstr
    end if

    mix=numerical_parameters%convergence%mixing

!
!   Calculation of the components of the Skyrme-Hartree-Fock potential:
!      * central potential
!      * inverse effective mass (hbar^2/2m*)
!      * spin-orbit potential
!      * spin potential S
!      * current : alpha
!      * tensor
!      * C
!

    if (debugmode=='high') write(nqp,'(13(A,F11.5,1X))') 'B1=',B1,' B2=',B2,' B3=',B3, &
         ' B4=',B4,' B7=',B7,' B8=',B8,' B9=',B9,' B10=',B10,' B11=',B11,' B12=',B12, &
         ' B13=',B13

    open(field_unit,file='spin_field.txt',status='replace',action='write')

    norm=fz*pi/(bz*bp2)
    integral_spin_field(:) = 0.0
    integral_s_contribution(:) = 0.0
    integral_T_contribution(:) = 0.0
    integral_so_contribution(:) = 0.0
    integral_Delta_s_contribution(:) = 0.0

    do isospin=1,2

       do iHerm=iHmin,basis%NGz/2

          if (iHerm==0) cycle
          z=basis%xHerm(iHerm)/bz
          Gz=basis%G_Herm(iHerm)

          do iLag=1,basis%NGr

             weight=Gz*basis%G_Lag(iLag)

             r=sqrt(basis%xLag(iLag))/bp

             rho_q=densities(isospin)%rho(iHerm,iLag)
             rho_n=densities(1)%rho(iHerm,iLag)
             rho_p=densities(2)%rho(iHerm,iLag)
             if (rho_p < 0.0) then
                rho_p=0.0
                call warning(modulename, subroutinename, &
                     'rho_p found < 0 ==> rho_p set to 0')
             end if
             rho_tot=rho_n+rho_p
             if (rho_tot < 0.0) then
                rho_a=0.0
                call warning(modulename, subroutinename, &
                     'rho_tot found < 0 ==> rho_tot set to 0')
             else
                rho_a=rho_tot**(alpha-1.0_rk)
             end if
             tau_q=densities(isospin)%tau(iHerm,iLag)
             tau_tot=densities(1)%tau(iHerm,iLag)+densities(2)%tau(iHerm,iLag)
             Delta_rho_q=densities(isospin)%Delta_rho(iHerm,iLag)
             Delta_rho_tot=densities(1)%Delta_rho(iHerm,iLag)+densities(2)%Delta_rho(iHerm,iLag)
             div_J_q=densities(isospin)%div_J(iHerm,iLag)
             div_J_tot=densities(1)%div_J(iHerm,iLag)+densities(2)%div_J(iHerm,iLag)
             s_tot_z=densities(1)%s%z(iHerm,iLag)+densities(2)%s%z(iHerm,iLag)
             s_tot_rho=densities(1)%s%rho(iHerm,iLag)+densities(2)%s%rho(iHerm,iLag)
             s_tot_phi=densities(1)%s%phi(iHerm,iLag)+densities(2)%s%phi(iHerm,iLag)

             rho_qbar=densities(3-isospin)%rho(iHerm,iLag)
             s_qbar=vector(densities(3-isospin)%s%z(iHerm,iLag),densities(3-isospin)%s%rho(iHerm,iLag), &
                  densities(3-isospin)%s%phi(iHerm,iLag))
             s2_qbar=real(scalar_product(s_qbar,s_qbar))

             Vcentr=2.0_rk*(B1*rho_tot+B2*rho_q)+B3*tau_tot+B4*tau_q+ &
                  2.0_rk*(B5*Delta_rho_tot+B6*Delta_rho_q)+ &
                  B7*(2.0_rk+alpha)*(rho_a*rho_tot**2)+ &
                  B8*rho_a*(alpha*(rho_n**2+rho_p**2)+2.0_rk*rho_tot*rho_q)+ &
                  B9*(div_J_tot+div_J_q)+alpha*rho_a*(B12*(s_tot_z**2+s_tot_rho**2+s_tot_phi**2)+ &
                  B13*(densities(1)%s%z(iHerm,iLag)**2+densities(1)%s%rho(iHerm,iLag)**2+ &
                  densities(1)%s%phi(iHerm,iLag)**2+densities(2)%s%z(iHerm,iLag)**2+ &
                  densities(2)%s%rho(iHerm,iLag)**2+densities(2)%s%phi(iHerm,iLag)**2))+ &
                  1.5_rk*u0*rho_q*rho_qbar+0.75_rk*(u0+2*v0*rho_q)*(rho_qbar**2-s2_qbar)

             Vconstr=0.0
             do k=1,Nconstr
                Vconstr=Vconstr+stiffness(k)*fconstr(k,z,r) ! Valid for linear as well as quadratic constraints
             end do

             Vcentr=Vcentr+Vconstr

             if (isospin==2) then
                Vcentr=Vcentr+HFfield(2)%Coulomb(iHerm,iLag)
                if (.not.exact_coul_exch) Vcentr=Vcentr-e2*(3.0_rk/pi*rho_p)**(1.0_rk/3.0_rk)*c_Slater*c_Coulomb
             end if

             HFfield(isospin)%central(iHerm,iLag)=mix*HFfield(isospin)%central(iHerm,iLag)+ &
                  (1.0_rk-mix)*Vcentr

             HFfield(isospin)%effective_mass(iHerm,iLag)=mix*HFfield(isospin)%effective_mass(iHerm,iLag)+ &
                  (1.0_rk-mix)*(hb0(isospin)*(1.0_rk-c_com*1.0_rk/nucleus%A)+B3*rho_tot+B4*rho_q)

             HFfield(isospin)%spin_orbit(iHerm,iLag)=mix*HFfield(isospin)%spin_orbit(iHerm,iLag)- &
                  (1.0_rk-mix)*B9*(rho_tot+rho_q)

             rot_j_tot_z=densities(1)%rot_j%z(iHerm,iLag)+densities(2)%rot_j%z(iHerm,iLag)
             rot_j_tot_rho=densities(1)%rot_j%rho(iHerm,iLag)+densities(2)%rot_j%rho(iHerm,iLag)
             rot_j_tot_phi=densities(1)%rot_j%phi(iHerm,iLag)+densities(2)%rot_j%phi(iHerm,iLag)
             t_tot_z=densities(1)%t%z(iHerm,iLag)+densities(2)%t%z(iHerm,iLag)
             t_tot_rho=densities(1)%t%rho(iHerm,iLag)+densities(2)%t%rho(iHerm,iLag)
             t_tot_phi=densities(1)%t%phi(iHerm,iLag)+densities(2)%t%phi(iHerm,iLag)
             Delta_s_tot_z=densities(1)%Delta_s%z(iHerm,iLag)+densities(2)%Delta_s%z(iHerm,iLag)
             Delta_s_tot_rho=densities(1)%Delta_s%rho(iHerm,iLag)+densities(2)%Delta_s%rho(iHerm,iLag)
             Delta_s_tot_phi=densities(1)%Delta_s%phi(iHerm,iLag)+densities(2)%Delta_s%phi(iHerm,iLag)
             F_tot_z=densities(1)%F%z(iHerm,iLag)+densities(2)%F%z(iHerm,iLag)
             F_tot_rho=densities(1)%F%rho(iHerm,iLag)+densities(2)%F%rho(iHerm,iLag)
             F_tot_phi=densities(1)%F%phi(iHerm,iLag)+densities(2)%F%phi(iHerm,iLag)

             HFfield(isospin)%S%z(iHerm,iLag)=mix*HFfield(isospin)%S%z(iHerm,iLag) + &
                  (1.0_rk-mix)*(B9*(rot_j_tot_z+densities(isospin)%rot_j%z(iHerm,iLag)) + &
                  2.0_rk*(B10+B12*rho_tot**alpha)*s_tot_z + &
                  (2.0_rk*(B11+B13*rho_tot**alpha) - 1.5_rk*(u0*rho_qbar+v0*(rho_qbar**2-s2_qbar)))* &
                  densities(isospin)%s%z(iHerm,iLag) - &
                  (B14*t_tot_z + B15*densities(isospin)%t%z(iHerm,iLag)) + &
                  2.0_rk*(B18*Delta_s_tot_z + B19*densities(isospin)%Delta_s%z(iHerm,iLag)) - &
                  2.0_rk*(B16*F_tot_z + B17*densities(isospin)%F%z(iHerm,iLag)))

             integral_spin_field(isospin) = integral_spin_field(isospin) + weight * &
                  HFfield(isospin)%S%z(iHerm,iLag) * rho_q
             integral_s_contribution(isospin) = integral_s_contribution(isospin) + weight * &
                  (2.0_rk*(B10+B12*rho_tot**alpha)*s_tot_z + 2.0_rk*(B11+B13*rho_tot**alpha) * &
                  densities(isospin)%s%z(iHerm,iLag)) * rho_q
             integral_T_contribution(isospin) = integral_T_contribution(isospin) - weight * &
                  (B14*t_tot_z + B15*densities(isospin)%t%z(iHerm,iLag)) * rho_q
             integral_so_contribution(isospin) = integral_so_contribution(isospin) + weight * &
                  B9 * (rot_j_tot_z + densities(isospin)%rot_j%z(iHerm,iLag))*rho_q
             integral_Delta_s_contribution(isospin) = integral_Delta_s_contribution(isospin) + &
                  weight*2.0_rk*(B18*Delta_s_tot_z + B19*densities(isospin)%Delta_s%z(iHerm,iLag))*rho_q

             write(field_unit,'(a1,6(1x,f10.6))') &
                  nucleon_type(isospin),z,r,HFfield(isospin)%S%z(iHerm,iLag), &
                  B9*(rot_j_tot_z+densities(isospin)%rot_j%z(iHerm,iLag)), &
                  2.0_rk*(B10+B12*rho_tot**alpha)*s_tot_z + &
                  2.0_rk*(B11+B13*rho_tot**alpha)*densities(isospin)%s%z(iHerm,iLag), &
                  -(B14*t_tot_z + B15*densities(isospin)%t%z(iHerm,iLag))

             HFfield(isospin)%S%rho(iHerm,iLag)=mix*HFfield(isospin)%S%rho(iHerm,iLag) + &
                  (1.0_rk-mix)*(B9*(rot_j_tot_rho+densities(isospin)%rot_j%rho(iHerm,iLag)) + &
                  2.0_rk*(B10+B12*rho_tot**alpha)*s_tot_rho + &
                  (2.0_rk*(B11+B13*rho_tot**alpha) - &
                  1.5_rk*(u0*rho_qbar+v0*(rho_qbar**2-s2_qbar)))*densities(isospin)%s%rho(iHerm,iLag) - &
                  B14*t_tot_rho - B15*densities(isospin)%t%rho(iHerm,iLag) + &
                  2.0_rk*(B18*Delta_s_tot_rho + B19*densities(isospin)%Delta_s%rho(iHerm,iLag)) - &
                  2.0_rk*(B16*F_tot_rho + B17*densities(isospin)%F%rho(iHerm,iLag)))

             HFfield(isospin)%S%phi(iHerm,iLag)=mix*HFfield(isospin)%S%phi(iHerm,iLag) + &
                  (1.0_rk-mix)*(B9*(rot_j_tot_phi+densities(isospin)%rot_j%phi(iHerm,iLag)) + &
                  2.0_rk*(B10+B12*rho_tot**alpha)*s_tot_phi + &
                  (2.0_rk*(B11+B13*rho_tot**alpha) - &
                  1.5_rk*(u0*rho_qbar+v0*(rho_qbar**2-s2_qbar)))*densities(isospin)%s%phi(iHerm,iLag) - &
                  B14*t_tot_phi - B15*densities(isospin)%t%phi(iHerm,iLag) + &
                  2.0_rk*(B18*Delta_s_tot_phi + B19*densities(isospin)%Delta_s%phi(iHerm,iLag)) - &
                  2.0_rk*(B16*F_tot_phi + B17*densities(isospin)%F%phi(iHerm,iLag)))

             HFfield(isospin)%div_s(iHerm,iLag)=mix*HFfield(isospin)%div_s(iHerm,iLag) + &
                  (1.0_rk-mix)*2.0_rk*(B21*densities(isospin)%div_s(iHerm,iLag) + &
                  B20*(densities(1)%div_s(iHerm,iLag)+densities(2)%div_s(iHerm,iLag)))

             j_tot_z=densities(1)%j%z(iHerm,iLag)+densities(2)%j%z(iHerm,iLag)
             j_tot_rho=densities(1)%j%rho(iHerm,iLag)+densities(2)%j%rho(iHerm,iLag)
             j_tot_phi=densities(1)%j%phi(iHerm,iLag)+densities(2)%j%phi(iHerm,iLag)
             rot_s_tot_z=densities(1)%rot_s%z(iHerm,iLag)+densities(2)%rot_s%z(iHerm,iLag)
             rot_s_tot_rho=densities(1)%rot_s%rho(iHerm,iLag)+densities(2)%rot_s%rho(iHerm,iLag)
             rot_s_tot_phi=densities(1)%rot_s%phi(iHerm,iLag)+densities(2)%rot_s%phi(iHerm,iLag)

             HFfield(isospin)%alpha%z(iHerm,iLag)=mix*HFfield(isospin)%alpha%z(iHerm,iLag) + &
                  (1.0_rk-mix)*(-2.0_rk*(B3*j_tot_z+B4*densities(isospin)%j%z(iHerm,iLag)) + &
                  B9*(rot_s_tot_z+densities(isospin)%rot_s%z(iHerm,iLag)))

             HFfield(isospin)%alpha%rho(iHerm,iLag)=mix*HFfield(isospin)%alpha%rho(iHerm,iLag) + &
                  (1.0_rk-mix)*(-2.0_rk*(B3*j_tot_rho+B4*densities(isospin)%j%rho(iHerm,iLag)) + &
                  B9*(rot_s_tot_rho+densities(isospin)%rot_s%rho(iHerm,iLag)))

             HFfield(isospin)%alpha%phi(iHerm,iLag)=mix*HFfield(isospin)%alpha%phi(iHerm,iLag) + &
                  (1.0_rk-mix)*(-2.0_rk*(B3*j_tot_phi+B4*densities(isospin)%j%phi(iHerm,iLag)) + &
                  B9*(rot_s_tot_phi+densities(isospin)%rot_s%phi(iHerm,iLag)))

             HFfield(isospin)%C%z(iHerm,iLag)=mix*HFfield(isospin)%C%z(iHerm,iLag) - &
                  (1.0_rk-mix)*(B14*s_tot_z+B15*densities(isospin)%s%z(iHerm,iLag))
             HFfield(isospin)%C%rho(iHerm,iLag)=mix*HFfield(isospin)%C%rho(iHerm,iLag) - &
                  (1.0_rk-mix)*(B14*s_tot_rho+B15*densities(isospin)%s%rho(iHerm,iLag))
             HFfield(isospin)%C%phi(iHerm,iLag)=mix*HFfield(isospin)%C%phi(iHerm,iLag) - &
                  (1.0_rk-mix)*(B14*s_tot_phi+B15*densities(isospin)%s%phi(iHerm,iLag))

             trace_J_q=densities(isospin)%tensor_J(1)%z(iHerm,iLag) + &
                  densities(isospin)%tensor_J(2)%rho(iHerm,iLag) + &
                  densities(isospin)%tensor_J(3)%phi(iHerm,iLag)
             trace_J_tot=trace_J_q+densities(3-isospin)%tensor_J(1)%z(iHerm,iLag) + &
                  densities(3-isospin)%tensor_J(2)%rho(iHerm,iLag) + &
                  densities(3-isospin)%tensor_J(3)%phi(iHerm,iLag)

!            W_q,z z
             HFfield(isospin)%tensor_so(1)%z(iHerm,iLag)= &
                     mix*HFfield(isospin)%tensor_so(1)%z(iHerm,iLag)+(1.0_rk-mix)*( &
                     B14*(densities(1)%tensor_J(1)%z(iHerm,iLag)+densities(2)%tensor_J(1)%z(iHerm,iLag)) + &
                     B15*densities(isospin)%tensor_J(1)%z(iHerm,iLag) + &
                     B16*(densities(1)%tensor_J(1)%z(iHerm,iLag)+densities(2)%tensor_J(1)%z(iHerm,iLag) + &
                     trace_J_tot) + &
                     B17*(densities(isospin)%tensor_J(1)%z(iHerm,iLag) + trace_J_q) )
!            W_q,z rho
             HFfield(isospin)%tensor_so(1)%rho(iHerm,iLag)= &
                     mix*HFfield(isospin)%tensor_so(1)%rho(iHerm,iLag)+(1.0_rk-mix)*( &
                     B14*(densities(1)%tensor_J(1)%rho(iHerm,iLag)+densities(2)%tensor_J(1)%rho(iHerm,iLag)) + &
                     B15*densities(isospin)%tensor_J(1)%rho(iHerm,iLag) + &
                     B16*(densities(1)%tensor_J(2)%z(iHerm,iLag) + &
                     densities(2)%tensor_J(2)%z(iHerm,iLag)) + &
                     B17*densities(isospin)%tensor_J(2)%z(iHerm,iLag) )
!            W_q,z phi
             HFfield(isospin)%tensor_so(1)%phi(iHerm,iLag)= &
                     mix*HFfield(isospin)%tensor_so(1)%phi(iHerm,iLag)+(1.0_rk-mix)*( &
                     B14*(densities(1)%tensor_J(1)%phi(iHerm,iLag)+densities(2)%tensor_J(1)%phi(iHerm,iLag)) + &
                     B15*densities(isospin)%tensor_J(1)%phi(iHerm,iLag) + &
                     B16*(densities(1)%tensor_J(3)%z(iHerm,iLag) + &
                     densities(2)%tensor_J(3)%z(iHerm,iLag)) + &
                     B17*densities(isospin)%tensor_J(3)%z(iHerm,iLag) )

!            W_q,rho z
             HFfield(isospin)%tensor_so(2)%z(iHerm,iLag)= &
                     mix*HFfield(isospin)%tensor_so(2)%z(iHerm,iLag)+(1.0_rk-mix)*( &
                     B14*(densities(1)%tensor_J(2)%z(iHerm,iLag)+densities(2)%tensor_J(2)%z(iHerm,iLag)) + &
                     B15*densities(isospin)%tensor_J(2)%z(iHerm,iLag) + &
                     B16*(densities(1)%tensor_J(1)%rho(iHerm,iLag)+densities(2)%tensor_J(1)%rho(iHerm,iLag)) + &
                     B17*densities(isospin)%tensor_J(1)%rho(iHerm,iLag) )
!            W_q,rho rho
             HFfield(isospin)%tensor_so(2)%rho(iHerm,iLag)= &
                     mix*HFfield(isospin)%tensor_so(2)%rho(iHerm,iLag)+(1.0_rk-mix)*( &
                     B14*(densities(1)%tensor_J(2)%rho(iHerm,iLag)+densities(2)%tensor_J(2)%rho(iHerm,iLag)) + &
                     B15*densities(isospin)%tensor_J(2)%rho(iHerm,iLag) + &
                     B16*(densities(1)%tensor_J(2)%rho(iHerm,iLag)+densities(2)%tensor_J(2)%rho(iHerm,iLag) + &
                     trace_J_tot) + &
                     B17*(densities(isospin)%tensor_J(2)%rho(iHerm,iLag) + trace_J_q) )
!            W_q,rho phi
             HFfield(isospin)%tensor_so(2)%phi(iHerm,iLag)= &
                     mix*HFfield(isospin)%tensor_so(2)%phi(iHerm,iLag)+(1.0_rk-mix)*( &
                     B14*(densities(1)%tensor_J(2)%phi(iHerm,iLag)+densities(2)%tensor_J(2)%phi(iHerm,iLag)) + &
                     B15*densities(isospin)%tensor_J(2)%phi(iHerm,iLag) + &
                     B16*(densities(1)%tensor_J(3)%rho(iHerm,iLag)+densities(2)%tensor_J(3)%rho(iHerm,iLag)) + &
                     B17*densities(isospin)%tensor_J(3)%rho(iHerm,iLag) )

!            W_q,phi z
             HFfield(isospin)%tensor_so(3)%z(iHerm,iLag)= &
                     mix*HFfield(isospin)%tensor_so(3)%z(iHerm,iLag)+(1.0_rk-mix)*( &
                     B14*(densities(1)%tensor_J(3)%z(iHerm,iLag)+densities(2)%tensor_J(3)%z(iHerm,iLag)) + &
                     B15*densities(isospin)%tensor_J(3)%z(iHerm,iLag) + &
                     B16*(densities(1)%tensor_J(1)%phi(iHerm,iLag)+densities(2)%tensor_J(1)%phi(iHerm,iLag)) + &
                     B17*densities(isospin)%tensor_J(1)%phi(iHerm,iLag) )
!            W_q,phi rho
             HFfield(isospin)%tensor_so(3)%rho(iHerm,iLag)= &
                     mix*HFfield(isospin)%tensor_so(3)%rho(iHerm,iLag)+(1.0_rk-mix)*( &
                     B14*(densities(1)%tensor_J(3)%rho(iHerm,iLag)+densities(2)%tensor_J(3)%rho(iHerm,iLag)) + &
                     B15*densities(isospin)%tensor_J(3)%rho(iHerm,iLag) + &
                     B16*(densities(1)%tensor_J(2)%phi(iHerm,iLag)+densities(2)%tensor_J(2)%phi(iHerm,iLag)) + &
                     B17*densities(isospin)%tensor_J(2)%phi(iHerm,iLag) )
!            W_q,phi phi
             HFfield(isospin)%tensor_so(3)%phi(iHerm,iLag)= &
                     mix*HFfield(isospin)%tensor_so(3)%phi(iHerm,iLag)+(1.0_rk-mix)*( &
                     B14*(densities(1)%tensor_J(3)%phi(iHerm,iLag)+densities(2)%tensor_J(3)%phi(iHerm,iLag)) + &
                     B15*densities(isospin)%tensor_J(3)%phi(iHerm,iLag) + &
                     B16*(densities(1)%tensor_J(3)%phi(iHerm,iLag)+densities(2)%tensor_J(3)%phi(iHerm,iLag) + &
                     trace_J_tot) + &
                     B17*(densities(isospin)%tensor_J(3)%phi(iHerm,iLag)+trace_J_q) )

             HFfield(isospin)%D%z(iHerm,iLag)=mix*HFfield(isospin)%D%z(iHerm,iLag) - &
                  (1.0_rk-mix)*2.0_rk*(B16*s_tot_z+B17*densities(isospin)%s%z(iHerm,iLag))
             HFfield(isospin)%D%rho(iHerm,iLag)=mix*HFfield(isospin)%D%rho(iHerm,iLag) - &
                  (1.0_rk-mix)*2.0_rk*(B16*s_tot_rho+B17*densities(isospin)%s%rho(iHerm,iLag))
             HFfield(isospin)%D%phi(iHerm,iLag)=mix*HFfield(isospin)%D%phi(iHerm,iLag) - &
                  (1.0_rk-mix)*2.0_rk*(B16*s_tot_phi+B17*densities(isospin)%s%phi(iHerm,iLag))

             if (debugmode=='high') then
                write(nqp,'(3(i3,2x),9(f12.6,2x))') isospin,iHerm,iLag, &
                     HFfield(isospin)%central(iHerm,iLag),HFfield(isospin)%effective_mass(iHerm,iLag), &
                     HFfield(isospin)%spin_orbit(iHerm,iLag),HFfield(isospin)%alpha%rho(iHerm,iLag), &
                     HFfield(isospin)%alpha%phi(iHerm,iLag),HFfield(isospin)%alpha%z(iHerm,iLag), &
                     HFfield(isospin)%S%rho(iHerm,iLag),HFfield(isospin)%S%phi(iHerm,iLag), &
                     HFfield(isospin)%S%z(iHerm,iLag)
             end if

          end do

       end do

       integral_spin_field(isospin) = norm * integral_spin_field(isospin)
       integral_s_contribution(isospin) = norm * integral_s_contribution(isospin)
       integral_T_contribution(isospin) = norm * integral_T_contribution(isospin)
       integral_so_contribution(isospin) = norm * integral_so_contribution(isospin)
       integral_Delta_s_contribution(isospin) = norm * integral_Delta_s_contribution(isospin)

    end do

    write(field_unit,'(/a)') &
         'Nucleon Int(S_q,z rho_q)/2   s contribution    T contribution    so contribution   Delta s contribution'
    do isospin = 1, 2
       write(field_unit,'(a1,5(8x,f10.3))') nucleon_type(isospin), &
            0.5_rk*integral_spin_field(isospin),0.5_rk*integral_s_contribution(isospin), &
            0.5_rk*integral_T_contribution(isospin),0.5_rk*integral_so_contribution(isospin), &
            0.5_rk*integral_Delta_s_contribution(isospin)
    end do

    close(field_unit)

    if (debugmode=='high') write(nqp,*)

    return
 
  end subroutine potential


  !=====================================================================
  !                       FUNCTION FCONSTR
  !---------------------------------------------------------------------
  ! Function entering the constraint field.
  !=====================================================================

  function fconstr(k,z,r)

!
!   Arguments
!
    integer,intent(in) :: k
    real(kind=rk),intent(in) :: z,r
!
!   Result
!
    real(kind=rk) :: fconstr
!
!   Local variables
!
    character(len=7),parameter :: subroutinename='fconstr'

    select case(k)
    case(1)                                             ! Center-of-mass position
       fconstr=z
    case(2)                                             ! Quadrupole moment
       fconstr=2.0_rk*z**2-r**2
    case(3)                                             ! Octupole moment
       fconstr=z*(z**2-1.5_rk*r**2)
    case(4)                                             ! Hexadecapole moment
       fconstr=z**4-3.0_rk*z**2*r**2+0.375_rk*r**4
    case(5)                                             ! Neck coordinate
       fconstr=exp(-(z/a_QN)**2) ! Temporary : z_min=0
    case(6)                                             ! Square radius
       fconstr=z**2+r**2
    case default
       call warning(modulename, subroutinename, 'k must be between 1 and Nconstr', &
            'set fconstr to 0')
       fconstr=0.0
    end select

    return

  end function fconstr


  !=====================================================================
  !                 SUBROUTINE SINGLE_PARTICLE_SUM_RULES
  !---------------------------------------------------------------------
  ! Calculates the sum of single-particle energies and the corresponding 
  ! integral form.
  !=====================================================================

  subroutine single_particle_sum_rules(nucleus, model, basis, HFfield, &
       densities, sp_states, Coulomb_exchange_potential)

!
!   Arguments
!
    type(nucleus_properties),intent(in)            :: nucleus
    type(modelization),intent(in)                  :: model
    type(cylindrical_basis),intent(in)             :: basis
    type(field),dimension(2),intent(in)            :: HFfield
    type(local_density),dimension(2),intent(in)    :: densities
    type(sp_state),dimension(basis%size,2),intent(in) :: sp_states
    type(operator),dimension(basis%Nblocks),optional,intent(in) :: &
         Coulomb_exchange_potential

!
!   Local variables
!
    character(len=25),parameter :: subroutinename='single_particle_sum_rules'
    integer :: isospin,ih,il,i,iblock,j,k,n,shift
    integer,dimension(2) :: c_coul
    real(kind=rk) :: Gz,Gr,weight,norm
    real(kind=rk) :: int_effective_mass,int_spin_kinetic,int_coul,int_central
    real(kind=rk) :: int_spin,int_current,int_so,int_spin_current,int_W_J
    real(kind=rk) :: int_Sz_rho,int_U_sz,int_Cz_tau,int_Tz,int_A_J,int_Wso_jphi,int_W_jphi
    real(kind=rk) :: Ci, Cj, Ecoul_exch
    real(kind=rk) :: B9,rho
    real(kind=rk) :: total,sum_e,sum_sz_e,avg_sz
    real(kind=rk) :: int_grad_rho
    real(kind=rk),dimension(3,3) :: int_J2

    c_coul(1) = 0
    c_coul(2) = 1
    norm=fz*pi/(bz*bp2)
    B9 = model%parameters%B9

    do isospin = 1, 2

       select case (isospin)
       case (1)
          write(out_unit,'(a/)') 'NEUTRONS'
       case (2)
          write(out_unit,'(/a/)') 'PROTONS'
       case default
          call fatal(modulename,subroutinename,'isospin should be 1 or 2')
       end select

       int_effective_mass = 0.0
       int_spin_kinetic = 0.0
       int_central = 0.0
       int_coul = 0.0
       int_spin = 0.0
       int_current = 0.0
       int_so = 0.0
       int_W_J = 0.0
       int_spin_current = 0.0
       int_Sz_rho = 0.0 
       int_U_sz = 0.0
       int_Cz_tau = 0.0
       int_Tz = 0.0
       int_A_J = 0.0
       int_Wso_jphi = 0.0
       int_W_jphi = 0.0
       int_J2 = 0.0
       int_grad_rho = 0.0

       do iH = iHmin,basis%NGz/2

          if (iH==0) cycle
          Gz=basis%G_Herm(iH)

          do iL = 1, basis%NGr

             Gr=basis%G_Lag(iL)
             rho = sqrt(basis%xLag(iL))/bp
             weight=Gz*Gr

             int_effective_mass = int_effective_mass + weight * &
                  HFfield(isospin)%effective_mass(iH,iL) * densities(isospin)%tau(iH,iL)

             int_spin_kinetic = int_spin_kinetic + weight * ( &
                  HFfield(isospin)%C%z(iH,iL) * densities(isospin)%T%z(iH,iL) + &
                  HFfield(isospin)%C%rho(iH,iL) * densities(isospin)%T%rho(iH,iL) + &
                  HFfield(isospin)%C%phi(iH,iL) * densities(isospin)%T%phi(iH,iL) )

             int_central = int_central + weight * (HFfield(isospin)%central(iH,iL) - &
                  c_coul(isospin) * HFfield(2)%Coulomb(iH,iL)) * densities(isospin)%rho(iH,iL)

             int_coul = int_coul + weight * c_coul(isospin) * HFfield(2)%Coulomb(iH,iL) * &
                  densities(2)%rho(iH,iL)

             int_spin = int_spin + weight * ( &
                  HFfield(isospin)%S%z(iH,iL) * densities(isospin)%s%z(iH,iL) + &
                  HFfield(isospin)%S%rho(iH,iL) * densities(isospin)%s%rho(iH,iL) + &
                  HFfield(isospin)%S%phi(iH,iL) * densities(isospin)%s%phi(iH,iL))

             int_current = int_current + weight * ( &
                  HFfield(isospin)%alpha%z(iH,iL) * densities(isospin)%j%z(iH,iL) + &
                  HFfield(isospin)%alpha%rho(iH,iL) * densities(isospin)%j%rho(iH,iL) + &
                  HFfield(isospin)%alpha%phi(iH,iL) * densities(isospin)%j%phi(iH,iL))

             int_so = int_so + weight * &
                  (-1)*HFfield(isospin)%spin_orbit(iH,iL) * densities(isospin)%div_J(iH,iL)

             do i = 1, 3
                int_spin_current = int_spin_current + weight * 2.0_rk * &
                     (HFfield(isospin)%tensor_so(i)%z(iH,iL) * densities(isospin)%tensor_J(i)%z(iH,iL) + &
                     HFfield(isospin)%tensor_so(i)%rho(iH,iL) * densities(isospin)%tensor_J(i)%rho(iH,iL) + &
                     HFfield(isospin)%tensor_so(i)%phi(iH,iL) * densities(isospin)%tensor_J(i)%phi(iH,iL) )
             end do

             do i = 1, 3
                int_J2(i,1) = int_J2(i,1) + weight * densities(isospin)%tensor_J(i)%z(iH,iL)**2
                int_J2(i,2) = int_J2(i,2) + weight * densities(isospin)%tensor_J(i)%rho(iH,iL)**2
                int_J2(i,3) = int_J2(i,3) + weight * densities(isospin)%tensor_J(i)%phi(iH,iL)**2
             end do

             int_grad_rho = int_grad_rho + weight * (densities(isospin)%grad_rho%z(iH,iL) * &
                  basis%xHerm(iH)/bz + densities(isospin)%grad_rho%rho(iH,iL) * rho)

             int_W_J = int_W_J + weight * (-B9) * &
                  ( (densities(1)%grad_rho%z(iH,iL) + densities(2)%grad_rho%z(iH,iL) + &
                  densities(isospin)%grad_rho%z(iH,iL)) * (densities(isospin)%tensor_J(2)%phi(iH,iL) - &
                  densities(isospin)%tensor_J(3)%rho(iH,iL)) + &
                  (densities(1)%grad_rho%rho(iH,iL) + densities(2)%grad_rho%rho(iH,iL) + &
                  densities(isospin)%grad_rho%rho(iH,iL)) * (densities(isospin)%tensor_J(3)%z(iH,iL) - &
                  densities(isospin)%tensor_J(1)%phi(iH,iL)))

             int_Sz_rho = int_Sz_rho + weight * &
                  HFfield(isospin)%S%z(iH,iL) * densities(isospin)%rho(iH,iL)
             int_U_sz = int_U_sz + weight * &
                  HFfield(isospin)%central(iH,iL) * densities(isospin)%s%z(iH,iL)
             int_Cz_tau = int_Cz_tau + weight * &
                  HFfield(isospin)%C%z(iH,iL) * densities(isospin)%tau(iH,iL)
             int_Tz = int_Tz + weight * &
                  HFfield(isospin)%effective_mass(iH,iL) * densities(isospin)%T%z(iH,iL)
             int_A_J = int_A_J + weight * HFfield(isospin)%alpha%phi(iH,iL) * &
                  densities(isospin)%tensor_J(3)%z(iH,iL)
             int_Wso_jphi = int_Wso_jphi + weight * &
                  (- B9) * (densities(1)%grad_rho%rho(iH,iL) + densities(2)%grad_rho%rho(iH,iL) + &
                  densities(isospin)%grad_rho%rho(iH,iL)) * densities(isospin)%j%phi(iH,iL)
             int_W_jphi = int_W_jphi + weight * &
                  2 * HFfield(isospin)%tensor_so(3)%z(iH,iL) * densities(isospin)%j%phi(iH,iL)

          end do

       end do

       int_effective_mass = norm * int_effective_mass
       int_spin_kinetic = norm * int_spin_kinetic
       int_central = norm * int_central
       int_coul = norm * int_coul
       int_spin = norm * int_spin
       int_current = norm * int_current
       int_so = norm * int_so
       int_spin_current = norm * int_spin_current
       int_W_J = norm * int_W_J

       int_J2 = int_J2 * norm
       int_grad_rho = int_grad_rho * norm

       int_Sz_rho = 0.5_rk * norm * int_Sz_rho
       int_U_sz = 0.5_rk * norm * int_U_sz
       int_Cz_tau = 0.5_rk * norm * int_Cz_tau
       int_Tz = 0.5_rk * norm * int_Tz
       int_A_J = 0.5_rk * norm * int_A_J
       int_Wso_jphi = 0.5_rk * norm * int_Wso_jphi
       int_W_jphi = 0.5_rk * norm * int_W_jphi

       Ecoul_exch = 0.0

       if (isospin == 2 .and. present(Coulomb_exchange_potential)) then

          do k = 1, basis%size
             if (abs(sp_states(k,2)%occupation) < 1e-6_rk) cycle
             iblock = sp_states(k,2)%block
             n = basis%block_size(iblock)
             shift = 0
             if (iblock>1) shift = sum(basis%block_size(1:iblock-1))
             do i = 1, n
                Ci = sp_states(k,2)%coef%cylHO(i+shift)
                do j = 1, n
                   Cj = sp_states(k,2)%coef%cylHO(j+shift)
                   Ecoul_exch = Ecoul_exch + 0.5_rk * sp_states(k,2)%occupation * &
                        Coulomb_exchange_potential(iblock)%HO(i,j) * Ci * Cj
                end do
             end do
          end do

       end if

       int_coul = int_coul + 2.0_rk * Ecoul_exch

       write(out_unit,'(a/)') 'Expectation value of h_HF'

       total = int_effective_mass + int_spin_kinetic + int_central + int_coul + &
            int_spin + int_current + int_so + int_spin_current

       write(out_unit,'(2(a,f15.6,2x))') 'int_effective_mass =',int_effective_mass, &
            'int_spin_kinetic =',int_spin_kinetic
       write(out_unit,'(2(a,f15.6,2x))') 'int_central        =',int_central, &
            'int_spin         =',int_spin
       if (isospin == 2) write(out_unit,'(a,f15.6,2xa,f13.6)') &
            'int_Coulomb        =',int_coul,'including exchange =',Ecoul_exch
       write(out_unit,'(a,f15.6)') 'int_current        =',int_current
       write(out_unit,'(2(a,f15.6,2x))') 'int_spin_orbit     =',int_so, &
            'int_spin_current   =',int_spin_current
       write(out_unit,'(a,f15.6)') 'Total     =',total

       if (debugmode /= 'low') then
          write(out_unit,'(a,f15.6)') 'int_W_J            =',int_W_J
          write(out_unit,'(a,f15.6)') 'int_grad_rho   =',int_grad_rho
          write(out_unit,'(/3(a,f15.6,2x))') &
               'int_J2_z z=',int_J2(1,1),'int_J2_z rho=',int_J2(1,2),'int_J2_z phi=',int_J2(1,3)
          write(out_unit,'(3(a,f15.6,2x))') &
               'int_J2_rho z=',int_J2(2,1),'int_J2_rho rho=',int_J2(2,2),'int_J2_rho phi=',int_J2(2,3)
          write(out_unit,'(3(a,f15.6,2x)/)') &
               'int_J2_phi z=',int_J2(3,1),'int_J2_phi rho=',int_J2(3,2),'int_J2_phi phi=',int_J2(3,3)
       end if

       sum_e = 0.0
       do k = 1, basis%size
          sum_e = sum_e + sp_states(k,isospin)%occupation * sp_states(k,isospin)%energy
       end do
       write(out_unit,'(/a,f15.6/)') 'sum_e     =',sum_e

       write(out_unit,'(a/)') 'Expectation value of sz x h_HF'

       total = int_Sz_rho + int_U_sz + int_Cz_tau + int_Tz + int_A_J + int_Wso_jphi + int_W_jphi

       write(out_unit,'(2(a,f11.6,2x))') 'int_Sz_rho =',int_Sz_rho,'int_U_sz     =',int_U_sz
       write(out_unit,'(2(a,f11.6,2x))') 'int_Cz_tau =',int_Cz_tau,'int_m*_Tz    =',int_Tz
       write(out_unit,'(3(a,f11.6,2x))') 'int_A_J    =',int_A_J,'int_Wso_jphi =',int_Wso_jphi, &
            'int_W_jphi =',int_W_jphi
       if (isospin == 1) then
          write(out_unit,'(a,f11.6)') 'Total      =',total
       else
          write(out_unit,'(a,f11.6,2x,a)') 'Total      =',total,'(without exact exchange Coulomb)'
       end if

       sum_sz_e = 0.0
       do k = 1, basis%size
          iblock = sp_states(k,isospin)%block
          n = basis%block_size(iblock)
          shift = 0
          if (iblock>1) shift = sum(basis%block_size(1:iblock-1))
          avg_sz = 0.0
          do i = 1, n
             avg_sz = avg_sz + sp_states(k,isospin)%coef%cylHO(i+shift)**2*basis%ns(i+shift)
          end do
          sum_sz_e = sum_sz_e + sp_states(k,isospin)%occupation * 0.5_rk * avg_sz * &
               sp_states(k,isospin)%energy
       end do
       write(out_unit,'(a,f11.6)') 'sum_sz_e   =',sum_sz_e

    end do

  end subroutine single_particle_sum_rules


  !=====================================================================
  !                       SUBROUTINE NILSSON_EXPANSION
  !---------------------------------------------------------------------
  ! Determines dominant Nilsson configurations in single-particle states
  !=====================================================================

  subroutine Nilsson_expansion(nucleus, model, symmetries, basis, sp_states)

    implicit none
    !
    ! Arguments
    !
    type(nucleus_properties),intent(in)                 :: nucleus
    type(modelization),intent(in)                       :: model
    type(symmetry_properties),intent(in)                :: symmetries
    type(cylindrical_basis),intent(in)                  :: basis
    type(sp_state),dimension(basis%size,2),intent(inout):: sp_states
    !
    ! Local variables
    !
    integer :: isospin, iHF, i, j
    integer,dimension(2) :: Npart
    integer,dimension(basis%size) :: permut
    real(kind=rk),dimension(basis%size) :: tab

    Npart(:)=(/ nucleus%N, nucleus%Z /)
    
    do isospin=1,2
       do iHF=1,basis%size
          sp_states(iHF,isospin)%Nilsson(:)%N = 0
          sp_states(iHF,isospin)%Nilsson(:)%nz = 0
          sp_states(iHF,isospin)%Nilsson(:)%Lambda = 0
          sp_states(iHF,isospin)%Nilsson(:)%Sigma = 0
          sp_states(iHF,isospin)%Nilsson(:)%amplitude = 0
          tab=abs(sp_states(iHF,isospin)%coef%cylHO)
          if (sp_states(iHF,isospin)%energy > sp_states(Npart(isospin),isospin)%energy + &
               2*model%truncmax) exit
          call sort(tab,basis%size,'decreasing',permut)
          do i = 1, nb_Nilsson
             j = permut(i)
             sp_states(iHF,isospin)%Nilsson(i)%amplitude = &
                  sp_states(iHF,isospin)%coef%cylHO(j)
             sp_states(iHF,isospin)%Nilsson(i)%N = basis%nz(j) + &
                  2*basis%nr(j)+abs(basis%Lambda(j))
             sp_states(iHF,isospin)%Nilsson(i)%nz = basis%nz(j)
             sp_states(iHF,isospin)%Nilsson(i)%Lambda = basis%lambda(j)
             sp_states(iHF,isospin)%Nilsson(i)%Sigma = basis%ns(j)
          end do
       end do
    end do

  end subroutine Nilsson_expansion


  !=====================================================================
  !                       SUBROUTINE DISPLAY_ITERATION_RESULTS
  !---------------------------------------------------------------------
  ! Display nuclear properties (moments, pairing and energy) at the end 
  ! of the current iteration
  !=====================================================================

  subroutine display_iteration_results(nucleus, model, symmetries, iter, &
       basis_size, sp_states)

!
!   Arguments
!
    type(nucleus_properties),intent(in)  :: nucleus
    type(modelization),intent(in)        :: model
    type(symmetry_properties),intent(in) :: symmetries
    integer,intent(in)                   :: iter, basis_size
    type(sp_state),dimension(basis_size,2),intent(in) :: sp_states

!
!   Local variables
!
    integer :: N,Z
    integer :: i, isospin, j
    integer :: pn, pp
    integer,dimension(2) :: Npart
    character(len=100) :: Nilsson_config, tmp
!
!   Formats
!
600 format(/,'<Sz>(n)=',f10.6,1x,'hbar',2x,'<Sz>(p)=',f10.6,1x,'hbar',2x,'<Sz>(tot)=',f10.6,1x,'hbar',/ &
         '<Jz>(n)=',f10.6,1x,'hbar',2x,'<Jz>(p)=',f10.6,1x,'hbar',2x,'<Jz>(tot)=',f10.6,1x,'hbar',/ &
         'r2=',f11.3,1x,'fm^2',/ &
         'Q20(n)=',f11.3,1x,'fm^2',2x,'Q20(p)=',f11.3,1x,'fm^2',2x,/ &
         'Q20=',f11.3,1x,'fm^2',2x,'Q40=',f10.1,1x,'fm^4',2x,'QN=',f7.3,/ &
         'Etot=',f14.6,1x,'MeV'/)
601 format(/,'lambda_n=',f9.4,' MeV  lambda_p=',f9.4,' MeV',/ &
         '<Delta>(n) =',f9.4,' MeV  <Delta>(p) =',f9.4,' MeV',/ &
!         'Delta_n =',f9.4,' MeV  Delta_p =',f9.4,' MeV (at Fermi level)',/ &
         '<Sz>(n)=',f10.6,1x,'hbar',2x,'<Sz>(p)=',f10.6,1x,'hbar',2x,'<Sz>(tot)=',f10.6,1x,'hbar',/ &
         '<Jz>(n)=',f10.6,1x,'hbar',2x,'<Jz>(p)=',f10.6,1x,'hbar',2x,'<Jz>(tot)=',f10.6,1x,'hbar',/ &
         'pi(n)=',i3,1x,'pi(p)=',i3,/ &
         '<r2>(n)=',f11.3,1x,'fm^2',2x,'<r2>(p)=',f11.3,1x,'fm^2',/ &
         '<r2>=',f11.3,1x,'fm^2',/ &
         'Q20(n)=',f11.3,1x,'fm^2',2x,'Q20(p)=',f11.3,1x,'fm^2',2x,/ &
         'Q20=',f11.3,1x,'fm^2',2x,'Q40=',f10.1,1x,'fm^4',2x,'QN=',f7.3,/ &
         'beta2=',f8.4,/ &
         'Etot=',f14.6,1x,'MeV'/)
602 format('Iteration ',i3,/ &
         /'<Sz>(n)=',f10.6,1x,'hbar',2x,'<Sz>(p)=',f10.6,1x,'hbar',2x,'<Sz>(tot)=',f10.6,1x,'hbar',/ &
         '<Jz>(n)=',f10.6,1x,'hbar',2x,'<Jz>(p)=',f10.6,1x,'hbar',2x,'<Jz>(tot)=',f10.6,1x,'hbar',/ &
         'zcm=',f8.3,1x,'fm',2x,'r2=',f11.3,1x,'fm^2',/ &
         'Q20(n)=',f11.3,1x,'fm^2',2x,'Q20(p)=',f11.3,1x,'fm^2',2x,/ &
         'Q20=',f11.3,1x,'fm^2',2x,'Q30=',f10.1,1x,'fm^3',2x, &
         'Q40=',f10.1,1x,'fm^4',2x,'QN=',f7.3,/'Etot=',f14.6,1x,'MeV'/)
603 format('lambda_n=',f9.4,' MeV  lambda_p=',f9.4,' MeV',/ &
         '<Delta>(n) =',f9.4,' MeV  <Delta>(p) =',f9.4,' MeV',/ &
!         'Delta_n =',f9.4,' MeV  Delta_p =',f9.4,' MeV (at Fermi level)',/ &
         '<Sz>(n)=',f10.6,1x,'hbar',2x,'<Sz>(p)=',f10.6,1x,'hbar',2x,'<Sz>(tot)=',f10.6,1x,'hbar',/ &
         '<Jz>(n)=',f10.6,1x,'hbar',2x,'<Jz>(p)=',f10.6,1x,'hbar',2x,'<Jz>(tot)=',f10.6,1x,'hbar',/ &
         'zcm=',f8.3,1x,'fm',2x,'r2=',f11.3,1x,'fm^2',/ &
         '<pi>(n)=',f6.3,1x,'<pi>(p)=',f6.3,1x,'<pi>(tot)=',f6.3,/ &
         'Q20(n)=',f11.3,1x,'fm^2',2x,'Q20(p)=',f11.3,1x,'fm^2',2x,/ &
         'Q20=',f11.3,1x,'fm^2',2x,'Q30=',f10.1,1x,'fm^3',2x, &
         'Q40=',f10.1,1x,'fm^4',2x,'QN=',f7.3,/ &
         'beta2=',f8.4,2x,'beta3=',f8.4,/ &
         'Etot=',f14.6,1x,'MeV'/)
610 format('Ekin=',f14.6,1x,'MeV',2x,'Evol =',f14.6,1x,'MeV',2x,'Esurf=',f14.6,1x,'MeV',/ &
         'Eso =',f14.6,1x,'MeV',2x,'Ecoul=',f14.6,1x,'MeV',2x,'Epair=',f14.6,1x,'MeV'/,&
         'Eodd=',f14.6,1x,'MeV'/)


    if (symmetries%parity) then

!       if (iter==1 .and. model%pairing=='DDDI') then
!          write(out_unit,600) iter, &
!               nucleus%Sz(1),nucleus%Sz(2),nucleus%Sz(0), &
!               nucleus%Jz(1),nucleus%Jz(2),nucleus%Jz(0), &
!               nucleus%r2(0),nucleus%Q20(1),nucleus%Q20(2),&
!               nucleus%Q20(0),nucleus%Q40(0),nucleus%QN(0), &
!               nucleus%energy%Etot          
!       else
          N=nucleus%N
          if (nucleus%N == 0) N=1
          Z=nucleus%Z
          if (nucleus%Z == 0) Z=1
          if (model%blocking /= 'no') then
             pn = 1
             do i = 1, size(nucleus%pi(1,:))
                if (nucleus%K(1,i) /= 0) pn = pn * nucleus%pi(1,i)
             end do
             pp = 1
             do i = 1, size(nucleus%pi(2,:))
                if (nucleus%K(2,i) /= 0) pp = pp * nucleus%pi(2,i)
             end do
          else
             pn = 1
             pp = 1
          end if
          write(out_unit,601) &
               nucleus%pairing(:)%chemical_potential, &
               nucleus%pairing(:)%average_gap, &
!               nucleus%pairing(1)%gap(nucleus%N/2),nucleus%pairing(2)%gap(nucleus%Z/2), &
               nucleus%Sz(1),nucleus%Sz(2),nucleus%Sz(0), &
               nucleus%Jz(1),nucleus%Jz(2),nucleus%Jz(0), &
               pn,pp, &
               nucleus%r2(1)/N,nucleus%r2(2)/Z,nucleus%r2(0)/nucleus%A,&
               nucleus%Q20(1),nucleus%Q20(2), &
               nucleus%Q20(0),nucleus%Q40(0),nucleus%QN(0), &
               nucleus%beta2,nucleus%energy%Etot
!       end if

    else

!       if (iter==1 .and. model%pairing=='DDDI') then
!          write(out_unit,602) iter, &
!               nucleus%Sz(1),nucleus%Sz(2),nucleus%Sz(0), &
!               nucleus%Jz(1),nucleus%Jz(2),nucleus%Jz(0), &
!               nucleus%zcm(0),nucleus%r2(0),nucleus%Q20(1),nucleus%Q20(2), &
!               nucleus%Q20(0),nucleus%Q30(0),nucleus%Q40(0),nucleus%QN(0), &
!               nucleus%energy%Etot
!       else
          write(out_unit,603) &
               nucleus%pairing(:)%chemical_potential, &
               nucleus%pairing(:)%average_gap, &
!               nucleus%pairing(1)%gap(nucleus%N/2),nucleus%pairing(2)%gap(nucleus%Z/2), &
               nucleus%Sz(1),nucleus%Sz(2),nucleus%Sz(0), &
               nucleus%Jz(1),nucleus%Jz(2),nucleus%Jz(0), &
               nucleus%zcm(0),nucleus%r2(0), &
               nucleus%avg_pi(1),nucleus%avg_pi(2),nucleus%avg_pi(0), &
               nucleus%Q20(1),nucleus%Q20(2), &
               nucleus%Q20(0),nucleus%Q30(0),nucleus%Q40(0),nucleus%QN(0), &
               nucleus%beta2,nucleus%beta3,nucleus%energy%Etot
 !      end if

    end if

    write(out_unit,610) nucleus%energy%Ekin, nucleus%energy%Evol, &
    nucleus%energy%Esurf, nucleus%energy%Eso, nucleus%energy%Ecoul, &
    nucleus%energy%Epair, nucleus%energy%Eodd

    write(out_unit,'(4(a,f13.6,2x))') 'E_B1 =', nucleus%energy%E_B(1), &
         'E_B2 =', nucleus%energy%E_B(2), 'E_B3 =', nucleus%energy%E_B(3), &
         'E_B4 =', nucleus%energy%E_B(4)
    write(out_unit,'(4(a,f13.6,2x))') 'E_B5 =', nucleus%energy%E_B(5), &
         'E_B6 =', nucleus%energy%E_B(6), 'E_B7 =', nucleus%energy%E_B(7), &
         'E_B8 =', nucleus%energy%E_B(8)
    write(out_unit,'(a,f13.6)') 'E_B9 =', nucleus%energy%E_B(9)
    write(out_unit,'(4(a,f13.6,2x))') 'E_B10=', nucleus%energy%E_B(10), &
         'E_B11=', nucleus%energy%E_B(11), 'E_B12=', nucleus%energy%E_B(12), &
         'E_B13=', nucleus%energy%E_B(13)
    write(out_unit,'(4(a,f13.6,2x))') 'E_B14=', nucleus%energy%E_B(14), &
         'E_B15=', nucleus%energy%E_B(15), 'E_B16=', nucleus%energy%E_B(16), &
         'E_B17=', nucleus%energy%E_B(17)
    write(out_unit,'(4(a,f13.6,2x))') 'E_B18=', nucleus%energy%E_B(18), &
         'E_B19=', nucleus%energy%E_B(19), 'E_B20=', nucleus%energy%E_B(20), &
         'E_B21=', nucleus%energy%E_B(21)
    write(out_unit,'(2(a,f13.6,2x))') 'E_u0 =', nucleus%energy%E_B(22), &
         'E_v0 =', nucleus%energy%E_B(23)
    write(out_unit,'(a,f13.6)') 'Total binding energy: Ekin + sum_i E_B(i) + Ecoul + Epair =',&
         nucleus%energy%Ekin + sum(nucleus%energy%E_B(:)) + nucleus%energy%Ecoul + nucleus%energy%Epair
    write(out_unit,'(a,f13.6)') 'Total NN potential energy: sum_i E_B(i) =',&
         sum(nucleus%energy%E_B(:))
    write(out_unit,'(/a,f12.6,a/)') 'Epot = <BCS|V-bar|BCS>/2 =',nucleus%energy%Epot,' MeV'

    if (model%blocking=='no' .and. (mod(nucleus%N,2) == 1 .or. mod(nucleus%Z,2) == 1)) then
       write(out_unit,'(/2(a,f12.6,2x)/)') &
            'False vacuum average particle number: <N>(n) =',nucleus%avg_particle_number(1), &
            '<N>(p) =',nucleus%avg_particle_number(2)
    end if

    write(out_unit,'(a,f8.2,1x,a)') 'Single-particle spectrum in window [Fermi-X;Fermi+X] ' // &
         'with X =',model%truncmax,'MeV:'
    if (symmetries%parity) then
       write(out_unit,'(/70a1/,a/,a/,70a1)') ('-',i=1,70), &
            'isospin label Omega pi rank    energy  occupation  Nilsson config.', &
            '                                                   [N nz /\ Sigma]', &
            ('-',i=1,70)
    else
       write(out_unit,'(/70a1/,a/,a/,70a1)') ('-',i=1,70), &
            'isospin label Omega <pi> rank    energy  occupation  Nilsson config.', &
            '                                                     [N nz /\ Sigma]', &
            ('-',i=1,70)
    end if
    Npart(:)=(/ nucleus%N , nucleus%Z /)
    do isospin = 1, 2
       do i = 1, basis_size
          if (sp_states(i,isospin)%energy < sp_states(Npart(isospin),isospin)%energy - &
               model%truncmax) cycle
          if (sp_states(i,isospin)%energy > sp_states(Npart(isospin),isospin)%energy + &
               model%truncmax) exit
          Nilsson_config=''
          do j = 1, min(2,nb_Nilsson)
             write(tmp,'(f5.1,a,2(i2,1x),i3,1x,a)') &
                  sp_states(i,isospin)%Nilsson(j)%amplitude**2*100,'% [', &
                  sp_states(i,isospin)%Nilsson(j)%N, sp_states(i,isospin)%Nilsson(j)%nz, &
                  sp_states(i,isospin)%Nilsson(j)%Lambda, &
                  parity_sign(sp_states(i,isospin)%Nilsson(j)%Sigma) // ']'
             Nilsson_config=trim(Nilsson_config)//trim(tmp)
          end do
          if (symmetries%parity) then
             write(out_unit,'(2x,a1,i9,i5,a2,1x,a1,1x,i4,1x,f11.6,f10.4,1x,a)') &
                  nucleon_type(isospin),i,sp_states(i,isospin)%Omega,'/2', &
                  parity_sign(sp_states(i,isospin)%pi),sp_states(i,isospin)%block_rank, &
                  sp_states(i,isospin)%energy,sp_states(i,isospin)%occupation, &
                  trim(Nilsson_config)
          else
             write(out_unit,'(2x,a1,i9,i5,a2,1x,f7.3,1x,i4,1x,f10.3,f10.4,1x,a)') &
                  nucleon_type(isospin),i,sp_states(i,isospin)%Omega,'/2', &
                  sp_states(i,isospin)%avg_pi,sp_states(i,isospin)%block_rank, &
                  sp_states(i,isospin)%energy,sp_states(i,isospin)%occupation, &
                  trim(Nilsson_config)
          end if
       end do
       write(out_unit,'(70a1)') ('-',i=1,70)
    end do
    return

  end subroutine display_iteration_results


  !=====================================================================
  !                       SUBROUTINE OBSERVABLES
  !---------------------------------------------------------------------
  ! Calculates physical quantities of interest associated to the SHFBCS
  ! solution passed as argument 5 'sp_states'.
  !=====================================================================

  subroutine observables(nucleus, model, symmetries, basis, sp_states, &
       densities, HFfield, Coulomb_exchange_potential)

!
!   Arguments
!
    type(nucleus_properties),intent(inout)              :: nucleus
    type(modelization),intent(in)                       :: model
    type(symmetry_properties),intent(in)                :: symmetries
    type(cylindrical_basis),intent(inout)               :: basis
    type(sp_state),dimension(basis%size,2),intent(inout):: sp_states
    type(local_density),dimension(2),intent(in)         :: densities
    type(field),dimension(2),intent(in)                 :: HFfield
    type(operator),dimension(basis%Nblocks),optional,intent(in) :: &
         Coulomb_exchange_potential

!
!   Local variables
!
    character(len=11),parameter :: subroutinename='observables'
    type(spherical_basis) :: sph_basis
    type(nucleus_properties),dimension(2) :: parent_nucleus
    character(len=10) :: operator_name
    real(kind=rk) :: avg_sz,avg_lz,quenching_factor
    real(kind=rk),dimension(nqp) :: s_odd
    real(kind=rk),dimension(2) :: ewsz
    integer :: iunit,i,j,k,l,iHFn,iHFn_max,iHFp,signature_n,signature_p
    integer :: iHF,jHF,isospin,iblock,shift
    integer,dimension(basis%size) :: neutron_analog, proton_analog
    real(kind=rk) :: sign_convention,ovlp,ovlp_max,norm,det
    real(kind=rk) :: tmp, alpha
    real(kind=rk),dimension(:,:),allocatable :: overlap_matrix
    integer,dimension(2) :: Npart
    integer :: count_hole,count_particle
    real(kind=rk) :: fermi,Emin,Emax,mean_particle,mean_hole,rms_particle,rms_hole
    real(kind=rk) :: Erot
    real(kind=rk),dimension(1:2) :: Q20core,Q20tot
    real(kind=rk),dimension(2,basis%size) :: Q20sp,r2sp
    real(kind=rk),dimension(basis%size,basis%size) :: j2matrix_HFbasis, eigenvec_j2
    real(kind=rk),dimension(basis%size) :: eigenval_j2
    real(kind=rk) :: z,r,Gz,Gr
    integer :: iHerm, iLag, Block_num, i_min, i_max
    real(kind=rk) :: v2, p2, r2, z2, vbar, C1, C2
    real(kind=rk),dimension(2) :: Ekin1, Ekin2, Ekin2dir, Ekin2exch, r2MB, Epot
    complex(kind=rk) :: px, pz
    real(kind=rk) :: py
    complex(kind=rk),dimension(3) :: PHF
    type(sp_state), dimension(2) :: HF_states
    character(len=1000) :: Nilsson_config, sph_config, str
    character(len=5),dimension(2*l_max+1) :: str_l_j
    real(kind=rk),dimension(2*l_max+1) :: weight_l_j
    integer,dimension(2*l_max+1) :: permut
    real(kind=rk),parameter :: weight_threshold = 2.0_rk ! in %
    
!
!   Header message
!
    
    write(out_unit,'(/80a1)') ('=',i=1,80)
    do i = 1, 3
       write(out_unit,'(a1,78x,a1)') '|','|'
    end do
    write(out_unit,'(a1,5x,a,52x,a1)') '|','O B S E R V A B L E S','|'
    do i = 1, 3
       write(out_unit,'(a1,78x,a1)') '|','|'
    end do
    write(out_unit,'(80a1/)') ('=',i=1,80)

    Npart(1) = nucleus%N
    Npart(2) = nucleus%Z

    if (model%pairing == 'HTDA') then
       nucleus%pairing(1)%chemical_potential = 0.5_rk*( & ! By default
            sp_states(nucleus%N,1)%energy + sp_states(nucleus%N+1,1)%energy)
       nucleus%pairing(2)%chemical_potential = 0.5_rk*( &
            sp_states(nucleus%Z,2)%energy + sp_states(nucleus%Z+1,2)%energy)
    end if

    !-----------------------------------------------------------------------
    !
    ! Decomposition of s.p. energies into Skyrme-HF fields contributions
    !
    !-----------------------------------------------------------------------

    if (model%h_tensor_perturbative) then
       open(unit=37,file='h_tensor_perturbative.txt',status='replace')
       call h_tensor_perturbative(nucleus, model, symmetries, basis, &
            HFfield, densities, 9999, 37, 'exact')
       close(37)
    end if
    
    !-----------------------------------------------------------------------
    !
    ! 1) Display single-particle spectra
    !
    !-----------------------------------------------------------------------

    if (nucleus%single_particle_spectra) then

       write(out_unit,'(/77a1/,/10x,a/,/77a1)') ('=',i=1,77), &
            'S I N G L E - P A R T I C L E   S P E C T R A', &
            ('-',i=1,77)

       write(out_unit,'(/a,f7.2,a)') 'Weight display threshold =',weight_threshold,'%'

       do isospin = 1, 2
          if (symmetries%parity) then 
             write(out_unit,'(/77a1/,a/,a/,77a1)') ('-',i=1,77), &
                  'iso. label Omega pi  rank  energy   occup.  pair   Nilsson config.',&
                  '                                          partner  [N nz /\ Sigma]', &
                  ('-',i=1,77)
          else
             write(out_unit,'(/77a1/,a/,a/,77a1)') ('-',i=1,77), &
                  'iso. label Omega <pi> rank  energy   occup.  pair  Nilsson config.', &
                  '                                           partner [N nz /\ Sigma]', &
                  ('-',i=1,77)
          end if
          do iHF = 1, basis%size
             if (sp_states(iHF,isospin)%energy > 2*model%truncmax + &
                  sp_states(Npart(isospin),isospin)%energy) exit
             Nilsson_config=''
             do j = 1, nb_Nilsson
                if (sp_states(iHF,isospin)%Nilsson(j)%amplitude**2*100 < weight_threshold) exit
                write(str,'(f5.1,a,2(i2,1x),i3,1x,a)') &
                     sp_states(iHF,isospin)%Nilsson(j)%amplitude**2*100,'% [', &
                     sp_states(iHF,isospin)%Nilsson(j)%N, &
                     sp_states(iHF,isospin)%Nilsson(j)%nz, &
                     sp_states(iHF,isospin)%Nilsson(j)%Lambda, &
                     parity_sign(sp_states(iHF,isospin)%Nilsson(j)%Sigma) // ']'
                Nilsson_config=trim(Nilsson_config)//trim(str)
             end do
             if (symmetries%parity) then
                write(out_unit,'(a1,2(1x,i3),a2,1x,a1,1x,i4,1x,f10.4,1x,f7.4,2x,i3,2x,a)') &
                     nucleon_type(isospin),iHF,sp_states(iHF,isospin)%Omega,'/2', &
                     parity_sign(sp_states(iHF,isospin)%pi),sp_states(iHF,isospin)%block_rank, &
                     sp_states(iHF,isospin)%energy,sp_states(iHF,isospin)%occupation, &
                     sp_states(iHF,isospin)%pair_partner(isospin),trim(Nilsson_config)
             else
                write(out_unit,'(a1,2(1x,i3),a2,1x,f7.3,1x,i4,1x,f10.4,1x,f7.4,2x,i3,2x,a)') &
                     nucleon_type(isospin),iHF,sp_states(iHF,isospin)%Omega,'/2', &
                     sp_states(iHF,isospin)%avg_pi,sp_states(iHF,isospin)%block_rank,&
                     sp_states(iHF,isospin)%energy, sp_states(iHF,isospin)%occupation, &
                     sp_states(iHF,isospin)%pair_partner(isospin),trim(Nilsson_config)
             end if
          end do
          write(out_unit,'(77a1)') ('-',i=1,77)
       end do

       write(out_unit,'(/77a1/)') ('=',i=1,77)

    end if

    !-----------------------------------------------------------------------
    !
    ! 2) Spherical expansion in isotropic HO basis
    !
    !-----------------------------------------------------------------------

    if (nucleus%spherical_expansion) then

       do isospin = 1, 2
          call cyl2sph(nucleus,symmetries,basis,sp_states(:,isospin), &
               sph_basis,isospin)
       end do

       write(out_unit,'(/77a1/,/10x,a/,/77a1)') ('=',i=1,77), &
            'S P H E R I C A L   E X P A N S I O N', &
            ('-',i=1,77)

       write(out_unit,'(/a,f7.2,a)') 'Weight display threshold =',weight_threshold,'%'

       do isospin = 1, 2
          if (symmetries%parity) then 
             write(out_unit,'(/77a1/,a/,77a1)') ('-',i=1,77), &
                  'iso. label Omega pi rank   energy  occup.   weights (%)', &
                  ('-',i=1,77)
          else
             write(out_unit,'(/77a1/,a/,77a1)') ('-',i=1,77), &
                  'iso. label Omega <pi> rank   energy  occup.   weights (%)', &
                  ('-',i=1,77)
          end if
          do iHF = 1, basis%size
             if (sp_states(iHF,isospin)%energy > 2*model%truncmax + &
                  sp_states(Npart(isospin),isospin)%energy) exit
             i = 0
             do l = 0, l_max
                do j = abs(2*l-1), 2*l+1, 2
                   i = i + 1
                   weight_l_j(i) = sp_states(iHF,isospin)%sph_weight(l,j)
                   str_l_j(i)=trim(spectroscopic_notation(l,j))
                end do
             end do
             if (i /= 2*l_max+1) stop 'Error of counting of (l,j) combinations'
             call sort(weight_l_j,2*l_max+1,'decreasing',permut)
             sph_config=''
             do i = 1, 2*l_max+1
                if (weight_l_j(i) < weight_threshold) exit ! Keep weights > 1%
                write(str,'(f5.1,a)') weight_l_j(i),'% ' // trim(str_l_j(permut(i)))
                sph_config = trim(sph_config) // trim(str)
             end do
             if (symmetries%parity) then
                write(out_unit,'(a1,2(1x,i3),a2,1x,a1,1x,i4,1x,f10.4,1x,f7.4,2x,a)') &
                     nucleon_type(isospin),iHF,sp_states(iHF,isospin)%Omega,'/2', &
                     parity_sign(sp_states(iHF,isospin)%pi),sp_states(iHF,isospin)%block_rank, &
                     sp_states(iHF,isospin)%energy,sp_states(iHF,isospin)%occupation, &
                     trim(sph_config)
             else
                write(out_unit,'(a1,2(1x,i3),a2,1x,f7.3,1x,i4,1x,f10.4,1x,f7.4,2x,a)') &
                     nucleon_type(isospin),iHF,sp_states(iHF,isospin)%Omega,'/2', &
                     sp_states(iHF,isospin)%avg_pi,sp_states(iHF,isospin)%block_rank,&
                     sp_states(iHF,isospin)%energy,sp_states(iHF,isospin)%occupation, &
                     trim(sph_config)
             end if
          end do
          write(out_unit,'(77a1)') ('-',i=1,77)
       end do

       write(out_unit,'(/77a1/)') ('=',i=1,77)

    end if

    !-----------------------------------------------------------------------
    !
    ! 3) Nuclear shape
    !
    !-----------------------------------------------------------------------

    if (nucleus%nuclear_shape) then
       
       write(out_unit,'(/70a1/,/a/,/70a1)') ('=',i=1,70), &
            '             N U C L E A R   S H A P E', &
            ('-',i=1,70)

       norm=fz/(2*bz*bp2)

       !   Q20 of all s.p. states
       open(unit=30,file='sp_quadrupole_moments.txt',status='replace')
       write(30,'(70a1)') ('-',i=1,70)
       write(30,'(a)') '# n/p    iHF  Omega pi   energy  occupation  q20(fm^2)'
       write(30,'(70a1)') ('-',i=1,70)
       do isospin = 1, 2
          Q20sp(isospin,:) = 0.0
          r2sp(isospin,:) = 0.0
          do iHF = 1, basis%size
             do iHerm=iHmin,basis%NGz/2
                if (iHerm==0) cycle
                z=basis%xHerm(iHerm)/bz
                Gz=basis%G_Herm(iHerm)
                do iLag=1,basis%NGr
                   r=sqrt(basis%xLag(iLag))/bp
                   Gr=basis%G_Lag(iLag)
                   Q20sp(isospin,iHF)=Q20sp(isospin,iHF)+Gz*Gr*(2.0_rk*z**2-r**2) * &
                        scalar_product(sp_states(iHF,isospin)%wf(iHerm,iLag), &
                        sp_states(iHF,isospin)%wf(iHerm,iLag))
                   r2sp(isospin,iHF)=r2sp(isospin,iHF)+Gz*Gr*(z**2+r**2) * &
                        scalar_product(sp_states(iHF,isospin)%wf(iHerm,iLag), &
                        sp_states(iHF,isospin)%wf(iHerm,iLag))
                end do
             end do
             Q20sp(isospin,iHF)=norm*Q20sp(isospin,iHF)
             r2sp(isospin,iHF)=norm*r2sp(isospin,iHF)
             if (sp_states(iHF,isospin)%energy > 2*model%truncmax + &
                  sp_states(Npart(isospin),isospin)%energy) cycle
             if (symmetries%parity) then
                write(30,'(2x,a1,i9,i5,a2,1x,a1,f10.3,f10.4,1(1x,f9.3))') &
                     nucleon_type(isospin),iHF,sp_states(iHF,isospin)%Omega,'/2', &
                     parity_sign(sp_states(iHF,isospin)%pi),sp_states(iHF,isospin)%energy, &
                     sp_states(iHF,isospin)%occupation,Q20sp(isospin,iHF)
             else
                write(30,'(2x,a1,i9,i5,a2,1x,f6.3,f10.3,f10.4,1(1x,f9.3))') &
                     nucleon_type(isospin),iHF,sp_states(iHF,isospin)%Omega,'/2', &
                     sp_states(iHF,isospin)%avg_pi,sp_states(iHF,isospin)%energy, &
                     sp_states(iHF,isospin)%occupation,Q20sp(isospin,iHF)
             end if
          end do
          write(30,'(70a1)') ('-',i=1,70)
       end do
       close(30)

       !   Q20 of the core for neutrons and protons
       write(out_unit,'(/a)') 'Intrinsic axial quadrupole moment breakdown:'
       do isospin = 1, 2
          write(out_unit,'(a)') '* '//trim(particle_type(isospin))//'s:'
          Q20core(isospin)=0.0_rk
          Q20tot(isospin)=0.0_rk
          loop_iHF: do iHF = 1, basis%size
             do i = 1, nqp
                if (iHF == i_odd(isospin,i)) then
                   Q20tot(isospin)=Q20tot(isospin)+Q20sp(isospin,iHF)
                   write(out_unit,'(2x,a,1x,a,i4,1x,a,f8.3,1x,i2,a,1x,2(a,f10.3,1x,a,1x))') &
                        'blocked '//nucleon_type(isospin),'iHF=',iHF,'e=',sp_states(iHF,isospin)%energy, &
                        sp_states(iHF,isospin)%Omega,'/2'//parity_sign(sp_states(iHF,isospin)%pi), &
                        'q20=',Q20sp(isospin,iHF),'fm^2','r2=',r2sp(isospin,iHF),'fm^2'
                   cycle loop_iHF
                end if
             end do
             Q20core(isospin)=Q20core(isospin)+Q20sp(isospin,iHF)*sp_states(iHF,isospin)%occupation
          end do loop_iHF
          Q20tot(isospin)=Q20tot(isospin)+Q20core(isospin)
          write(out_unit,'(2x,a,f12.3,1x,a)') 'core:',Q20core(isospin),'fm^2'
          write(out_unit,'(2x,a,f12.3,1x,a)') '=> total for '// &
               nucleon_type(isospin)//':',Q20tot(isospin),'fm^2'
       end do
       write(out_unit,'(a,f12.3,1x,a)') '=> Total Q20 =',Q20tot(1)+Q20tot(2),'fm^2'

       !
       !   beta parameters of an isodensity surface (spherical expansion of the 
       !   nuclear radius as a function of theta)
       !

       write(out_unit,'(/a,/a,f8.3)') 'Bohr parameters from multipole moments:', &
         'beta2 =',nucleus%beta2

!       call nuclear_shape(sp_states, basis, symmetries, nucleus)

       write(out_unit,*)
       write(out_unit,'(70a1/)') ('=',i=1,70)

    end if

    !-----------------------------------------------------------------------
    !
    ! 4) Moment of inertia of the core
    !
    !-----------------------------------------------------------------------

    if (nucleus%moment_of_inertia) then


      write(out_unit,'(/70a1/,/a//,a/,a/,/a/,a/,/70a1/)') ('=',i=1,70), &
            '                 A N G U L A R   M O M E N T U M', &
            '                M O M E N T   O F   I N E R T I A', &
            '                  (Inglis-Belyaev approximation)', &
            '    C O L L E C T I V E   G Y R O M A G N E T I C   R A T I O', &
            '                    O F   T H E   C O R E', &
            ('-',i=1,70)

      if (.not. nucleus%spherical_expansion) call warning(modulename, subroutinename, &
            "Spherical expansion not requested => <J^2>blocked  = 0")

      operator_name='J2'
       call expectation_value(nucleus, symmetries, basis, sp_states, &
            trim(operator_name), model, nucleus%J2_core, result2=nucleus%inertia, &
            result3=nucleus%gR_core)
       if (1.0_rk+4.0_rk*nucleus%J2_core<=0.0) call fatal(modulename, &
            subroutinename, '<J^2>core is not positive')
       if (trim(model%pairing) /= 'HTDA') then
          Erot = 0.5_rk/(nucleus%inertia*1.32_rk)*(2*nucleus%Jz(0)-nucleus%J2_core)
          write(out_unit,'(a,f9.3,1x,a/,/70a1/)') &
               'Rotational energy correction for bandhead hbar^2/(2*I*1.32)*(2*K - <J^2>core) =', &
               Erot,'MeV',('=',i=1,70)
       end if
       
    end if
    
    !-----------------------------------------------------------------------
    !
    ! 5) Magnetic moments
    !
    !-----------------------------------------------------------------------

    if (nucleus%magnetic_moment) then

       iunit=10
       open(unit=iunit,file='magnetic_moments',status='unknown')
       write(iunit,'(50a1//,20x,a//,50a1)') ('-',i=1,50),'Magnetic moments',('-',i=1,50)

       !
       !   Expectation values of Sz, Lz and mu_z for each nucleon type
       !

       write(out_unit,'(/70a1/,/a/,/70a1)') ('=',i=1,70), &
            '                M A G N E T I C   M O M E N T', &
            ('-',i=1,70)

       do isospin=1,2

          write(iunit,'(/6(a,2x)/)') 'Nucleon','state #','energy',' Omega pi',' <i|s_z|i>', &
               '<i|s_z|i> + <i~|s_z|i~>'

          avg_sz=0.0
          mean_particle = 0.0
          mean_hole = 0.0
          count_hole = 0
          count_particle = 0
          if (Npart(isospin) >= 2) then
             fermi = sp_states(Npart(isospin) - mod(Npart(isospin),2),isospin)%energy
          else
             fermi = sp_states(1,isospin)%energy
          end if
          Emin =  fermi - 0.5_rk * model%UE_cutoff(isospin)
          Emax =  fermi + 0.5_rk * model%UE_cutoff(isospin)

          do iHF=1,basis%size
             iblock=sp_states(iHF,isospin)%block
             shift=sum(basis%block_size(1:iblock))-basis%block_size(iblock)
             tmp=0.0
             do i=1,basis%block_size(iblock)
                tmp = tmp + sp_states(iHF,isospin)%coef%cylHO(i+shift)**2*basis%ns(i+shift)
             end do
             sp_states(iHF,isospin)%avg_sz=0.5_rk*tmp
             if (sp_states(iHF,isospin)%Omega > 0) then
                jHF = sp_states(iHF,isospin)%pair_partner(isospin)
                if (jHF == 0) stop 'partner not found'
                iblock=sp_states(jHF,isospin)%block
                shift=sum(basis%block_size(1:iblock))-basis%block_size(iblock)
                tmp=0.0
                do i=1,basis%block_size(iblock)
                   tmp = tmp + sp_states(jHF,isospin)%coef%cylHO(i+shift)**2* &
                        basis%ns(i+shift)
                end do
                tmp = 0.5_rk*tmp
                if (sp_states(iHF,isospin)%energy <= fermi .and. &
                     sp_states(iHF,isospin)%energy >= Emin) then
                   mean_hole = mean_hole + sp_states(iHF,isospin)%avg_sz+tmp
                   count_hole = count_hole + 1
                else if (sp_states(iHF,isospin)%energy > fermi .and. &
                     sp_states(iHF,isospin)%energy <= Emax) then
                   mean_particle = mean_particle + sp_states(iHF,isospin)%avg_sz+tmp
                   count_particle = count_particle + 1
                end if
                write(iunit,'(3x,a1,4x,i4,2x,f12.6,2x,i3,a2,1x,a1,1x,f10.6,5x,f10.6)') &
                     nucleon_type(isospin),iHF,sp_states(iHF,isospin)%energy, &
                     sp_states(iHF,isospin)%Omega,'/2',signparity(sp_states(iHF,isospin)%pi), &
                     sp_states(iHF,isospin)%avg_sz,sp_states(iHF,isospin)%avg_sz+tmp
             else
                write(iunit,'(3x,a1,4x,i4,2x,f12.6,2x,i3,a2,1x,a1,1x,f10.6)') &
                     nucleon_type(isospin),iHF,sp_states(iHF,isospin)%energy, &
                     sp_states(iHF,isospin)%Omega,'/2',signparity(sp_states(iHF,isospin)%pi), &
                     sp_states(iHF,isospin)%avg_sz
             end if
             avg_sz=avg_sz+sp_states(iHF,isospin)%occupation*sp_states(iHF,isospin)%avg_sz
          end do

          if (count_hole /= 0) mean_hole = mean_hole/count_hole
          if (count_particle /= 0) mean_particle = mean_particle/count_particle

          rms_particle = 0.0
          rms_hole = 0.0

          do iHF=1,basis%size
             if (sp_states(iHF,isospin)%Omega > 0) then
                jHF = sp_states(iHF,isospin)%pair_partner(isospin)
                iblock=sp_states(jHF,isospin)%block
                shift=sum(basis%block_size(1:iblock))-basis%block_size(iblock)
                tmp=0.0
                do i=1,basis%block_size(iblock)
                   tmp = tmp + sp_states(jHF,isospin)%coef%cylHO(i+shift)**2* &
                        basis%ns(i+shift)
                end do
                tmp = 0.5_rk*tmp
                if (sp_states(iHF,isospin)%energy <= fermi .and. &
                     sp_states(iHF,isospin)%energy >= Emin) then
                   rms_hole = rms_hole + (sp_states(iHF,isospin)%avg_sz+tmp - mean_hole)**2
                else if (sp_states(iHF,isospin)%energy > fermi .and. &
                     sp_states(iHF,isospin)%energy <= Emax) then
                   rms_particle = rms_particle + &
                        (sp_states(iHF,isospin)%avg_sz+tmp - mean_particle)**2
                end if
             end if
          end do

          write(iunit,'(4(a,3x))') 'Quasi-pairs','Number', &
               'Mean value of <i|s_z|i> + <i~|s_z|i~>','rms deviation'
          if (count_hole /= 0 .and. rms_hole >= 0.0) then
             rms_hole = sqrt(rms_hole/count_hole)
             write(iunit,'(a,10x,i4,33x,f10.6,2x,f10.6)') &
                  'hole',count_hole,mean_hole,rms_hole
          end if
          if (count_particle /= 0.and. rms_particle >= 0.0) then
             rms_particle = sqrt(rms_particle/count_particle)
             write(iunit,'(a,6x,i4,33x,f10.6,2x,f10.6)') &
                  'particle',count_particle,mean_particle,rms_particle
          end if

          s_odd(:)=0.0
          do i=1,nqp
             if (i_odd(isospin,i)/=0) s_odd(i)=sp_states(i_odd(isospin,i),isospin)%avg_sz
          end do
          avg_lz=0.5_rk*sum(nucleus%K(isospin,:))-avg_sz
          nucleus%Lz(isospin)=avg_lz
          nucleus%Jz(isospin)=0.5_rk*sum(nucleus%K(isospin,:))
          nucleus%mu_z(isospin)=gs(isospin)*avg_sz+gl(isospin)*avg_lz
          write(iunit,'(/a1)') nucleon_type(isospin)
          write(iunit,'(a,f10.6,2x,a,4(f10.6,1x),1x,a,f10.6)') &
               '   <Sz>=',avg_sz,'odd particles:',s_odd(1:nqp),'core:',avg_sz-sum(s_odd(:))
          write(iunit,'(a,f10.6,2x,a,4(f10.6,1x),1x,a,f10.6)') &
               '   <Lz>=',avg_lz,'odd particles:',0.5_rk*nucleus%K(isospin,1:nqp)-s_odd(1:nqp), &
               'core:',avg_lz-(nucleus%Jz(isospin)-sum(s_odd(:)))
          write(iunit,'(a,f10.6,2x,a,4(f10.6,1x),1x,a,f10.6)') &
               '   <Jz>=',nucleus%Jz(isospin),'odd particles:',0.5_rk*nucleus%K(isospin,1:nqp),&
               'core:',0.0_rk
          write(iunit,'(a,f10.6)') &
               '   <mu_z>=',nucleus%mu_z(isospin)
          write(out_unit,'(/a1)') nucleon_type(isospin)
          write(out_unit,'(a,f10.6,2x,a,4(f10.6,1x),1x,a,f10.6)') &
               '   <Sz>=',avg_sz,'odd particles:',s_odd(1:nqp),'core:',avg_sz-sum(s_odd(:))
          write(out_unit,'(a,f10.6,2x,a,4(f10.6,1x),1x,a,f10.6)') &
               '   <Lz>=',avg_lz,'odd particles:',0.5_rk*nucleus%K(isospin,1:nqp)-s_odd(1:nqp), &
               'core:',avg_lz-(nucleus%Jz(isospin)-sum(s_odd(:)))
          write(out_unit,'(a,f10.6,2x,a,4(f10.6,1x),1x,a,f10.6)') &
               '   <Jz>=',nucleus%Jz(isospin),'odd particles:',0.5_rk*nucleus%K(isospin,1:nqp),&
               'core:',0.0_rk
          write(out_unit,'(a,f10.6)') &
               ' <mu_z>=',nucleus%mu_z(isospin)
       end do

       !
       !   Total expectation values of Lz, Jz and mu_z
       !

       nucleus%Lz(0)=nucleus%Lz(1)+nucleus%Lz(2)
       nucleus%Jz(0)=nucleus%Jz(1)+nucleus%Jz(2)
       nucleus%mu_z(0)=nucleus%mu_z(1)+nucleus%mu_z(2)
       write(iunit,'(a,2x,4(a,f10.6,2x)/)') 'Total:','<Sz>=',nucleus%Sz(0),'<Lz>=',nucleus%Lz(0), &
            '<Jz>=',nucleus%Jz(0),'<mu_z>=',nucleus%mu_z(0)
       write(out_unit,'(/1(a,f10.6,2x)/)') 'Total <mu_z> =',nucleus%mu_z(0)

       !
       !   Global quenching factor of the spin gyromagnetic ratio
       !

       avg_sz = 0.0
       avg_lz = 0.0
       do isospin = 1,2
          do i=1,nqp
             if (i_odd(isospin,i) /= 0) then
                avg_sz = avg_sz + gs(isospin) * sp_states(i_odd(isospin,i),isospin)%avg_sz
                avg_lz = avg_lz + gl(isospin) * (0.5_rk*nucleus%K(isospin,i) - &
                     sp_states(i_odd(isospin,i),isospin)%avg_sz)
             end if
          end do
       end do
       if (abs(avg_sz) > 1e-6_rk) then
          quenching_factor = (nucleus%mu_z(0) - avg_lz)/avg_sz
       else
          quenching_factor = 1.0_rk
       end if
       write(iunit,'(a,f10.6/)') 'Global gs quenching factor=',quenching_factor
       close(iunit)
       write(out_unit,'(a,f10.6/)') 'Global gs quenching factor=',quenching_factor

       write(out_unit,'(a,f10.6)') 'Sum_blocked g_l^(q)*<Lz>blocked =', avg_lz
       write(out_unit,'(a,f10.6)') 'Sum_blocked g_s^(q)*<Sz>blocked =', avg_sz
       write(out_unit,'(a,f10.6/)') 'Sum_blocked [g_l^(q)*<Lz>blocked + ' // &
            'g_s^(q)*<Sz>blocked] =', avg_lz + avg_sz

       write(out_unit,'(a,f8.4)') 'Intrinsic magnetic moment:  mu_intr = K/(K+1) * <mu_z>  =', &
            nucleus%mu_z(0)*nucleus%Jz(0)/(nucleus%Jz(0)+1)
       write(out_unit,'(a,f8.4)') 'Collective magnetic moment: mu_coll = K/(K+1) * gR_core =', &
            nucleus%gR_core*nucleus%Jz(0)/(nucleus%Jz(0)+1)
       write(out_unit,'(a,f8.4)') 'Total magnetic moment:      mu      = mu_intr + mu_coll =', &
            (nucleus%mu_z(0)+nucleus%gR_core)*nucleus%Jz(0)/(nucleus%Jz(0)+1)

       write(out_unit,'(/4(a/))') &
            'REMARK: gR_core calculated as follows', &
            '        * Inglis-Belyaev formula excluding blocked states', &
            '        * assuming i-tilde = i-bar', &
            '        * assuming Eqp(k)=Eqp(k-tilde)'

       write(out_unit,'(70a1/)') ('=',i=1,70)

    end if

    !-----------------------------------------------------------------------
    !
    ! 6) Single-particle angular-momentum quantities
    !
    !-----------------------------------------------------------------------

    if (nucleus%sp_angular_momentum) then

       !
       !   Matrix elements of j^2 in HF basis for each isospin and eigenvalues
       !

       write(out_unit,'(/70a1/,/a//,2(/a/),/70a1/)') ('=',i=1,70), &
            '                A N G U L A R   M O M E N T U M', &
            '              S I N G L E - P A R T I C L E   j^2', &
            '  M A T R I X   E L E M E N T S   A N D   E I G E N V A L U E S', &
            ('-',i=1,70)

       do isospin = 1, 2

          j2matrix_HFbasis = 0
          do iHF = 1, basis%size
             if (sp_states(iHF,isospin)%Omega < 0) then
                j2matrix_HFbasis(iHF,iHF) = 1e6_rk
             end if
             do jHF = 1, iHF
                if (sp_states(jHF,isospin)%Omega < 0) cycle
                tmp = 0.0
                do i = 1, sph_basis%size
                   tmp = tmp + sp_states(iHF,isospin)%coef%sphHO(i) * &
                        sp_states(jHF,isospin)%coef%sphHO(i) &
                        * 0.5_rk*sph_basis%j2(i) * (0.5_rk*sph_basis%j2(i)+1)
                end do
                j2matrix_HFbasis(iHF,jHF) = tmp
                j2matrix_HFbasis(jHF,iHF) = j2matrix_HFbasis(iHF,jHF)
             end do
          end do

          call diagon(j2matrix_HFbasis,eigenval_j2,eigenvec_j2,basis%size)

          if (isospin == 1) then
             write(out_unit,'(/a/)') 'Eigenvalues of j^2 in HF basis for neutrons:'
          else
             write(out_unit,'(/a/)') 'Eigenvalues of j^2 in HF basis for protons:'
          end if
          do i = 1, basis%size
             if (eigenval_j2(i) > 100.0_rk) cycle
             write(out_unit,'(2x,a,i4,1x,2(a,f9.3,1x),i4)') 'i=',i,'eigenvalue=',eigenval_j2(i), &
                  'j=',0.5_rk*(sqrt(1+4*eigenval_j2(i))-1),maxloc(abs(eigenvec_j2(:,i)))
          end do

       end do

       write(out_unit,*)
       write(out_unit,'(70a1/)') ('=',i=1,70)

       !
       !   Matrix elements of j^+ in HF basis for each isospin
       !

       write(out_unit,'(/70a1/,/a//,/a/,/70a1/)') ('=',i=1,70), &
            '                 A N G U L A R   M O M E N T U M', &
            'S I N G L E - P A R T I C L E   j+  M A T R I X  E L E M E N T S', &
            ('-',i=1,70)

       do isospin = 1, 2
          write(out_unit,'(/a1,1x,a)') nucleon_type(isospin),'matrix elements of j+'
          do iHF = 1, basis%size
             if (sp_states(iHF,isospin)%Omega < 1 .or. &
                  abs(sp_states(Npart(isospin),isospin)%energy - sp_states(iHF,isospin)%energy) > 5.0_rk) cycle
             do jHF = 1, basis%size
                if (sp_states(iHF,isospin)%Omega /= sp_states(jHF,isospin)%Omega + 2 .or. &
                     sp_states(iHF,isospin)%pi /= sp_states(jHF,isospin)%pi .or. &
                     abs(sp_states(Npart(isospin),isospin)%energy - sp_states(jHF,isospin)%energy) > 5.0_rk) cycle
                write(out_unit,'(2(a,i4,2x,f12.6,2x,i3,a2,1x,a1,2x),a,f12.6)') &
                     'iHF=',iHF,sp_states(iHF,isospin)%energy,sp_states(iHF,isospin)%Omega,'/2', &
                     signparity(sp_states(iHF,isospin)%pi), &
                     'jHF=',jHF,sp_states(jHF,isospin)%energy,sp_states(jHF,isospin)%Omega,'/2', &
                     signparity(sp_states(jHF,isospin)%pi), &
                     '< iHF | j+ | jHF > =',HF_Jplus([sp_states(iHF,isospin),sp_states(jHF,isospin)], &
                     basis)
             end do
          end do
       end do

       write(out_unit,*)
       write(out_unit,'(70a1/)') ('=',i=1,70)

    end if

    !-------------------------------------------------
    !
    ! 7) Isospin mixing
    !
    !-------------------------------------------------

    if (nucleus%isospin_mixing) then

       !
       !   Expectation value of T^2
       !
       
       write(out_unit,'(/70a1/,/a//,3(a/),/70a1/)') ('=',i=1,70), &
            '         I S O S P I N - M I X I N G   P A R A M E T E R', &
            '         <BCS|T^2|BCS> - T0(T0+1)', &
            ' alpha = ------------------------  with T0 = |Tz|', &
            '                  2(T0+1)', ('-',i=1,70)
       operator_name='T2'
       call expectation_value(nucleus, symmetries, basis, sp_states, &
            trim(operator_name), model, nucleus%T2)

       if (1.0_rk+4.0_rk*nucleus%T2<=0.0) call warning(modulename, &
            subroutinename, '<T^2> is not positive')
       nucleus%T = 0.5_rk * abs(nucleus%N-nucleus%Z)
       alpha = (nucleus%T2 - nucleus%T*(nucleus%T+1))/(2*(nucleus%T+1))

       write(out_unit,'(a,1x,f11.6,/a,1x,f11.6,/a,1x,f11.6,1x,a1//,70a1/)') &
            'T0    =',nucleus%T, &
            '<T^2> =',nucleus%T2, &
            'alpha =',alpha*100,'%', ('=',i=1,70)

       !
       !   Neutron-proton overlaps and analog partners
       !
       !   Sign convention: n and p corresponding states 
       !   (same rank in symmetry block) should have a positive 
       !   (and large enough) overlap
       !

       neutron_analog(:) = 0
       proton_analog(:) = 0
       allocate(overlap_matrix(basis%size,basis%size))
       overlap_matrix(:,:)=0.0

       open(unit=10,file='overlaps')

       do isospin=1,2
          do i=1,basis%size
             sp_states(i,isospin)%pair_partner(3-isospin) = 0
          end do
       end do
       
       do iHFp=1,basis%size

          iHFn_max=iHFp
          ovlp_max=1e-9_rk
          do iHFn=1,basis%size
             if (sp_states(iHFn,1)%block /= sp_states(iHFp,2)%block) cycle
             ovlp=overlap(sp_states(iHFn,1),sp_states(iHFp,2))
             if (abs(ovlp)>abs(ovlp_max)) then
                iHFn_max=iHFn
                ovlp_max=ovlp
             end if
          end do
          iHFn=iHFn_max
          ovlp=ovlp_max
          neutron_analog(iHFp) = iHFn
          proton_analog(iHFn) = iHFp
          if (sp_states(iHFp,2)%block_rank /= sp_states(iHFn,1)%block_rank .and. &
               abs(ovlp)<0.5_rk) then
             write(out_unit,*) 'iHFn=',iHFn,sp_states(iHFn,1)%Omega,'/2',&
                  signparity(sp_states(iHFn,1)%pi),' rank=',sp_states(iHFn,1)%block_rank
             write(out_unit,*) 'iHFp=',iHFp,sp_states(iHFp,2)%Omega,'/2',&
                  signparity(sp_states(iHFp,2)%pi),' rank=',sp_states(iHFp,2)%block_rank
             write(out_unit,*) 'overlap=',ovlp
             close(10)
          end if
          sp_states(iHFp,2)%pair_partner(1) = sp_states(neutron_analog(iHFp),1)%pair_partner(1)
          sp_states(iHFn,1)%pair_partner(2) = sp_states(proton_analog(iHFn),2)%pair_partner(2)

          sign_convention = sign(1.0_rk,ovlp)
          sp_states(iHFp,2)%coef%cylHO(:) = sp_states(iHFp,2)%coef%cylHO(:) &
               * sign_convention
          sp_states(iHFp,2)%coef%sphHO(:) = sp_states(iHFp,2)%coef%sphHO(:) &
               * sign_convention

       end do

       do iHFn=1,basis%size
          norm=0.0
          do iHFp=1,basis%size
             ovlp = overlap(sp_states(iHFn,1),sp_states(iHFp,2))
             overlap_matrix(iHFn,iHFp)=ovlp
             norm=norm+ovlp**2
             if (abs(ovlp) > 1e-3) then
                if (symmetries%parity) then
                   write(10,'(2(a,i4,a,i3,3(a)),f7.4)') &
                        '<',iHFn,'(',sp_states(iHFn,1)%Omega,'/2', &
                        signparity(sp_states(iHFn,1)%pi),')n ','|', &
                        iHFp,'(',sp_states(iHFp,2)%Omega,'/2', &
                        signparity(sp_states(iHFp,2)%pi),')p > = ',ovlp
                else
                   write(10,'(2(a,i4,a,i3,2(a)),e13.6)') &
                        '<',iHFn,'(',sp_states(iHFn,1)%Omega,'/2',')n ', &
                        '|',iHFp,'(',sp_states(iHFp,2)%Omega,'/2',')p > = ',ovlp
                end if
             end if
          end do
          if (symmetries%parity) then
             write(10,'(a,i4,a,i3,3(a),f7.4,a,1x,i4/)') &
                  '||',iHFn,'(',sp_states(iHFn,1)%Omega,'/2',&
                  signparity(sp_states(iHFn,1)%pi),')n ||^2 = ',norm,&
                  ' proton analog=',proton_analog(iHFn)
          else
             write(10,'(a,i4,a,i3,2(a),f7.4,a,1x,i4/)') &
                  '||',iHFn,'(',sp_states(iHFn,1)%Omega,'/2',')n ||^2 = ',norm,&
                  ' proton analog=',proton_analog(iHFn)
          end if
          if (abs(norm-1.0_rk) > 1e-3_rk) stop 'Norm should be 1'
       end do

       !
       !   det = -/+ 1 when n and p states are identical, 
       !   depending on the order in which the n and p 
       !   states appear in the matrix
       !

       if (debugmode /= 'low') then
          write(out_unit,'(a/)') 'Overlaps of n and p wave functions:'
          det=determinant(overlap_matrix,basis%size)
          write(out_unit,'(a,f9.6)') 'det(overlap_matrix)=',det

          signature_n=permutation_signature(neutron_analog,basis%size)
          signature_p=permutation_signature(proton_analog,basis%size)  
          if (signature_n /= signature_p) then
             write(out_unit,*) 'Signature(proton_analog)=',signature_p
             write(out_unit,*) 'Signature(neutron_analog)=',signature_n
          else if (det*signature_n < 0) then
             write(out_unit,*) 'Approximate value of the overlap determinant = ',signature_n
          end if
          write(out_unit,'(a,i3)') 'Approximate value of the overlap determinant =',signature_n
       end if
       
       deallocate(overlap_matrix)
       close(10)

    end if
    
    !-----------------------------------------------------------------------
    !
    ! 8) Center-of-mass correction
    !
    !-----------------------------------------------------------------------

    if (nucleus%two_body_center_of_mass_correction) then

       open(unit=20,file='center_of_mass_correction.txt',status='replace')

       write(20,'(/70a1/,/a/,/70a1)') ('=',i=1,70), &
            '           C E N T E R   O F   M A S S   C O R R E C T I O N', &
            ('-',i=1,70)

       write(20,'(/a/,a/)') 'WARNING: 2-body center of mass correction implemented', &
            'for a pure HF solution only'

       write(out_unit,'(/70a1/,/a/,/70a1)') ('=',i=1,70), &
            '           C E N T E R   O F   M A S S   C O R R E C T I O N', &
            ('-',i=1,70)

       write(out_unit,'(/a/,a/)') 'WARNING: 2-body center of mass correction implemented', &
            'for a pure HF solution only'

       Ekin1 = 0.0_rk
       Ekin2 = 0.0_rk
       Ekin2dir = 0.0_rk
       Ekin2exch = 0.0_rk
       r2MB = 0.0_rk
       do isospin = 1, 2
          PHF(:) = cmplx(0.0)
          do iHF = 1, basis%size
             v2 = sp_states(iHF,isospin)%occupation
             if (v2 < 1e-6) cycle
             HF_states(1) = sp_states(iHF,isospin)
             HF_states(2) = HF_states(1)
             p2 = HF_psqr(HF_states,basis)
             r2 = HF_xsqr_plus_ysqr(HF_states,basis)
             z2 = HF_zsqr(HF_states,basis)
             z = HF_z(HF_states,basis)
             Ekin1(isospin) = Ekin1(isospin) + v2 * p2*hb0(isospin)
             PHF(1) = PHF(1) + v2 * HF_px(HF_states,basis)
             PHF(2) = PHF(2) + v2 * HF_py(HF_states,basis)
             PHF(3) = PHF(3) + v2 * HF_pz(HF_states,basis)
             r2MB(isospin) = r2MB(isospin) + v2 * (r2+z2)
             write(20,'(i3,1x,a,f8.3,1x,i3,a2,a1,2(1x,a,f8.3))') &
                  iHF,'e =',sp_states(iHF,isospin)%energy,sp_states(iHF,isospin)%Omega,'/2', &
                  parity_sign(sp_states(iHF,isospin)%pi),'v2 =',v2,'<i|p^2/(2m)|i>=',p2*hb0(isospin)
             do jHF = 1, basis%size
                if (sp_states(jHF,isospin)%occupation < 1e-6) cycle
                HF_states(2) = sp_states(jHF,isospin)
                px = HF_px(HF_states,basis)
                py = HF_py(HF_states,basis)
                pz = HF_pz(HF_states,basis)
                if (abs(px)**2 + abs(py)**2 + abs(pz)**2 >= 1e-3_rk) &
                     write(20,'(2x,a,i3,1x,a,f8.3,1x,i3,a2,a1,3(1x,a,f8.3))') 'jHF=',jHF, &
                     'e =',sp_states(jHF,isospin)%energy,sp_states(jHF,isospin)%Omega,'/2', &
                     parity_sign(sp_states(jHF,isospin)%pi),'Im(px) =',aimag(px), &
                     'py=',py,'Im(pz)=',aimag(pz)
                Ekin2exch(isospin) = Ekin2exch(isospin) - v2 * sp_states(jHF,isospin)%occupation* &
                     (abs(px)**2 + abs(py)**2 + abs(pz)**2)
             end do
          end do
          Ekin2dir(isospin) = (abs(PHF(1))**2 + abs(PHF(2))**2 + abs(PHF(3))**2)* &
               hb0(isospin)/nucleus%A
          Ekin2exch(isospin) = Ekin2exch(isospin)*hb0(isospin)/nucleus%A
          Ekin2(isospin) = Ekin2dir(isospin) + Ekin2exch(isospin)
          write(20,'(a,2(1x,a,f9.3),1x,a)') nucleon_type(isospin), &
               '1-body kinetic energy (uncorrected) =',Ekin1(isospin),' MeV  <r^2>/nucleon=', &
               r2MB(isospin)/Npart(isospin),'fm^2'
          write(20,'(a,f10.3,1x,a)') nucleon_type(isospin) // &
               ' two-body center-of-mass kinetic energy =',Ekin2(isospin),'MeV'
          write(20,'(a,f10.3,1x,a)') '* direct contribution =',Ekin2dir(isospin),'MeV'
          write(20,'(a,f10.3,1x,a)') '* exchange contribution =',Ekin2exch(isospin),'MeV'
          write(out_unit,'(a,2(1x,a,f9.3),1x,a)') nucleon_type(isospin), &
               '1-body kinetic energy (uncorrected) =',Ekin1(isospin),' MeV  <r^2>/nucleon=', &
               r2MB(isospin)/Npart(isospin),'fm^2'
          write(out_unit,'(a,f10.3,1x,a)') nucleon_type(isospin) // &
               ' two-body center-of-mass kinetic energy =',Ekin2(isospin),'MeV'
          write(out_unit,'(a,f10.3,1x,a)') '* direct contribution =',Ekin2dir(isospin),'MeV'
          write(out_unit,'(a,f10.3,1x,a)') '* exchange contribution =',Ekin2exch(isospin),'MeV'
       end do

       write(20,'(/a,f10.3,1x,a)') 'Total kinetic energy of the nucleus: Ekin_tot =', &
            Ekin1(1)+Ekin1(2),'MeV'

       write(20,'(/a,f10.3,1x,a)') 'Total one-body center-of-mass kinetic energy: Ekin_cm(1) =', &
            (Ekin1(1)+Ekin1(2))/nucleus%A,'MeV'

       write(20,'(/a,f10.3,1x,a)') 'Total two-body center-of-mass kinetic energy: Ekin_cm(2) =', &
            Ekin2(1)+Ekin2(2),'MeV'
       write(20,'(a,f10.3,1x,a)') '* direct contribution =',Ekin2dir(1)+Ekin2dir(2),'MeV'
       write(20,'(a,f10.3,1x,a)') '* exchange contribution =',Ekin2exch(1)+Ekin2exch(2),'MeV'

       write(20,'(/a,/a,f10.3,1x,a)') 'Total kinetic energy of the nucleus ' // &
            'with 1-body c.o.m. correction:','Ekin_tot-Ekin_cm(1) =', &
            (Ekin1(1)+Ekin1(2))*(1.0_rk-1.0_rk/nucleus%A),'MeV'

       write(20,'(/a,/a,f10.3,1x,a)') 'Total kinetic energy of the nucleus ' // &
            'with 1-body+2-body c.o.m. corrections:','Ekin_tot-Ekin_cm(1)-Ekin_cm(2) =', &
            (Ekin1(1)+Ekin1(2))*(1.0_rk-1.0_rk/nucleus%A)-(Ekin2(1)+Ekin2(2)),'MeV'

       write(out_unit,'(/a,f10.3,1x,a)') 'Total kinetic energy of the nucleus: Ekin_tot =', &
            Ekin1(1)+Ekin1(2),'MeV'

       write(out_unit,'(/a,f10.3,1x,a)') 'Total one-body center-of-mass kinetic energy: Ekin_cm(1) =', &
            (Ekin1(1)+Ekin1(2))/nucleus%A,'MeV'

       write(out_unit,'(/a,f10.3,1x,a)') 'Total two-body center-of-mass kinetic energy: Ekin_cm(2) =', &
            Ekin2(1)+Ekin2(2),'MeV'
       write(out_unit,'(a,f10.3,1x,a)') '* direct contribution =',Ekin2dir(1)+Ekin2dir(2),'MeV'
       write(out_unit,'(a,f10.3,1x,a)') '* exchange contribution =',Ekin2exch(1)+Ekin2exch(2),'MeV'

       write(out_unit,'(/a,/a,f10.3,1x,a)') 'Total kinetic energy of the nucleus ' // &
            'with 1-body c.o.m. correction:','Ekin_tot-Ekin_cm(1) =', &
            (Ekin1(1)+Ekin1(2))*(1.0_rk-1.0_rk/nucleus%A),'MeV'

       write(out_unit,'(/a,/a,f10.3,1x,a)') 'Total kinetic energy of the nucleus ' // &
            'with 1-body+2-body c.o.m. corrections:','Ekin_tot-Ekin_cm(1)-Ekin_cm(2) =', &
            (Ekin1(1)+Ekin1(2))*(1.0_rk-1.0_rk/nucleus%A)-(Ekin2(1)+Ekin2(2)),'MeV'

       close(20)

       write(out_unit,*)
       write(out_unit,'(70a1/)') ('=',i=1,70)

    end if

    !-----------------------------------------------------------------------
    !
    ! 9) Spectroscopic factors
    !
    !-----------------------------------------------------------------------
    !
    if (nucleus%spectroscopic_factors) then

       !
       !   Set up particle number of parent nuclei
       !      * 1: neutron parent nucleus (n-capture equivalent reaction on even-even target)
       !      * 2: proton parent nucleus (p-capture equivalent reaction on even-even target)
       !

       parent_nucleus(1)%N=nucleus%N+1
       parent_nucleus(1)%Z=nucleus%Z
       parent_nucleus(1)%A=nucleus%A+1
       parent_nucleus(1)%name=Zsymbol(parent_nucleus(1)%Z) // '-' // &
            int2char(parent_nucleus(1)%A)

       parent_nucleus(2)%N=nucleus%N
       parent_nucleus(2)%Z=nucleus%Z+1
       parent_nucleus(2)%A=nucleus%A+1
       parent_nucleus(2)%name=Zsymbol(parent_nucleus(2)%Z) // '-' // &
            int2char(parent_nucleus(2)%A)

       !
       !   Calculating GS spectroscopic factors for present nucleus
       !   and parent nuclei (1 neutron more or 1 proton more, mass A+1)
       !

       if (nucleus%A >= 2) then
          call spectroscopic_factors(nucleus, basis, sp_states, parent_nucleus)
          write(out_unit,'(/70a1/,/a/,/70a1/)') ('=',i=1,70), &
               'G R O U N D - S T A T E   S P E C T R O S C O P I C   F A C T O R S', &
               ('-',i=1,70)
          write(out_unit,'(2(3a,f12.6,3x)//)') &
               'SF_n(',trim(nucleus%name),')=',nucleus%GS_spectroscopic_factor(1), &
               'SF_p(',trim(nucleus%name),')=',nucleus%GS_spectroscopic_factor(2)
          write(out_unit,'(2(3a,f12.6,3x)//)') &
               'SF_n(',trim(parent_nucleus(1)%name),')=',parent_nucleus(1)%GS_spectroscopic_factor(1),&
               'SF_p(',trim(parent_nucleus(1)%name),')=',parent_nucleus(1)%GS_spectroscopic_factor(2)
          write(out_unit,'(2(3a,f12.6,3x))') &
               'SF_n(',trim(parent_nucleus(2)%name),')=',parent_nucleus(2)%GS_spectroscopic_factor(1),&
               'SF_p(',trim(parent_nucleus(2)%name),')=',parent_nucleus(2)%GS_spectroscopic_factor(2)
          write(out_unit,'(/70a1/)') ('=',i=1,70)
          
       end if

!       call spherical_expansion_DSD(basis, sp_states, sph_basis)

    end if

    !-------------------------------
    !
    ! 10)  Single-particle sum rules
    !
    !-------------------------------

    if (nucleus%sp_sum_rules) then

       write(out_unit,'(/70a1/,/10x,a/,/70a1/)') ('=',i=1,70), &
            'S I N G L E - P A R T I C L E   S U M   R U L E S', &
            ('-',i=1,70)

       if (model%exact_coul_exch .and. .not. present(Coulomb_exchange_potential)) &
            stop 'model%exact_coul_exch = .true. but present(Coulomb_exchange_potential)=.false.'
       if (present(Coulomb_exchange_potential)) then
          call single_particle_sum_rules(nucleus, model, basis, HFfield, densities, sp_states, &
               Coulomb_exchange_potential)
       else
          call single_particle_sum_rules(nucleus, model, basis, HFfield, densities, sp_states)
       end if

       write(out_unit,'(/70a1/)') ('=',i=1,70)

       open(unit=15,file='spin_field.txt',status='old',action='write',access='append')
       write(15,*)
       do isospin = 1, 2
          ewsz(isospin) = 0
          loop_k: do k = 1, basis%size
             ewsz(isospin) = ewsz(isospin) + sp_states(k,isospin)%occupation * &
                  sp_states(k,isospin)%energy * sp_states(k,isospin)%avg_sz
          end do loop_k
          write(15,'(a1,2x,a,f8.3)') nucleon_type(isospin),'sum_k(v_k^2 e_k <s_z>_k) =', &
               ewsz(isospin)
       end do
       close(15)

    end if

    !-------------------------------
    !
    ! 11) Charge density
    !
    !-------------------------------

    if (nucleus%charge_density) then
       
       write(out_unit,'(/70a1/,/a/,/70a1)') ('=',i=1,70), &
            '           C H A R G E   D E N S I T Y', &
            ('-',i=1,70)

       call charge_density(densities(2),symmetries,basis,nucleus)

       write(out_unit,'(/,70a1/)') ('=',i=1,70)

    end if

    !------------------------------------
    !
    ! 12) Integral of current densities
    !
    !------------------------------------

    if (nucleus%integral_currents) then
       
       write(out_unit,'(/70a1/,2(/a/),/70a1)') ('=',i=1,70), &
            '   I N T E G R A L   O V E R   S P A C E', &
            '  O F   C U R R E N T   D E N S I T I E S', &
            ('-',i=1,70)

       call integral_currents(densities,basis,nucleus)

       write(out_unit,'(/,70a1/)') ('=',i=1,70)

    end if

  end subroutine observables


  !=====================================================================
  !                       SUBROUTINE SPECTROSCOPIC_FACTORS
  !---------------------------------------------------------------------
  ! Calculates the spectroscopic factors in the BCS approximation
  ! associated to the single-particle states passed as argument 3
  ! 'sp_states'.
  !=====================================================================

  subroutine spectroscopic_factors(nucleus, basis, sp_states, parent_nucleus)

!
!   Arguments
!
    type(nucleus_properties),intent(inout)                    :: nucleus
    type(cylindrical_basis),intent(in)                        :: basis
    type(sp_state),dimension(basis%size,2),intent(in) :: sp_states
    type(nucleus_properties),dimension(2),intent(inout)       :: parent_nucleus
!
!   Local variables
!
    character(len=21),parameter :: subroutinename='spectroscopic_factors'
    integer,dimension(2) :: Npart
    integer :: isospin,iHF,i
    real(kind=rk) :: j,v2,u2
    character(len=150) :: msg1,msg2

!
!   Debugging note
!

    if (debugmode=='medium'.or.debugmode=='high') call debugnote(modulename, &
         subroutinename)

!
!   Auxiliary variables
!

    Npart(:)=(/ nucleus%N , nucleus%Z /)

!
!   GS neutron and proton spectroscopic factors
!

    nucleus%GS_spectroscopic_factor(:)=0.0

    do isospin=1,2

       if (mod(Npart(isospin),2)==0) then

!
!         GS spectroscopic factor for the present even-even nucleus
!

          iHF=Npart(isospin)/2
          j=sp_states(iHF,isospin)%avg_j
          if (abs(j)<1e-3_rk) then
             msg1='Need single-particle <j> to calculate spectroscopic factors (set to 0)'
             msg2='Call cyl2sph first'
             call warning(modulename, subroutinename, trim(msg1), trim(msg2))
          end if
          v2=sp_states(iHF,isospin)%occupation
          nucleus%GS_spectroscopic_factor(isospin)=(2.0_rk*j+1.0_rk)*v2

!
!         GS spectroscopic factor for the neutron and proton parent odd nuclei
!         Approximation: GS determined by the lowest qp created on the even-even core
!         Proton spectroscopic factor of neutron parent nucleus = that of even-even 
!         nucleus, and conversely
!

          iHF=nucleus%pairing(isospin)%min_qp_index
          if (iHF < 0 .or. iHF > ubound(sp_states,dim=2)-1) then
             iHF = 0
             call warning(modulename, subroutinename,'Unfound min_qp_index')
          end if
          v2=sp_states(iHF+1,isospin)%occupation
          u2=1.0_rk-v2
          do i=1,2
             parent_nucleus(i)%GS_spectroscopic_factor(isospin)= &
                  kronecker(isospin,i)*u2 + &
                  (1.0_rk-kronecker(isospin,i))*nucleus%GS_spectroscopic_factor(isospin)
          end do

       else

          msg1='Odd ' // ' number'
          msg2=particle_type(isospin) // 'spectroscopic factor for this nucleus ' // &
               'require to calculate the GS solution for the even core'
          call warning(modulename, subroutinename, trim(msg1), trim(msg2))
          cycle

       end if

    end do

    return

  end subroutine spectroscopic_factors


  !=====================================================================
  !                       SUBROUTINE EXPECTATION_VALUE
  !---------------------------------------------------------------------
  ! Calculates the expectation value of the operator given as argument 4.
  !=====================================================================

  subroutine expectation_value(nucleus, symmetries, basis, sp_states, &
       operator_name, model, expval, result2, result3)

    implicit none

!
!   Arguments
!
    type(nucleus_properties),intent(inout)               :: nucleus
    type(symmetry_properties),intent(in)              :: symmetries
    type(cylindrical_basis),intent(in)                :: basis
    type(sp_state),dimension(basis%size,2),intent(in) :: sp_states
    character(len=*),intent(in)                       :: operator_name
    type(modelization),intent(in)                     :: model
    real(kind=rk),intent(out)                         :: expval
    real(kind=rk),intent(out),optional                :: result2, result3

!
!   Local variables
!
    character(len=17),parameter :: subroutinename='expectation_value'
    real(kind=rk),parameter :: tol = 1.0e-4_rk
    integer :: isospin,iHF,jHF,i,k,l,Omega_k,Omega_l
    integer :: count
    integer,dimension(2) :: blocked_state
    real(kind=rk) :: sum_v4,sum_overlap,v2i,v2j,ovlp
    real(kind=rk) :: e,Eqp_k,Eqp_k_tilde,Eqp_l,Eqp_l_tilde,uk,vk,ul,vl,ui,vi,tmp,tmp_W
    real(kind=rk) :: jp,jpb,jpb1,jpb2,jp1,jp2,jp3,jp4,sp,sp1,sp2,sp3,sp4,spb,spb1,spb2
    real(kind=rk) :: sx1,sx2,sx3,sx4,sxt1,sxt2,stx3,stx4,spt1,spt2,spt3,spt4
    real(kind=rk) :: jx1,jx2,jx3,jx4,jxt1,jxt2,jtx3,jtx4,jpt1,jpt2,jpt3,jpt4
    real(kind=rk) :: J2,J2core,J2blocked,J2crossterm,J2core_crossterm,gR
    real(kind=rk) :: J2intr,J2polblock,J2z,J2p
    type(expansion) :: psi_bar
    real(kind=rk),dimension(2) :: W,Wbel,Wcr
    type(sp_state), dimension(2) :: HF_states
    real(kind=rk),dimension(2) :: Icr,Ibel
    logical :: blocked_k, blocked_l


    select case(operator_name)
    case('T2')

!-----------------------------------
!
!      Isospin: <BCS|T^2|BCS>
!
!-----------------------------------

!      Sum of v^4

       sum_v4=0.0
       do isospin=1,2
          do iHF=1,basis%size
             sum_v4=sum_v4+sp_states(iHF,isospin)%occupation**2
          end do
       end do
       sum_v4=0.5_rk*sum_v4

!      Overlap between neutron and proton sp wave functions

       sum_overlap=0.0
       do iHF=1,basis%size
          v2i=sp_states(iHF,1)%occupation
          do jHF=1,basis%size
             v2j=sp_states(jHF,2)%occupation
             if (sp_states(jHF,2)%Omega * sp_states(iHF,1)%Omega < 0 .or. &
                  sp_states(jHF,2)%Pi * sp_states(iHF,1)%Pi == -1) cycle
             sum_overlap=sum_overlap+v2i*v2j* &
                  dot_product(sp_states(iHF,1)%coef%cylHO(:), &
                  sp_states(jHF,2)%coef%cylHO(:))**2
          end do
       end do

!      <BCS|T^2|BCS>

       expval=real(nucleus%N+nucleus%Z)+0.25_rk*real(nucleus%N-nucleus%Z)**2 &
            -(sum_v4+sum_overlap)

    case ('J2')

!-----------------------------------
!
!  Angular momentum: <BCS|J^2|BCS>
!
!-----------------------------------

       do isospin = 1, 2

          count = 0
          blocked_state(isospin) = 0
          do i = 1, nqp
             if (i_odd(isospin,i) /= 0) then
                count = count + 1
                blocked_state(isospin) = i_odd(isospin,i)
             end if
          end do

          if (count > 1) then
             write(out_unit,'(a)') 'WARNING: Blocked configuration of seniority > 1 for ' &
                  // trim(adjustl(particle_type(isospin))) // 's ==> <J^2> incomplete'
          end if

       end do

!
!      1) Calculation assuming i-tilde = i-bar
!

       if (trim(model%approx) == 'min_pol' .or. &
            trim(model%approx) == 'all') then

       write(out_unit,'(/70a1)') ((/'+'/),i=1,70)
       write(out_unit,'(a)') '1) Expectation value of J^2 assuming i-tilde = i-bar'
       write(out_unit,'(a)') 'WARNING: the blocked state is assumed to have Omega > 0'
       write(out_unit,'(70a1)') ((/'+'/),i=1,70)
       gR = 0.0_rk
       Icr(:)=0.0
       Ibel(:)=0.0
       Wcr(:) = 0.0
       Wbel(:) = 0.0
       J2 = 0.0_rk
       J2core = 0.0_rk
       J2crossterm = 0.0_rk

       do isospin = 1, 2
          
          do k = 1, basis%size

             Omega_k = sp_states(k,isospin)%Omega
             blocked_k = .false.
             do i = 1, nqp
                if (i_odd(isospin,i) == 0) cycle
                if (k == i_odd(isospin,i)) blocked_k = .true.
             end do
             if (Omega_k < 0 .or. blocked_k) cycle

             e = sp_states(k,isospin)%energy
             vk = sqrt(sp_states(k,isospin)%occupation)
             if (sp_states(k,isospin)%occupation > 1.0_rk+tol) then
                stop 'v_k^2 > 1'
             end if
             uk = sqrt(abs(1.0_rk-sp_states(k,isospin)%occupation))
             if (sp_states(k,isospin)%occupation < 1e-6_rk .or. &
                  1-sp_states(k,isospin)%occupation < 1e-6_rk) then
                Eqp_k = sqrt((e-nucleus%pairing(isospin)%chemical_potential)**2 + &
                     nucleus%pairing(isospin)%gap(k)**2)
             else
                Eqp_k = nucleus%pairing(isospin)%gap(k)/(2*uk*vk)
             end if
!             Eqp_k = sqrt((e-nucleus%pairing(isospin)%chemical_potential)**2 + &
!                  nucleus%pairing(isospin)%average_gap**2)

             i = sp_states(k,isospin)%pair_partner(isospin)
             if (e <= nucleus%pairing(isospin)%chemical_potential+6 .and. &
                     debugmode=='high') then
             write(out_unit,'(a,i2,1x,a,i4,1x,a,f9.3,1x,i3,a,1x,5(a,f6.3,1x),a,i4,2(1x,a,f9.3))') &
                  'isospin=',isospin,'iHF=',k, &
                  'e_i=',e,sp_states(k,isospin)%Omega,'/2' // &
                  parity_sign(sp_states(k,isospin)%pi),'v_i=',vk,&
                  'Delta_i=',nucleus%pairing(isospin)%gap(k), &
                  '<Delta>=',nucleus%pairing(isospin)%average_gap,&
                  'Eqp_i=',sqrt((e-nucleus%pairing(isospin)%chemical_potential)**2 + &
                  nucleus%pairing(isospin)%gap(k)**2), &
                  'Delta_i/(2*u_i*v_i)=',nucleus%pairing(isospin)%gap(k)/(2*uk*vk), &
                  'i-tilde=',i,'e_i-tilde=',sp_states(i,isospin)%energy, &
                  'Delta_i-tilde=',nucleus%pairing(isospin)%gap(i)
             end if
             e = sp_states(i,isospin)%energy
!             Eqp_k_tilde = sqrt((e-nucleus%pairing(isospin)%chemical_potential)**2 + &
!                  nucleus%pairing(isospin)%average_gap**2)
             if (sp_states(i,isospin)%occupation < 1e-6_rk .or. &
                  1-sp_states(i,isospin)%occupation < 1e-6_rk) then
                Eqp_k_tilde = sqrt((e-nucleus%pairing(isospin)%chemical_potential)**2 + &
                     nucleus%pairing(isospin)%gap(i)**2)
             else
                vi = sqrt(sp_states(i,isospin)%occupation)
                if (sp_states(i,isospin)%occupation > 1.0_rk+tol) then
                   stop 'v_{k-tilde}^2 > 1'
                end if
                ui = sqrt(abs(1-sp_states(i,isospin)%occupation))
                Eqp_k_tilde = nucleus%pairing(isospin)%gap(i)/(2*ui*vi)
             end if

             if (count == 1 .and. blocked_state(isospin) /= 0) then ! implemented for seniority-1 only

                HF_states(1) = sp_states(k,isospin)
                HF_states(2) = sp_states(blocked_state(isospin),isospin)
                jp1 = HF_Jplus(HF_states,basis)

                HF_states(1) = sp_states(blocked_state(isospin),isospin)
                HF_states(2) = sp_states(k,isospin)
                jp2 = HF_Jplus(HF_states,basis)

                jpb = 0
                if (Omega_k == 1 .and. sp_states(blocked_state(isospin),isospin)%Omega == 1) &
                     jpb = HF_Jplus_time_reversal(HF_states,basis)

                J2crossterm = J2crossterm + vk**2 * (jp1**2 + jp2**2 + jpb**2)

             end if

             do l = 1, basis%size

                Omega_l = sp_states(l,isospin)%Omega
                blocked_l = .false.
                do i = 1, nqp
                   if (i_odd(isospin,i) == 0) cycle
                   if (l == i_odd(isospin,i)) blocked_l = .true.
               end do
               if (Omega_l < 0 .or. blocked_l) cycle

                e = sp_states(l,isospin)%energy
                vl = sqrt(sp_states(l,isospin)%occupation)
                if (sp_states(l,isospin)%occupation > 1.0_rk+tol) then
                   write(6,'(a,i4,1x,a,f9.3,1x,a,f9.6,1x,a,e15.7)') &
                        'l=',l,'e=',e,'occ=',sp_states(l,isospin)%occupation, &
                        'v_l=',vl
                   stop 'v_l^2 > 1'
                end if
                ul = sqrt(abs(1.0_rk-sp_states(l,isospin)%occupation))
                if (sp_states(l,isospin)%occupation < 1e-6_rk .or. &
                     1-sp_states(l,isospin)%occupation < 1e-6_rk) then
                   Eqp_l = sqrt((e-nucleus%pairing(isospin)%chemical_potential)**2 + &
                        nucleus%pairing(isospin)%gap(l)**2)
                else
                   Eqp_l = nucleus%pairing(isospin)%gap(l)/(2*ul*vl)
                end if
!                Eqp_l = sqrt((e-nucleus%pairing(isospin)%chemical_potential)**2 + &
!                     nucleus%pairing(isospin)%average_gap**2)

                i = sp_states(l,isospin)%pair_partner(isospin)
                e = sp_states(i,isospin)%energy
!                Eqp_l_tilde = sqrt((e-nucleus%pairing(isospin)%chemical_potential)**2 + &
!                     nucleus%pairing(isospin)%average_gap**2)
                if (sp_states(i,isospin)%occupation < 1e-6_rk .or. &
                     1-sp_states(i,isospin)%occupation < 1e-6_rk) then
                   Eqp_l_tilde = sqrt((e-nucleus%pairing(isospin)%chemical_potential)**2 + &
                        nucleus%pairing(isospin)%gap(i)**2)
                else
                   vi = sqrt(sp_states(i,isospin)%occupation)
                   if (sp_states(i,isospin)%occupation > 1.0_rk+tol) then
                      stop 'v_{l-tilde}^2 > 1'
                   end if
                   ui = sqrt(abs(1-sp_states(i,isospin)%occupation))
                   Eqp_l_tilde = nucleus%pairing(isospin)%gap(i)/(2*ui*vi)
                end if

                HF_states(1) = sp_states(k,isospin)
                HF_states(2) = sp_states(l,isospin)
                jp = HF_Jplus(HF_states,basis)
                sp = HF_Splus(HF_states,basis)
                tmp = jp**2
                tmp_W = sp * jp
                if (abs(tmp_W*(uk*vl-ul*vk)**2/(Eqp_k + Eqp_l)) > 0.1_rk .and. &
                     debugmode=='high') then
                   write(out_unit,'(a,i2,1x,2(a,i4,1x,a,f12.6,1x,i3,a,1x),a,f12.6)') &
                        'isospin=',isospin,'k=',k,'e=',sp_states(k,isospin)%energy, &
                        sp_states(k,isospin)%Omega,'/2'//parity_sign(sp_states(k,isospin)%pi), &
                        'l=',l,'e=',sp_states(l,isospin)%energy, &
                        sp_states(l,isospin)%Omega,'/2'//parity_sign(sp_states(l,isospin)%pi), &
                        '<l|s-|k>*<k|j+|l>/(Ek+El)*(uk*vl-ul*vk)**2=', &
                        tmp_W*(uk*vl-ul*vk)**2/(Eqp_k + Eqp_l)
                end if

                jpb = 0.0
                spb = 0.0
                if (Omega_k == 1 .and. Omega_l == 1) then
                   jpb = HF_Jplus_time_reversal(HF_states,basis)
                   spb = HF_Splus_time_reversal(HF_states,basis)
                   tmp = tmp + 0.5_rk * jpb**2
                   tmp_W = tmp_W + 0.5_rk * spb * jpb
                end if

                tmp = tmp * (uk*vl-ul*vk)**2
                tmp_W = tmp_W * (uk*vl-ul*vk)**2

                J2core = J2core + tmp
                Ibel(isospin)=Ibel(isospin)+tmp/(Eqp_k + Eqp_l)
                Wbel(isospin) = Wbel(isospin) + tmp_W/(Eqp_k + Eqp_l)

                Icr(isospin) = Icr(isospin) + (uk*vl-ul*vk)**2*(jp**2/(Eqp_k + Eqp_l_tilde) + &
                     0.25_rk * jpb**2 * (1.0_rk/(Eqp_k+Eqp_l)+1.0_rk/(Eqp_k_tilde+Eqp_l_tilde)))

                HF_states(1) = sp_states(l,isospin)
                HF_states(2) = sp_states(k,isospin)
                jp1 = HF_Jplus(HF_states,basis) ! <l|j+|k>
                sp1 = HF_Splus(HF_states,basis) ! <l|s+|k>
                spb1 = HF_Sminus_time_reversal(HF_states,basis)! <l-bar|s+|k> = <k|s-|l-bar>
                tmp_W = 0.25_rk*(uk*vl - ul*vk)**2 * ((sp + sp1) * (jp + jp1) / (Eqp_k+Eqp_l_tilde) + &
                     0.5_rk*(spb + spb1) * jpb * (1.0_rk/(Eqp_k+Eqp_l)+1.0_rk/(Eqp_k_tilde+Eqp_l_tilde)))
                Wcr(isospin) = Wcr(isospin) + 2 * tmp_W

             end do

          end do

       end do

       J2core = J2core + J2crossterm
       write(out_unit,'(/a,f9.3)') '<J^2>core      =',J2core
       if (count <= 1) write(out_unit,'(a,f9.3)') '<J^2>crossterm =',-J2crossterm
       j2blocked = 0
       do isospin = 1, 2
          do i = 1, nqp
             if (i_odd(isospin,i)>0) then
                iHF = i_odd(isospin,i)
                j2blocked = j2blocked + sp_states(iHF,isospin)%avg_j * (sp_states(iHF,isospin)%avg_j+1.0_rk)
!                write(out_unit,'(a,i3,1x,i3,a2,a1,2x,a,f9.3,1x,a)') 'blocked state #',i,&
!                     sp_states(iHF,isospin)%Omega,'/2',parity_sign(sp_states(iHF,isospin)%pi), &
!                     '<J^2>=',sp_states(iHF,isospin)%avg_j * (sp_states(iHF,isospin)%avg_j+1.0_rk)
             end if
          end do
       end do
       write(out_unit,'(a,f9.3)') '<J^2>blocked   =',j2blocked
       J2 = J2core - J2crossterm + j2blocked
       if (count <= 1) write(out_unit,'(a,f9.3)') '==> <J^2>      =',J2

       write(out_unit,'(/70a1)') ((/'-'/),i=1,70)
       write(out_unit,'(/a)') 'Assuming Eqp(k)=Eqp(k-tilde):'
       write(out_unit,'(/70a1)') ((/'-'/),i=1,70)
       write(out_unit,'(/a,/a)') 'Moment of inertia of the core in the cranking approximation with BCS pairing:', &
            '(excluding blocked state and its conjugate)'
       write(out_unit,'(3(a,f10.6,2x))') 'Icr(n)=',Ibel(1),'Icr(p)=',Ibel(2),'Icr(tot)=',Ibel(1)+Ibel(2)
       write(out_unit,'(a,f10.3,1x,a)') '==> hbar^2/(2*Icr(tot))=',0.5_rk/(Ibel(1)+Ibel(2))*1000.0_rk,'keV'
       write(out_unit,'(/a,/a)') 'Moment of inertia in the cranking approximation with BCS pairing and ' &
            // 'Thouless-Valatin correction:','(excluding blocked state and its conjugate)'
       write(out_unit,'(a,f10.3,1x,a)') 'hbar^2/(2*Icr(tot)*1.32)=',0.5_rk/(Ibel(1)+Ibel(2))*1000.0_rk/1.32_rk,'keV'
       write(out_unit,'(/a,/a)') 'Excitation energy of the first 2+ state in the cranking approximation with BCS pairing and ' &
            // 'Thouless-Valatin correction:','(excluding blocked state and its conjugate)'
       write(out_unit,'(a,f10.3,1x,a)') 'E(2+)=',3.0_rk/(Ibel(1)+Ibel(2))*1000.0_rk/1.32_rk,'keV'

       gR = (Ibel(2) + (gs(2)-1)*Wbel(2) + gs(1)*Wbel(1))/(Ibel(1)+Ibel(2))
       write(out_unit,'(/a,/a)') 'Collective gyromagnetic ratio of the core in the cranking approximation with BCS pairing:', &
            '(excluding blocked state and its conjugate)'
       write(out_unit,'(2(a,f10.6,2x))') 'W(n)=',Wbel(1),'W(p)=',Wbel(2)
       write(out_unit,'(a,f10.3)') 'gR = ',gR
       result3 = gR

       if (trim(model%pairing) == 'HTDA') then
          expval = J2
          result2 = Icr(1) + Icr(2)
          return
       end if

       write(out_unit,'(/70a1)') ((/'-'/),i=1,70)
       write(out_unit,'(/a)') 'Without assuming Eqp(k)=Eqp(k-tilde):'
       write(out_unit,'(/70a1)') ((/'-'/),i=1,70)
       write(out_unit,'(/a,/a)') 'Moment of inertia of the core in the cranking approximation with BCS pairing:', &
            '(excluding blocked state and its conjugate)'
       write(out_unit,'(3(a,f10.6,2x))') 'Icr(n)=',Icr(1),'Icr(p)=',Icr(2),'Icr(tot)=',Icr(1)+Icr(2)
       write(out_unit,'(a,f10.3,1x,a)') '==> hbar^2/(2*Icr(tot))=',0.5_rk/(Icr(1)+Icr(2))*1000.0_rk,'keV'
       write(out_unit,'(/a,/a)') 'Moment of inertia in the cranking approximation with BCS pairing and ' &
            // 'Thouless-Valatin correction:','(excluding blocked state and its conjugate)'
       write(out_unit,'(a,f10.3,1x,a)') 'hbar^2/(2*Icr(tot)*1.32)=',0.5_rk/(Icr(1)+Icr(2))*1000.0_rk/1.32_rk,'keV'

       gR = (Icr(2) + (gs(2)-1)*Wcr(2) + gs(1)*Wcr(1))/(Icr(1)+Icr(2))
       write(out_unit,'(/a,/a)') 'Collective gyromagnetic ratio of the core in the cranking approximation with BCS pairing:', &
            '(excluding blocked state and its conjugate)'
       write(out_unit,'(2(a,f10.6,2x))') 'W(n)=',Wcr(1),'W(p)=',Wcr(2)
       write(out_unit,'(a,f10.3)') 'gR = ',gR

       end if

!
!      2) Calculation without assuming i-tilde = i-bar
!
       if (trim(model%approx) == 'full_pol' .or. &
            trim(model%approx) == 'all') then

       write(out_unit,'(/70a1)') ((/'+'/),i=1,70)
       write(out_unit,'(a)') '2) Expectation value of J^2 without assuming i-tilde = i-bar'
       write(out_unit,'(a)') 'WARNING: the blocked state is assumed to have Omega > 0 (1qp only)'
       write(out_unit,'(70a1)') ((/'+'/),i=1,70)
       Icr(:)=0.0
       Wcr(:) = 0.0
       J2 = 0.0_rk
       J2core = 0.0_rk
       J2crossterm = 0.0_rk
       J2core_crossterm = 0.0_rk

       do isospin = 1, 2

          do k = 1, basis%size

             Omega_k = sp_states(k,isospin)%Omega
             blocked_k = .false.
             do i = 1, nqp
                if (i_odd(isospin,i) == 0) cycle
                if (k == i_odd(isospin,i)) blocked_k = .true.
             end do
             if (Omega_k < 0 .or. blocked_k) cycle

             e = sp_states(k,isospin)%energy
!             Eqp_k = sqrt((e-nucleus%pairing(isospin)%chemical_potential)**2 + &
!                  nucleus%pairing(isospin)%average_gap**2)
             Eqp_k = sqrt((e-nucleus%pairing(isospin)%chemical_potential)**2 + &
                  nucleus%pairing(isospin)%gap(k)**2)
             vk = sqrt(sp_states(k,isospin)%occupation)
             uk = sqrt(1.0_rk-sp_states(k,isospin)%occupation)

             if (count == 1 .and. blocked_state(isospin) /= 0) then

                HF_states(1) = sp_states(blocked_state(isospin),isospin)
                HF_states(2) = sp_states(k,isospin)
                jp1 = HF_Jplus(HF_states,basis) ! <alpha|j+|k>

                i = sp_states(blocked_state(isospin),isospin)%pair_partner(isospin)
                HF_states(1) = sp_states(i,isospin)
                i = sp_states(k,isospin)%pair_partner(isospin)
                HF_states(2) = sp_states(i,isospin)
                jp2 = HF_Jplus(HF_states,basis) ! <alpha-tilde|j+|k-tilde>

                jpb = 0.0
                if (Omega_k == 1 .and. sp_states(blocked_state(isospin),isospin)%Omega == 1) then
                   HF_states(1) = sp_states(blocked_state(isospin),isospin)
                   i = sp_states(k,isospin)%pair_partner(isospin)
                   HF_states(2) = sp_states(i,isospin)
                   jpb = HF_Jplus(HF_states,basis) ! <alpha|j+|k-tilde>
                end if

                J2core_crossterm = J2core_crossterm + vk**2 * (jp1**2 + jp2**2 + jpb**2)

                HF_states(1) = sp_states(k,isospin)
                HF_states(2) = sp_states(blocked_state(isospin),isospin)
                jp3 = HF_Jplus(HF_states,basis) ! <k|j+|alpha>

                J2crossterm = J2crossterm - vk**2 * (jp1**2 + jp3**2 + jpb**2)

                if (vk**2 * (jp1**2 + jp3**2 + jpb**2) > 1e-6_rk .and. debugmode /= 'low') then
                   write(out_unit,'(a1,1x,i4,1x,2(a,f10.3,1x),i3,a,1x,2(a,f10.6,1x),2(a,f9.3,1x))') &
                        nucleon_type(isospin),k,'e=',e,'v^2=',vk**2, &
                        Omega_k,'/2'//parity_sign(sp_states(k,isospin)%pi), &
                        '<alpha-tilde|j+|k-tilde>**2=',jp2**2,'<k|j+|alpha>**2=',jp3**2,&
                        'J2core_crossterm=',J2core_crossterm,'J2crossterm=',J2crossterm
                end if

             end if

             do l = 1, basis%size

                Omega_l = sp_states(l,isospin)%Omega
                blocked_l = .false.
                do i = 1, nqp
                   if (i_odd(isospin,i) == 0) cycle
                   if (l == i_odd(isospin,i)) blocked_l = .true.
                end do
                if (Omega_l < 0 .or. blocked_l) cycle

                e = sp_states(l,isospin)%energy
!                Eqp_l = sqrt((e-nucleus%pairing(isospin)%chemical_potential)**2 + &
!                     nucleus%pairing(isospin)%average_gap**2)
                Eqp_l = sqrt((e-nucleus%pairing(isospin)%chemical_potential)**2 + &
                     nucleus%pairing(isospin)%gap(l)**2)
                vl = sqrt(sp_states(l,isospin)%occupation)
                ul = sqrt(1.0_rk-sp_states(l,isospin)%occupation)

                HF_states(1) = sp_states(l,isospin)
                HF_states(2) = sp_states(k,isospin)
                jp1 = HF_Jplus(HF_states,basis) ! <l|j+|k>
                sp1 = HF_Splus(HF_states,basis) ! <l|s+|k>

                i = sp_states(l,isospin)%pair_partner(isospin)
                HF_states(1) = sp_states(i,isospin)
                e = sp_states(i,isospin)%energy
!                Eqp_l_tilde = sqrt((e-nucleus%pairing(isospin)%chemical_potential)**2 + &
!                     nucleus%pairing(isospin)%average_gap**2)
                Eqp_l_tilde = sqrt((e-nucleus%pairing(isospin)%chemical_potential)**2 + &
                     nucleus%pairing(isospin)%gap(i)**2)
                i = sp_states(k,isospin)%pair_partner(isospin)
                HF_states(2) = sp_states(i,isospin)
                e = sp_states(i,isospin)%energy
!                Eqp_k_tilde = sqrt((e-nucleus%pairing(isospin)%chemical_potential)**2 + &
!                     nucleus%pairing(isospin)%average_gap**2)
                Eqp_k_tilde = sqrt((e-nucleus%pairing(isospin)%chemical_potential)**2 + &
                     nucleus%pairing(isospin)%gap(i)**2)
                jp2 = HF_Jplus(HF_states,basis) ! <l-tilde|j+|k-tilde>
                sp2 = HF_Splus(HF_states,basis) ! <l-tilde|s+|k-tilde>

                i = sp_states(k,isospin)%pair_partner(isospin)
                HF_states(1) = sp_states(i,isospin)
                i = sp_states(l,isospin)%pair_partner(isospin)
                HF_states(2) = sp_states(i,isospin)
                jp3 = HF_Jplus(HF_states,basis) ! <k-tilde|j+|l-tilde>
                sp3 = HF_Splus(HF_states,basis) ! <k-tilde|s+|l-tilde>

                HF_states(1) = sp_states(k,isospin)
                HF_states(2) = sp_states(l,isospin)
                jp4 = HF_Jplus(HF_states,basis) ! <k|j+|l>
                sp4 = HF_Splus(HF_states,basis) ! <k|s+|l>

                jpt1 = 0.0
                jpt2 = 0.0
                jpt3 = 0.0 ! <k-tilde|j+|l> = 0 (axial symmetry)
                jpt4 = 0.0 ! <l-tilde|j+|k> = 0 (axial symmetry)
                spt1 = 0.0
                spt2 = 0.0
                spt3 = 0.0
                spt4 = 0.0
                if (Omega_k == 1 .and. Omega_l == 1) then

                   HF_states(1) = sp_states(l,isospin)
                   i = sp_states(k,isospin)%pair_partner(isospin)
                   HF_states(2) = sp_states(i,isospin)
                   jpt1 = HF_Jplus(HF_states,basis) ! <l|j+|k-tilde>
                   spt1 = HF_Splus(HF_states,basis) ! <l|s+|k-tilde>

                   HF_states(1) = sp_states(k,isospin)
                   i = sp_states(l,isospin)%pair_partner(isospin)
                   HF_states(2) = sp_states(i,isospin)
                   jpt2 = HF_Jplus(HF_states,basis) ! <k|j+|l-tilde>
                   spt2 = HF_Splus(HF_states,basis) ! <k|s+|l-tilde>

                   i = sp_states(k,isospin)%pair_partner(isospin)
                   HF_states(1) = sp_states(i,isospin)
                   HF_states(2) = sp_states(l,isospin)
                   spt3 = HF_Splus(HF_states,basis) ! <k-tilde|s+|l>

                   i = sp_states(l,isospin)%pair_partner(isospin)
                   HF_states(1) = sp_states(i,isospin)
                   HF_states(2) = sp_states(k,isospin)
                   spt4 = HF_Splus(HF_states,basis) ! <l-tilde|s+|k>

                end if

                tmp = ul**2 * vk**2 * (jp1**2 + jp2**2 + jpt1**2) + &
                     uk*vk*ul*vl * (2 * jp1 * jp3 - jpt1 * jpt2)

                J2core = J2core + tmp

                tmp = (uk*vl*jp4 + ul*vk*jp2)**2/(Eqp_k+Eqp_l_tilde) + &
                     0.25_rk * (uk*vl*jpt2 - ul*vk*jpt1)**2 * &
                     (1.0_rk/(Eqp_k+Eqp_l_tilde) + 1.0_rk/(Eqp_k+Eqp_l_tilde))

                Icr(isospin) = Icr(isospin) + tmp

                jx1 = 0.5_rk * (jp1 + jp4) ! <l|jx|k> = (<l|j+|k> + <k|j+|l>)/2 = <k|jx|l>
                sx1 = 0.5_rk * (sp1 + sp4) ! <l|sx|k> = (<l|s+|k> + <k|s+|l>)/2 = <k|sx|l>

                jx2 = 0.5_rk * (jp2 + jp3) ! <l-tilde|jx|k-tilde> = (<l-tilde|j+|k-tilde> + <k-tilde|j+|l-tilde>)/2
                sx2 = 0.5_rk * (sp2 + sp3) ! <l-tilde|sx|k-tilde> = (<l-tilde|s+|k-tilde> + <k-tilde|s+|l-tilde>)/2

                jx3 = 0.5_rk * jpt2        ! <l-tilde|jx|k> = (<l-tilde|j+|k> + <k|j+|l-tilde>)/2 = <k|jx|l-tilde>
                sx3 = 0.5_rk * (spt4+spt2) ! <l-tilde|sx|k> = (<l-tilde|s+|k> + <k|s+|l-tilde>)/2 = <k|sx|l-tilde>

                jx4 = 0.5_rk * jpt1        ! <k-tilde|jx|l> = (<k-tilde|j+|l> + <l|j+|k-tilde>)/2 = <l|jx|k-tilde>
                sx4 = 0.5_rk * (spt3+spt1) ! <k-tilde|sx|l> = (<k-tilde|s+|l> + <l|s+|k-tilde>)/2 = <l|sx|k-tilde>

                tmp = (uk*vl*sx1+ul*vk*sx2)*(uk*vl*jx1+ul*vk*jx2)/(Eqp_k+Eqp_l_tilde) + &
                     0.5_rk*(uk*vl*sx4-ul*vk*sx3)*(uk*vl*jx4-ul*vk*jx3)*(1.0_rk/(Eqp_k_tilde+Eqp_l_tilde) + &
                     1.0_rk/(Eqp_k+Eqp_l))

                Wcr(isospin) = Wcr(isospin) + 2 * tmp

             end do

          end do

       end do

       write(out_unit,'(/a,f9.3)') '<J^2>core(bulk)      =',J2core
       if (count <= 1) write(out_unit,'(a,f9.3)') '<J^2>core(crossterm) =',J2core_crossterm
       J2core = J2core + J2core_crossterm
       if (count <= 1) then
          write(out_unit,'(a,f9.3)') '<J^2>core(tot)       =',J2core
          write(out_unit,'(a,f9.3)') '<J^2>crossterm       =',J2crossterm
       end if
       j2blocked = 0
       do isospin = 1, 2
          do i = 1, nqp
             if (i_odd(isospin,i)>0) then
                iHF = i_odd(isospin,i)
                j2blocked = j2blocked + sp_states(iHF,isospin)%avg_j * (sp_states(iHF,isospin)%avg_j+1.0_rk)
!                write(out_unit,'(a,i3,1x,i3,a2,a1,2x,a,f9.3,1x,a)') 'blocked state #',i,&
!                     sp_states(iHF,isospin)%Omega,'/2',parity_sign(sp_states(iHF,isospin)%pi), &
!                     '<J^2>=',sp_states(iHF,isospin)%avg_j * (sp_states(iHF,isospin)%avg_j+1.0_rk)
             end if
          end do
       end do
       write(out_unit,'(a,f9.3)') '<J^2>blocked         =',j2blocked
       J2 = J2core + J2crossterm + j2blocked
       if (count <= 1) write(out_unit,'(a,f9.3)') '==> <J^2>            =',J2
          
       write(out_unit,'(/a,/a)') 'Moment of inertia of the core in the cranking approximation with BCS pairing:', &
            '(excluding blocked state and its conjugate)'
       write(out_unit,'(3(a,f10.6,2x))') 'Icr(n)=',Icr(1),'Icr(p)=',Icr(2),'Icr(tot)=',Icr(1)+Icr(2)
       write(out_unit,'(a,f10.3,1x,a)') '==> hbar^2/(2*Icr(tot))=',0.5_rk/(Icr(1)+Icr(2))*1000.0_rk,'keV'
       write(out_unit,'(/a,/a)') 'Moment of inertia in the cranking approximation with BCS pairing and ' &
            // 'Thouless-Valatin correction:','(excluding blocked state and its conjugate)'
       write(out_unit,'(a,f10.3,1x,a)') 'hbar^2/(2*Icr(tot)*1.32)=',0.5_rk/(Icr(1)+Icr(2))*1000.0_rk/1.32_rk,'keV'

       gR = (Icr(2) + (gs(2)-1)*Wcr(2) + gs(1)*Wcr(1))/(Icr(1)+Icr(2))
       write(out_unit,'(/a,/a)') 'Collective gyromagnetic ratio of the core in the cranking approximation ' &
            // 'with BCS pairing:','(excluding blocked state and its conjugate)'
       write(out_unit,'(2(a,f10.6,2x))') 'W(n)=',Wcr(1),'W(p)=',Wcr(2)
       write(out_unit,'(a,f10.3)') 'gR = ',gR

       write(out_unit,*)


       nucleus%moi_n = Icr(1)
       nucleus%moi_p = Icr(2)

       end if

!
!      3) General calculation
!
       if (trim(model%approx) == 'general' .or. &
            trim(model%approx) == 'all') then

       write(out_unit,'(/70a1)') ((/'+'/),i=1,70)
       write(out_unit,'(a)') '3) Expectation value of J^2 in the general case'
       write(out_unit,'(70a1)') ((/'+'/),i=1,70)
       Icr(:)=0.0
       Wcr(:) = 0.0
       J2 = 0.0_rk
       J2core = 0.0_rk
       J2intr = 0.0_rk
       J2polblock = 0.0_rk
       J2blocked = 0.0_rk
       J2z = 0.0_rk
       J2p = 0.0_rk

       do isospin = 1, 2

          do k = 1, basis%size

             Omega_k = sp_states(k,isospin)%Omega
             blocked_k = .false.
             do i = 1, nqp
                if (i_odd(isospin,i) == 0) cycle
                if (k == i_odd(isospin,i) .or. &
                     k == sp_states(i_odd(isospin,i),isospin)%pair_partner(isospin)) &
                     blocked_k = .true.
             end do
             if (Omega_k < 0 .or. blocked_k) cycle

             e = sp_states(k,isospin)%energy
             Eqp_k = sqrt((e-nucleus%pairing(isospin)%chemical_potential)**2 + &
                  nucleus%pairing(isospin)%gap(k)**2)
             vk = sqrt(sp_states(k,isospin)%occupation)
             uk = sqrt(1.0_rk-sp_states(k,isospin)%occupation)

             do i = 1, nqp

                if (i_odd(isospin,i) == 0) cycle

                iHF = sp_states(i_odd(isospin,i),isospin)%pair_partner(isospin)
                HF_states(1) = sp_states(iHF,isospin)
                
                HF_states(2) = sp_states(k,isospin)
                jp1 = HF_Jplus(HF_states,basis) ! <alpha-tilde|j+|k>

                iHF = sp_states(k,isospin)%pair_partner(isospin)
                HF_states(2) = sp_states(iHF,isospin)
                jp2 = HF_Jplus(HF_states,basis) ! <alpha-tilde|j+|k-tilde>

                HF_states(1) = sp_states(k,isospin)
                HF_states(2) = sp_states(i_odd(isospin,i),isospin)
                jp3 = HF_Jplus(HF_states,basis) ! <k|j+|alpha>

                iHF = sp_states(k,isospin)%pair_partner(isospin)
                HF_states(1) = sp_states(iHF,isospin)
                HF_states(2) = sp_states(i_odd(isospin,i),isospin)
                jp4 = HF_Jplus(HF_states,basis) ! <k-tilde|j+|alpha>

                J2polblock = J2polblock + vk**2 * (jp1**2 + jp2**2 - (jp3**2 + jp4**2))

             end do

             do l = 1, basis%size

                Omega_l = sp_states(l,isospin)%Omega
                blocked_l = .false.
                do i = 1, nqp
                   if (i_odd(isospin,i) == 0) cycle
                   if (l == i_odd(isospin,i) .or. &
                        l == sp_states(i_odd(isospin,i),isospin)%pair_partner(isospin)) &
                        blocked_l = .true.
                end do
                if (Omega_l < 0 .or. blocked_l) cycle

                e = sp_states(l,isospin)%energy
                Eqp_l = sqrt((e-nucleus%pairing(isospin)%chemical_potential)**2 + &
                     nucleus%pairing(isospin)%gap(l)**2)
                vl = sqrt(sp_states(l,isospin)%occupation)
                ul = sqrt(1.0_rk-sp_states(l,isospin)%occupation)

                HF_states(1) = sp_states(l,isospin)
                HF_states(2) = sp_states(k,isospin)
                jp1 = HF_Jplus(HF_states,basis) ! <l|j+|k>
                sp1 = HF_Splus(HF_states,basis) ! <l|s+|k>

                i = sp_states(l,isospin)%pair_partner(isospin)
                HF_states(1) = sp_states(i,isospin)
                e = sp_states(i,isospin)%energy
                Eqp_l_tilde = sqrt((e-nucleus%pairing(isospin)%chemical_potential)**2 + &
                     nucleus%pairing(isospin)%gap(i)**2)
                i = sp_states(k,isospin)%pair_partner(isospin)
                HF_states(2) = sp_states(i,isospin)
                e = sp_states(i,isospin)%energy
                Eqp_k_tilde = sqrt((e-nucleus%pairing(isospin)%chemical_potential)**2 + &
                     nucleus%pairing(isospin)%gap(i)**2)
                jp2 = HF_Jplus(HF_states,basis) ! <l-tilde|j+|k-tilde>
                sp2 = HF_Splus(HF_states,basis) ! <l-tilde|s+|k-tilde>

                i = sp_states(k,isospin)%pair_partner(isospin)
                HF_states(1) = sp_states(i,isospin)
                i = sp_states(l,isospin)%pair_partner(isospin)
                HF_states(2) = sp_states(i,isospin)
                jp3 = HF_Jplus(HF_states,basis) ! <k-tilde|j+|l-tilde>
                sp3 = HF_Splus(HF_states,basis) ! <k-tilde|s+|l-tilde>

                HF_states(1) = sp_states(k,isospin)
                HF_states(2) = sp_states(l,isospin)
                jp4 = HF_Jplus(HF_states,basis) ! <k|j+|l>
                sp4 = HF_Splus(HF_states,basis) ! <k|s+|l>

                jpt1 = 0.0
                jpt2 = 0.0
                jpt3 = 0.0 ! <k-tilde|j+|l> = 0 (axial symmetry)
                jpt4 = 0.0 ! <l-tilde|j+|k> = 0 (axial symmetry)
                spt1 = 0.0
                spt2 = 0.0
                spt3 = 0.0
                spt4 = 0.0
                if (Omega_k == 1 .and. Omega_l == 1) then

                   HF_states(1) = sp_states(l,isospin)
                   i = sp_states(k,isospin)%pair_partner(isospin)
                   HF_states(2) = sp_states(i,isospin)
                   jpt1 = HF_Jplus(HF_states,basis) ! <l|j+|k-tilde>
                   spt1 = HF_Splus(HF_states,basis) ! <l|s+|k-tilde>

                   HF_states(1) = sp_states(k,isospin)
                   i = sp_states(l,isospin)%pair_partner(isospin)
                   HF_states(2) = sp_states(i,isospin)
                   jpt2 = HF_Jplus(HF_states,basis) ! <k|j+|l-tilde>
                   spt2 = HF_Splus(HF_states,basis) ! <k|s+|l-tilde>

                   i = sp_states(k,isospin)%pair_partner(isospin)
                   HF_states(1) = sp_states(i,isospin)
                   HF_states(2) = sp_states(l,isospin)
                   spt3 = HF_Splus(HF_states,basis) ! <k-tilde|s+|l>

                   i = sp_states(l,isospin)%pair_partner(isospin)
                   HF_states(1) = sp_states(i,isospin)
                   HF_states(2) = sp_states(k,isospin)
                   spt4 = HF_Splus(HF_states,basis) ! <l-tilde|s+|k>

                end if

                tmp = ul**2 * vk**2 * (jp1**2 + jp2**2 + jpt1**2) + &
                     uk*vk*ul*vl * (2 * jp1 * jp3 - jpt1 * jpt2)

                J2core = J2core + tmp

                tmp = (uk*vl*jp4 + ul*vk*jp2)**2/(Eqp_k+Eqp_l_tilde) + &
                     0.25_rk * (uk*vl*jpt2 - ul*vk*jpt1)**2 * &
                     (1.0_rk/(Eqp_k+Eqp_l_tilde) + 1.0_rk/(Eqp_k+Eqp_l_tilde))

                Icr(isospin) = Icr(isospin) + tmp

                jx1 = 0.5_rk * (jp1 + jp4) ! <l|jx|k> = (<l|j+|k> + <k|j+|l>)/2 = <k|jx|l>
                sx1 = 0.5_rk * (sp1 + sp4) ! <l|sx|k> = (<l|s+|k> + <k|s+|l>)/2 = <k|sx|l>

                jx2 = 0.5_rk * (jp2 + jp3) ! <l-tilde|jx|k-tilde> = (<l-tilde|j+|k-tilde> + <k-tilde|j+|l-tilde>)/2
                sx2 = 0.5_rk * (sp2 + sp3) ! <l-tilde|sx|k-tilde> = (<l-tilde|s+|k-tilde> + <k-tilde|s+|l-tilde>)/2

                jx3 = 0.5_rk * jpt2        ! <l-tilde|jx|k> = (<l-tilde|j+|k> + <k|j+|l-tilde>)/2 = <k|jx|l-tilde>
                sx3 = 0.5_rk * (spt4+spt2) ! <l-tilde|sx|k> = (<l-tilde|s+|k> + <k|s+|l-tilde>)/2 = <k|sx|l-tilde>

                jx4 = 0.5_rk * jpt1        ! <k-tilde|jx|l> = (<k-tilde|j+|l> + <l|j+|k-tilde>)/2 = <l|jx|k-tilde>
                sx4 = 0.5_rk * (spt3+spt1) ! <k-tilde|sx|l> = (<k-tilde|s+|l> + <l|s+|k-tilde>)/2 = <l|sx|k-tilde>

                tmp = (uk*vl*sx1+ul*vk*sx2)*(uk*vl*jx1+ul*vk*jx2)/(Eqp_k+Eqp_l_tilde) + &
                     0.5_rk*(uk*vl*sx4-ul*vk*sx3)*(uk*vl*jx4-ul*vk*jx3)*(1.0_rk/(Eqp_k_tilde+Eqp_l_tilde) + &
                     1.0_rk/(Eqp_k+Eqp_l))

                Wcr(isospin) = Wcr(isospin) + 2 * tmp

             end do

          end do

          do i = 1, nqp
             if (i_odd(isospin,i) == 0) cycle
             J2blocked = J2blocked + sp_states(i_odd(isospin,i),isospin)%avg_j * &
                  (sp_states(i_odd(isospin,i),isospin)%avg_j + 1.0_rk)
             J2z = J2z - (0.5_rk*sp_states(i_odd(isospin,i),isospin)%Omega)**2
             HF_states(2) = sp_states(i_odd(isospin,i),isospin)
             do k = 1, nqp
                if (i_odd(isospin,k) == 0) cycle
                HF_states(1) = sp_states(i_odd(isospin,k),isospin)
                jp1 = HF_Jplus(HF_states,basis)
                J2p = J2p - jp1**2
             end do
          end do
 
       end do

       J2intr = J2blocked + nucleus%Jz(0)**2 + J2z + J2p

       write(out_unit,'(/4(a,f9.3/))') &
            'sum_blocked <a|J^2|a>      =',J2blocked, &
            'K^2                        =',nucleus%Jz(0)**2, &
            '-sum_blocked Omega^2       =',J2z, &
            '-sum_blocked |<a''|J+|a>|^2 =',J2p
       
       write(out_unit,'(a,f9.3)')'<J^2>intr       =',J2intr
       write(out_unit,'(a,f9.3)') '<J^2>core       =',J2core
       write(out_unit,'(a,f9.3)') '<J^2>pol.block. =',J2polblock
       J2 = J2intr + J2core + J2polblock
       write(out_unit,'(a,f9.3)') '==> <J^2>       =',J2
          
       write(out_unit,'(/a,/a)') 'Moment of inertia of the core in the cranking approximation with BCS pairing:', &
            '(excluding blocked state and its conjugate)'
       write(out_unit,'(3(a,f10.6,2x))') 'Icr(n)=',Icr(1),'Icr(p)=',Icr(2),'Icr(tot)=',Icr(1)+Icr(2)
       write(out_unit,'(a,f10.3,1x,a)') '==> hbar^2/(2*Icr(tot))=',0.5_rk/(Icr(1)+Icr(2))*1000.0_rk,'keV'
       write(out_unit,'(/a,/a)') 'Moment of inertia in the cranking approximation with BCS pairing and ' &
            // 'Thouless-Valatin correction:','(excluding blocked state and its conjugate)'
       write(out_unit,'(a,f10.3,1x,a)') 'hbar^2/(2*Icr(tot)*1.32)=',0.5_rk/(Icr(1)+Icr(2))*1000.0_rk/1.32_rk,'keV'

       gR = (Icr(2) + (gs(2)-1)*Wcr(2) + gs(1)*Wcr(1))/(Icr(1)+Icr(2))
       write(out_unit,'(/a,/a)') 'Collective gyromagnetic ratio of the core in the cranking approximation ' &
            // 'with BCS pairing:','(excluding blocked state and its conjugate)'
       write(out_unit,'(2(a,f10.6,2x))') 'W(n)=',Wcr(1),'W(p)=',Wcr(2)
       write(out_unit,'(a,f10.3)') 'gR = ',gR

       write(out_unit,*)

!
!      Store results
!

       expval = J2core
       result2 = Icr(1) + Icr(2)

       end if

    case default

       call warning(modulename, subroutinename, 'Wrong operator name.')

    end select

    return

  end subroutine expectation_value


  !=====================================================================
  !                       FUNCTION DETERMINANT
  !---------------------------------------------------------------------
  ! Calculates the determinant of a real square matrix of size n x n.
  !=====================================================================

  function determinant(matrix, n)

    implicit none
!
!   Arguments
!
    integer,intent(in) :: n
    real(kind=rk),dimension(n,n) :: matrix
    real(kind=rk) :: determinant
!
!   Local variables
!
    real(kind=rk) :: m,tmp
    integer :: i, j, k, l
    logical :: DetExists = .true.


    l = 1

!
!   Conversion to upper triangular form
!

    do k = 1, n-1
       if (matrix(k,k) == 0) then
          DetExists = .false.
          do i = k+1, n
             if (matrix(i,k) /= 0) then
                do j = 1, n
                   tmp = matrix(i,j)
                   matrix(i,j)= matrix(k,j)
                   matrix(k,j) = tmp
                end do
                DetExists = .true.
                l=-l
                EXIT
             endif
          end do
          if (.not.DetExists) then
             determinant = 0
             return
          end if
       endif
       do j = k+1, n
          m = matrix(j,k)/matrix(k,k)
          do i = k+1, n
             matrix(j,i) = matrix(j,i) - m*matrix(k,i)
          end do
       end do
    end do

!    
!   Calculation of determinant by finding product of diagonal elements
!

    determinant = l
    do i = 1, n
       determinant = determinant * matrix(i,i)
    end do

    return
    
  end function determinant


  !=====================================================================
  !                       FUNCTION PERMUTATION_SIGNATURE
  !---------------------------------------------------------------------
  ! Calculates the signature of a permutation of n integer elements.
  !=====================================================================

  function permutation_signature(P,n)

    implicit none

!
!   Arguments
!
    integer :: n
    integer,dimension(n) :: P
    integer :: permutation_signature

!
!   Local variables
!
    integer :: j,k
    integer,dimension(n) :: tmp

    if (minval(P) < 1 .or. maxval(P) > n) then
       permutation_signature = 0
       return
    end if
    tmp(:) = P(:)
    permutation_signature=1
    do j=1,n
       k=tmp(j)
       if (k /= j) then
          permutation_signature = -permutation_signature
          tmp(j)=k
          tmp(k)=j
       end if
    end do

    return

  end function permutation_signature


  !====================================================================!
  !               Subroutine Exact_Coulomb_potential                   !
  !--------------------------------------------------------------------!
  ! Calculates exactly the matrix elements in the cylindrical HO basis !
  ! of the Coulomb contribution to the Hartree--Fock potential, i.e    !
  ! the one-body reduction of the Coulomb interaction for the HF solu- !
  ! tion.                                                              !
  !====================================================================!

  subroutine Exact_Coulomb_potential(term,Coul_op,basis,sp_states,Z)

    implicit none

!
!   Arguments
!
    character(len=3), intent(in) :: term
    type(cylindrical_basis),intent(in) :: basis
    type(operator), dimension(basis%Nblocks), intent(inout) :: Coul_op
    type(sp_state), dimension(basis%size,2), intent(in) :: sp_states
    integer, intent(in) :: Z

!
!   Local variables
!
    character(len=23),parameter :: subroutinename='Exact_Coulomb_potential'
    integer :: sgn, np, i, j, iHF, a, b
    integer :: iBlock, k_block
    integer :: first1, last1, first2, last2, shift
    real(kind=rk) :: C_ak, C_bk, HO_ME, HF_ME
    integer, dimension(4) :: nz, lz, alpha, beta, sigma
    integer, dimension(2) :: ref_index


!
!   Debugging note
!

    if (debugmode=='high') call debugnote(modulename, subroutinename)

    if (term == 'dir') then
       ref_index = (/3,4/)
       sgn=1
    else if (term == 'exc') then
       ref_index = (/4,3/)
       sgn=-1
    else
       STOP "Unrecognized 'term' in Exact_Coulomb_potential."
    end if
    
    call setup_openmp_parameters(1, basis%size, 16, OMP_parameters%set_bounds, &
         OMP_parameters%nthreads)
    call omp_set_num_threads(OMP_parameters%nthreads)
    !$OMP PARALLEL DO &
    !$OMP PRIVATE(iBlock,first1,last1,shift,i,j) &
    !$OMP PRIVATE(nz,lz,alpha,beta,sigma,np) &
    !$OMP PRIVATE(HF_ME,iHF,k_block,first2,last2,a,b,C_ak,C_bk) &
    !$OMP SHARED(basis)
    
    do iBlock = 1, basis%Nblocks

       Coul_op(iBlock)%HO(:,:) = 0.0_rk

       first1 = MERGE(1,SUM(basis%block_size(1:iBlock-1))+1,(iBlock == 1))
       last1 = first1 + basis%block_size(iBlock) - 1
       shift = first1 - 1

       do i = first1, last1

          nz(1) = basis%nz(i)
          lz(1) = basis%Lambda(i)
          np=2*basis%nr(i)+abs(lz(1))
          beta(1) = (np-lz(1))/2
          alpha(1) = beta(1) + lz(1)
          sigma(1) = basis%ns(i)

          do j = i, last1

             nz(ref_index(1)) = basis%nz(j)
             lz(ref_index(1)) = basis%Lambda(j)
             np=2*basis%nr(j)+abs(lz(ref_index(1)))
             beta(ref_index(1)) = (np-lz(ref_index(1)))/2
             alpha(ref_index(1)) = beta(ref_index(1)) + lz(ref_index(1))
             sigma(ref_index(1)) = basis%ns(j)

             if ((term == 'dir').AND.(sigma(1) /= sigma(3))) CYCLE

             HF_ME = 0.0

             do iHF = 1, basis%size

                if (abs(sp_states(iHF,2)%occupation)<1e-06_rk) cycle

                k_block = sp_states(iHF,2)%block
                first2 = MERGE(1,SUM(basis%block_size(1:k_Block-1))+1,(k_Block == 1))
                last2 = first2 + basis%block_size(k_Block) - 1
                   
                do a = first2, last2

                   nz(2) = basis%nz(a)
                   lz(2) = basis%Lambda(a)
                   np=2*basis%nr(a)+abs(lz(2))
                   beta(2) = (np-lz(2))/2
                   alpha(2) = beta(2) + lz(2)
                   sigma(2) = basis%ns(a)

                   if ((term == 'exc').AND.(sigma(2) /= sigma(4))) CYCLE

                   C_ak = sp_states(iHF,2)%coef%cylHO(a)
                   if (abs(C_ak)<1e-06_rk) cycle
                      
                   do b = first2, last2
                         
                      nz(ref_index(2)) = basis%nz(b)
                      lz(ref_index(2)) = basis%Lambda(b)
                      np=2*basis%nr(b)+abs(lz(ref_index(2)))
                      beta(ref_index(2)) = (np-lz(ref_index(2)))/2
                      alpha(ref_index(2))= beta(ref_index(2)) + lz(ref_index(2))
                      sigma(ref_index(2)) = basis%ns(b)

                      if (sigma(Ref_index(2)) /= sigma(Ref_index(2)-2)) CYCLE
                      if (MOD(nz(1)+nz(2),2) /= MOD(nz(3)+nz(4),2)) CYCLE

                      C_bk = sp_states(iHF,2)%coef%cylHO(b)
                      if (abs(C_bk)<1e-06_rk) cycle

                      HO_ME = HO_Coulomb(nz,alpha,beta)
                         
                      HF_ME = HF_ME + HO_ME*C_ak*C_bk*sp_states(iHF,2)%occupation

                   end do! b loop

                end do ! a loop

             end do ! i_HF loop

             Coul_op(iBlock)%HO(i-shift,j-shift) = sgn*HF_ME
             Coul_op(iBlock)%HO(j-shift,i-shift) = Coul_op(iBlock)%HO(i-shift,j-shift)

          end do ! j loop

       end do! i loop

    end do ! iBlock loop

    !$OMP END PARALLEL DO

    return

  end subroutine Exact_Coulomb_potential


  !=======================================================!
  !            subroutine Display_Coulomb                 !
  !-------------------------------------------------------!
  ! Displays the matrix elements of the Coulomb potential !
  ! in the cylindrical HO basis.                          !
  !=======================================================!

  subroutine Display_Coulomb(Coul_dir,Coul_exch,HFfield,basis,Nblocks,unit)

    implicit none

!
!   Arguments
!
    integer, intent(in) :: Nblocks, unit
    type(operator), dimension(Nblocks), intent(in) :: Coul_dir, Coul_exch
    type(field),intent(in) :: HFfield
    type(cylindrical_basis),intent(in) :: basis

!
!   Local variables
!
    character(len=15),parameter :: subroutinename='Display_Coulomb'
    integer :: N, shift, iBlock, i, j
    integer :: ih,il,nzi,nri,li,si,Delta_i,nzj,nrj,lj,sj,Delta_j
    real(kind=rk) :: tmp,xi,eta
    type(spinor) :: chi_i,chi_j,phi_i,phi_j


!
!   Debugging note
!

    if (debugmode=='high') call debugnote(modulename, subroutinename)

    shift = 0

    do iBlock = 1, Nblocks

       N = Size(Coul_dir(iBlock)%HO,dim=1)
       write(unit,'(a,1x,I2,2(1x,a,1x,I3))') 'Block No',iBlock,'=> du',shift+1,'au',shift+N

       do i = 1, N
          nzi=basis%nz(i+shift)
          nri=basis%nr(i+shift)
          li=basis%lambda(i+shift)
          si=basis%ns(i+shift)
          Delta_i=(li-abs(li))/2
          chi_i=spinor(kronecker(si,1),kronecker(si,-1))
          do j = 1, i
             nzj=basis%nz(j+shift)
             nrj=basis%nr(j+shift)
             lj=basis%lambda(j+shift)
             sj=basis%ns(j+shift)
             Delta_j=(lj-abs(lj))/2
             chi_j=spinor(kronecker(sj,1),kronecker(sj,-1))
             tmp=0.0
             if (mod(nzi+nzj,2)==0) then
                do ih=iHmin,basis%NGz/2
                   if (ih==0) cycle
                   xi=basis%xHerm(ih)
                   do il=1,basis%NGr
                      eta=basis%xLag(il)
                      phi_i=sqrt(2.0*bz)*bp*basis%P_Herm(nzi,ih)*sqrt(eta**abs(li))*basis%P_Lag(nri,abs(li),il)* &
                           exp(-0.5_rk*(xi**2+eta))*(-1)**Delta_i*chi_i
                      phi_j=sqrt(2.0*bz)*bp*basis%P_Herm(nzj,ih)*sqrt(eta**abs(lj))*basis%P_Lag(nrj,abs(lj),il)* &
                           exp(-0.5_rk*(xi**2+eta))*(-1)**Delta_j*chi_j
                      tmp=tmp+basis%G_Herm(ih)*basis%G_Lag(il)*real(scalar_product(phi_i,phi_j))* &
                           HFfield%Coulomb(ih,il)
                   end do
                end do
                tmp=tmp*fz/(2.0_rk*bz*bp2)
             end if

             write(unit,'(1x,2(1x,a,1x,I3),4(3x,a,1x,e16.4))') 'iHF',i+shift,'jHF',j+shift,&
                  'Coul_dir =',Coul_dir(iBlock)%HO(i,j),'Coul_exch =',-Coul_exch(iBlock)%HO(i,j), &
                  'int_Coul_dir=',tmp,' diff=',abs(tmp-Coul_dir(iBlock)%HO(i,j))

          end do ! j loop
       end do ! i loop

       shift = shift + N

    end do ! iBlock loop

  end subroutine Display_Coulomb


  !=====================================================================
  !                       SUBROUTINE CHARGE_DENSITY
  !---------------------------------------------------------------------
  ! Calculates the charge density by convolution of proton density and
  ! gaussian proton form factor.
  !=====================================================================

  subroutine charge_density(proton_density,symmetries,basis,nucleus)

    implicit none
    ! Arguments
    type(local_density),intent(in) :: proton_density
    type(symmetry_properties),intent(in) :: symmetries
    type(cylindrical_basis),intent(in) :: basis
    type(nucleus_properties),intent(in) :: nucleus
    !   Local variables
    integer :: NGz, NGr, NGphi, iH, iL, iHprime, iLprime, iphi
    real(kind=rk),parameter :: sigma = 0.85_rk ! in fm
    real(kind=rk) :: z, rho, dz, drho
    real(kind=rk) :: zmin, zmax, rhomin,rhomax
    real(kind=rk) :: Gz, Gr, weight, zprime, rhoprime
    real(kind=rk) :: rho_ch, rc_sqr, rc, rms_p, rms_c
    real(kind=rk) :: int_rho_ch, int_rho_p, int_fp, int_exp_rho_rhoprime_cos
    real(kind=rk),dimension(:),allocatable :: roots, weights
    real(kind=rk),dimension(:,:),allocatable :: tab_rho_ch
    
    rhomin = 0
    rhomax = 10.0_rk ! in fm
    zmax = rhomax
    zmin = -zmax
    if (symmetries%parity) zmin = 0
    dz = 1.0_rk ! in fm
    drho = dz

    NGz = basis%NGz
    NGr = basis%NGr
    NGphi = 20

    allocate(tab_rho_ch(-NGz/2:NGz/2,NGr))
    allocate(roots(NGphi),weights(NGphi))

    call Gauss_Legendre(0.0_rk,2*pi,NGphi,roots,weights)

    !
    ! Check normalization of form factor
    !

    rho_ch = 0.0_rk
    do iH = -NGz/2, NGz/2
       if (iH==0) cycle
       Gz = basis%G_Herm(abs(iH))
       zprime = sign(1.0_rk,iH*1.0_rk) * basis%xHerm(abs(iH))/bz
       do iL = 1, NGr
          weight = Gz * basis%G_Lag(iL)
          rhoprime = sqrt(basis%xLag(iL))/bp
          rho_ch = rho_ch + weight * form_factor(rhoprime,zprime,sigma)
       end do
    end do
    rho_ch = rho_ch * pi/(bz*bp2)
    write(out_unit,'(/a)') &
         'Proton form factor: fp(rho,z) = exp(-(rho^2+z^2)/sigma^2)/' // &
         '(sqrt(pi)*sigma)^3'
    write(out_unit,'(a,f8.3,1x,a,2x,a,f12.9,1x,a)') &
         'sigma =',sigma,'fm','2*pi*integral(fp(rho,z) rho drho dz) =',rho_ch,'e'

    !
    ! Integral of form factor centered at Gauss-Hermite and Gauss-Laguerre points
    !

    write(out_unit,'(/a)') 'Integral of form factor centered at ' // &
         'Gauss-Hermite and Gauss-Laguerre points:'
    
    do iH = iHmin, NGz/2
       if (iH==0) cycle
       z = sign(1.0_rk,iH*1.0_rk) * basis%xHerm(abs(iH))/bz
       do iL = 1, NGr
          rho = sqrt(basis%xLag(iL))/bp
          int_fp = 0
          do iHprime = -NGz/2, NGz/2
             if (iHprime==0) cycle
             zprime = sign(1.0_rk,iHprime*1.0_rk) * basis%xHerm(abs(iHprime))/bz
             Gz = basis%G_Herm(abs(iHprime))
             do iLprime = 1, NGr
                weight = Gz * basis%G_Lag(iLprime)
                rhoprime = sqrt(basis%xLag(iLprime))/bp
                int_exp_rho_rhoprime_cos = 0.0
                do iphi = 1, NGphi
                   int_exp_rho_rhoprime_cos = int_exp_rho_rhoprime_cos + &
                        weights(iphi) * exp(2*rho*rhoprime*cos(roots(iphi))/sigma**2)
                end do
                int_fp = int_fp + weight * &
                     exp(-zprime**2/sigma**2+2*z*zprime/sigma**2) * &
                     exp(-rhoprime**2/sigma**2) * &
                     int_exp_rho_rhoprime_cos
             end do
          end do
          int_fp = form_factor(rho,z,sigma) * int_fp / (bz*2*bp2)
          write(out_unit,'(2(a,f9.3,1x),1x,a,f12.6)') &
               'z =',z,'rho =',rho,'int(fp(r''-r) d^3r'') =',int_fp
       end do
    end do

    !
    ! Tabulate charge density at Gauss-Hermite and Gauss-Laguerre points
    ! to calculate charge radius
    !

    write(out_unit,'(/a)') 'Tabulating charge density at ' // &
         'Gauss-Hermite and Gauss-Laguerre points...'
    
    tab_rho_ch = 0
    do iH = -NGz/2, NGz/2
       if (iH==0) cycle
       z = sign(1.0_rk,iH*1.0_rk) * basis%xHerm(abs(iH))/bz
       do iL = 1, NGr
          rho = sqrt(basis%xLag(iL))/bp
          rho_ch = 0
          do iHprime = -NGz/2, NGz/2
             if (iHprime==0) cycle
             zprime = sign(1.0_rk,iHprime*1.0_rk) * basis%xHerm(abs(iHprime))/bz
             Gz = basis%G_Herm(abs(iHprime))
             do iLprime = 1, NGr
                weight = Gz * basis%G_Lag(iLprime)
                rhoprime = sqrt(basis%xLag(iLprime))/bp
                int_exp_rho_rhoprime_cos = 0.0
                do iphi = 1, NGphi
                   int_exp_rho_rhoprime_cos = int_exp_rho_rhoprime_cos + &
                        weights(iphi) * exp(2*rho*rhoprime*cos(roots(iphi))/sigma**2)
                end do
                rho_ch = rho_ch + weight * &
                     proton_density%rho(abs(iHprime),iLprime) * &
                     exp(-zprime**2/sigma**2+2*z*zprime/sigma**2) * &
                     exp(-rhoprime**2/sigma**2) * &
                     int_exp_rho_rhoprime_cos
             end do
          end do
          rho_ch = form_factor(rho,z,sigma) * rho_ch / (bz*2*bp2)
          tab_rho_ch(iH,iL) = rho_ch
       end do
    end do

    !
    ! Charge radius from root mean square radius of charge density
    !

    rc_sqr = 0
    int_rho_ch = 0
    int_rho_p = 0
    do iHprime = -NGz/2, NGz/2
       if (iHprime==0) cycle
       Gz = basis%G_Herm(abs(iHprime))
       zprime = sign(1.0_rk,iHprime*1.0_rk) * basis%xHerm(abs(iHprime))/bz
       do iLprime = 1, NGr
          weight = Gz * basis%G_Lag(iLprime)
          rhoprime = sqrt(basis%xLag(iLprime))/bp
          rc_sqr = rc_sqr + weight * (zprime**2+rhoprime**2) * &
               tab_rho_ch(iHprime,iLprime)
          int_rho_ch = int_rho_ch + weight * tab_rho_ch(iHprime,iLprime)
          int_rho_p = int_rho_p + weight * proton_density%rho(abs(iHprime),iLprime)
       end do
    end do
    rc_sqr = rc_sqr * pi/(bz*bp2)
    int_rho_ch = int_rho_ch * pi/(bz*bp2)
    int_rho_p = int_rho_p * pi/(bz*bp2)
    rc = sqrt(rc_sqr/int_rho_ch)
    rms_p = sqrt(nucleus%r2(2)/nucleus%Z)
    rms_c = sqrt(nucleus%r2(2)/nucleus%Z + sigma**2)
    
    write(out_unit,'(/a,f9.3,1x,a)') 'Integral(rho_ch) =',int_rho_ch,'e'
    write(out_unit,'(a,f9.3)')       'Integral(rho_p)  =',int_rho_p
    write(out_unit,'(/a)') 'Charge radius:'
    write(out_unit,'(a/,a,f8.3,1x,a)') &
         'Root mean square charge-density radius: rc = sqrt(int(r^2 rho_ch)/int(rho_ch))', &
         '                                           =',rc,'fm'
    ! rc = sqrt(rc_sqr/nucleus%Z)
    ! write(out_unit,'(a/,a,f8.3,1x,a)') &
    !      '                                        rc = sqrt(int(r^2 rho_ch)/Z)', &
    !      '                                           =',rc,'fm'
    write(out_unit,'(a/,a,f8.3,1x,a)') &
         'Root mean square proton-density radius: rp = sqrt(r^2_p/Z)', &
         '                                           =',rms_p,'fm'
    write(out_unit,'(a,f8.3,1x,a)') &
         'Width of proton form factor:         sigma =',sigma,'fm'
    write(out_unit,'(a/,a,f8.3,1x,a)') &
         'Approximate rms charge radius:       rms_c = sqrt(rp^2+sigma^2)', &
         '                                           =',rms_c,'fm'

    deallocate(tab_rho_ch)
    
    return
    
    !
    ! Print rho_ch on regular mesh for plotting
    !

    open(unit=20,file='charge_density.txt',status='replace')

    write(20,'(a)') '#  z (fm)   rho (fm)   charge density (fm^{-3})'
    
    z = zmin
    do while (z <= zmax)
       rho = rhomin
       do while (rho <= rhomax)
          rho_ch = 0.0_rk
          do iH = iHmin, NGz/2
             if (iH==0) cycle
             Gz = basis%G_Herm(iH)
             zprime = basis%xHerm(iH)/bz
             do iL = 1, NGr
                weight = Gz * basis%G_Lag(iL)
                rhoprime = sqrt(basis%xLag(iL))/bp
                rho_ch = rho_ch + weight * proton_density%rho(iH,iL) * &
                     form_factor(rhoprime-rho,zprime-z,sigma)
             end do
          end do
          rho_ch = rho_ch * fz*pi/(bz*bp2)
          write(20,'(2(f8.3,2x),f12.9)') z,rho,rho_ch
          rho = rho + drho
       end do
       write(20,*)
       z = z + dz
    end do

    close(20)    

  contains

    function form_factor(rho,z,sigma) result(f)

      implicit none
      real(kind=rk) :: rho,z,sigma,f

      f = exp(-(rho**2+z**2)/sigma**2)/(sqrt(pi)*sigma)**3

    end function form_factor

  end subroutine charge_density


  !=====================================================================
  !                       SUBROUTINE INTEGRAL_CURRENTS
  !---------------------------------------------------------------------
  ! Calculates the integral over space of the current densities.
  !=====================================================================

  subroutine integral_currents(densities,basis,nucleus)!, &
!       int_j_n,int_j_p,int_rot_s_n,int_rot_s_p)

    implicit none
    ! Arguments
    type(local_density),dimension(2),intent(in) :: densities
!    type(symmetry_properties),intent(in) :: symmetries
    type(cylindrical_basis),intent(in) :: basis
    type(nucleus_properties),intent(in) :: nucleus
    !   Local variables
    integer :: NGz, NGr, NGphi, iH, iL, iHprime, iLprime, iphi
    real(kind=rk) :: z, r
    real(kind=rk) :: Gz, Gr, norm
    real(kind=rk) :: int_rho_n, int_rho_p
    type(vector) :: j_n,j_p,s_n,s_p,rot_s_n,rot_s_p
    type(vector) :: int_j_n,int_j_p,int_s_n,int_s_p,int_rot_s_n,int_rot_s_p    

    
    NGz = basis%NGz
    NGr = basis%NGr

    norm = fz*pi/(bz*bp**2)

    int_j_n = vector(cmplx(0.0_rk),cmplx(0.0_rk),cmplx(0.0_rk))
    int_j_p = int_j_n
    int_s_n = int_j_n
    int_s_p = int_j_n
    int_rot_s_n = int_j_n
    int_rot_s_p = int_j_n
    int_rho_n = 0
    int_rho_p = 0

    do iH = iHmin, NGz/2
       if (iH==0) cycle
       Gz = basis%G_Herm(iH)
       z = sign(1.0_rk,iH*1.0_rk) * basis%xHerm(abs(iH))/bz
       do iL = 1, NGr
          Gr = basis%G_Lag(iL)
          r = sqrt(basis%xLag(iL))/bp
          j_n=vector(densities(1)%j%z(iH,iL),densities(1)%j%rho(iH,iL), &
               densities(1)%j%phi(iH,iL))
          j_p=vector(densities(2)%j%z(iH,iL),densities(2)%j%rho(iH,iL), &
               densities(2)%j%phi(iH,iL))
          s_n=vector(densities(1)%s%z(iH,iL),densities(1)%s%rho(iH,iL), &
               densities(1)%s%phi(iH,iL))
          s_p=vector(densities(2)%s%z(iH,iL),densities(2)%s%rho(iH,iL), &
               densities(2)%s%phi(iH,iL))
          rot_s_n=vector(densities(1)%rot_s%z(iH,iL),densities(1)%rot_s%rho(iH,iL), &
               densities(1)%rot_s%phi(iH,iL))
          rot_s_p=vector(densities(2)%rot_s%z(iH,iL),densities(2)%rot_s%rho(iH,iL), &
               densities(2)%rot_s%phi(iH,iL))
          int_rho_n = int_rho_n + norm*Gz*Gr*densities(1)%rho(iH,iL)
          int_rho_p = int_rho_p + norm*Gz*Gr*densities(2)%rho(iH,iL)
          int_j_n = int_j_n + norm*Gz*Gr*j_n
          int_j_p = int_j_p + norm*Gz*Gr*j_p
          int_s_n = int_s_n + norm*Gz*Gr*0.5_rk*s_n
          int_s_p = int_s_p + norm*Gz*Gr*0.5_rk*s_p
          int_rot_s_n = int_rot_s_n + norm*Gz*Gr*rot_s_n
          int_rot_s_p = int_rot_s_p + norm*Gz*Gr*rot_s_p
       end do
    end do

    write(out_unit,'(/2(a,f8.3,2x))') 'int_rho_n =',int_rho_n,'int_rho_p =',int_rho_p
    write(out_unit,'(2(a,i4,2x))')    '        N =',nucleus%N,'            Z =',nucleus%Z

    write(out_unit,'(/2(a,f10.6,1x))') 'int_s_n%rho =',real(int_s_n%rho),'int_s_p%rho =',real(int_s_p%rho)
    write(out_unit,'(2(a,f10.6,1x))') 'int_s_n%z   =',real(int_s_n%z),'int_s_p%z   =',real(int_s_p%z)
    write(out_unit,'(2(a,f10.6,1x))') '<Sz>_n      =',nucleus%Sz(1), &
         '<Sz>_p      =',nucleus%Sz(2)

    write(out_unit,'(/2(a,f10.6,1x))') 'int_j_n%phi =',real(int_j_n%phi), &
         'int_j_p%phi =',real(int_j_p%phi)
    write(out_unit,'(2(a,f10.6,1x))')  '<Lz>_n      =',nucleus%Jz(1)-nucleus%Sz(1), &
         '<Lz>_p      =',nucleus%Jz(2)-nucleus%Sz(2)

    
  end subroutine integral_currents

  
end module hfbcs
