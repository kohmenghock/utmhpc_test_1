module pairing_fit

!
! Required modules
!
  use database
  use init
  use io
  use hfbcs
  use tools
  use mb_operators

  implicit none
  character(len=*),parameter :: modulename='pairing_fit'
  type(nucleus_properties)                  :: nucleus              ! Nucleus under study
  type(modelization)                        :: model                ! Model used (Skyrme and pairing)
  type(symmetry_properties)                 :: symmetries           ! Symmetries imposed to the solution
  type(cylindrical_basis)                   :: HO_basis             ! HO expansion basis for sp states
  type(numeric)				    :: numerical_parameters ! Numerical parameters
  type(field),dimension(2)		    :: HFfield              ! One-body potential
  type(local_density),dimension(2)          :: densities            ! Local densities (for n and p)
  type(sp_state),dimension(:,:),allocatable :: sp_states            ! Single-particle states
  type(many_body_basis)                     :: HTDA_basis   ! Many Body Basis
  integer :: Nitermax, HO_basis_size

  private
  public :: fit_pairing_strength, fit_seniority

contains

  !=====================================================================
  !                SUBROUTINE FIT_PAIRING_STRENGTH
  !---------------------------------------------------------------------
  ! Determine the optimal pairing strength through a fit based on 
  ! average pairing gap at the Fermi level.
  ! Restrictions:
  ! * applied only for seniority force
  ! * pairing windows around 6 MeV within Fermi level
  !=====================================================================

  subroutine fit_pairing_strength(nucleus, model, symmetries, HO_basis, &
             numerical_parameters, HFfield, sp_states, densities)

!
!   Arguments
!
    type(nucleus_properties),intent(inout)	:: nucleus              ! Nucleus under study
    type(modelization),intent(inout)		:: model                ! Model used (Skyrme and pairing)
    type(symmetry_properties),intent(in)	:: symmetries           ! Symmetries imposed to the solution
    type(cylindrical_basis)			:: HO_basis             ! HO expansion basis fo sp states
    type(numeric)				:: numerical_parameters ! Numerical parameters
    type(field),dimension(2),intent(inout)	:: HFfield              ! One-body potential
    type(local_density),dimension(2),intent(inout):: densities            ! Local densities (for n and p)
    type(sp_state),dimension(:,:),allocatable,intent(inout):: sp_states            ! Single-particle states
    type(operator),dimension(:),allocatable	:: Coul_exch              ! Exchange Coulomb field

!
!   Local variables
!
    character(len=28),parameter :: subroutinename='fit_pairing_strength'
    integer 			:: i,j,k,l,m,n, isospin,iter,array_size,iHF,jHF
    real(kind=rk)		:: pair_gap, np_ratio
    integer,dimension(2) 	:: Npart
    real(kind=rk),dimension(2) 		:: strength, keep_strength_dddi,max_ave_kHFB_sen,max_ave_kHFB_dddi, pair_gap_keep
    real(kind=rk)		:: pairing_matrix
    character(10)		:: mean_field, keep_model_pairing, blocking_keep
    real(kind=rk)		:: target_Q20, target_Q40, Gn_keep, Gp_keep
    real(kind=rk)		:: E_precision, Q20_precision, Nitermax_keep
    real(kind=rk),dimension(:),allocatable :: density,ener_dens,ave_dens,g_HFB,ave_gHFB,ener_gHFB,k_HFB,ave_kHFB,k_HFB_sen
    real(kind=rk),dimension(:),allocatable :: ave_sp_gap,ave_pair_gap,ave_pair_gap_old,ave_kHFB_old,ave_dens_old
    real(kind=rk), parameter		:: gamma_coeff=1.2,delta_ener=0.01,epsilon=5.0*0.1*delta_ener
    real(kind=rk)			:: gamma,e_fermi,integral,ener_min,ener_max,Epair_dddi,ave_gap
    integer, parameter   	:: order=2
!    integer			:: count
    real(kind=rk)		:: count
    real(kind=rk)		:: test_strength,test_strength_1,test_strength_2,Epair_dddi_1,Epair_dddi_2
    real(kind=rk),dimension(1:2)	:: pme_old, pme_new
    real(kind=rk),dimension(2)	:: Epair_sen, Fermi_energy
    real(kind=rk)		:: ratio
    logical			:: return_fermi
    real(kind=rk)		:: ave_pair, sum_ave_khfb
    real(kind=rk)		:: mu, center, width, fj, fi
    real(kind=rk)		:: Vpair_keep_1, Vpair_keep_2, V_incr, sum_tilde_k2, sum_ui2_vi2, ovlp, tmp
    integer,dimension(2,4)	:: i_odd
    real(kind=rk),dimension(HO_basis%size):: overlap
    integer,dimension(HO_basis%size):: partner
    type(expansion) :: psi_bar

!
!   In the case of DDDI force, first round of calculation is done
!   with seniority force to obtain the estimated constant pairing strength first.
!

    if (model%pairing .eq. 'seniority') then
       write(out_unit,'(a85)') '===================================================================================='
       write(out_unit,'(a85)') '                                                                                    '
       write(out_unit,'(a85)') '                    ENTERING FIT OF BCS PAIRING STRENGTHS                           '
       write(out_unit,'(a85)') '                             SENIORITY FORCE                                        '
       write(out_unit,'(a85)') '                                                                                    '
       write(out_unit,'(a85)') '===================================================================================='
    elseif (model%pairing .eq. 'DDDI') then
       print*, 'Entering fit of pairing strength for DDDI'
       keep_model_pairing = model%pairing
!       keep_strength_dddi(:) =  model%strength(:)
       if (model%eta .eq. 0.0) then 
          keep_strength_dddi(:) = -350.0
       elseif (model%eta .eq. 0.5) then
          keep_strength_dddi(:) = -500.0
       elseif (model%eta .eq. 1.0) then
          keep_strength_dddi(:) = -1000.0
       else
          print*, 'Warning: Using Vq strength in the input file'
          keep_strength_dddi(:) =  model%strength(:)
       end if

       write(out_unit,'(a85)') '===================================================================================='
       write(out_unit,'(a85)') '                                                                                    '
       write(out_unit,'(a85)') '                    ENTERING FIT OF BCS PAIRING STRENGTHS                           '
       write(out_unit,'(a85)') '                     DENSITY DEPENDENT DELTA INTERACTION                            '
       write(out_unit,*)
       write(out_unit,'(a85)') '                     (PERFORM CONSTANT PAIRING FIT FIRST)                           '
       write(out_unit,'(a85)') '                                                                                    '
       write(out_unit,'(a85)') '===================================================================================='

       model%pairing = 'seniority'
       if ((model%truncmin .eq. 6.0) .and. (model%truncmax .eq. 6.0)) then
          model%strength(1) = 19.0
          model%strength(2) = 19.0
       elseif ((model%truncmin .eq. 999.0) .and. (model%truncmax .eq. 6.0)) then
          model%strength(1) = 16.0
          model%strength(2) = 16.0
       else
          model%strength(1) = 16.0
          model%strength(2) = 16.0
          !write(out_unit,*) ' First HFBCS calculations will use the pairing strengths given in input file'
       end if

       mu=model%diffuseness
       center=model%truncmax-model%truncmin
       width=0.5_rk*(model%truncmax+model%truncmin)

    end if

    if (model%fit_pairing .eq. 'jensen') then
       write(out_unit,*) 'Pairing gap type for the fit is JENSEN'
    elseif (model%fit_pairing .eq. 'madland') then
       write(out_unit,*) 'Pairing gap type for the fit is MADLAND'
    elseif (model%fit_pairing .eq. 'moller') then
       write(out_unit,*) 'Pairing gap type for the fit is MOLLER-NIX'
    elseif (model%fit_pairing .eq. 'rmn') then
       write(out_unit,*) 'Pairing gap type for the fit is REDUCED MOLLER-NIX'
    else
       write(out_unit,*) 'Type of pairing gap is incorrect. Please check again.'
       stop
    end if
    write(out_unit,*)


!
!   Some limitations as implemented until May 2024
!

    if ((model%pairing .ne. 'seniority') .and. (model%pairing .ne. 'DDDI'))then
       write(out_unit,*) 'FIT PROCEDURE ONLY LIMITED TO SENIORITY FORCE'
       write(out_unit,*) 'FURTHER DEVELOPMENT MAY BE CONSIDERED'
       write(out_unit,*)
       stop ' ### FIT PROCEDURE ONLY LIMITED TO SENIORITY AND DDDI ONLY ###'
    elseif (model%pairing .eq. 'DDDI') then
       write(out_unit,*) '### FIT PROCEDURE FOR DDDI IS IN PROGRESS ###'
       !stop
    end if

    if ((mod(nucleus%Z,2).ne. 0) .or. (mod(nucleus%N,2) .ne. 0) ) then
!       write(out_unit,*) 'Fitting procedure for odd-mass nucleus not implemented'
!       write(out_unit,*) 'Consider taking the average strength of two neighbouring nuclei'
!       stop 'End fit of pairing strength'
       blocking_keep = model%blocking
       model%blocking = 'no'
    end if


!
!   The pairing fit is limited to the case of parity symmetric shape
!   Fit for exact coulomb is not implemented 
!   Convergence criteria relaxed to speed up calculations (one order lesser)
!
    E_precision = numerical_parameters%convergence%E_precision
    Q20_precision = numerical_parameters%convergence%Q20_precision
    Nitermax_keep = 300	!numerical_parameters%convergence%Nitermax
    iter = 0

    numerical_parameters%convergence%E_precision = E_precision*5.0
    numerical_parameters%convergence%Q20_precision = Q20_precision * 5.0


!
!  Implementing loop calculations
!
   if (symmetries%parity) then
      call cyl_HO_basis(HO_basis, symmetries)
      HO_basis_size = HO_basis%size
      call memory_allocation(nucleus, model, HO_basis, HFfield, sp_states, densities)

      if ((model%exact_coul_exch) .and. (model%fit_pairing .ne. 'no')) then
         write(out_unit,*) 'Pairing fit not implemented for exact Coulomb'
         stop
      else
         call initialize_potential(nucleus, model, HO_basis, HFfield, densities)
      end if
   else
      call fatal(modulename,subroutinename,'Pairing fit only for parity symmetric case')
   end if

   Npart(:)=(/ nucleus%N, nucleus%Z /)
   pme_old(1) = model%strength(1)/(11.0+Npart(1))
   pme_old(2) = model%strength(2)/(11.0+Npart(2))
   
   pme_new(1) = 1.2*pme_old(1)
   pme_new(2) = 1.2*pme_old(2)


   numerical_parameters%convergence%Nitermax = 300

   do while ((abs(pme_new(1) - pme_old(1)) .gt. 0.001) .or. (abs(pme_new(2) - pme_old(2)) .gt. 0.001))
   
      !
      !   First round of HFBCS calculations
      !
      iter = 0
      
      write(out_unit,'(a37)') '=== Iterative HF+BCS calculations ==='
      
      call solve_HFBCS_equations(nucleus, model, symmetries, HO_basis, &
           numerical_parameters, HFfield, sp_states, densities, iter, &
           print_header_message=.true.)

      !
      !   Calculate the pairing gap used for the fit of pairing strength
      !   ala Moller and Nix Nuclear Physics A536 (1992).
      !   The proton pairing gap will be updated 
      ! 

      do isospin = 1, 2

         pme_old(isospin) = pme_new(isospin)

         if (model%fit_pairing .eq. 'jensen') then
            if (isospin .eq. 1) then
               pair_gap = 7.36*(1-(8.15*np_ratio**(2.0)))/(Npart(1)+Npart(2))**(1.0/3.0)
            elseif (isospin .eq. 2) then
               pair_gap = 7.55*(1-(6.07*np_ratio**(2.0)))/(Npart(1)+Npart(2))**(1.0/3.0)
            end if
         elseif (model%fit_pairing .eq. 'sqrtA') then
            if (isospin .eq. 1) then
               pair_gap = 12.5/sqrt(real(Npart(1)+Npart(2)))
            elseif (isospin .eq. 2) then
               pair_gap = 12.5/sqrt(real(Npart(1)+Npart(2)))
            end if
         elseif (model%fit_pairing .eq. 'moller') then
            if (isospin .eq. 1) then
               pair_gap = 4.80/real(Npart(1))**(1.0/3.0)
            elseif (isospin .eq. 2) then
               pair_gap = 4.80/real(Npart(2))**(1.0/3.0)
            end if
         elseif (model%fit_pairing .eq. 'rmn') then
            if (isospin .eq. 1) then
               pair_gap = 4.80/real(Npart(1))**(1.0/3.0)
            elseif (isospin .eq. 2) then
               pair_gap = 4.8 * (0.0181*abs(nucleus%pairing(2)%pairing_energy) + 0.781) &
			/ real(Npart(2))**(1.0/3.0)
            end if
         elseif (model%fit_pairing .eq. 'no') then
            exit
         end if
         pair_gap_keep(isospin) = pair_gap


         !
         !   Fix the energy intervals for the Struntinsky energy averaging method
         !
         if (allocated(ener_dens)) deallocate(ener_dens)
         ener_min = sp_states(1,isospin)%energy - 15.0
         ener_max = sp_states(HO_basis%size,isospin)%energy + 15.0
         array_size = int((ener_max - ener_min)/delta_ener) + 1
         allocate(ener_dens(array_size))


         ! assign ener_dens with values of continuous energy e
         ener_dens(:) = 0.0
         ener_dens(1) = ener_min
         do i = 2,array_size
            ener_dens(i) = ener_dens(i-1) + delta_ener
         end do


         !
         !   Calculate the average level density using the Strutinsky energy averaging method
         !

         gamma = gamma_coeff * 41.0/((Npart(1)+Npart(2))**(1.0/3.0))

         !if (nucleus%pairing(isospin)%average_gap .eq. 0.0) then
         !   !
         !   !previous version (determine level density using discrete levels only)
         !   !
         !   print*, 'Estimation of seniority pairing matrix element using g(HF)'
         !
         !   if (allocated(density)) deallocate(density)
         !   allocate(density(array_size))
         !   call evaluate_level_density(nucleus,sp_states,HO_basis,isospin,ener_dens,array_size,delta_ener,epsilon,density)
         !
         !   if (allocated(ave_dens)) deallocate(ave_dens)
         !   allocate(ave_dens(array_size))
         !   return_fermi = .true.
         !   call ave_density_discrete(nucleus,model,HO_basis,sp_states,isospin,gamma,order,density,&
         !        ener_dens,array_size,delta_ener,ave_dens,epsilon,e_fermi,return_fermi)
         !   write(out_unit,'(a37,i1,a3,f7.3,a4)') 'Fermi energy (seniority) for isopsin ', &
         !               isospin, ' = ', e_fermi, ' MeV'
         !else	!to close for normal use
            print*, 'Estimation of seniority pairing matrix element using g(HFB)'

            if (allocated(density)) deallocate(density)
            allocate(density(array_size))
            call evaluate_g_HFB(nucleus,sp_states,HO_basis,isospin,ener_dens,array_size,delta_ener,epsilon,density)
            
            if (allocated(ave_dens)) deallocate(ave_dens)
            allocate(ave_dens(array_size))
            !return_fermi = .true.
            call ave_density(nucleus,model,HO_basis,sp_states,isospin,gamma,order,density,ener_dens,&
                 array_size,delta_ener,ave_dens,epsilon,e_fermi)

            write(out_unit,'(a37,i1,a3,f7.3,a4)') 'Fermi energy (seniority) for isospin ', &
                           isospin, ' = ', e_fermi, ' MeV'

            !!
            !! Recalculate g_HFB using previous g_HFB
            !!
            !if (allocated(ave_dens_old)) deallocate(ave_dens_old)
            !allocate(ave_dens_old(array_size))
            !ave_dens_old(:) = ave_dens(:)
            !
            !if (allocated(ave_dens)) deallocate(ave_dens)
            !allocate(ave_dens(array_size))
            !
            !return_fermi = .false.
            !call ave_density(nucleus,model,HO_basis,sp_states,isospin,gamma,order,ave_dens_old,ener_dens,&
            !     array_size,delta_ener,ave_dens,epsilon,e_fermi,return_fermi)
            !
            !write(out_unit,'(a50,i1,a3,f7.3,a4)') 'Fermi energy (seniority) recalc g_hfb for isospin ', &
            !               isospin, ' = ', e_fermi, ' MeV'

         !end if


         !
         ! Calculate pairing matrix element by integrating over a pairing window
         !
         integral = 0.0
         do i = 1, array_size
            if (ener_dens(i) .lt. (e_fermi-model%truncmin)) cycle
            if (ener_dens(i) .gt. (e_fermi+model%truncmax)) exit
            integral = integral + ave_dens(i) / sqrt((ener_dens(i) - e_fermi)**2.0 + pair_gap**2.0) * delta_ener
         end do


         pairing_matrix = 2.0 / (integral/2.0)
         strength(isospin) = pairing_matrix*(11.0 + Npart(isospin))

         !
         !   Store the seniority pairing strength and the target Q20
         !
         model%strength(:) = strength(:)
         Gn_keep = strength(1)
         Gp_keep = strength(2)
         target_Q20 = nucleus%Q20(0)
         target_Q40 = nucleus%Q40(0)
         
         pme_new(isospin) = strength(isospin) / (11.0 + Npart(isospin)) 
		 Fermi_energy(isospin) = e_fermi

        !
        ! Set the maximum iteration to 1 (only perform BCS)
        !
        numerical_parameters%convergence%Nitermax = 1
   
      end do   !isospin loop

      write(out_unit,'(a64,f12.6,2x,f12.6,3x,4(f7.4,3x))') &
		'Nuclear Q20, Nuclear Q40, V_n fitted, Vp_fitted >>> ', nucleus%Q20(0), nucleus%Q40(0), &
		pme_new(1), pme_new(2), strength(1), strength(2)
      write(out_unit,*) 'Fermi energy neutron = ',Fermi_energy(1),' and proton = ', Fermi_energy(2)

   end do !while loop to correct the pairing matrix elements   
   

   !
   ! HF+BCS calculations with estimated constant matrix element 
   !
   ! Allow to explore minimum with the estimated pairing strength
   iter = 0
   numerical_parameters%convergence%Nitermax = Nitermax_keep	!1
   
!   model%fit_pairing = 'no'
   call solve_HFBCS_equations(nucleus, model, symmetries, HO_basis, &
           numerical_parameters, HFfield, sp_states, densities, iter, &
           print_header_message=.true.)
             
   call observables(nucleus, model, symmetries, HO_basis, sp_states, &
        densities, HFfield)

   call output(nucleus, model, symmetries, HO_basis, numerical_parameters, &
          HFfield, densities)


   Epair_sen(:) = 0.0
   do isospin = 1,2

         !
         !   Fix the energy intervals for the Struntinsky energy averaging method
         !
         if (allocated(ener_dens)) deallocate(ener_dens)
         ener_min = sp_states(1,isospin)%energy - 15.0
         ener_max = sp_states(HO_basis%size,isospin)%energy + 15.0
         array_size = int((ener_max - ener_min)/delta_ener) + 1
         allocate(ener_dens(array_size))


         ! assign ener_dens with values of continuous energy e
         ener_dens(:) = 0.0
         ener_dens(1) = ener_min
         do i = 2,array_size
            ener_dens(i) = ener_dens(i-1) + delta_ener
         end do
         
         !
         ! Determine Fermi level from g_HFB
         !
         if (allocated(density)) deallocate(density)
         allocate(density(array_size))
         call evaluate_g_HFB(nucleus,sp_states,HO_basis,isospin,ener_dens,array_size,delta_ener,epsilon,density)
         
         do i = 1, array_size
            write(65,*) isospin, ener_dens(i), density(i)
         end do
         write(65,*)

         if (allocated(ave_dens)) deallocate(ave_dens)
         allocate(ave_dens(array_size))
         !return_fermi = .true.
         call ave_density(nucleus,model,HO_basis,sp_states,isospin,gamma,order,density,ener_dens,&
                array_size,delta_ener,ave_dens,epsilon,e_fermi)
         write(out_unit,'(a33,i1,a3,f7.3,a4)') 'Fermi energy (sen) for isopsin ', &
                        isospin, ' = ', e_fermi, ' MeV'

         Fermi_energy(isospin) = e_fermi

         print*,'Seniority Fermi energy for isospin (',isospin,') = ',e_fermi,' MeV'


         !
         ! Determine k_HFB at energy intervals 
         !
         if (allocated(k_HFB)) deallocate(k_HFB)
         allocate(k_HFB(array_size))
         k_HFB(:) = 0.0

         if (allocated(ave_kHFB)) deallocate(ave_kHFB)
	 allocate(ave_kHFB(HO_basis_size))


         write(779,*)  
         write(779,*) '=== Seniority average level density with isospin =', isospin,' ==='
         write(779,*) ' with strength = ', model%strength(isospin)
         ave_kHFB(:) = 0.0
         !return_fermi = .false.
         call ave_density_discrete(nucleus,model,HO_basis,sp_states,isospin,gamma,order,&
              k_HFB,ener_dens,array_size,delta_ener,ave_kHFB,epsilon,Fermi_energy(isospin))


         !
         ! Calculate the condensation energy using average gap and average k_HFB
         !
         Epair_sen(isospin) = 0.0
         do i = 1, HO_basis_size-1
            if (sp_states(i,isospin)%omega<0) cycle
            Epair_sen(isospin) = Epair_sen(isospin)+ ave_kHFB(i) * 0.50_rk
         end do      

         print*,'sum_ave_khfb seniority force for isospin ',isospin,' = ', Epair_sen(isospin)
         Epair_sen(isospin) = Epair_sen(isospin)**2.0 * model%strength(isospin) / (Npart(isospin) + 11.0)
!         Epair_sen(isospin) = Epair_sen(isospin)*0.5_rk

         print*, 'Epair_condensation averaged seniority (',isospin,')= ', Epair_sen(isospin),' MeV'
         print*, 'Epair_condensation from Delta^2/G = ', &
                  nucleus%pairing(isospin)%average_gap**2.0*(Npart(isospin) + 11.0)/model%strength(isospin)
         print*


         Epair_sen(isospin) = nucleus%pairing(isospin)%average_gap**2.0*(Npart(isospin) + 11.0)/model%strength(isospin)


         write(out_unit,'(a39,i1,a4,f8.3,a4)') 'Epair_condensation averaged seniority (',isospin,')= ', &
                 Epair_sen(isospin),' MeV'

   end do
 
	
   if (keep_model_pairing .eq. 'DDDI') then    

      !
      ! HF+BCS calc with delta force residual interaction
      !
   
      model%pairing = trim(keep_model_pairing)
      model%strength(:) = keep_strength_dddi(:)
      numerical_parameters%convergence%Nitermax = Nitermax_keep	!1
      iter = 0

      write(out_unit,*)
      write(out_unit,'(a45,f8.3,2x,f8.3)') &
              'HF+BCS CALC FOR DDDI WITH INITIAl STRENGTHS ',model%strength(1), model%strength(2)
      write(out_unit,*)

      call solve_HFBCS_equations(nucleus, model, symmetries, HO_basis, &
            numerical_parameters, HFfield, sp_states, densities, iter, &
            print_header_message=.true.)
   
      do isospin = 1, 2  

         !
         ! Allocate arrays
         !
         if (allocated(ener_dens)) deallocate(ener_dens)
         ener_min = sp_states(1,isospin)%energy - 15.0
         ener_max = sp_states(HO_basis%size,isospin)%energy + 15.0
         array_size = int((ener_max - ener_min)/delta_ener) + 1
         allocate(ener_dens(array_size))

         ener_dens(:) = 0.0
         ener_dens(1) = ener_min
         do i = 2,array_size
            ener_dens(i) = ener_dens(i-1) + delta_ener
         end do

         if (allocated(density)) deallocate(density)
         allocate(density(array_size))
         density(:) = 0.0

         if (allocated(ave_dens)) deallocate(ave_dens)
         allocate(ave_dens(array_size))
         ave_dens(:) = 0.0

         if (allocated(k_HFB)) deallocate(k_HFB)
         allocate(k_HFB(array_size))
         k_HFB(:) = 0.0

         if (allocated(ave_kHFB)) deallocate(ave_kHFB)
         allocate(ave_kHFB(HO_basis_size))
         ave_kHFB(:) = 0.0

         if (allocated(ave_pair_gap)) deallocate(ave_pair_gap)
         allocate(ave_pair_gap(HO_basis_size))
         ave_pair_gap(:) = 0.0

         !
         ! Loop to estimate delta strength
         !
         Epair_dddi_1 = 0.0
         Epair_dddi_2 = 0.0

         do while ((abs(Epair_dddi_1 - Epair_sen(isospin)) .gt. 0.05) .and. &
             (abs(Epair_dddi_2 - Epair_sen(isospin)) .gt. 0.05))


         !
         ! Determine Fermi level from g_HFB
         !
         call evaluate_g_HFB(nucleus,sp_states,HO_basis,isospin,ener_dens,array_size,delta_ener,epsilon,density)

         call ave_density(nucleus,model,HO_basis,sp_states,isospin,gamma,order,density,ener_dens,&
                array_size,delta_ener,ave_dens,epsilon,e_fermi)
         Fermi_energy(isospin) = e_fermi

         write(out_unit,'(a33,i1,a3,f7.3,a4)') 'Fermi energy (DDDI) for isospin ', &
                    isospin, ' = ', Fermi_energy(isospin), ' MeV'



         !
         ! Determine average k_HFB at discrete eigenenergy
         !
         k_HFB(:) = 0.0
         ave_kHFB(:) = 0.0

         call ave_density_discrete(nucleus,model,HO_basis,sp_states,isospin,gamma,order,&
              k_HFB,ener_dens,array_size,delta_ener,ave_kHFB,epsilon,Fermi_energy(isospin))

         !
         ! Calculate the average pairing gap
         !
         Epair_dddi = 0.0
         ave_pair_gap(:) = 0.0
		
         do i = 1, HO_basis_size
            do j = 1, HO_basis_size
               ave_pair_gap(i) = ave_pair_gap(i)-0.50_rk*ave_kHFB(j)*nucleus%pairing(isospin)%Vpair(i,j)

!               !hockt testing 17 Sept 2025
!               ave_pair_gap(i) = ave_pair_gap(i)-ave_kHFB(j)*nucleus%pairing(isospin)%Vpair(i,j)
            end do
         end do

         !
         ! Calculate the condensation energy using average gap and average k_HFB
         !
         Epair_dddi = 0.0
         do i = 1, HO_basis_size
            Epair_dddi = Epair_dddi + ave_kHFB(i)*ave_pair_gap(i) 
	 end do

		
         Epair_dddi_1 = Epair_dddi_2
         Epair_dddi_2 = Epair_dddi
	
         print*, isospin, model%strength(isospin), &
                 Epair_dddi_1, Epair_dddi_2, Epair_sen(isospin)
 
         write(out_unit,'(a15,2x,i1,2x,4(f8.3,2x))') 'Epair_cond dddi',isospin, model%strength(isospin), &
                 Epair_dddi_1, Epair_dddi_2, Epair_sen(isospin)

	
         !
         ! Adjust pairing strength based on Epair_dddi with Epair_sen
         !
         if ((abs(Epair_dddi_1 - Epair_sen(isospin)) .gt. 0.05) .and. &
                (abs(Epair_dddi_2 - Epair_sen(isospin)) .gt. 0.05)) then
	    if ((Epair_sen(isospin) .gt. Epair_dddi_1).and.(Epair_sen(isospin) .gt. Epair_dddi_2)&
	              .and. (abs(Epair_sen(isospin)-Epair_dddi_2) .gt. 2.0)) then
               model%strength(isospin) = model%strength(isospin) - 10.0 
            elseif ((Epair_sen(isospin) .lt. Epair_dddi_1).and.(Epair_sen(isospin) .lt. Epair_dddi_2)&
                      .and. (abs(Epair_sen(isospin) - Epair_dddi_2) .gt. 2.0)) then
               model%strength(isospin) = model%strength(isospin) + 10.0 
            elseif (abs(Epair_sen(isospin) - Epair_dddi_2) .lt. 2.0) then
               ratio = Epair_sen(isospin) - Epair_dddi_2
               model%strength(isospin) = model%strength(isospin)*(1.0 + ratio*0.02)
            end if
         end if


         !
         !search for sp state with opposite Omega and largest overlap with psi_bar
         !
         partner(:) = 0
         do iHF =1,HO_basis%size
            if (sp_states(iHF,isospin)%Omega<0) cycle
            !definition of time-reversed state
            psi_bar=time_reversal(sp_states(iHF,isospin)%coef,HO_basis)
            ovlp=1e-09_rk
            !partner(iHF)=0

            do jHF=1,HO_basis%size
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

         !
         ! Call pairing matrix element subroutine
         !
         call pairing_matrix_elements(sp_states(:,isospin), densities, isospin, &
              nucleus, model, HO_basis, nucleus%pairing(isospin)%Vpair(:,:), partner)

         !
         ! Call BCS subroutine to compute the pairing gaps (needed to calculate k_HFB)
         !
         ! CAUTION!!!!!
         ! Temporary set no blocked states 
         !
         i_odd(isospin,:) = 0
         numerical_parameters%convergence%Nitermax = 1

         call BCS(nucleus, model, HO_basis, sp_states(:,isospin), nucleus%pairing(isospin)%Vpair(:,:), &
                 isospin, i_odd(isospin,:))

!      numerical_parameters%convergence%Nitermax = 1
!      iter = 0
!      call solve_HFBCS_equations(nucleus, model, symmetries, HO_basis, &
!            numerical_parameters, HFfield, sp_states, densities, iter, &
!            print_header_message=.false.)


      end do 
  
      end do
      write(out_unit,'(a26,2(f12.3,2x,a4))') &
		'Estimated DDDI strength = ', model%strength(1),' MeV', model%strength(2),' MeV'
	  
   end if	! end if statement for DDDI

   call observables(nucleus, model, symmetries, HO_basis, sp_states, &
        densities, HFfield)

   call output(nucleus, model, symmetries, HO_basis, numerical_parameters, &
          HFfield, densities)

   write(out_unit,*)
   write(out_unit,'(a85)') '===================================================================================='
   write(out_unit,*)
   write(out_unit,'(a85)') '                                                                                    '
   write(out_unit,'(a85)') '                       END OF PAIRING STRENGTHS FIT                                 '
   write(out_unit,*)
   write(out_unit,'(a85)') '                                                                                    '
   write(out_unit,'(a85)') '===================================================================================='


   ! Allow to explore minimum with the estimated pairing strength
   model%blocking = blocking_keep
   iter = 0
   numerical_parameters%convergence%Nitermax = Nitermax_keep  
!   model%fit_pairing = 'no'

   call memory_deallocation(nucleus, model, HO_basis, HFfield, sp_states, densities)


   return      !comment the statement if do not want to do second calculations
               !not necessary now since another calculations will be performed upon returning to main.90


  end subroutine fit_pairing_strength



   !=====================================================================
   !                SUBROUTINE FIT_SENIORITY
   !---------------------------------------------------------------------
   ! This subroutine is to determine the SENIORITY FORCE pairing strengths
   ! based on the average pair condensation energy obtained using
   ! DELTA force pairing.
   !=====================================================================

   subroutine fit_seniority(nucleus, model, symmetries, HO_basis, &
             numerical_parameters, HFfield, sp_states, densities)

!
!  Arguments
!
   type(nucleus_properties),intent(inout)       :: nucleus              ! Nucleus under study
   type(modelization),intent(inout)             :: model                ! Model used (Skyrme and pairing)
   type(symmetry_properties),intent(in) :: symmetries           ! Symmetries imposed to the solution
   type(cylindrical_basis)                      :: HO_basis             ! HO expansion basis fo sp states
   type(numeric)                                :: numerical_parameters ! Numerical parameters
   type(field),dimension(2),intent(inout)       :: HFfield              ! One-body potential
   type(local_density),dimension(2),intent(inout):: densities            ! Local densities (for n and p)
   type(sp_state),dimension(:,:),allocatable,intent(inout):: sp_states            ! Single-particle states
   type(operator),dimension(:),allocatable      :: Coul_exch              ! Exchange Coulomb field


!
!   Local variables
!
    character(len=28),parameter :: subroutinename='fit_seniority'
    integer                     :: i,j,k,l,m,n, isospin,iter,array_size,iHF,jHF
    integer,dimension(2)        :: Npart
    real(kind=rk),dimension(2)  :: strength, keep_strength_dddi,max_ave_kHFB_sen,max_ave_kHFB_dddi, pair_gap_keep
    real(kind=rk)               :: pairing_matrix
    character(10)               :: mean_field, keep_model_pairing, blocking_keep
    real(kind=rk)               :: target_Q20, target_Q40, Gn_keep, Gp_keep
    real(kind=rk)               :: E_precision, Q20_precision, Nitermax_keep
    real(kind=rk),dimension(:),allocatable :: density,ener_dens,ave_dens,g_HFB,ave_gHFB,ener_gHFB,k_HFB,ave_kHFB,k_HFB_sen
    real(kind=rk),dimension(:),allocatable :: ave_sp_gap,ave_pair_gap,ave_pair_gap_old,ave_kHFB_old,ave_dens_old
    real(kind=rk), parameter    :: gamma_coeff=1.2,delta_ener=0.01,epsilon=5.0*0.1*delta_ener
    real(kind=rk)               :: gamma,e_fermi,integral,ener_min,ener_max,ave_gap
    integer, parameter          :: order=2
    real(kind=rk)               :: count
    real(kind=rk)               :: test_strength,test_strength_1,test_strength_2,Epair_dddi_1,Epair_dddi_2
    real(kind=rk),dimension(1:2):: pme_old, pme_new
    real(kind=rk),dimension(2)  :: Fermi_energy, Epair_dddi
    real(kind=rk)               :: ratio
    logical                     :: return_fermi
    real(kind=rk)               :: ave_pair, sum_ave_khfb

    real(kind=rk)               :: mu, center, width, fj, fi
    real(kind=rk)               :: Vpair_keep_1, Vpair_keep_2, V_incr, sum_tilde_k2, sum_ui2_vi2
    integer,dimension(2,4)      :: i_odd
    real(kind=rk)               :: Epair_sen_1, Epair_sen_2, seniority_strength, ovlp, tmp
    real(kind=rk),dimension(HO_basis%size):: overlap
    integer,dimension(HO_basis%size):: partner
    type(expansion)             :: psi_bar
    real                        :: sum_dens,diff_kHFB,diff_E,sum_n,diff_n,n_exact,kHFB_exact
    real                        :: E_pairing_exact, E_pairing_ave, E_particle_exact, E_particle_ave, sum_ave_n
    real                        :: E_total_exact, E_total_ave, E_total_exact_sum, E_total_ave_sum
    real(kind=rk),dimension(:),allocatable:: ni, ave_n
    integer                     :: nval



   !
   ! Check if the initial pairing mode is DDDI
   !
   if (model%pairing .ne. 'DDDI') then
      write(out_unit,'(a90)') '========================================================================================='
      write(out_unit,'(a90)') '     Error: fit of seniority force based on average pair condensation energy activated   '
      write(out_unit,'(a90)') '     But initial pairing is not DDDI type                                                '
      write(out_unit,'(a90)') '     Change pairing to DDDI                                                              '
      write(out_unit,'(a90)') '========================================================================================='
      stop
   else
      write(out_unit,'(a90)') '========================================================================================='
      write(out_unit,'(a90)') '     ONE HFBCS CALCULATION USING DDDI STRENGTHS IN INPUT DATA FILE                       '
      write(out_unit,'(a90)') '========================================================================================='
   end if
   

   !
   ! Call cylindrical HO subroutine
   !
   if (symmetries%parity) then
      call cyl_HO_basis(HO_basis, symmetries)
      call memory_allocation(nucleus, model, HO_basis, HFfield, sp_states, densities)

      if ((model%exact_coul_exch) .and. (model%fit_pairing .ne. 'no')) then
         write(out_unit,*) 'This subroutine is not implemented for exact Coulomb'
         stop
      else
         call initialize_potential(nucleus, model, HO_basis, HFfield, densities)
      end if
   else
      call fatal(modulename,subroutinename,'This is limited to parity symmetric case for now')
   end if
   
   !
   ! Perform one HF+BCS iteration using the DDDI strength
   !
   numerical_parameters%convergence%Nitermax = 1
   model%fit_pairing = 'rmn'
   iter = 0
   call solve_HFBCS_equations(nucleus, model, symmetries, HO_basis, &
      numerical_parameters, HFfield, sp_states, densities, iter, &
      print_header_message=.true.)


   !
   ! Estimate the average pair condensation energy from DDDI calculation
   !
   mu=model%diffuseness
   center=model%truncmax-model%truncmin
   width=0.5_rk*(model%truncmax+model%truncmin)

   Epair_dddi(:) = 0.0
   do isospin = 1, 2
   
      !
      ! Allocate arrays
      !
      if (allocated(ener_dens)) deallocate(ener_dens)
      ener_min = sp_states(1,isospin)%energy - 15.0
      ener_max = sp_states(HO_basis%size,isospin)%energy + 15.0
      array_size = int((ener_max - ener_min)/delta_ener) + 1
      allocate(ener_dens(array_size))

      ener_dens(:) = 0.0
      ener_dens(1) = ener_min
      do i = 2,array_size
         ener_dens(i) = ener_dens(i-1) + delta_ener
      end do

      if (allocated(density)) deallocate(density)
      allocate(density(array_size))
      density(:) = 0.0

      if (allocated(ave_dens)) deallocate(ave_dens)
      allocate(ave_dens(array_size))
      ave_dens(:) = 0.0

      if (allocated(k_HFB)) deallocate(k_HFB)
!      allocate(k_HFB(HO_basis%size))
      allocate(k_HFB(array_size))

      k_HFB(:) = 0.0

      if (allocated(ave_kHFB)) deallocate(ave_kHFB)
      allocate(ave_kHFB(HO_basis%size))
!         allocate(ave_kHFB(array_size))

      ave_kHFB(:) = 0.0

      if (allocated(ave_pair_gap)) deallocate(ave_pair_gap)
      allocate(ave_pair_gap(HO_basis%size))
      ave_pair_gap(:) = 0.0


      !
      !   Calculate the average level density using the Strutinsky energy averaging method
      !
      Npart(:)=(/ nucleus%N, nucleus%Z /)
      gamma = gamma_coeff * 41.0/((Npart(1)+Npart(2))**(1.0/3.0))

      !         
      ! Determine Fermi level from g_HFB
      !
      density(:) = 0.0
      ave_dens(:) = 0.0
      call evaluate_g_HFB(nucleus,sp_states,HO_basis,isospin,ener_dens,array_size,delta_ener,epsilon,density)


      !return_fermi = .true.
      call ave_density(nucleus,model,HO_basis,sp_states,isospin,2.0*gamma,order,density,ener_dens,&
             array_size,delta_ener,ave_dens,epsilon,e_fermi)
      write(out_unit,'(a33,i1,a3,f7.3,a4)') 'Fermi energy (DDDI) for isospin ', isospin, ' = ', e_fermi, ' MeV'
      Fermi_energy(isospin) = e_fermi 


      !
      ! Determine average k_HFB at discrete eigenenergy
      !
      k_HFB(:) = 0.0
      ave_kHFB(:) = 0.0

      call ave_density_discrete(nucleus,model,HO_basis,sp_states,isospin,gamma,order,&
          k_HFB,ener_dens,array_size,delta_ener,ave_kHFB,epsilon,Fermi_energy(isospin))

      !
      ! Calculate the average pairing gap
      !
      ave_pair_gap(:) = 0.0
		
      do i = 1, HO_basis%size
         do j = 1, HO_basis%size
            ave_pair_gap(i) = ave_pair_gap(i)-0.50_rk*ave_kHFB(j)*nucleus%pairing(isospin)%Vpair(i,j)
         end do			  
      end do
		
      !
      ! Calculate the condensation energy using average gap and average k_HFB
      !
      do i = 1, HO_basis%size
         if (sp_states(i,isospin)%omega<0) cycle

         Epair_dddi(isospin) = Epair_dddi(isospin) + ave_kHFB(i)*ave_pair_gap(i)
      end do

      write(out_unit,'(a38,i1,a3,f8.3,a4)') 'Epair_condensation (DDDI) for isospin ', &
                isospin,' = ', Epair_dddi(isospin), ' MeV'
      write(out_unit,*)

   end do 	!isospin loop

   !====================================
   ! CALCULATIONS FOR SENIORITY CASE
   !====================================

   !
   ! Perform one HF+BCS iteration using the seniority force
   !
   model%pairing = 'seniority'
   model%strength(:) = -19.0    !starting value for Gn and Gp
   model%constraints(:)%readj = 'n'
   model%constraints(:)%stiffness = 0.0_rk
   numerical_parameters%convergence%Nitermax = 1
!   iter = 0
!   call solve_HFBCS_equations(nucleus, model, symmetries, HO_basis, &
!        numerical_parameters, HFfield, sp_states, densities, iter, &
!        print_header_message=.true.)
!
!   write(out_unit,*) '========================================='
!   write(out_unit,*) '======== COMPLETE HF+BCS ================'
!   write(out_unit,*) '========================================='

   do isospin = 1, 2

      !
      ! Allocate ener_dens array
      !
      if (allocated(ener_dens)) deallocate(ener_dens)
      ener_min = sp_states(1,isospin)%energy - 15.0
      ener_max = sp_states(HO_basis%size,isospin)%energy + 15.0
      array_size = int((ener_max - ener_min)/delta_ener) + 1
      allocate(ener_dens(array_size))

      ener_dens(:) = 0.0
      ener_dens(1) = ener_min
      do i = 2,array_size
         ener_dens(i) = ener_dens(i-1) + delta_ener
      end do

      Npart(:)=(/ nucleus%N, nucleus%Z /)
      gamma = gamma_coeff * 41.0/((Npart(1)+Npart(2))**(1.0/3.0))


      do while (abs(Epair_dddi(isospin)-Epair_sen_2).gt.0.01)

        Epair_sen_2 = (-1.0)*nucleus%pairing(isospin)%average_gap**2.0*(Npart(isospin) + 11.0)/model%strength(isospin)

        ratio = Epair_dddi(isospin)-Epair_sen_2
        model%strength(isospin) = model%strength(isospin)*(1.0+ratio*0.01)

        write(out_unit,*)
        write(out_unit,'(a48,i1)') 'Fitting seniority pairing strength for isospin ', isospin
        write(out_unit,'(i2,3(f8.3,2x))') isospin, Epair_dddi(isospin), Epair_sen_2, model%strength(isospin)
        write(out_unit,*)
        print*, isospin, Epair_dddi(isospin), Epair_sen_2, model%strength(isospin)

        !
        ! Call pairing matrix element subroutine
        !
        call pairing_matrix_elements(sp_states(:,isospin), densities, isospin, &
             nucleus, model, HO_basis, nucleus%pairing(isospin)%Vpair(:,:), partner)

        !
        ! Call BCS subroutine to compute the pairing gaps (needed to calculate k_HFB)
        !
        ! CAUTION!!!!!
        ! Temporary set no blocked states
        !
        i_odd(isospin,:) = 0
        numerical_parameters%convergence%Nitermax = 1

        call BCS(nucleus, model, HO_basis, sp_states(:,isospin), nucleus%pairing(isospin)%Vpair(:,:), &
                isospin, i_odd(isospin,:))

     end do

     write(out_unit,*) 'Epair_sen_2 = ', Epair_sen_2
     write(out_unit,'(a35,i1,a3,f8.3,a4)') 'Seniority strength for isosopin ', &
           isospin,' = ', model%strength(isospin),' MeV'
     write(out_unit,*)

   end do !isospin

   
   !
   ! Perform one HF+BCS iteration with new occupation probabilities
   !
   numerical_parameters%convergence%Nitermax = 1
   model%fit_pairing = 'rmn'
   iter = 0
   write(out_unit,*)
   write(out_unit,*) '==============================================='
   write(out_unit,*) '              STRUTINSKY ENERGY AVE            '
   write(out_unit,*) '==============================================='


   do isospin = 1,2

      !
      ! Determine tilde_n for exact sp energies
      !
      if (allocated(ni)) deallocate(ni)
      allocate(ni(array_size))

      if (allocated(ave_n)) deallocate(ave_n)
      allocate(ave_n(HO_basis%size))

      ni(:) = 0.0_rk
      ave_n(:) = 0.0_rk

      call tilde_n(nucleus,model,HO_basis,sp_states,isospin,gamma,order,ni,ener_dens,nucleus%pairing(isospin)%gap(:),&
                array_size,delta_ener,ave_n,epsilon,Fermi_energy(isospin))

      !
      ! Calculate average energy for exact sp energies
      !
      E_total_exact = 0.0_rk

      do i = 1, HO_basis%size
         if (sp_states(i,isospin)%omega<0) cycle
         E_total_exact = E_total_exact + sp_states(i,isospin)%energy * ave_n(i) &
                           - ave_pair_gap(i) * ave_kHFB(i)
      end do

      write(out_unit,*) 'Strutinsky: Total energy exact sp for isospin', isospin,' = ', E_total_exact

      !
      ! Sum up total average energy
      !
      E_total_exact_sum = E_total_exact_sum + E_total_exact
      

      !
      ! Replace occupation probabilities with ave_n
      !
      do i = 1, HO_basis%size
         sp_states(i,isospin)%occupation = ave_n(i)
      end do

      !
      ! Update one body densities
      !
      call local_one_body_densities(sp_states(:,isospin), HO_basis, densities(isospin), &
            isospin, iter)

   end do    !isospin loop

   !
   ! Update HF potential
   !
   call potential(nucleus, model, symmetries, HO_basis, numerical_parameters, &
        HFfield, densities, model%exact_coul_exch)

   !
   ! Perform HF+BCS calculations with average n
   !
   call solve_HFBCS_equations(nucleus, model, symmetries, HO_basis, &
        numerical_parameters, HFfield, sp_states, densities, iter, &
        print_header_message=.true.)



   do isospin = 1, 2
      !
      ! Determine tilde_n for new (average) sp energies
      !
      if (allocated(ni)) deallocate(ni)
      allocate(ni(array_size))

      if (allocated(ave_n)) deallocate(ave_n)
      allocate(ave_n(HO_basis%size))

      ni(:) = 0.0_rk
      ave_n(:) = 0.0_rk

      call tilde_n(nucleus,model,HO_basis,sp_states,isospin,gamma,order,ni,ener_dens,nucleus%pairing(isospin)%gap(:),&
                array_size,delta_ener,ave_n,epsilon,Fermi_energy(isospin))

      !
      ! Determine tilde_kHFB for new (average) sp energies
      !
      k_HFB(:) = 0.0
      ave_kHFB(:) = 0.0
      return_fermi = .false.
      call ave_density_discrete(nucleus,model,HO_basis,sp_states,isospin,gamma,order,&
           k_HFB,ener_dens,array_size,delta_ener,ave_kHFB,epsilon,Fermi_energy(isospin))

      !
      ! Calculate average energy for new (average) sp energies
      !
      E_total_ave = 0.0_rk

      do i = 1, HO_basis%size
         if (sp_states(i,isospin)%omega<0) cycle
         E_total_ave = E_total_ave + sp_states(i,isospin)%energy * ave_n(i) &
                           - ave_pair_gap(i) * ave_kHFB(i)
      end do

      write(out_unit,*) 'Strutinsky: Total energy ave sp for isospin', isospin,' = ', E_total_ave

      !
      ! Sum up total average energy
      !
      E_total_ave_sum = E_total_ave_sum + E_total_ave

   end do !isospin dddi

   write(out_unit,*) 'Strutinsky: Total exact energy = ', E_total_exact_sum
   write(out_unit,*) 'Strutinsky: Total average energy = ', E_total_ave_sum

   write(out_unit,*)
   write(out_unit,*) 'Strutinsky: Average Epair_cond dddi = ', Epair_dddi(1), Epair_dddi(2)

!   !
!   ! Redo one HF+BCS to restore back original sp solutions
!   !
!   numerical_parameters%convergence%Nitermax = 1
!   model%fit_pairing = 'rmn'
!   iter = 0
!   write(out_unit,*)
!   write(out_unit,*) '==============================================='
!   write(out_unit,*) '            REDO ONE HF+BCS ITERATION          '
!   write(out_unit,*) '==============================================='
!   call solve_HFBCS_equations(nucleus, model, symmetries, HO_basis, &
!       numerical_parameters, HFfield, sp_states, densities, iter, &
!       print_header_message=.true.)


   
   numerical_parameters%convergence%Nitermax = Nitermax_keep
   call memory_deallocation(nucleus, model, HO_basis, HFfield, sp_states, densities)

   stop 'Fit of seniority strength completed (without HF+BCS calculation)'


   end subroutine fit_seniority





!=======================================================================================
!
!              Subroutine to determine average level density
!
!  This subroutine should be (can be) used to evaluate the average quantitites 
!  e.g. average level density (first used for the constant pairing strength)
!  or the density of HFB 
!
!  This subroutine takes in discrete quantities in the form of
!  >>> dens and ener_dens with an array size called array
!
!  to obtain an average density >>> ave_dens.
!  
!  Values of gamma, order, delta(pair_gap) are transferred from main part of the code
!  to ensure no inconsistency if any changes is made to the code.
!  
!  DENS IS REPRESENTING LEVEL DENSITY WITH CONTINUOUS ENERGIES
!  ENER_DENS IS NOW REPRESENTING THE CONTINUOUS ENERGY (NOT THE DISCRETE SP ENERGIES)
!
!=======================================================================================
   subroutine ave_density_discrete(nucleus,model,HO_basis,sp_states,isospin,gamma,order,dens,ener_dens,&
		array,delta_ener,ave_dens,epsilon,e_fermi)

   type(nucleus_properties)                  :: nucleus              ! Nucleus under study
   type(modelization)                        :: model                ! Model used (Skyrme and pairing)
   type(sp_state),dimension(:,:),allocatable :: sp_states            ! Single-particle states
   type(cylindrical_basis)                   :: HO_basis             ! HO expansion basis for sp states

   integer			::	i, j, k, array, isospin
   real(kind=rk)		::	ener_min, ener_max, e_test, x, weight, &
					 sum_particle,gamma,delta_ener,epsilon
   real(kind=rk)		::	e_fermi,poly,ratio,dens_fermi
   real(kind=rk), parameter	::	pi = 3.1415926536	!, epsilon=5.0*0.1*delta_ener
   integer			::	particle_num, mass_num, order
   real(kind=rk), dimension(array)	::	ener_dens
   real(kind=rk), dimension(array)	::     dens
   real(kind=rk), dimension(HO_basis%size)::	ave_dens
   real(kind=rk), dimension(0:4,0:8)	::	a2n
   logical			::	return_fermi
   real(kind=rk)		::	density, ener_diff


   !
   ! Allocate values for coefficients in Laguerre polynomial according to order M
   !subroutine ave_
   a2n(:,:) = 0.0
   if (order .eq. 0) then
      a2n(0,0) = 1.0
   elseif (order .eq. 1) then
      a2n(1,0) = 3.0/2.0
      a2n(1,2) = -1.0
   elseif (order .eq. 2) then
      a2n(2,0) = 15.0/8.0
      a2n(2,2) = -5.0/2.0
      a2n(2,4) = 1.0/2.0
   elseif (order .eq. 3) then
      a2n(3,0) = 35.0/16.0
      a2n(3,2) = -35.0/8.0
      a2n(3,4) = 7.0/4.0
      a2n(3,6) = -1.0/6.0
   elseif (order .eq. 4) then
      a2n(4,0) = 945.0/384.0
      a2n(4,2) = -315.0/48.0
      a2n(4,4) = 63.0/16.0
      a2n(4,6) = -9.0/12.0
      a2n(4,8) = 1.0/24.0
   end if

   !
   ! calculate the average density
   ! ener(:) = energy of the single-particle file
   ! ener_dens(:) = energy starting from ener(1) - 2*hbaromega
   ! save results in ave_dens(:)
   !
   ave_dens(:) = 0.0
   do i = 1, HO_basis%size
      poly = 0.0

      dens(:) = 0.0
      do j = 1, array

         x = (sp_states(i,isospin)%energy - ener_dens(j)) / gamma
         weight = exp(-(x)**2.0)/sqrt(pi)
         do k = 0, 2*order, 2
            poly = poly + a2n(order,k)*(x**real(k))
         end do

         dens(j) = dens(j) + 0.50_rk * nucleus%pairing(isospin)%gap(i)/ &
             sqrt((ener_dens(j) - e_fermi)**2.0 + (nucleus%pairing(isospin)%gap(i))**2.0)

         
!         ave_dens(i) = ave_dens(i) + 0.5*(1.0/gamma)*(dens(j)+dens(j+1))*poly*weight*&
!                        abs((sp_states(j,isospin)%energy-sp_states(j+1,isospin)%energy))

         ave_dens(i) = ave_dens(i) + (1.0/gamma)*dens(j)*poly*weight*delta_ener

         poly = 0.0
      end do
      !if (isospin .eq. 1) write(555,'(i1,2x,3(f12.6,2x))') isospin, sp_states(i,isospin)%energy, dens(i),ave_dens(i)
   end do

   
   end subroutine



!=======================================================================================
!
!              Subroutine to determine average level density
!
!  This subroutine should be (can be) used to evaluate the average quantitites 
!  e.g. average level density (first used for the constant pairing strength)
!  or the density of HFB 
!
!  This subroutine takes in discrete quantities in the form of
!  >>> dens and ener_dens with an array size called array
!
!  to obtain an average density >>> ave_dens.
!  
!  Values of gamma, order, delta(pair_gap) are transferred from main part of the code
!  to ensure no inconsistency if any changes is made to the code.
!  
!  DENS IS REPRESENTING LEVEL DENSITY WITH CONTINUOUS ENERGIES
!  ENER_DENS IS NOW REPRESENTING THE CONTINUOUS ENERGY (NOT THE DISCRETE SP ENERGIES)
!
!=======================================================================================
   subroutine ave_density(nucleus,model,HO_basis,sp_states,isospin,gamma,order,dens,ener_dens,&
		array,delta_ener,ave_dens,epsilon,e_fermi)

   type(nucleus_properties)                  :: nucleus              ! Nucleus under study
   type(modelization)                        :: model                ! Model used (Skyrme and pairing)
   type(sp_state),dimension(:,:),allocatable :: sp_states            ! Single-particle states
   type(cylindrical_basis)                   :: HO_basis             ! HO expansion basis for sp states

   integer			::	i, j, k, array, isospin
   real(kind=rk)		::	ener_min, ener_max, e_test, x, weight, &
					 sum_particle,gamma,delta_ener,epsilon
   real(kind=rk)		::	e_fermi,poly,ratio,dens_fermi
   real(kind=rk), parameter	::	pi = 3.1415926536	!, epsilon=5.0*0.1*delta_ener
   integer			::	particle_num, mass_num, order
   real(kind=rk), dimension(array)	::	ener_dens, ave_dens
   real(kind=rk), dimension(array)	::     dens
!   real, dimension(HO_basis%size)::	dens
   real(kind=rk), dimension(0:4,0:8)	::	a2n
   logical			::	return_fermi
   real(kind=rk)		::	density


   !
   ! Allocate values for coefficients in Laguerre polynomial according to order M
   !
   a2n(:,:) = 0.0
   if (order .eq. 0) then
      a2n(0,0) = 1.0
   elseif (order .eq. 1) then
      a2n(1,0) = 3.0/2.0
      a2n(1,2) = -1.0
   elseif (order .eq. 2) then
      a2n(2,0) = 15.0/8.0
      a2n(2,2) = -5.0/2.0
      a2n(2,4) = 1.0/2.0
   elseif (order .eq. 3) then
      a2n(3,0) = 35.0/16.0
      a2n(3,2) = -35.0/8.0
      a2n(3,4) = 7.0/4.0
      a2n(3,6) = -1.0/6.0
   elseif (order .eq. 4) then
      a2n(4,0) = 945.0/384.0
      a2n(4,2) = -315.0/48.0
      a2n(4,4) = 63.0/16.0
      a2n(4,6) = -9.0/12.0
      a2n(4,8) = 1.0/24.0
   end if

   !
   ! calculate the average density
   ! ener(:) = energy of the single-particle file
   ! ener_dens(:) = energy starting from ener(1) - 2*hbaromega
   ! save results in ave_dens(:)
   !
   ave_dens(:) = 0.0_rk
   do i = 1, array
      poly = 0.0_rk
      do j = 1, array	!HO_basis%size
         x = (ener_dens(i) - ener_dens(j)) / gamma
         weight = exp(-(x)**2.0_rk)/sqrt(pi)
         do k = 0, 2*order, 2
            poly = poly + a2n(order,k)*(x**real(k))
         end do

         ave_dens(i) = ave_dens(i) + (1.0_rk/gamma)*dens(j)*poly*weight*delta_ener
         poly = 0.0_rk
      end do
      !if (isospin .eq. 1) write(779,*) isospin, ener_dens(i), ave_dens(i)
   end do

!   if (return_fermi) then
      !
      ! Calculate Fermi levels from average density
      !
      
      if (isospin .eq. 1) then
         particle_num = nucleus%N
      elseif (isospin .eq. 2) then
         particle_num = nucleus%Z
      end if
      mass_num = nucleus%Z + nucleus%N

      sum_particle = 0.0
      e_fermi = 0.0
      do i = 1, array
         sum_particle = sum_particle + ave_dens(i)*delta_ener
         if (isospin .eq. 1) write(779,*) isospin, ener_dens(i), ave_dens(i), sum_particle

         if (((sum_particle - real(particle_num)) .gt. epsilon) .and. (i .gt. 1)) then
            e_fermi = 0.5*(ener_dens(i-1)+ener_dens(i))
            !write(556,*) 'Fermi energy for isospin(',isospin,') = ', e_fermi,' MeV'
            exit
         end if

      end do
 

      !
      ! Determine average density at the Fermi level
      !
      dens_fermi = 0.0
      do i =1, array-1
         if ((e_fermi .gt. ener_dens(i)) .and. (e_fermi .lt. ener_dens(i+1)) ) then
            ratio = (e_fermi - ener_dens(i)) / (ener_dens(i+1) - ener_dens(i))
            dens_fermi = ratio*(ave_dens(i+1) - ave_dens(i)) + ave_dens(i)
         end if
      end do
      write(out_unit,'(a33,2x,f7.3)') 'Average density at Fermi level = ', dens_fermi

!   end if
   
   end subroutine




!=======================================================================================
!
!   TO COMPUTE THE DISCRETE LEVEL DENSITY
!
!=======================================================================================
   subroutine evaluate_level_density(nucleus,sp_states,HO_basis,isospin,ener,array,delta_ener,epsilon,dens)

    implicit none
    type(nucleus_properties)                    :: nucleus              ! Nucleus under study
    type(sp_state),dimension(:,:), allocatable  :: sp_states            ! Single-particle states
    type(cylindrical_basis)                     :: HO_basis             ! HO expansion basis fo sp states


    integer                             :: i,j,isospin,array
    real(kind=rk)			:: delta_ener, epsilon
    !real, parameter                     :: epsilon=5.0*0.1*delta_ener
    real(kind=rk),dimension(:),allocatable	:: dens,ener
 


   !
   ! sum discrete level density for at discrete sp energy
   !
   !allocate(dens(array))


   dens(:) = 0.0
   do j = 1, array	!HO_basis%size
      do i = 1, HO_basis%size
         if (abs(ener(j) - sp_states(i,isospin)%energy) .le. epsilon) then
            dens(j) = dens(j) + 1.0
         end if
      end do
   end do

   end subroutine


!=======================================================================================
!
!   TO COMPUTE THE AVERAGE LEVEL DENSITY g_HFB AS IN TABLE 1 OF 
!   M. BRACK & P. QUENTIN NPA 361 (1981)
!
!=======================================================================================
   subroutine evaluate_g_HFB(nucleus,sp_states,HO_basis,isospin,ener,array,delta_ener,epsilon,g_HFB)

    implicit none
    type(nucleus_properties)                    :: nucleus              ! Nucleus under study
    type(sp_state),dimension(:,:), allocatable  :: sp_states            ! Single-particle states
    type(cylindrical_basis)                     :: HO_basis             ! HO expansion basis fo sp states


    integer				:: i,j,isospin,array,iHF,jHF
    real(kind=rk)			:: ener_min,ener_max,sum_g,delta_ener,epsilon
    !real, parameter                     :: delta_ener=0.001,epsilon=5.0*0.1*delta_ener
    real(kind=rk),dimension(:)		:: g_HFB,ener


   g_HFB(:) = 0.0
   do i = 1, array		!HO_basis%size
      sum_g = 0.0
      do iHF = 1, HO_basis%size

         sum_g = sum_g + ((nucleus%pairing(isospin)%gap(iHF))**2.0_rk / &
                 (2.0*((ener(i) - sp_states(iHF,isospin)%energy)**2.0_rk &
                  + (nucleus%pairing(isospin)%gap(iHF))**2.0_rk)**1.5_rk))

      end do

      g_HFB(i) = sum_g !* 2.0


   end do
   
   end subroutine




!=======================================================================================
!
!   TO COMPUTE THE DISCRETE ABNORMAL DENSITY k_HFB AS IN TABLE 1 OF
!   M. BRACK & P. QUENTIN NPA 361 (1981)
!
!=======================================================================================
   subroutine evaluate_k_HFB(nucleus,model,sp_states,HO_basis,isospin,e_fermi,ener,array,delta_ener,epsilon,k_HFB)

    implicit none
    type(nucleus_properties)                    :: nucleus              ! Nucleus under study
    type(sp_state),dimension(:,:), allocatable  :: sp_states            ! Single-particle states
    type(cylindrical_basis)                     :: HO_basis             ! HO expansion basis fo sp states
    type(modelization),intent(inout)            :: model                ! Model used (Skyrme and pairing)

    integer                             :: i,j,isospin,array,iHF,jHF,count
    real(kind=rk)                       :: ener_min,ener_max,sum_k,e_fermi,epsilon,delta_ener,sp_gap
    !real, parameter                    :: delta_ener=0.001,epsilon=5.0*0.1*delta_ener
    !real,dimension(array)				:: k_HFB,ener
    real(kind=rk),dimension(array)	:: ener
    real(kind=rk),dimension(HO_basis%size):: k_HFB
    real(kind=rk)			:: mu, center, width, fi



    mu=model%diffuseness
    center=model%truncmax-model%truncmin
    width=0.5_rk*(model%truncmax+model%truncmin)


    k_HFB(:) = 0.0
    do i = 1, HO_basis%size
       do j = 1, array

          k_HFB(i) = k_HFB(i) + 0.50_rk * nucleus%pairing(isospin)%gap(i)/ &
             sqrt((ener(j) - e_fermi)**2.0 + (nucleus%pairing(isospin)%gap(i))**2.0)

       end do
       write(55,*) isospin, k_HFB(i)
    end do
    write(55,*)

   end subroutine

  

!=============================================================================================
!  COMPUTE AVERAGE PARTICLE NUMBER
!  SIMILAR STRUCTURE TO AVERAGE K_HFB
!=============================================================================================
   subroutine tilde_n(nucleus,model,HO_basis,sp_states,isospin,gamma,order,n,ener_dens,ave_pair_gap, &
                array,delta_ener,ave_n,epsilon,e_fermi)

   type(nucleus_properties)                  :: nucleus              ! Nucleus under study
   type(modelization)                        :: model                ! Model used (Skyrme and pairing)
   type(sp_state),dimension(:,:),allocatable :: sp_states            ! Single-particle states
   type(cylindrical_basis)                   :: HO_basis             ! HO expansion basis for sp states

   integer                      ::      i, j, k, array, isospin
   real(kind=rk)                ::      ener_min, ener_max, e_test, x, weight, &
                                         sum_particle,gamma,delta_ener,epsilon
   real(kind=rk)                ::      e_fermi,poly,ratio,dens_fermi
   real(kind=rk), parameter     ::      pi = 3.1415926536       !, epsilon=5.0*0.1*delta_ener
   real(kind=rk), dimension(0:4,0:8)    ::      a2n
   integer                              ::      particle_num, mass_num, order
   real(kind=rk), dimension(array)      ::     ener_dens
   real(kind=rk), dimension(array)      ::     n
   real(kind=rk), dimension(HO_basis%size)      ::     ave_n, ave_pair_gap

   !
   ! Allocate values for coefficients in Laguerre polynomial according to order M
   !
   a2n(:,:) = 0.0
   if (order .eq. 0) then
      a2n(0,0) = 1.0
   elseif (order .eq. 1) then
      a2n(1,0) = 3.0/2.0
      a2n(1,2) = -1.0
   elseif (order .eq. 2) then
      a2n(2,0) = 15.0/8.0
      a2n(2,2) = -5.0/2.0
      a2n(2,4) = 1.0/2.0
   elseif (order .eq. 3) then
      a2n(3,0) = 35.0/16.0
      a2n(3,2) = -35.0/8.0
      a2n(3,4) = 7.0/4.0
      a2n(3,6) = -1.0/6.0
   elseif (order .eq. 4) then
      a2n(4,0) = 945.0/384.0
      a2n(4,2) = -315.0/48.0
      a2n(4,4) = 63.0/16.0
      a2n(4,6) = -9.0/12.0
      a2n(4,8) = 1.0/24.0
   end if


   !
   ! calculate the average particle number
   ! ener_dens(:) = energy starting from ener(1) - 2*hbaromega
   ! save results in ave_n(:)
   !
   ave_n(:) = 0.0
   do i = 1, HO_basis%size
      poly = 0.0

      n(:) = 0.0
      do j = 1, array

         x = (sp_states(i,isospin)%energy - ener_dens(j)) / gamma
         weight = exp(-(x)**2.0)/sqrt(pi)
         do k = 0, 2*order, 2
            poly = poly + a2n(order,k)*(x**real(k))
         end do

!         n(j) =  0.50_rk * (1.0_rk + (e_fermi-ener_dens(j)) / &
!                  sqrt((ener_dens(j) - e_fermi)**2.0 + (nucleus%pairing(isospin)%gap(i))**2.0))

         n(j) =  0.50_rk * (1.0_rk + (e_fermi-ener_dens(j)) / &
                  sqrt((ener_dens(j) - e_fermi)**2.0 + (ave_pair_gap(i))**2.0))


         ave_n(i) = ave_n(i) + (1.0/gamma)*n(j)*poly*weight*delta_ener

         poly = 0.0
      end do

   end do

   end subroutine



end module pairing_fit
