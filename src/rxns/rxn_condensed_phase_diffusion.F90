! Copyright (C) 2021 Barcelona Supercomputing Center and University of
! Illinois at Urbana-Champaign
! SPDX-License-Identifier: MIT

!> \file
!> The camp_rxn_condensed_phase_diffusion module.

!> \page camp_rxn_condensed_phase_diffusion CAMP: Condensed Phase Diffusion Reaction
!!
!! Condensed phase diffusion reactions are based on Fick's Law of 
!! diffusion and align with kinetic modeling (e.g. \cite Shiraiwa et al. (2010))
!! diffusion representation. 
!!
!!
!! Input data for condensed phase diffusion reactions have the following format :
!! \code{.json}
!!     {
!!    "name" : "condensed phase diffusion",
!!    "type" : "MECHANISM",
!!    "reactions" : [
!!      {
!!        "type" : "DIFFUSION_LAYERS",
!!        "species": [{
!!          "phase": "aqueous",
!!          "name": "H2O_aq"
!!      },
!!      {
!!        "phase": "organic",
!!        "name": "H2O_org"
!!      }]
!!      }
!!    ]}
!! \endcode
!! The key-value pairs \b condensed phase and associated \b species 
!! of the diffusing species are required. The species indicated must 
!! have an associated diffusion coefficient [m2 s-1] listed as a
!! property assocaited with the phase it exists in.
!!
!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> The rxn_condensed_phase_diffusion_t type and associated functions.
module camp_rxn_condensed_phase_diffusion

  use camp_aero_phase_data
  use camp_aero_rep_data
  use camp_chem_spec_data
  use camp_constants,                        only: const
  use camp_camp_state
  use camp_property
  use camp_rxn_data
  use camp_util,                             only: i_kind, dp, to_string, &
                                                  assert, assert_msg, &
                                                  die_msg, string_t

  implicit none
  private

#define NUM_AERO_PHASE_ this%condensed_data_int(1)
#define NUM_AERO_SPECIES_ this%condensed_data_int(2)
#define NUM_INT_PROP_ 2
#define NUM_REAL_PROP_ 0
#define NUM_ENV_PARAM_ 4 
#define DIFF_COEFF_(x) this%condensed_data_int(NUM_AERO_SPECIES_+x)

#define AERO_SPEC_(x) this%condensed_data_int(NUM_INT_PROP_+x)
#define AERO_PHASE_ID_(x) this%condensed_data_int(NUM_INT_PROP_+2*NUM_AERO_PHASE_+x)
#define AERO_REP_ID_(x) this%condensed_data_int(NUM_INT_PROP_+3*NUM_AERO_PHASE_+x)
#define DERIV_ID_(x) this%condensed_data_int(NUM_INT_PROP_+4*NUM_AERO_PHASE_+x)
#define JAC_ID_(x) this%condensed_data_int(NUM_INT_PROP_+1+7*NUM_AERO_PHASE_+x)
#define PHASE_INT_LOC_(x) this%condensed_data_int(NUM_INT_PROP_+2+10*NUM_AERO_PHASE_+x)
#define PHASE_REAL_LOC_(x) this%condensed_data_int(NUM_INT_PROP_+2+11*NUM_AERO_PHASE_+x)
#define NUM_AERO_PHASE_JAC_ELEM_(x) this%condensed_data_int(PHASE_INT_LOC_(x))
#define PHASE_JAC_ID_(x,s,e) this%condensed_data_int(PHASE_INT_LOC_(x)+(s-1)*NUM_AERO_PHASE_JAC_ELEM_(x)+e)
#define EFF_RAD_JAC_ELEM_(x,e) this%condensed_data_real(PHASE_REAL_LOC_(x)-1+e)
#define NUM_CONC_JAC_ELEM_(x,e) this%condensed_data_real(PHASE_REAL_LOC_(x)-1+NUM_AERO_PHASE_JAC_ELEM_(x)+e)
#define MASS_JAC_ELEM_(x,e) this%condensed_data_real(PHASE_REAL_LOC_(x)-1+2*NUM_AERO_PHASE_JAC_ELEM_(x)+e)
#define MW_JAC_ELEM_(x,e) this%condensed_data_real(PHASE_REAL_LOC_(x)-1+3*NUM_AERO_PHASE_JAC_ELEM_(x)+e)

  public :: rxn_condensed_phase_diffusion_t

  !> Generic test reaction data type
  type, extends(rxn_data_t) :: rxn_condensed_phase_diffusion_t
  contains
    !> Reaction initialization
    procedure :: initialize
    !> Finalize the reaction
    final :: finalize, finalize_array
  end type rxn_condensed_phase_diffusion_t

  !> Constructor for rxn_condensed_phase_diffusion_t
  interface rxn_condensed_phase_diffusion_t
    procedure :: constructor
  end interface rxn_condensed_phase_diffusion_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor for Phase transfer reaction
  function constructor() result(new_obj)

    !> A new reaction instance
    type(rxn_condensed_phase_diffusion_t), pointer :: new_obj

    allocate(new_obj)
    new_obj%rxn_phase = AERO_RXN

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize the reaction data, validating component data and loading
  !! any required information into the condensed data arrays for use during
  !! solving
  subroutine initialize(this, chem_spec_data, aero_phase, aero_rep, n_cells)

    !> Reaction data
    class(rxn_condensed_phase_diffusion_t), intent(inout) :: this
    !> Chemical species data
    type(chem_spec_data_t), intent(in) :: chem_spec_data
    !> Aerosol phase data
    type(aero_phase_data_ptr), intent(in) :: aero_phase(:)
    !> Aerosol representations
    type(aero_rep_data_ptr), pointer, intent(in) :: aero_rep(:)
    !> Number of grid cells to solve simultaneously
    integer(kind=i_kind), intent(in) :: n_cells

    type(string_t), allocatable :: diffusion_phase_names(:)
    type(string_t), allocatable :: diffusion_species_names(:)
    !! QQQ: how to i treat this variable?
    type(aero_phase_data_t), pointer :: aero_phase_data

    type(property_t), pointer :: species, spec_props
    character(len=:), allocatable :: key_name, aero_spec_name
    character(len=:), allocatable :: phase_name, species_name, error_msg
    integer(kind=i_kind) :: i_spec, i_aero_rep, n_aero_ids, i_aero_id
    integer(kind=i_kind) :: i_phase, i_species, n_aero_jac_elem, tmp_size
    type(string_t), allocatable :: unique_spec_names(:), unique_act_names(:)
    integer(kind=i_kind), allocatable :: phase_ids(:)
    type(index_pair_t), allocatable :: adjacent_phases(:)
    real(kind=dp) :: temp_real

    ! Get the property set
    if (.not. associated(this%property_set)) call die_msg(300992470, &
            "Missing property set needed to initialize reaction")

    ! Get the species involved in diffusion
    key_name = "species"
    call assert_msg(712818751, &
                    this%property_set%get_property_t(key_name, species), &
                    "Missing species for condensed phase diffusion "// &
                    "reaction")
                    call assert_msg(340551815, species%size() .gt. 0, &
                    "No species specified for the condensed phase "// &
                    "diffusion reaction.")
                    call assert_msg(023007260, species%size() .lt. 3, &
                    "Too many species specified for condensed phase "// &
                    "diffusion reaction (two species maximum).")


    ! Allocate space for phases and species involved in reaction
    allocate(diffusion_phase_names(species%size()))
    allocate(diffusion_species_names(species%size()))
    
    call species%iter_reset()
    do i_species = 1, species%size()

      ! Get the species properties
      call assert_msg(815257799, species%get_property_t(val=species_props), &
              "Invalid structure for species '"// &
              diffusion_species_names(i_layer)%string// &
              "' in condensed phase diffusion reaction.")

      ! Get the phase names
      key_name = "phase"
      call assert_msg(354574496, species_props%get_string(key_name, phase_name), &
              "Missing phase name in condensed phase diffusion reaction.")
      diffusion_phase_names(i_species)%string = phase_name

      ! Get the associated species names
      key_name = "name"
      call assert_msg(629919883, species_props%get_string(key_name, species_name), &
              "Missing species name in condensed phase reaction.")
      diffusion_species_names(i_species)%string = species_name

      ! Load phase dataset
      ! Get the species specific diffusion coefficient
      ! QQQ: this needs work, how do i load the phase property set?
      key_name = "diffusion coefficient [m2 s-1]"
      call assert_msg(690036421, aero_phase_set%val%get_spec_property_set(species_name)% &
              get_real(key, temp_real), "Missing property 'diffusion coefficient [m2 s-1]' &
              for '" species_name//error_msg)

      call species%iter_next()
    end do

    ! Check that the species exist in adjacent layers. 
    ! For the modal/binned aerosol represetnation (no layers) the adjacent_phases array
    ! is alwasys 0. 
    adjacent_phases = aero_rep%adjacent_phases(diffusion_phase_names(1), &
       diffusion_phase_names(SIZE(diffusion_phase_names)))
    call assert_msg(051987857, size(adjacent_phases) .gt. 0, &
       "No adjacent phases found condensed phase diffusion reaction.")

    ! Set up a general error message
    error_msg = " for condensed phase diffusion of aerosol species '"// &
                diffusion_species_names(1)//"' to aerosol species '"// &
                diffusion_species_names(SIZE(diffusion_species_name))

    ! Check for aerosol representations
    call assert_msg(161043212, associated(aero_rep), &
            "Missing aerosol representation"//error_msg)
    call assert_msg(411220610, size(aero_rep).gt.0, &
            "Missing aerosol representation"//error_msg)

    ! Count the instances of this phase/species pair
    n_aero_ids = 0
    n_aero_jac_elem = 0
    do i_aero_rep = 1, size(aero_rep)
      do i_species = 1, species%size()

        ! Get the unique names in this aerosol representation for the
        ! partitioning species
        unique_spec_names = aero_rep(i_aero_rep)%val%unique_names( &
                phase_name = diffusion_phase_names(i_species), &
                spec_name = diffusion_species_names(i_species))

        ! Skip aerosol representations that do not contain this phase
        if (.not.allocated(unique_spec_names)) cycle

        ! Add these instances to the list
        n_aero_ids = n_aero_ids + size(unique_spec_names)

        ! Get the number of Jacobian elements for calculations of mass, volume,
        ! number, etc. for this partitioning into this phase
        ! QQQ: this needs work to save the phase ids for each species
        phase_ids = aero_rep(i_aero_rep)%val%phase_ids(diffusion_phase_names(i_species))
        do i_phase = 1, size(phase_ids)
          n_aero_jac_elem = n_aero_jac_elem + &
                  aero_rep(i_aero_rep)%val%num_jac_elem(phase_ids(i_phase))
        end do

        deallocate(unique_spec_names)
      end do

    end do

    call assert_msg(314000134, n_aero_ids.gt.0, &
                    "Aerosol species not found"//error_msg)

    ! Allocate space in the condensed data arrays
    allocate(this%condensed_data_int(NUM_INT_PROP_ + 2 + n_aero_ids * 13 + &
                                     n_aero_jac_elem * 2))
    allocate(this%condensed_data_real(NUM_REAL_PROP_ + n_aero_jac_elem * 4))
    this%condensed_data_int(:) = int(0, kind=i_kind)
    this%condensed_data_real(:) = real(0.0, kind=dp)

    ! Save space for the environment-dependent parameters
    this%num_env_params = NUM_ENV_PARAM_

    ! Set the number of aerosol-species instances
    NUM_AERO_PHASE_ = n_aero_ids

    ! Set the ids of each aerosol-phase species instance
    i_aero_id = 1
    PHASE_INT_LOC_(i_aero_id)  = NUM_INT_PROP_+12*NUM_AERO_PHASE_+3
    PHASE_REAL_LOC_(i_aero_id) = NUM_REAL_PROP_+1
    do i_aero_rep = 1, size(aero_rep)
      do i_species = 1, species%size()
      ! Get the unique names in this aerosol representation for the
      ! partitioning species
      unique_spec_names = aero_rep(i_aero_rep)%val%unique_names( &
              phase_name = diffusion_phase_names(i_species), &
              spec_name = diffusion_species_names(i_species))

      ! Get the phase ids for this aerosol phase
      phase_ids = aero_rep(i_aero_rep)%val%phase_ids(diffusion_phase_names(i_species))
      ! Add the species concentration and activity coefficient ids to
      ! the condensed data, and set the number of Jacobian elements for
      ! the aerosol representations and the locations of the real data
        do i_spec = 1, size(unique_spec_names)
          NUM_AERO_PHASE_JAC_ELEM_(i_aero_id) = &
                aero_rep(i_aero_rep)%val%num_jac_elem(phase_ids(i_spec))
          AERO_SPEC_(i_aero_id) = &
                aero_rep(i_aero_rep)%val%spec_state_id( &
                unique_spec_names(i_spec)%string)
          AERO_PHASE_ID_(i_aero_id) = phase_ids(i_spec)
          AERO_REP_ID_(i_aero_id) = i_aero_rep
          i_aero_id = i_aero_id + 1
          if (i_aero_id .le. NUM_AERO_PHASE_) then
            PHASE_INT_LOC_(i_aero_id)  = PHASE_INT_LOC_(i_aero_id - 1) + 1 + &
                                       2*NUM_AERO_PHASE_JAC_ELEM_(i_aero_id - 1)
            PHASE_REAL_LOC_(i_aero_id) = PHASE_REAL_LOC_(i_aero_id - 1) + &
                                       4*NUM_AERO_PHASE_JAC_ELEM_(i_aero_id - 1)
          end if
        end do
      end do

      deallocate(unique_spec_names)

    end do

    ! Check the sizes of the data arrays
    tmp_size = PHASE_INT_LOC_(i_aero_id - 1) + 1 + &
               2*NUM_AERO_PHASE_JAC_ELEM_(i_aero_id - 1) - 1
    call assert_msg(625802519, size(this%condensed_data_int) .eq. tmp_size, &
                    "int array size mismatch"//error_msg)
    tmp_size = PHASE_REAL_LOC_(i_aero_id - 1) + &
               4*NUM_AERO_PHASE_JAC_ELEM_(i_aero_id - 1) - 1
    call assert_msg(391089510, size(this%condensed_data_real) .eq. tmp_size, &
                    "real array size mismatch"//error_msg)

  end subroutine initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalize the reaction
  subroutine finalize(this)

    !> Reaction data
    type(rxn_condensed_phase_diffusion_t), intent(inout) :: this

    if (associated(this%property_set)) &
            deallocate(this%property_set)
    if (allocated(this%condensed_data_real)) &
            deallocate(this%condensed_data_real)
    if (allocated(this%condensed_data_int)) &
            deallocate(this%condensed_data_int)

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalize an array of reactions
  subroutine finalize_array(this)
  
    !> Array of reaction data
    type(rxn_condensed_phase_diffusion_t), intent(inout) :: this(:)

    integer(kind=i_kind) :: i

    do i = 1, size(this)
      call finalize(this(i))
    end do

  end subroutine finalize_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module camp_rxn_condensed_phase_diffusion
