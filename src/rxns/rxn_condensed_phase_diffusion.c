/* Copyright (C) 2021 Barcelona Supercomputing Center and University of
 * Illinois at Urbana-Champaign
 * SPDX-License-Identifier: MIT
 *
 * Condensed phase diffusion reaction solver functions
 *
 */
/** \file
 * \brief Condensed phase diffusion reaction solver functions
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <camp/aero_rep_solver.h>
#include <camp/aero_phase_solver.h>
#include <camp/rxns.h>
#include <camp/sub_model_solver.h>
#include <camp/util.h>

// TODO Lookup environmental indices during initialization
//#define TEMPERATURE_K_ env_data[0]
//#define PRESSURE_PA_ env_data[1]

#define NUM_ADJACENT_PAIRS_ int_data[0]

#define NUM_INT_PROP_ 1
#define NUM_FLOAT_PROP_ 0
#define NUM_ENV_PARAM_ 0
#define JAC_INNER_ 0
#define JAC_OUTER_ 1

#define DIFF_COEFF_INNER_(x) (float_data[(NUM_FLOAT_PROP_) + (x)])
#define DIFF_COEFF_OUTER_(x) (float_data[(NUM_FLOAT_PROP_) + (NUM_ADJACENT_PAIRS_) + (x)])
#define PHASE_ID_INNER_(x) (int_data[(NUM_INT_PROP_) + (x)]-1)
#define PHASE_ID_OUTER_(x) (int_data[(NUM_INT_PROP_) + (NUM_ADJACENT_PAIRS_) + (x)]-1)
#define NUM_AERO_PHASE_JAC_ELEM_INNER_(x) (int_data[(NUM_INT_PROP_) + (2*NUM_ADJACENT_PAIRS_) + (x)])
#define NUM_AERO_PHASE_JAC_ELEM_OUTER_(x) (int_data[(NUM_INT_PROP_) + (3*NUM_ADJACENT_PAIRS_) + (x)])
#define NUM_AERO_PHASE_JAC_ELEM_(x) (int_data[(NUM_INT_PROP_) + (4*NUM_ADJACENT_PAIRS_) + (x)])
#define NUM_JAC_ELEM_INNER_TOTAL_(x) (int_data[(NUM_INT_PROP_) + (5*NUM_ADJACENT_PAIRS_) + (x)])
#define AERO_SPEC_INNER_(x) (int_data[(NUM_INT_PROP_) + (6*NUM_ADJACENT_PAIRS_) + (x)]-1)
#define AERO_SPEC_OUTER_(x) (int_data[(NUM_INT_PROP_) + (7*NUM_ADJACENT_PAIRS_) + (x)]-1)
#define AERO_REP_ID_(x) (int_data[(NUM_INT_PROP_) + (8*NUM_ADJACENT_PAIRS_) + (x)]-1)

#define DERIV_ID_INNER_(x) (int_data[(NUM_INT_PROP_) + (9*NUM_ADJACENT_PAIRS_) + (x)])
#define DERIV_ID_OUTER_(x) (int_data[(NUM_INT_PROP_) + (10*NUM_ADJACENT_PAIRS_) + (x)])
#define JAC_ID_INNER_INNER_(x) \
  int_data[(NUM_INT_PROP_) + 11 * (NUM_ADJACENT_PAIRS_) + (x)]
#define JAC_ID_INNER_OUTER_(x) \
  int_data[(NUM_INT_PROP_) + 12 * (NUM_ADJACENT_PAIRS_) + (x)]
#define JAC_ID_OUTER_OUTER_(x) \
  int_data[(NUM_INT_PROP_) + 13 * (NUM_ADJACENT_PAIRS_) + (x)]
#define JAC_ID_OUTER_INNER_(x) \
  int_data[(NUM_INT_PROP_) + 14 * (NUM_ADJACENT_PAIRS_) + (x)]

#define PHASE_JAC_ID_INNER_(x, s, e) \
  int_data[(NUM_INT_PROP_) + 15 * (NUM_ADJACENT_PAIRS_) + NUM_AERO_PHASE_JAC_ELEM_INNER_(x) + (s) + (e)]
#define PHASE_JAC_ID_OUTER_(x, s, e) \
  int_data[(NUM_INT_PROP_) + 15 * (NUM_ADJACENT_PAIRS_) + NUM_JAC_ELEM_INNER_TOTAL_(x) + (s) * NUM_AERO_PHASE_JAC_ELEM_OUTER_(x) + (e)]

#define LAYER_THICKNESS_JAC_ELEM_INNER_(x,e) \
  (float_data[(NUM_FLOAT_PROP_) + 2 * (NUM_ADJACENT_PAIRS_) + \
              NUM_AERO_PHASE_JAC_ELEM_(x) + (e)])
#define LAYER_THICKNESS_JAC_ELEM_OUTER_(x,e) \
  (float_data[(NUM_FLOAT_PROP_) + 2 * (NUM_ADJACENT_PAIRS_) + NUM_AERO_PHASE_JAC_ELEM_(x) + \
              NUM_AERO_PHASE_JAC_ELEM_INNER_(x) + (e)])
#define PHASE_VOLUME_JAC_ELEM_INNER_(x,e) \
  (float_data[(NUM_FLOAT_PROP_) + 2 * (NUM_ADJACENT_PAIRS_) + \
              NUM_AERO_PHASE_JAC_ELEM_(x) + (e)])
#define PHASE_VOLUME_JAC_ELEM_OUTER_(x,e) \
  (float_data[(NUM_FLOAT_PROP_) + 2 * (NUM_ADJACENT_PAIRS_) + NUM_AERO_PHASE_JAC_ELEM_(x) + \
              NUM_AERO_PHASE_JAC_ELEM_INNER_(x) + (e)])
#define INTERFACE_SURFACE_AREA_JAC_ELEM_(x,e) \
  (float_data[(NUM_FLOAT_PROP_) + 2 * (NUM_ADJACENT_PAIRS_) + \
              NUM_AERO_PHASE_JAC_ELEM_(x) + NUM_AERO_PHASE_JAC_ELEM_INNER_(x) + \
              NUM_AERO_PHASE_JAC_ELEM_OUTER_(x) + (e)])


/** \brief Flag Jacobian elements used by this reaction
 *
 * \param model_data Pointer to the model data
 * \param rxn_int_data Pointer to the reaction integer data
 * \param rxn_float_data Pointer to the reaction floating-point data
 * \param jac Jacobian
 */
void rxn_condensed_phase_diffusion_get_used_jac_elem(ModelData *model_data,
                                                 int *rxn_int_data,
                                                 double *rxn_float_data,
                                                 Jacobian *jac) {
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;

  bool *aero_jac_elem =
      (bool *)malloc(sizeof(bool) * model_data->n_per_cell_state_var);
  if (aero_jac_elem == NULL) {
    printf(
        "\n\nERROR allocating space for 1D Jacobian structure array for "
        "condensed phase diffusion reaction\n\n");
    exit(1);
  }

  for (int i_adj_pairs = 0; i_adj_pairs < NUM_ADJACENT_PAIRS_; ++i_adj_pairs) {
    // Direct dependence on adjacent-pair species concentrations
    jacobian_register_element(jac, AERO_SPEC_INNER_(i_adj_pairs),
                              AERO_SPEC_INNER_(i_adj_pairs));
    jacobian_register_element(jac, AERO_SPEC_INNER_(i_adj_pairs),
                              AERO_SPEC_OUTER_(i_adj_pairs));
    jacobian_register_element(jac, AERO_SPEC_OUTER_(i_adj_pairs),
                              AERO_SPEC_INNER_(i_adj_pairs));
    jacobian_register_element(jac, AERO_SPEC_OUTER_(i_adj_pairs),
                              AERO_SPEC_OUTER_(i_adj_pairs));

    // Discover state dependencies contributed by inner-phase geometry terms.
    int n_inner_jac_elem_total = 0;
    for (int i_elem = 0; i_elem < model_data->n_per_cell_state_var; ++i_elem) {
      aero_jac_elem[i_elem] = false;
    }
    int n_inner_jac_elem =
        aero_rep_get_used_jac_elem(model_data, AERO_REP_ID_(i_adj_pairs),
                                   PHASE_ID_INNER_(i_adj_pairs), aero_jac_elem);
    if (n_inner_jac_elem > NUM_AERO_PHASE_JAC_ELEM_INNER_(i_adj_pairs)) {
      printf(
          "\n\nERROR Received more inner phase Jacobian elements than expected "
          "for condensed phase diffusion reaction. Got %d, expected <= %d",
          n_inner_jac_elem, NUM_AERO_PHASE_JAC_ELEM_INNER_(i_adj_pairs));
      exit(1);
    }
    n_inner_jac_elem_total += n_inner_jac_elem;

    // Any inner-phase dependency affects both inner and outer tendencies.
    int i_used_elem = 0;
    for (int i_elem = 0; i_elem < model_data->n_per_cell_state_var; ++i_elem) {
      if (aero_jac_elem[i_elem] == true) {
        jacobian_register_element(jac, AERO_SPEC_INNER_(i_adj_pairs), i_elem);
        jacobian_register_element(jac, AERO_SPEC_OUTER_(i_adj_pairs), i_elem);
        PHASE_JAC_ID_INNER_(i_adj_pairs, JAC_INNER_, i_used_elem) = i_elem;
        printf("Jacobian element for inner phase %d, elem %d: %d\n", i_adj_pairs,
               i_used_elem, PHASE_JAC_ID_INNER_(i_adj_pairs, JAC_INNER_, i_used_elem));
        ++i_used_elem;
      }
    }

    for (; i_used_elem < NUM_AERO_PHASE_JAC_ELEM_INNER_(i_adj_pairs);
         ++i_used_elem) {
      PHASE_JAC_ID_INNER_(i_adj_pairs, JAC_INNER_, i_used_elem) = -1;
    }
    if (i_used_elem != n_inner_jac_elem) {
      printf(
          "\n\nERROR Error setting used Jacobian elements in condensed phase "
          "diffusion reaction %d %d\n\n",
          i_used_elem, n_inner_jac_elem);
      exit(1);
    }

    // Discover state dependencies contributed by outer-phase geometry terms.
    int n_outer_jac_elem_total = 0;
    for (int i_elem = 0; i_elem < model_data->n_per_cell_state_var; ++i_elem) {
      aero_jac_elem[i_elem] = false;
    }
    int n_outer_jac_elem =
        aero_rep_get_used_jac_elem(model_data, AERO_REP_ID_(i_adj_pairs),
                                   PHASE_ID_OUTER_(i_adj_pairs), aero_jac_elem);
    if (n_outer_jac_elem > NUM_AERO_PHASE_JAC_ELEM_OUTER_(i_adj_pairs)) {
      printf(
          "\n\nERROR Received more outer phase Jacobian elements than expected "
          "for condensed phase diffusion reaction. Got %d, expected <= %d",
          n_outer_jac_elem, NUM_AERO_PHASE_JAC_ELEM_OUTER_(i_adj_pairs));
      exit(1);
    }
    n_outer_jac_elem_total += n_outer_jac_elem;

    // Any outer-phase dependency also affects both tendency equations.
    i_used_elem = 0;
    for (int i_elem = 0; i_elem < model_data->n_per_cell_state_var; ++i_elem) {
      if (aero_jac_elem[i_elem] == true) {
        jacobian_register_element(jac, AERO_SPEC_INNER_(i_adj_pairs), i_elem);
        jacobian_register_element(jac, AERO_SPEC_OUTER_(i_adj_pairs), i_elem);
        PHASE_JAC_ID_OUTER_(i_adj_pairs, JAC_OUTER_, i_used_elem) = i_elem;
        printf("Jacobian element for outer phase %d, elem %d: %d\n", i_adj_pairs,
               i_elem, PHASE_JAC_ID_OUTER_(i_adj_pairs, JAC_OUTER_, i_used_elem));
        ++i_used_elem;
      }
    }

    for (; i_used_elem < NUM_AERO_PHASE_JAC_ELEM_OUTER_(i_adj_pairs);
         ++i_used_elem) {
      PHASE_JAC_ID_OUTER_(i_adj_pairs, JAC_OUTER_, i_used_elem) = -1;
    }
    if (i_used_elem != n_outer_jac_elem) {
      printf(
          "\n\nERROR Error setting used Jacobian elements in condensed phase "
          "diffusion reaction %d %d\n\n",
          i_used_elem, n_outer_jac_elem);
      exit(1);
    }

    // Validate the packed scratch-space capacity used by downstream Jacobian code.
    if (n_inner_jac_elem_total + n_outer_jac_elem_total >
        NUM_AERO_PHASE_JAC_ELEM_(i_adj_pairs)) {
      printf(
          "\n\nERROR Received more total Jacobian elements than expected for "
          "condensed phase diffusion reaction. Got %d, expected <= %d",
          n_inner_jac_elem_total + n_outer_jac_elem_total,
          NUM_AERO_PHASE_JAC_ELEM_(i_adj_pairs));
      exit(1);
    }
  }
  free(aero_jac_elem);

}

/** \brief Update the time derivative and Jacbobian array indices
 *
 * \param model_data Pointer to the model data for finding sub model ids
 * \param deriv_ids Id of each state variable in the derivative array
 * \param jac Jacobian
 * \param rxn_int_data Pointer to the reaction integer data
 * \param rxn_float_data Pointer to the reaction floating-point data
 */
void rxn_condensed_phase_diffusion_update_ids(ModelData *model_data, int *deriv_ids,
                                          Jacobian jac, int *rxn_int_data,
                                          double *rxn_float_data) {
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;

  for (int i_adj_pairs = 0; i_adj_pairs < NUM_ADJACENT_PAIRS_; ++i_adj_pairs) {
     DERIV_ID_INNER_(i_adj_pairs) = deriv_ids[AERO_SPEC_INNER_(i_adj_pairs)];
     DERIV_ID_OUTER_(i_adj_pairs) = deriv_ids[AERO_SPEC_OUTER_(i_adj_pairs)];
  }

  // Update the Jacobian element ids
  int i_jac = 0;
  for (int i_adj_pairs = 0; i_adj_pairs < NUM_ADJACENT_PAIRS_; ++i_adj_pairs) {
    JAC_ID_INNER_INNER_(i_jac++) = jacobian_get_element_id(jac, AERO_SPEC_INNER_(i_adj_pairs), AERO_SPEC_INNER_(i_adj_pairs));
    JAC_ID_INNER_OUTER_(i_jac++) = jacobian_get_element_id(jac, AERO_SPEC_INNER_(i_adj_pairs), AERO_SPEC_OUTER_(i_adj_pairs));
    JAC_ID_OUTER_OUTER_(i_jac++) = jacobian_get_element_id(jac, AERO_SPEC_OUTER_(i_adj_pairs), AERO_SPEC_OUTER_(i_adj_pairs));
    JAC_ID_OUTER_INNER_(i_jac++) = jacobian_get_element_id(jac, AERO_SPEC_OUTER_(i_adj_pairs), AERO_SPEC_INNER_(i_adj_pairs));
    printf("AERO_SPEC_INNER_(%d) = %d, AERO_SPEC_OUTER_(%d) = %d\n", i_adj_pairs, AERO_SPEC_INNER_(i_adj_pairs), i_adj_pairs, AERO_SPEC_OUTER_(i_adj_pairs));
  }

    // Save non-zero Jacobian element indices for aerosol representation
    // function dependencies. We use the state-variable indices stored
    // previously in PHASE_JAC_ID_INNER_ and PHASE_JAC_ID_OUTER_ to look up the 
    // corresponding Jacobian element index in the flattened sparse matrix.
    // We do this for both the dependence of the INNER species and the
    // OUTER species on each independent variable used by the aerosol
    // representation functions.
    for (int i_adj_pairs = 0; i_adj_pairs < NUM_ADJACENT_PAIRS_; ++i_adj_pairs) {
      for (int i_elem = 0; i_elem < NUM_AERO_PHASE_JAC_ELEM_INNER_(i_adj_pairs); ++i_elem) {
        if (PHASE_JAC_ID_INNER_(i_adj_pairs, JAC_INNER_, i_elem) > 0) {
          PHASE_JAC_ID_INNER_(i_adj_pairs, JAC_INNER_, i_elem) =
              jacobian_get_element_id(jac, AERO_SPEC_INNER_(i_adj_pairs),
                                      PHASE_JAC_ID_INNER_(i_adj_pairs, JAC_INNER_, i_elem));
        }
      }
      for (int i_elem = 0; i_elem < NUM_AERO_PHASE_JAC_ELEM_OUTER_(i_adj_pairs); ++i_elem) {
        if (PHASE_JAC_ID_OUTER_(i_adj_pairs, JAC_OUTER_, i_elem) > 0) {
          PHASE_JAC_ID_OUTER_(i_adj_pairs, JAC_OUTER_, i_elem) =
              jacobian_get_element_id(jac, AERO_SPEC_OUTER_(i_adj_pairs),
                                      PHASE_JAC_ID_OUTER_(i_adj_pairs, JAC_OUTER_, i_elem));
        }
      }
    }
}

/** \brief Update reaction data for new environmental conditions
 *
 *
 * \param model_data Pointer to the model data
 * \param rxn_int_data Pointer to the reaction integer data
 * \param rxn_float_data Pointer to the reaction floating-point data
 * \param rxn_env_data Pointer to the environment-dependent parameters
 */
void rxn_condensed_phase_diffusion_update_env_state(ModelData *model_data,
                                                int *rxn_int_data,
                                                double *rxn_float_data,
                                                double *rxn_env_data) {
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;
  double *env_data = model_data->grid_cell_env;

  return;
}

/** \brief Calculate contributions to the time derivative \f$f(t,y)\f$ from
 * this reaction.
 *
 * \param model_data Pointer to the model data, including the state array
 * \param time_deriv TimeDerivative object
 * \param rxn_int_data Pointer to the reaction integer data
 * \param rxn_float_data Pointer to the reaction floating-point data
 * \param rxn_env_data Pointer to the environment-dependent parameters
 * \param time_step Current time step being computed (s)
 */
#ifdef CAMP_USE_SUNDIALS
void rxn_condensed_phase_diffusion_calc_deriv_contrib(
    ModelData *model_data, TimeDerivative time_deriv, int *rxn_int_data,
    double *rxn_float_data, double *rxn_env_data, realtype time_step) {
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;
  double *state = model_data->grid_cell_state;
  double *env_data = model_data->grid_cell_env;

  // Calculate derivative contributions for each aerosol phase
  for (int i_adj_pairs = 0, i_deriv = 0; i_adj_pairs < NUM_ADJACENT_PAIRS_; i_adj_pairs++) {

    /* Get the layer thickness for inner phase id (m) */
    realtype layer_thickness_inner;
    aero_rep_get_layer_thickness__m(
      model_data, //model data 
      AERO_REP_ID_(i_adj_pairs), // aerosol representation index
      PHASE_ID_INNER_(i_adj_pairs), // inner phase id
      &layer_thickness_inner, // layer thickness 
      NULL); // partial derivative

    // Get the layer thickness for outer phase id (m)
    realtype layer_thickness_outer;
    aero_rep_get_layer_thickness__m(
      model_data, //model data 
      AERO_REP_ID_(i_adj_pairs), // aerosol representation index
      PHASE_ID_OUTER_(i_adj_pairs), // outer phase id
      &layer_thickness_outer, // layer thickness 
      NULL); // partial derivative

    // Get the volume of the inner phase
    realtype volume_phase_inner;
    aero_rep_get_phase_volume__m3_m3(
        model_data, //model data
        AERO_REP_ID_(i_adj_pairs), // aerosol representation index
        PHASE_ID_INNER_(i_adj_pairs), // inner phase id
        &volume_phase_inner, // volume of inner phase
        NULL); // partial derivative

    // Get the volume of the outer phase
    realtype volume_phase_outer;
    aero_rep_get_phase_volume__m3_m3(
        model_data, //model data
        AERO_REP_ID_(i_adj_pairs), // aerosol representation index
        PHASE_ID_OUTER_(i_adj_pairs), // outer phase id
        &volume_phase_outer, // volume of outer phase
        NULL); // partial derivative

    // Get the interface surface area (m2)
    realtype eff_sa;
    aero_rep_get_interface_surface_area__m2(
        model_data, //model data 
        AERO_REP_ID_(i_adj_pairs), // aerosol representation index
        PHASE_ID_INNER_(i_adj_pairs), // inner phase id
        PHASE_ID_OUTER_(i_adj_pairs), // outer phase id
        &eff_sa, // interface surface area 
        NULL); // partial derivative

    // Calculate the rate constant for diffusion limited mass transfer between
    // particle layers
    double rate_inner_loss = (double)(eff_sa / volume_phase_inner);
    double rate_inner_prod = (double)(eff_sa / volume_phase_inner);
    double rate_outer_loss = (double)(eff_sa / volume_phase_outer);
    double rate_outer_prod = (double)(eff_sa / volume_phase_outer);

    rate_inner_loss *= ((DIFF_COEFF_INNER_(i_adj_pairs) / layer_thickness_inner) 
                    * state[AERO_SPEC_INNER_(i_adj_pairs)]);

    rate_inner_prod *= ((DIFF_COEFF_OUTER_(i_adj_pairs) / layer_thickness_outer) 
                    * state[AERO_SPEC_OUTER_(i_adj_pairs)]);

    rate_outer_loss *= ((DIFF_COEFF_OUTER_(i_adj_pairs) / layer_thickness_outer) 
                    * state[AERO_SPEC_OUTER_(i_adj_pairs)]);

    rate_outer_prod *= ((DIFF_COEFF_INNER_(i_adj_pairs) / layer_thickness_inner)
                    * state[AERO_SPEC_INNER_(i_adj_pairs)]);
    
    if (DERIV_ID_INNER_(i_adj_pairs) >= 0) {
      time_derivative_add_value(time_deriv, DERIV_ID_INNER_(i_adj_pairs),
                                  -rate_inner_loss);
    }
    if (DERIV_ID_INNER_(i_adj_pairs) >= 0) {
      time_derivative_add_value(time_deriv, DERIV_ID_INNER_(i_adj_pairs),
                                  rate_inner_prod);
    }
    if (DERIV_ID_OUTER_(i_adj_pairs) >= 0) {
      time_derivative_add_value(time_deriv, DERIV_ID_OUTER_(i_adj_pairs), -rate_outer_loss);
    }
    if (DERIV_ID_OUTER_(i_adj_pairs) >= 0) {
      time_derivative_add_value(time_deriv, DERIV_ID_OUTER_(i_adj_pairs), rate_outer_prod);
    }

  }
  return;
}
#endif

/** \brief Calculate contributions to the Jacobian from this reaction
 *
 * \param model_data Pointer to the model data
 * \param jac Reaction Jacobian
 * \param rxn_int_data Pointer to the reaction integer data
 * \param rxn_float_data Pointer to the reaction floating-point data
 * \param rxn_env_data Pointer to the environment-dependent parameters
 * \param time_step Current time step being calculated (s)
 */
#ifdef CAMP_USE_SUNDIALS
void rxn_condensed_phase_diffusion_calc_jac_contrib(ModelData *model_data,
                                                Jacobian jac, int *rxn_int_data,
                                                double *rxn_float_data,
                                                double *rxn_env_data,
                                                realtype time_step) {
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;
  double *state = model_data->grid_cell_state;
  double *env_data = model_data->grid_cell_env;

  // Calculate derivative contributions for each adjacent phase pair
  // Fill Jacobians from aerosol representation functions for each adjacent phase pair
  for (int i_adj_pairs = 0, i_deriv = 0; i_adj_pairs < NUM_ADJACENT_PAIRS_; i_adj_pairs++) {

    /* Get the layer thickness for inner phase id (m) */
    realtype layer_thickness_inner;
    aero_rep_get_layer_thickness__m(
      model_data, //model data 
      AERO_REP_ID_(i_adj_pairs), // aerosol representation index
      PHASE_ID_INNER_(i_adj_pairs), // inner phase id
      &layer_thickness_inner, // layer thickness 
      &(LAYER_THICKNESS_JAC_ELEM_INNER_(i_adj_pairs, 0))); // partial derivative

    // Get the layer thickness for outer phase id (m)
    realtype layer_thickness_outer;
    aero_rep_get_layer_thickness__m(
      model_data, //model data 
      AERO_REP_ID_(i_adj_pairs), // aerosol representation index
      PHASE_ID_OUTER_(i_adj_pairs), // outer phase id
      &layer_thickness_outer, // layer thickness 
      &(LAYER_THICKNESS_JAC_ELEM_OUTER_(i_adj_pairs, 0))); // partial derivative

    // Get the volume of the inner phase
    realtype volume_phase_inner;
    aero_rep_get_phase_volume__m3_m3(
        model_data, //model data
        AERO_REP_ID_(i_adj_pairs), // aerosol representation index
        PHASE_ID_INNER_(i_adj_pairs), // inner phase id
        &volume_phase_inner, // volume of inner phase
        &(PHASE_VOLUME_JAC_ELEM_INNER_(i_adj_pairs, 0))); // partial derivative

    // Get the volume of the outer phase
    realtype volume_phase_outer;
    aero_rep_get_phase_volume__m3_m3(
        model_data, //model data
        AERO_REP_ID_(i_adj_pairs), // aerosol representation index
        PHASE_ID_OUTER_(i_adj_pairs), // outer phase id
        &volume_phase_outer, // volume of outer phase
        &(PHASE_VOLUME_JAC_ELEM_OUTER_(i_adj_pairs, 0))); // partial derivative

    // Get the interface surface area (m2)
    realtype eff_sa;
    aero_rep_get_interface_surface_area__m2(
        model_data, //model data 
        AERO_REP_ID_(i_adj_pairs), // aerosol representation index
        PHASE_ID_INNER_(i_adj_pairs), // inner phase id
        PHASE_ID_OUTER_(i_adj_pairs), // outer phase id
        &eff_sa, // interface surface area 
        &(INTERFACE_SURFACE_AREA_JAC_ELEM_(i_adj_pairs, 0))); // partial derivative

    // Calculate the rate constant for diffusion limited mass transfer between
    // particle layers
    double rate_inner_loss = (double)(eff_sa / volume_phase_inner);
    double rate_inner_prod = (double)(eff_sa / volume_phase_inner);
    double rate_outer_loss = (double)(eff_sa / volume_phase_outer);
    double rate_outer_prod = (double)(eff_sa / volume_phase_outer);

    rate_inner_loss *= ((DIFF_COEFF_INNER_(i_adj_pairs) / layer_thickness_inner));

    rate_inner_prod *= ((DIFF_COEFF_OUTER_(i_adj_pairs) / layer_thickness_outer));

    rate_outer_loss *= ((DIFF_COEFF_OUTER_(i_adj_pairs) / layer_thickness_outer));

    rate_outer_prod *= ((DIFF_COEFF_INNER_(i_adj_pairs) / layer_thickness_inner));

    if (JAC_ID_INNER_INNER_(i_adj_pairs) >= 0) {
      jacobian_add_value(jac, (unsigned int)JAC_ID_INNER_INNER_(i_adj_pairs), JACOBIAN_LOSS, -rate_inner_loss);
    }
    if (JAC_ID_INNER_OUTER_(i_adj_pairs) >= 0) {
      jacobian_add_value(jac, (unsigned int)JAC_ID_INNER_OUTER_(i_adj_pairs), JACOBIAN_PRODUCTION, rate_inner_prod);
    }
    if (JAC_ID_OUTER_OUTER_(i_adj_pairs) >= 0) {
      jacobian_add_value(jac, (unsigned int)JAC_ID_OUTER_OUTER_(i_adj_pairs), JACOBIAN_LOSS, -rate_outer_loss);
    }
    if (JAC_ID_OUTER_INNER_(i_adj_pairs) >= 0) {
      jacobian_add_value(jac, (unsigned int)JAC_ID_OUTER_INNER_(i_adj_pairs), JACOBIAN_PRODUCTION, rate_outer_prod);
    }

    realtype gamma = (eff_sa * DIFF_COEFF_INNER_(i_adj_pairs)) / (layer_thickness_inner * layer_thickness_inner);
    realtype alpha = (eff_sa * DIFF_COEFF_INNER_(i_adj_pairs)) / (layer_thickness_inner * volume_phase_outer);
    realtype epsilon = (eff_sa * DIFF_COEFF_OUTER_(i_adj_pairs)) / (layer_thickness_outer * volume_phase_outer);
    realtype beta = (eff_sa * DIFF_COEFF_OUTER_(i_adj_pairs)) / (layer_thickness_outer * volume_phase_inner);

    realtype deriv_inner_loss = -(gamma * state[AERO_SPEC_INNER_(i_adj_pairs)]);
    realtype deriv_inner_prod = (beta * state[AERO_SPEC_OUTER_(i_adj_pairs)]);
    realtype deriv_outer_loss = -(epsilon * state[AERO_SPEC_OUTER_(i_adj_pairs)]);
    realtype deriv_outer_prod = (alpha * state[AERO_SPEC_INNER_(i_adj_pairs)]);

    for (int i_elem = 0; i_elem < NUM_AERO_PHASE_JAC_ELEM_INNER_(i_adj_pairs); ++i_elem) {
      deriv_inner_loss += ( - gamma / eff_sa) * state[AERO_SPEC_INNER_(i_adj_pairs)] *
           INTERFACE_SURFACE_AREA_JAC_ELEM_(i_adj_pairs, i_elem) + 
           (gamma / volume_phase_inner) * state[AERO_SPEC_INNER_(i_adj_pairs)] * 
           PHASE_VOLUME_JAC_ELEM_INNER_(i_adj_pairs, i_elem) + ( gamma / layer_thickness_inner ) * 
           state[AERO_SPEC_INNER_(i_adj_pairs)] * LAYER_THICKNESS_JAC_ELEM_INNER_(i_adj_pairs, i_elem);
      deriv_inner_prod += ( beta / eff_sa) * state[AERO_SPEC_OUTER_(i_adj_pairs)] *
           INTERFACE_SURFACE_AREA_JAC_ELEM_(i_adj_pairs, i_elem) - (beta / volume_phase_inner) * 
           state[AERO_SPEC_OUTER_(i_adj_pairs)] * PHASE_VOLUME_JAC_ELEM_INNER_(i_adj_pairs, i_elem) - 
           ( beta / layer_thickness_inner ) * state[AERO_SPEC_OUTER_(i_adj_pairs)] * 
           LAYER_THICKNESS_JAC_ELEM_INNER_(i_adj_pairs, i_elem);

      if (PHASE_JAC_ID_INNER_(i_adj_pairs, JAC_INNER_, i_elem) >= 0) {
        jacobian_add_value(jac, (unsigned int)PHASE_JAC_ID_INNER_(i_adj_pairs, JAC_INNER_, i_elem), JACOBIAN_LOSS, deriv_inner_loss);
        jacobian_add_value(jac, (unsigned int)PHASE_JAC_ID_INNER_(i_adj_pairs, JAC_INNER_, i_elem), JACOBIAN_PRODUCTION, deriv_inner_prod);
      }

    }
    for (int i_elem = 0; i_elem < NUM_AERO_PHASE_JAC_ELEM_OUTER_(i_adj_pairs); ++i_elem) {
      deriv_outer_loss += ( - epsilon / eff_sa) * state[AERO_SPEC_OUTER_(i_adj_pairs)] *
           INTERFACE_SURFACE_AREA_JAC_ELEM_(i_adj_pairs, i_elem) + (epsilon / volume_phase_outer) * 
           state[AERO_SPEC_OUTER_(i_adj_pairs)] * PHASE_VOLUME_JAC_ELEM_OUTER_(i_adj_pairs, i_elem) + 
           (epsilon / layer_thickness_outer) * state[AERO_SPEC_OUTER_(i_adj_pairs)] *
           LAYER_THICKNESS_JAC_ELEM_OUTER_(i_adj_pairs, i_elem);
      deriv_outer_prod += ( alpha / eff_sa) * state[AERO_SPEC_INNER_(i_adj_pairs)] *
           INTERFACE_SURFACE_AREA_JAC_ELEM_(i_adj_pairs, i_elem) - (alpha / volume_phase_outer) * 
           state[AERO_SPEC_INNER_(i_adj_pairs)] * PHASE_VOLUME_JAC_ELEM_OUTER_(i_adj_pairs, i_elem) - 
           (alpha / layer_thickness_outer) * state[AERO_SPEC_INNER_(i_adj_pairs)] *
           LAYER_THICKNESS_JAC_ELEM_OUTER_(i_adj_pairs, i_elem);

      if (PHASE_JAC_ID_OUTER_(i_adj_pairs, JAC_OUTER_, i_elem) >= 0) {
        jacobian_add_value(jac, (unsigned int)PHASE_JAC_ID_OUTER_(i_adj_pairs, JAC_OUTER_, i_elem), JACOBIAN_LOSS, deriv_outer_loss);
        jacobian_add_value(jac, (unsigned int)PHASE_JAC_ID_OUTER_(i_adj_pairs, JAC_OUTER_, i_elem), JACOBIAN_PRODUCTION, deriv_outer_prod);
      }
    }
  }
    return;
}
#endif

/** \brief Print the Phase Transfer reaction parameters
 *
 * \param rxn_int_data Pointer to the reaction integer data
 * \param rxn_float_data Pointer to the reaction floating-point data
 */
void rxn_condensed_phase_diffusion_print(int *rxn_int_data,
                                     double *rxn_float_data) {
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;

  printf("\n\nCondensed Phase Diffusion reaction\n");
  for (int i = 0; i < NUM_ADJACENT_PAIRS_; ++i) {
    printf("\n  Diffusion coefficient inner: %g", DIFF_COEFF_INNER_(i));
    printf("\n  Diffusion coefficient outer: %g", DIFF_COEFF_OUTER_(i));
    printf("\n  Aerosol phase id inner: %d", PHASE_ID_INNER_(i));
    printf("\n  Aerosol phase id outer: %d", PHASE_ID_OUTER_(i));
    printf("\n  Aerosol species id inner: %d", AERO_SPEC_INNER_(i));
    printf("\n  Aerosol species id outer: %d", AERO_SPEC_OUTER_(i));
    printf("\n  Aerosol representation id: %d", AERO_REP_ID_(i));
}

  return;
}
