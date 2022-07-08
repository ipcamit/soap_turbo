from typing import List
import numpy as np
from soap_turbo_radial import get_ortho_mat_poly3

def get_goap(n_species:int,
             central_weight:np.ndarray,
             rcut_hard:np.ndarray,
             alpha_max:np.ndarray,
             l_max:int,
             n_atom_pairs:int,
             n_sites:int
             ):
    
    # ! Check if we need to expand the central atom
    do_central = np.where(central_weight !=0, True, False)
    rcut_max = np.max(rcut_hard)
    
    # ! This assigns begin and end indices to the components of the radial basis
    i_end = np.cumsum(alpha_max)
    i_beg = np.zeros_like(i_end)
    i_beg[-(n_species - 1):] = i_end[:-1]

    k_max = 1 + l_max * (l_max + 1)/2 + l_max

    # ! This is to build the radial basis
    n_max = np.sum(alpha_max)

    # Assuming recompute to be always true

    W = np.zeros((n_max, n_max))
    S = np.zeros((n_max, n_max))

    # !   This is done per species. Each loop iteration modifies the slice of the W and S matrices that
    # !   corresponds to the species in question. <b>The radial basis functions for species A are always
    # !   assumed orthogonal to the basis functions for species B. W and S are therefore block diagonal.</b>

    for i in range(n_species):
        # Only assuming basis to be poly3
        S_temp, W_temp = get_ortho_mat_poly3(alpha_max[i])
        S[i_beg[i]:i_end[i],i_beg[i]:i_end[i]] = S_temp
        W[i_beg[i]:i_end[i],i_beg[i]:i_end[i]] = W_temp

    # ! Uncompressed SOAP dimensions
    n_soap_uncompressed = n_max*(n_max+1)/2 * (l_max+1)

    radial_exp_coeff = np.zeros((n_max, n_atom_pairs) )
    angular_exp_coeff = np.zeros((k_max, n_atom_pairs) )

    cnk = np.zeros((k_max, n_max, n_sites))

    # skipping compressing part for now

    



#   real*8, intent(in) :: rjs(:), thetas(:), phis(:)
#   real*8, intent(in) :: amplitude_scaling(:), atom_sigma_r_scaling(:), atom_sigma_t(:), atom_sigma_t_scaling(:)
#   real*8, intent(in) :: central_weight(:), atom_sigma_r(:), global_scaling(:)
#   real*8, intent(in) :: nf(:), rcut_hard(:), rcut_soft(:)
#   real*8, intent(in) :: compress_P_el(:)

#   integer, intent(in) :: n_species, radial_enhancement, species(:,:), species_multiplicity(:)
#   integer, intent(in) :: n_sites, n_neigh(:), , n_atom_pairs, alpha_max(:), compress_P_nonzero, &
#                          compress_P_i(:), compress_P_j(:)

#   logical, intent(in) :: do_derivatives, do_timing, mask(:,:), compress_soap

#   character(*), intent(in) :: basis, scaling_mode

# ! Output variables
#   real*8, intent(inout) :: soap(:,:), soap_cart_der(:,:,:)
# !-------------------


# !-------------------
# ! Internal variables
#   complex*16, allocatable :: angular_exp_coeff(:,:), cnk(:,:,:)
#   complex*16, allocatable :: angular_exp_coeff_rad_der(:,:), angular_exp_coeff_azi_der(:,:), cnk_rad_der(:,:,:)
#   complex*16, allocatable :: cnk_azi_der(:,:,:), angular_exp_coeff_pol_der(:,:), cnk_pol_der(:,:,:)
#   complex*16, allocatable :: eimphi(:), prefm(:), eimphi_rad_der(:)

#   real*8, allocatable, save :: W(:,:), S(:,:), multiplicity_array(:)
#   real*8, allocatable :: soap_rad_der(:,:), sqrt_dot_p(:), soap_azi_der(:,:)
#   real*8, allocatable :: W_temp(:,:), S_temp(:,:)
#   real*8, allocatable :: radial_exp_coeff(:,:), soap_pol_der(:,:)
#   real*8, allocatable :: preflm(:), plm_array(:), prefl(:), fact_array(:), prefl_rad_der(:)
#   real*8, allocatable :: radial_exp_coeff_der(:,:)
#   real*8, allocatable :: this_soap(:), this_soap_rad_der(:), this_soap_azi_der(:), this_soap_pol_der(:)
#   real*8 :: amplitude, multiplicity, pi, rcut_max
#   real*8 :: radial_time, angular_time, coeff_time, time3, total_time, soap_time, time1, time2, compress_time, &
#             memory_time, basis_time

#   integer, allocatable :: i_beg(:), i_end(:)
#   integer, save :: n_max_prev, l_max_prev
#   integer :: k_max, n_max
#   integer :: i, counter, j, k, n_soap, k2, k3, n, l, m, np, counter2, n_soap_uncompressed
#   logical, allocatable :: do_central(:), skip_soap_component(:,:,:), skip_soap_component_flattened(:)
#   logical, allocatable, save :: skip_soap_component_flattened_prev(:)
#   logical, save :: recompute_basis = .true., recompute_multiplicity_array = .true.