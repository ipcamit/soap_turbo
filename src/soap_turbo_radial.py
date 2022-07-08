from signal import Sigmasks
import numpy as np
from scipy.special import erf

def N_a(rcut, a):
    return np.sqrt(rcut/(2 * a + 5))


def get_ortho_mat_poly3(alpha_max):
    """
    ! This subroutine returns the overlap matrix S and the orthonormalization
    ! matrix W = S^{-1/2} needed to construct the orthonormal basis from the
    ! polynomial basis set. It is also needed to transform between original
    ! basis and orthonormal basis. It does not requires blas/lapack to work.
    """
    S = np.zeros((alpha_max, alpha_max))
    W = np.zeros((alpha_max, alpha_max))

    for i in range(alpha_max):
        S[i, i] = 1.
        for j in range(i + 1):
            S[i, j] = np.sqrt((5 + 2*i) * (5 + 2*j))/(5 + i + j) # why?
            S[j, i] = S[i, j]

    U,W,V = np.linalg.svd(S)
    W = U @ np.diag(np.sqrt(W)) @ V
    W = np.linalg.inv(W)

    return S, W


def get_rad_exp_coeff_poly3(alpha_max:int,
                            n_sites:int,
                            n_neigh:np.ndarray,
                            r_cut_soft_in:float, r_cut_hard_in:float,
                            atom_sigma_in:float,
                            do_central:bool,
                            rjs_in:np.ndarray,
                            mask:np.ndarray,
                            atom_sigma_scaling:float,
                            scaling_mode:str
                            ):
    r_cut_soft = r_cut_soft_in/r_cut_hard_in
    r_cut_hard = 1.0
    atom_sigma = atom_sigma_in/r_cut_hard_in
    dr = 1.0 - r_cut_soft

    k = -1
    for i in range(n_sites):
        for j in range(n_neigh[i]):
            k += 1
            if (j == 0) and not do_central:
                continue
            if (rjs_in[k] < r_cut_hard_in and mask[k]):
                exp_coeff_temp1 = 0.
                exp_coeff_temp2 = 0.
                exp_coeff_der_temp = 0.
                rj = rjs_in[k]/r_cut_hard_in
                atom_sigma_scaled = atom_sigma + atom_sigma_scaling*rj
                s2 = atom_sigma_scaled**2
                




#   integer, intent(in) :: alpha_max, n_neigh(:), n_sites, radial_enhancement
#     real*8, intent(in) :: rcut_soft_in, rcut_hard_in, rjs_in(:), atom_sigma_in, nf, atom_sigma_scaling
#     real*8, intent(in) :: amplitude_scaling, central_weight
#     real*8 :: rcut_soft, rcut_hard, atom_sigma, atom_sigma_scaled, amplitude
#     logical, intent(in) :: mask(:), do_derivatives, do_central
#     character(*), intent(in) :: scaling_mode