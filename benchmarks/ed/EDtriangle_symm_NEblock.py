#!/usr/bin/env python3

import numpy as np
import time
import logging
from quspin.basis.user import user_basis, next_state_sig_32, op_sig_32, map_sig_32
from quspin.operators import hamiltonian
from numba import cfunc, carray, uint32, int32, float64
import math

# Set up logging
logging.basicConfig(filename="calc.log", level=logging.INFO, format="%(asctime)s - %(message)s")

def read_params(filename="params.txt"):
    try:
        with open(filename, 'r') as f:
            params = {}
            for line in f:
                if line.strip() and not line.startswith("#"):
                    key, value = line.split("=")
                    params[key.strip()] = float(value.strip())
            return params
    except FileNotFoundError:
        logging.error(f"Error: {filename} not found.")
        raise

import numpy as np
from quspin.basis.user import user_basis, next_state_sig_32, op_sig_32, map_sig_32
from quspin.operators import hamiltonian
from numba import cfunc, carray, uint32, int32, float64
import math

#############################################
# 1. Mixed Radix External Basis Generation Function (Two Species)
#############################################
def generate_two_species_basis(Lx, Ly, NE_b, NE_c):
    """
    For a two-dimensional lattice, assume the first Lq = Lx*Ly sites represent b particles,
    and the next Lq sites represent c particles.
    For b particles, the local dimension is sps_b = NE_b + 1, and for c particles, sps_c = NE_c + 1.
    Only states with a fixed total particle number (NE_b for b and NE_c for c) are generated using mixed radix encoding.

    Returns:
      - An array of external basis states (dtype=uint32)
      - The local dimension array sps_arr
      - The full_factor array for each digit (mixed radix factors)
    """
    Lq = Lx * Ly
    L = 2 * Lq
    sps_b = NE_b + 1
    sps_c = NE_c + 1
    # Construct the local dimension array for the whole system
    sps_arr = np.array([sps_b] * Lq + [sps_c] * Lq, dtype=np.uint32)
    # Calculate the mixed radix factors: full_factor[0]=1, full_factor[i] = full_factor[i-1] * sps_arr[i-1]
    full_factor = np.empty(L, dtype=np.uint32)
    full_factor[0] = 1
    for i in range(1, L):
        full_factor[i] = full_factor[i - 1] * sps_arr[i - 1]

    # Generate b and c configurations separately
    def gen_configs(n_sites, NE, d):
        configs = []
        config = [0] * n_sites
        def rec(i, rem):
            if i == n_sites - 1:
                if rem < d:
                    config[i] = rem
                    configs.append(config.copy())
                return
            for n in range(min(d, rem + 1)):
                config[i] = n
                rec(i + 1, rem - n)
        rec(0, NE)
        return configs

    b_configs = gen_configs(Lq, NE_b, sps_b)
    c_configs = gen_configs(Lq, NE_c, sps_c)
    valid_states = []
    for b_conf in b_configs:
        # Encode the b particle part
        state_b = 0
        for i in range(Lq):
            state_b += b_conf[i] * full_factor[i]
        for c_conf in c_configs:
            state = state_b
            for i in range(Lq):
                state += c_conf[i] * full_factor[Lq + i]
            valid_states.append(state)
    return np.array(valid_states, dtype=np.uint32), sps_arr, full_factor

#############################################
# 2. Local Operator Function for Two Species (Mixed Op)
#############################################
@cfunc(op_sig_32, locals=dict(n=int32, new_state=uint32, fac=float64))
def op_two_species(op_struct_ptr, op_str, ind, N, args):
    """
    Apply a local operator on site 'ind':
      - Uses the mixed local dimension array and factors provided in args (each of length N)
      - op_str is one of '+' (creation), '-' (annihilation), or 'n' (number operator)
    """
    arr = carray(args, 2 * N)  # First N: local dimensions; next N: full_factor
    sps = arr[:N]
    full_fac = arr[N:2 * N]
    op_struct = carray(op_struct_ptr, 1)[0]
    local_sps = sps[ind]
    base = full_fac[ind]
    n = (op_struct.state // base) % local_sps
    new_state = op_struct.state
    fac = 1.0
    if op_str == ord('+'):
        if n < local_sps - 1:
            fac = math.sqrt(n + 1)
            new_state = op_struct.state + base
        else:
            fac = 0.0
    elif op_str == ord('-'):
        if n > 0:
            fac = math.sqrt(n)
            new_state = op_struct.state - base
        else:
            fac = 0.0
    elif op_str == ord('n'):
        fac = float(n)
    else:
        fac = 0.0
    op_struct.state = new_state
    op_struct.matrix_ele *= fac
    return 0

#############################################
# 3. Translation Mapping Function (Dummy Implementation)
#############################################
# --- translation_map_2species function ---
@cfunc(map_sig_32, locals=dict(L_total=uint32, offset1=uint32, offset2=uint32, offset3=uint32, i=uint32, inv=uint32, base=uint32, digit=uint32, sps_val=uint32, Lq=uint32, new_state=uint32))
def translation_map_2species(state, N_dummy, sign_ptr, args):
    """
    Apply translation mapping on a state for two species (using mixed radix encoding).

    Parameters:
      - state: the integer encoding of the original state
      - N_dummy: number of lattice sites in the translation direction (not directly used here)
      - sign_ptr: pointer to the sign factor (can be ignored if no antisymmetry is present)
      - args: a uint32 array storing in order:
            [ L_total,  sps_arr (length L_total), full_factor (length L_total), inv_array (length L_total) ]
            where L_total = 2 * Lq; the first Lq entries correspond to b particles, and the next Lq to c particles.

    Algorithm:
      For each site i (handled separately for b and c parts):
        - Get the mapped original position: inv = inv_array[i]
        - Extract the digit from the original state: digit = (state // full_factor[inv]) % sps_arr[inv]
        - Place the digit in the i-th position of the new state (weighted by full_factor[i])
    """
    # Read the total number of sites
    L_total = args[0]
    # Define offsets for the segments:
    offset1 = 1                # sps_arr starts at args[1], length L_total
    offset2 = offset1 + L_total  # full_factor starts at args[offset2], length L_total
    offset3 = offset2 + L_total  # inv_array starts at args[offset3], length L_total

    L_total_int = int(L_total)
    # For this model, the system is divided into two parts, each of length Lq
    Lq = L_total_int // 2

    new_state = 0

    # Process the b particle part: indices 0 to Lq-1
    for i in range(Lq):
        inv = args[offset3 + i]
        base = args[offset2 + int(inv)]
        sps_val = args[offset1 + int(inv)]
        digit = (state // base) % sps_val
        new_state += digit * args[offset2 + i]

    # Process the c particle part: indices Lq to L_total_int-1
    for i in range(Lq, L_total_int):
        inv = args[offset3 + i]
        base = args[offset2 + int(inv)]
        sps_val = args[offset1 + int(inv)]
        digit = (state // base) % sps_val
        new_state += digit * args[offset2 + i]

    return new_state

# --- next_state function ---
@cfunc(next_state_sig_32)
def next_state(s, counter, N, args):
    """
    Returns the next state in the external basis (stored in args).
    """
    return args[counter + 1]

#############################################
# 4. Build user_basis Wrapper Function for Two Species
#############################################
def build_two_species_user_basis(Lx, Ly, NE_b, NE_c, kx, ky):
    """
    Construct a user_basis for two species (b and c) with fixed particle numbers (NE_b, NE_c),
    while block-diagonalizing the basis using lattice translation symmetry (kx, ky).
    """
    Lq = Lx * Ly
    L = 2 * Lq
    # Generate the external basis, local dimension array, and full_factor array
    valid_states, sps_arr, full_factor = generate_two_species_basis(Lx, Ly, NE_b, NE_c)

    # Build the translation mapping: first, construct the mapping for a single-layer 2D lattice
    def set_lattice(Lx, Ly):
        site_list = np.zeros((Lx * Ly, 2), dtype=np.uint32)
        inv_site_list = np.zeros((Lx, Ly), dtype=np.uint32)
        nc = 0
        for ny in range(Ly):
            for nx in range(Lx):
                site_list[nc, :] = [nx + 1, ny + 1]
                inv_site_list[nx, ny] = nc
                nc += 1
        return site_list, inv_site_list

    def set_translation(Lx, Ly):
        site_list, inv_site_list = set_lattice(Lx, Ly)
        Tx = np.empty(Lx * Ly, dtype=np.uint32)
        Ty = np.empty(Lx * Ly, dtype=np.uint32)
        for ny in range(Ly):
            for nx in range(Lx):
                site = inv_site_list[nx, ny]
                Tx[site] = inv_site_list[(nx + 1) % Lx, ny]
                Ty[site] = inv_site_list[nx, (ny + 1) % Ly]
        return Tx, Ty

    Tx_b, Ty_b = set_translation(Lx, Ly)
    # For the c layer, add an offset of Lq
    Tx_full = np.concatenate([Tx_b, Tx_b + Lq])
    Ty_full = np.concatenate([Ty_b, Ty_b + Lq])
    # Compute the inverse mapping
    invTx = np.empty(L, dtype=np.uint32)
    invTy = np.empty(L, dtype=np.uint32)
    for i in range(L):
        invTx[Tx_full[i]] = i
        invTy[Ty_full[i]] = i

    # Package the mapping parameters:
    # Packaging order: [L, sps_arr (length L), full_factor (length L), invTx] or [L, sps_arr, full_factor, invTy]
    args_x = np.concatenate([np.array([L], dtype=np.uint32), sps_arr, full_factor, invTx])
    args_y = np.concatenate([np.array([L], dtype=np.uint32), sps_arr, full_factor, invTy])

    # Build the pcon dictionary (using the external basis next_state interface)
    class function_wrapper:
        def __init__(self, basis_arr):
            self.basis = basis_arr
        def get_s0_pcon(self, N, Np):
            return self.basis[0]
        def get_Ns_pcon(self, N, Np):
            return self.basis.size
    FW = function_wrapper(valid_states)
    pcon_dict = dict(
        Np=(),
        next_state=next_state,
        next_state_args=valid_states,
        get_Ns_pcon=FW.get_Ns_pcon,
        get_s0_pcon=FW.get_s0_pcon,
    )

    # op_args: package the local dimension array and full_factor for the operator function (length 2 * L)
    op_args = np.concatenate([sps_arr, full_factor])

    # Construct the user_basis object (note that data type, number of sites N, and op_dict must be provided)
    basis = user_basis(
        np.uint32,
        L,
        op_dict=dict(op=op_two_species, op_args=op_args),
        allowed_ops=set("+-n"),
        sps=np.max(sps_arr),  # Use the maximum local dimension as a dummy value instead of 0
        pcon_dict=pcon_dict,
        Tx_block=(translation_map_2species, Lx, kx, args_x),
        Ty_block=(translation_map_2species, Ly, ky, args_y)
    )
    return basis

#############################################
# 5. Construct Hamiltonian for Two Species (Hopping and Interaction Terms)
#############################################
def set_hopping_two_species(Lx, Ly, t):
    """
    Construct nearest-neighbor hopping terms for both species on a two-dimensional triangular lattice.
    Note: b layer site indices are 0 to Lq-1, and c layer indices are Lq to 2*Lq-1.
    """
    Lq = Lx * Ly
    hopping = []
    def set_lattice(Lx, Ly):
        site_list = np.zeros((Lx * Ly, 2), dtype=np.uint32)
        inv_site_list = np.zeros((Lx, Ly), dtype=np.uint32)
        nc = 0
        for ny in range(Ly):
            for nx in range(Lx):
                site_list[nc, :] = [nx + 1, ny + 1]
                inv_site_list[nx, ny] = nc
                nc += 1
        return site_list, inv_site_list
    site_list, inv_site_list = set_lattice(Lx, Ly)
    # Hopping for the b layer
    for ny in range(Ly):
        for nx in range(Lx):
            i = inv_site_list[nx, ny]
            n1 = inv_site_list[(nx + 1) % Lx, ny]
            n2 = inv_site_list[nx, (ny + 1) % Ly]
            n3 = inv_site_list[(nx - 1) % Lx, (ny + 1) % Ly]
            n4 = inv_site_list[(nx - 1) % Lx, ny]
            n5 = inv_site_list[nx, (ny - 1) % Ly]
            n6 = inv_site_list[(nx + 1) % Lx, (ny - 1) % Ly]
            hopping += [[t, i, n1], [t, i, n2], [t, i, n3],
                        [t, i, n4], [t, i, n5], [t, i, n6]]
    # Hopping for the c layer (with offset Lq)
    for ny in range(Ly):
        for nx in range(Lx):
            i = inv_site_list[nx, ny] + Lq
            n1 = inv_site_list[(nx + 1) % Lx, ny] + Lq
            n2 = inv_site_list[nx, (ny + 1) % Ly] + Lq
            n3 = inv_site_list[(nx - 1) % Lx, (ny + 1) % Ly] + Lq
            n4 = inv_site_list[(nx - 1) % Lx, ny] + Lq
            n5 = inv_site_list[nx, (ny - 1) % Ly] + Lq
            n6 = inv_site_list[(nx + 1) % Lx, (ny - 1) % Ly] + Lq
            hopping += [[t, i, n1], [t, i, n2], [t, i, n3],
                        [t, i, n4], [t, i, n5], [t, i, n6]]
    return hopping

def get_Hamiltonian_static_two_species(basis, Lx, Ly, t, U1, U2):
    """
    Construct the Hamiltonian including:
      - Hopping terms for both b and c particles.
      - Cross terms: for each lattice site i, add 2*(U1 - U2)*n_b*n_c (with b and c corresponding to sites i and i + Lq respectively)
      - Self-energy terms: for all sites, add (U1 + U2)*n_i^2
    """
    Lq = Lx * Ly
    hopping = set_hopping_two_species(Lx, Ly, t)
    interaction = [[2 * (U1 - U2), i, i + Lq] for i in range(Lq)]
    self_energy = [[(U1 + U2), i, i] for i in range(2 * Lq)]
    static = [
        ["+-", hopping],
        ["nn", interaction],
        ["nn", self_energy],
    ]
    H = hamiltonian(static, [], basis=basis,
                    check_symm=False, check_herm=False, check_pcon=False)
    return H

def get_kinetic_two_species(basis, Lx, Ly, t):
    hopping = set_hopping_two_species(Lx, Ly, t)
    static = [
        ["+-", hopping],
    ]
    Op_kinetic = hamiltonian(static, [], basis=basis,
                             check_symm=False, check_herm=False, check_pcon=False)
    return Op_kinetic

def grand_ED_full(H, t, U1, U2):
    E, V = H.eigh()
    return E, V

def grand_fullED_observe(Op, E, V, beta, mu, NE):
    W = np.exp(-(E-mu*NE)*beta)
    O = Op.matrix_ele(V, V, diagonal=True)
    OZ = np.sum(W*O)
    Z  = np.sum(W)
    return OZ, Z

def grand_at_NE_block_2flavor(Lx, Ly, t, U1, U2, beta, mu, NE_b, NE_c, kx, ky):
    NE = NE_b + NE_c
    basis = build_two_species_user_basis(Lx, Ly, NE_b, NE_c, kx, ky)
    logging.info(f'        kx={kx}, ky={ky}, basis.Ns={basis.Ns}')
    H = get_Hamiltonian_static_two_species(basis, Lx, Ly, t, U1, U2)
    Op_kinetic = get_kinetic_two_species(basis, Lx, Ly, t)
    E, V = grand_ED_full(H, t, U1, U2)
    OZ, Z = grand_fullED_observe(Op_kinetic, E, V, beta, mu, NE)
    return np.real(OZ), Z

def grand_at_NE(Lx, Ly, t, U1, U2, beta, mu, NE_b, NE_c):
    OZ = 0.0
    Z = 0.0
    for kx in range(Lx):
        for ky in range(Ly):
            OZ_block, Z_block = grand_at_NE_block_2flavor(Lx, Ly, t, U1, U2, beta, mu, NE_b, NE_c, kx, ky)
            OZ += OZ_block
            Z  += Z_block
    return OZ, Z

def main():
    logging.info('######################### EDtriangle_symm BEGIN !!!!')

    # Read parameters from file
    logging.info(f'... Read parameters from file')
    try:
        params = read_params()
        Lx, Ly = int(params["Lx"]), int(params["Ly"])
        t, U1, U2 = params["t"], params["U1"], params["U2"]
        beta, mu = params["beta"], params["mu"]
    except Exception as e:
        logging.error(f"Parameter reading failed: {e}")
        return

    start_time = time.time()

    Ztot, NE_Ztot, kinetic_Ztot = 0.0, 0.0, 0.0


    logging.info(f'... Start')

    # NE = 0
    logging.info(f' ')
    logging.info(f'NE = {0} \t NE_b = {0} \t NE_c = {0}')
    Ztot += 1.0
    NE_Ztot += 0.0
    kinetic_Ztot += 0.0
    logging.info(f'==> KE * weight = {0} \t weight = {1} \t KE(NE) = {0}')

    # NE > 0
    for NE in range(1, Lx*Ly*Lx*Ly+1):
        OZ = np.zeros([NE+1,NE+1])
        Z  = np.zeros([NE+1,NE+1])
        for NE_b in range(NE+1):
            NE_c = NE - NE_b
            logging.info(f' ')
            logging.info(f'NE = {NE} \t NE_b = {NE_b} \t NE_c = {NE_c}')
            if NE_b <= NE_c:
                OZ[NE_b,NE_c], Z[NE_b,NE_c] = grand_at_NE(Lx, Ly, t, U1, U2, beta, mu, NE_b, NE_c)
            else:
                OZ[NE_b,NE_c], Z[NE_b,NE_c] = OZ[NE_c,NE_b], Z[NE_c,NE_b]
            logging.info(f'==> KE * weight = {OZ[NE_b,NE_c]}, weight = {Z[NE_b,NE_c]}, KE(NE) = {OZ[NE_b,NE_c]/Z[NE_b,NE_c]}')

        OZ_NE = np.sum(OZ)
        Z_NE  = np.sum(Z)
        Ztot    += Z_NE
        NE_Ztot += NE * Z_NE
        kinetic_Ztot += OZ_NE

        res_NE = NE_Ztot / Ztot
        res_kinetic = kinetic_Ztot / Ztot

        err_NE = NE * Z_NE / Ztot
        err_kinetic = OZ_NE / Ztot

        # Write temp results to results.txt
        with open("results.txt", "w") as result_file:
            result_file.write(f"{res_NE}\n")
            result_file.write(f"{res_kinetic}\n")
            result_file.write(f"{err_NE}\n")          # as temp error bar
            result_file.write(f"{err_kinetic}\n")     # as temp error bar
            logging.info(f' ')
            logging.info(f'>>.............. concluding at NE = {NE}')
            logging.info(f'>>.............. temporary results: <NE> = {res_NE}')
            logging.info(f'>>.............. temporary results: <NE> +- {err_NE}')
            logging.info(f'>>.............. temporary results: <KE> = {res_kinetic}')
            logging.info(f'>>.............. temporary results: <KE> +- {err_kinetic}')

        if abs(OZ_NE) < abs(0.001 * kinetic_Ztot):
            break

    res_NE = NE_Ztot / Ztot
    res_kinetic = kinetic_Ztot / Ztot

    err_NE = NE * Z_NE / Ztot
    err_kinetic = OZ_NE / Ztot

    # Write final results to results.txt
    with open("results.txt", "w") as result_file:
        result_file.write(f"{res_NE}\n")
        result_file.write(f"{res_kinetic}\n")
        result_file.write(f"{err_NE}\n")          # as temp error bar
        result_file.write(f"{err_kinetic}\n")     # as temp error bar

    execution_time = time.time() - start_time
    logging.info(f"Execution completed in {execution_time:.2f} seconds.")
    logging.info(f'Final results: <NE> = {res_NE}')
    logging.info(f'Final results: <NE> +- {err_NE}')
    logging.info(f'Final results: <KE> = {res_kinetic}')
    logging.info(f'Final results: <KE> +- {err_kinetic}')
    logging.info(f'\n\n')

if __name__ == "__main__":
    main()
