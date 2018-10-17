# Jackson Parsons
# Applied Cryptography, Spring 2018
# April 27th, 2018
import numpy
import os.path
import re
from ast import literal_eval

# Calculate round constants:
# 7 constants per round for 24 rounds for a total of 168 constants.
rc = [1, 0, 0, 0, 0, 0, 0, 0]
while len(rc) < 168:
    end = len(rc)
    rc.append(rc[end - 8] ^ rc[end - 4] ^ rc[end - 3] ^ rc[end - 2])


# SHA3-1
def expand(flat_list):
    """ Expand a 1d list into a SHA3-style 3d list according to the algorithm
        3-dim[i][j][k] = 1-dim[64(5j+i)+k].
    Params:
        flat_list (list)
    Returns:
        keccak (list of list of lists): A three-dimensional list whose values
            are mapped from flat_list
    """
    # Assert that the one-dimensional list has 1600 entries.
    assert len(flat_list) == 1600, "the message is not 1600"

    # Instantiate a blank 5x5x64 cube with -1 in all cells.
    # Then, scan each cubecell and follow the map to pull a value from the list.
    keccak = [[[-1 for k in range(64)] for j in range(5)] for i in range(5)]
    for i in range(5):  # 0 <= i <= 4
        for j in range(5):  # 0 <= j <= 4
            for k in range(64):  # 0 <= k <= 63
                keccak[i][j][k] = flat_list[64 * (5 * j + i) + k]
    return keccak


# SHA3-2
def collapse(keccak):
    """ Collapse a SHA3-style 3d list to 1d list according to the algorithm
        1-dim[64(5j+i)+k] = 3-dim[i][j][k].
    Params:
        keccak (list of list of lists): a 5 x 5 x 64 keccak
    Returns:
        flat_list (list): A list whose values are mapped from keccak
    """
    # Assert that the cube is 5x5x64.
    assert len(keccak) == 5, "the height of the 3d array is not always 5"
    for row in keccak:
        assert len(row) == 5, "the width of the array is not always is not 5"
        for drawer in row:
            assert len(drawer) == 64, "the depth of the array is not always 64"

    # Declare a blank, 1600-entry list to hold the values from the cube.
    # Then, scan each cubecell and follow the map to send a value to the list.
    flat_list = [-1] * 1600
    for i in range(5):
        for j in range(5):
            for k in range(64):
                flat_list[64 * (5 * j + i) + k] = keccak[i][j][k]
    return flat_list


# SHA3-3
def theta(theta_in):
    """ Apply theta diffusion to each entry in the SHA3 cube.
    Params:
        theta_in (list of list of lists): 5 x 5 x 64 keccak
    Returns:
        theta_out (list of list of lists): the post-diffusion keccak
    """
    # Assert that the cube is 5x5x65.
    # This way inner functions like rho and pi don't have to check.
    assert len(theta_in) == 5, "the height of the 3d array is not always 5"
    for row in theta_in:
        assert len(row) == 5, "the width of the array is not always is not 5"
        for drawer in row:
            assert len(drawer) == 64, "the depth of the array is not always 64"

    # Instantiate the return cube, then fill it with values from theta_in.
    theta_out = [[[-1 for k in range(64)] for j in range(5)] for i in range(5)]
    for i in range(5):
        for j in range(5):
            for k in range(64):
                interm_val = theta_in[i][j][k]
                # Diffuse with North row. Wrap i = -1 to 4.
                for j_prime in range(5):
                    interm_val = interm_val ^ theta_in[(i - 1) % 5][j_prime][k]
                # Diffuse with South-West row. Wrap i = 5 to 0.
                for j_prime in range(5):
                    interm_val = interm_val ^ theta_in[(i + 1) % 5][j_prime][(k - 1) % 64]
                theta_out[i][j][k] = interm_val
    return theta_out


# SHA3-4
def rho(rho_in):
    """ Apply rho diffusion to each 'drawer' (i, j independent, k dependent) of
        entries in the SHA3 cube after theta diffusion.
    Params: rho_in
        rho_in (list of lists of lists): 5 x 5 x 64 keccak
    Returns:
        rho_out (list of lists of lists): the post-diffusion keccak
    """
    # Instantiate the [0 1\n 2 3] square matrix and the [1\n 2] column matrix.
    # Default int type in Windows is in32, so int64 is explicitly specified.
    special_matrix = numpy.matrix([[0, 1], [2, 3]], dtype=numpy.int64)
    col_matrix = numpy.matrix([[1], [0]])

    # Instantiate the return cube and fill it with values from rho_in.
    rho_out = [[[-1 for k in range(64)] for j in range(5)] for i in range(5)]
    for t in range(24):  # Working in mod 5 space here
        ij_matrix = (special_matrix ** t % 5) * col_matrix
        i = int(ij_matrix[0, 0])
        j = int(ij_matrix[1, 0])
        for k in range(64):  # Working in mod64 space here.
            rho_out[i][j][k] = rho_in[i][j][(int(k - (t + 1) * (t + 2) / 2)) % 64]
    for k in range(64):
        rho_out[0][0][k] = rho_in[0][0][k]
    return rho_out


# SHA3-5
def pi(pi_in):
    """ Apply pi diffusion to each slice (k-dependent with i and j free to move)
        of entries in the SHA3 cube.
    Params:
        pi_in (list of lists of lists): 5 x 5 x 64 keccak
    Returns:
        pi_out (list of lists of lists): the post-diffusion keccak
    """
    # Instantiate the return cube and fill it with values from pi_in.
    pi_out = [[[-1 for k in range(64)] for j in range(5)] for i in range(5)]
    for i in range(5):
        for j in range(5):
            for k in range(64):
                pi_out[j][(2 * i + 3 * j) % 5][k] = pi_in[i][j][k]
    return pi_out


# SHA3-6
def chi(chi_in):
    """ Apply chi diffusion to each cell in pi_out.
    Params:
        chi_in (list of lists of lists): a 5 x 5 x 64 keccak
    Returns:
        chi_out (list of lists of lists): the post-diffusion keccak
    """
    # Instantiate the return cube and fill it with values from chi_in.
    chi_out = [[[-1 for k in range(64)] for j in range(5)] for i in range(5)]
    for i in range(5):
        for j in range(5):
            for k in range(64):
                chi_out[i][j][k] = chi_in[i][j][k] ^ ((chi_in[(i + 1) % 5][j][k] ^ 1) * chi_in[(i + 2) % 5][j][k])
    return chi_out


# SHA3-7
def iota(iota_in, constants, round_num):
    """ Apply iota diffusion to each cell in pi_out.
    Params:
        iota_in (list of lists of lists): 5 x 5 x 64 keccak
        constants (list of bits): A list containing the first 24 constants
            from x^round in the finite field F2[x] / (x^8 + x^6 + x^5 + x^4 + 1)
        round_num (int): The current round number.
    Returns:
        iota_out (list of lists of lists): the post-diffusion keccak
    """
    # Fill one drawer in a 5x5x65 cube round constants. The other cells are 0.
    bit = [[[0 for k in range(64)] for j in range(5)] for i in range(5)]
    for l in range(7):
        bit[0][0][2 ** l - 1] = constants[l + 7 * round_num]

    # Instantiate the return cube and fill it with values from iota_in ^ bit.
    iota_out = [[[-1 for k in range(64)] for j in range(5)] for i in range(5)]
    for i in range(5):
        for j in range(5):
            for k in range(64):
                iota_out[i][j][k] = iota_in[i][j][k] ^ bit[i][j][k]
    return iota_out


# SHA3-8
def sha3(msg):
    """ Wrapper function for the SHA3 algorithm.
    Params:
        msg (list): A list of bits (0 or 1) of arbitrary length. For simplification
        purposes, the length of msg is limited to a max of 1086.
    Returns: The 256-bit hash of the message.
    """
    num_bits = len(msg)
    assert num_bits <= 1086, "For simplification, please limit to 1086 entries."
    for msg_bit in msg:
        assert (msg_bit == 0 or 1), "Your message can only contain bits."

    # Create a pad whose length is the difference between the number of bits in
    # the message and 1086. Set the beginning and end entries of that pad to 1.
    pad = [0] * (1088 - num_bits)
    pad[0] = pad[-1] = 1
    for pad_bit in pad:
        msg.append(pad_bit)
    # Append 512 0s to the padded message - the new message is now 1600 bits.
    while len(msg) < 1600:
        msg.append(0)

    # Convert the msg to a keccak and apply 24 rounds of theta, rho, pi, chi,
    # and iota in sequence.
    keccak = expand(msg)
    for t in range(24):
        keccak = iota(chi(pi(rho(theta(keccak)))), rc, t)
    # Collapse the final keccak to a list and return the first 256 bits.
    return collapse(keccak)[0:256]


# Driver code
# Open SHA3inout.txt (assuming it's at root) and pull arrays for testing.
rootpath = os.path.expanduser("~")
filepath = os.path.join(rootpath, "SHA3inout.txt")
with open(filepath, 'r') as file:
    for line in file:
        if re.match("randomvec", line):
            randomvec = literal_eval(line.split("=")[1])
        elif re.match("thetaout", line):
            thetaout = literal_eval(line.split("=")[1])
        elif re.match("v3box", line):
            v3box = literal_eval(line.split("=")[1])
        elif re.match("rhoout", line):
            rhoout = literal_eval(line.split("=")[1])
        elif re.match("v4box", line):
            v4box = literal_eval(line.split("=")[1])
        elif re.match("piout", line):
            piout = literal_eval(line.split("=")[1])
        elif re.match("v5box", line):
            v5box = literal_eval(line.split("=")[1])
        elif re.match("chiout", line):
            chiout = literal_eval(line.split("=")[1])
        elif re.match("v6box", line):
            v6box = literal_eval(line.split("=")[1])
        elif re.match("iotainround20", line):
            iotainround20 = literal_eval(line.split(" = ")[1])
        elif re.match("iotaoutround20", line):
            iotaoutround20 = literal_eval(line.split(" = ")[1])
        elif re.match("SHA3in", line):
            sha3in = literal_eval(line.split("=")[1])
        elif re.match("SHA3out", line):
            sha3out = literal_eval(line.split("=")[1])
    rhoin = thetaout
    piin = rhoout
    chiin = piout

    # Function testing:
    # SHA3-4 Output:
    print("SHA3-4")
    if rho(rhoin) == rhoout:
        print(" rhoin is equal to rhooout!")
    print(" Bits 11-20 from rho(v3box):", rho(v3box)[2][3][11:21])
    # SHA3-5 Output:
    print("SHA3-5")
    if pi(piin) == piout:
        print(" piin is equal to piout!")
    print(" Bits 11-20 from pi(v4box):", pi(v4box)[2][3][11:21])
    # SHA3-6 Output:
    print("SHA3-6")
    if chi(chiin) == chiout:
        print(" chiin is equal to chiout!")
    print(" Bits 11-20 from chi(v5box):", chi(v5box)[2][3][11:21])
    # SHA3-7 Output:
    print("SHA3-7")
    if iota(iotainround20, rc, 20) == iotaoutround20:
        print(" iotainround20 is equal to iotaoutround20!")
    print(" Bits 1-10 from iota(v6box) at round 13:", iota(v6box, rc, 13)[2][3][11:21])
    # SHA3-8 Output:
    print("SHA3-8")
    if sha3(sha3in) == sha3out:
        print(" SHA3in is equal to SHA3out!")
    print(" Bits 0-15 of sha3(empty_string): ", sha3(list())[0:15])