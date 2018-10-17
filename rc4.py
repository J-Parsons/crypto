# RC4 4-2 Script for Applied Cryptography Spring 2018
# Jackson Parsons
# April 13th, 2018

from math import ceil


def decimal_to_binary(num_bits, decimal_value):
    """ Convert a decimal value to binary using the divide by 2 algorithm.
    Args:
        num_bits (int): The number of elements allowed in the list
            representation of the binary value to be returned.
            i.e. the number of bits in binary representation
        decimal_value (int): The decimal value to be converted to binary.
    Returns:
        bin_stack (list): A list of binary digits representing the binary value
            of decimal_value. The most significant / leftmost binary digit is
            the first element in the list. The least significant / rightmost
            binary digit is the last element in the list.
    """
    assert (type(num_bits) == int and num_bits > 0), \
        "num_bits must be a positive integer"
    assert (type(decimal_value) == int and decimal_value > 0), \
        "decimal_value must be a positive integer"
    bin_stack = []

    # Store the remainders of division by 2, updating decimal_value to the floor
    # of division by 2 afterwards.
    # I.e. perform division by 2, store the remainder, then update decimal_value
    while decimal_value > 0:
        remainder = decimal_value % 2
        bin_stack.append(remainder)
        decimal_value = decimal_value // 2

    # If the number of bits that you need at minimum to represent decimal_value
    # is different than the number of bits passed to decimal_to_binary,
    if len(bin_stack) > num_bits:  # if too many bits, truncate most significant
        while len(bin_stack) > num_bits:
            bin_stack.pop()
    elif len(bin_stack) < num_bits:  # if too few, sign-extend with 0s
        while len(bin_stack) < num_bits:
            bin_stack.append(0)

    # At this point, the list is represented from low to high, so reverse it so
    # that the most significant digits are written from left to right.
    bin_stack.reverse()
    return bin_stack


def binary_to_decimal(bit_array):  # Subroutine B
    """ Convert a list representing a binary value into a decimal value.
    Args:
        bit_array:  (list): A list of 1s and 0s representing a binary value.
            This uses the same formatting used by decimal_to_binary() function.
            The first element in the list is the most significant/leftmost
            binary digit. The last element in the list is the least significant/
            rightmost binary digit.
    Returns:
        decimal_value (int): A decimal integer whose value corresponds to the
            value represented by bin_array.
    """
    assert (len(bit_array) > 0), "binary_to_decimal: bit_array can't be empty"
    for digit in bit_array:
        assert (digit == 0 or digit == 1), \
            "binary_to_decimal: digits in bit_array must be binary"

    decimal_value = 0
    for i, val in enumerate(reversed(bit_array)):  # iterate backwards
        decimal_value += val * (2 ** i)  # bit_array[i] * i^2

    return decimal_value


def rc4(n, l, bit_array):
    """ Generate a string of pseudo-random bits using the RC4 algorithm.
    Args:
        n (int): A positive integer representing the number of bits that will be
            used to represent each binary digit in the key as a decimal values.
        l (int): The ceiling of the length of the BINARY ciphertext.
        bit_array (list): A list of binary integers that represents the key
            (a list of positive decimal integers)
    Returns:
        keystream (list): A list of pseudo-random decimal values.
    """
    assert (n > 0 and type(n) == int), "rc4: n must be a positive integer"
    assert (l > 0 and type(l) == int), "rc4: l must be a positive integer"
    for digit in bit_array:
        assert (digit == 0 or digit == 1), "rc4: digits in key must be binary"

    key = []
    for i in range(0, len(bit_array), n):
        num = binary_to_decimal(bit_array[i:(i + n)])  # decimalize n bits
        key.append(num)

    # Add copies of key entries until the key has 2^n elements.
    # Typically you would need to i to 0 when i == 2^n, but since the key is
    # always increasing in length with same-order values, we can just let i
    # increment without modding it.
    i = 0
    while len(key) != (2 ** n):
        key.append(key[i])
        i += 1

    # Initialize the state to 1, 2... 4
    state = []
    for i in range(2 ** n):  # i from 0 to 2^n - 1
        state.append(i)

    j = 0
    for i in range(2 ** n):
        j = ((j + state[i] + key[i]) % (2 ** n))  # j = j + Si + Ki
        state[i], state[j] = state[j], state[i]  # swap state[i] and state[j]

    key_stream = []
    i = j = 0
    t = 0
    for r in range(0, l):
        i = (i + 1) % (2 ** n)  # i = (i + 1) mod 2^n
        j = (j + state[i]) % (2 ** n)  # j = (j + Si) mod 2^n
        state[i], state[j] = state[j], state[i]
        t = (state[i] + state[j]) % (2 ** n)  # t = (Si + Sj) mod 2^n
        key_stream.append(state[t])

    return key_stream


# Driver code
# Part i, test binary_to_decimal aka SubroutineA
print("Part i test 1:\n", decimal_to_binary(32, 2181038082), "\n")

# Part i, test decimal_to_binary aka SubroutineB
print("Part i test 2:\n", binary_to_decimal(
    [1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 1, 0]), "\n")

# Part ii, RC4-1 comparison test
n = 3
plaintext = [1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0,
             1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0]
l = ceil(len(plaintext) / n)
shared_key = [1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1]
print("Part ii RC401 comparison test:\n", rc4(n, l, shared_key), "\n")

# Part iii, test 1
n = 8
l = 28
shared_key = [1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0,
              1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0,
              1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0,
              0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1]
decimal_key_stream = rc4(n, l, shared_key)
print("Part iii test 1:\n", decimal_key_stream, "\n")

# Part iii, test 2
ciphertext = [1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0,
              1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0,
              0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0,
              0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0,
              1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0,
              0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0,
              1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1,
              1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1,
              0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1,
              1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0,
              0, 1, 0, 0]

binary_key_stream = []
for entry in decimal_key_stream:
    binary_value = decimal_to_binary(n, entry)
    binary_key_stream.append(binary_value)
# Flatten the keystream, which was a list of lists of binary digits, into a
# single, flat list of binary digits.
binary_key_stream = [item for sublist in binary_key_stream for item in sublist]

# XOR the ciphertext and the keystream to get the plaintext in binary.
plaintext = []
for i in range(len(ciphertext)):
    plaintext.append(ciphertext[i] ^ binary_key_stream[i])  # bitwise xor

# Convert every byte in the plaintext to an ASCII character
ascii_plaintext = []
for i in range(0, len(plaintext), n):
    decimal = binary_to_decimal(plaintext[i:(i + n)])  # a pt byte in decimal
    ascii_plaintext.append(chr(decimal))  # convert 0-255 to ascii

print("Part iii test 2:\n", ''.join(ascii_plaintext), "\n")