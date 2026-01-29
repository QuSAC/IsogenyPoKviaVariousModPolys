# This file computes the number of constraints, variables and NNZ for the
# constraint systems given in the paper.
#
# The structure is as follows: variables are defined in a dictionary, which
# also stores how many variables are necessary to express it. For Fp^2, this
# will always be one, for Fp we keep track of how efficient it is to express
# the real part, the imaginary part, and their sum. We then use this knowledge
# when asserting equalities.
#
# The functions `assert_product` and `assert_square` are used every time a
# constraint is added. Over Fp^2, this translates to a row of the R1CS, for Fp
# this is then lifted to 3 and 2 rows, respectively.

var('k X j0 j1')

NO_VARS_USED = 'no_vars_used'
VARS = 'variables'

def phi_coefs_for_ell(ell):
    s = 12/gcd(12, ell - 1)
    phi = DedekindEtaModularPolynomialDatabase()[ell]

    var('x j')

    # Compute the coefficients
    coefs = []
    for [coef, power] in (phi(x=x,j=j) + j*x).coefficients():
        coefs.append(coef)

    # The alternative coefficients
    coefsp = []
    for i in range(0, ell + 2):
        coefsp.append(coefs[ell + 1 - i] * ell^(s * (ell - i)))
    return phi, s, coefs, coefsp

# Returns
# - The number of constraints used to compute the product
# - The number of variables used to express the product
# - The number of non-zero entries in the constraint matrices
def assert_product(in_fp, left_input, right_input, output, expected_output = 0):
    if in_fp:
        constraints = 3
        # We just have u2 as additional variable
        no_variables = 1

        # The number of non-zero entries in the matrix
        [x1, x2, sum_x] = left_input['no_vars_used']
        [y1, y2, sum_y] = right_input['no_vars_used']
        [z1, z2, sum_z] = output['no_vars_used']
        # u2 = x2 * y2
        non_zero_entries = 1 + x2 + y2
        # z1 - d*u2 = x1 * y1
        non_zero_entries += z1 + 1 + x1 + y1
        # zs + (1-d)u2 = xs + ys
        non_zero_entries += sum_z + 1 + sum_x + sum_y
    else:
        constraints = 1
        no_variables = 0

        # The number of non-zero entries in the matrix
        non_zero_entries = left_input['no_vars_used'] + right_input['no_vars_used'] + output['no_vars_used']

    # Check that the constraint is valid
    diff = left_input['vars'] * right_input['vars'] - output['vars'] - expected_output
    if diff != 0:
        print(expand(diff))
    assert diff == 0
    
    return [constraints, no_variables, non_zero_entries]

# Returns
# - The number of constraints used to compute the square
# - The number of variables used to express the square
#   of the output
# - The number of non-zero entries in the constraint matrices
def assert_square(in_fp, inp, outp, expected_output = 0):
    if in_fp:
        constraints = 2
        no_variables = 0

        # The number of non-zero entries in the matrix
        [x1, x2, sum_x] = inp['no_vars_used']
        [z1, z2, sum_z] = outp['no_vars_used']
        # z2 = 2 x1 * x2
        non_zero_entries = z2 + x1 + x2
        # z1 + (d+1)/2 z2 = sum_x * (x1 + d*x2)
        non_zero_entries += z1 + min(z2, sum_z) + sum_x + x1 + min(x2, sum_x)
    else:
        constraints = 1
        no_variables = 0

        # The number of non-zero entries in the matrix
        non_zero_entries = inp['no_vars_used'] + inp['no_vars_used'] + outp['no_vars_used']

    # Check that the constraint is valid
    diff = inp['vars'] * inp['vars'] - outp['vars'] - expected_output
    if diff != 0:
        print(expand(diff))
    assert diff == 0

    return [constraints, no_variables, non_zero_entries]

def scale_var(in_fp, coef, entry):
    if in_fp:
        return {
            'vars': entry['vars'] * coef,
            'no_vars_used': [x for x in entry['no_vars_used']]
        }
    else:
        return {
            'vars': entry['vars'] * coef,
            'no_vars_used': entry['no_vars_used']
        }

# Sum variable dictionaries
def sum_vars(in_fp, coefs_and_entries):
    if in_fp:
        result = {
            'vars': 0,
            'no_vars_used': [0, 0, 0]
        }
    else:
        result = {
            'vars': 0,
            'no_vars_used': 0
        }
    for [coef, entry] in coefs_and_entries:
        result['vars'] += coef * entry['vars']

        if in_fp:
            result['no_vars_used'] = [sum(x) for x in zip(result['no_vars_used'], entry['no_vars_used'])]
        else:
            result['no_vars_used'] += entry['no_vars_used']

    return result

basis_regular = [1, 1, 2]
basis_sum = [1, 2, 1]

VARIABLE_ONE = { 'vars': 1, 'no_vars_used': [1, 0, 1] }

def original_paper(in_fp):
    phi = ClassicalModularPolynomialDatabase()[2]
    constraints = 0
    non_zero_entries = 0

    # We don't count the 1 in the variable count.
    # We duplicate j0 and j1 in the dictionary, but this is merely for
    # asserting correctness of the constraints.
    if in_fp:
        variables = 2 * (k + 1) \
            + 2 * (k + 1) \
            + 2 * (k + 1) \
            + 2 * k
        var_dict = {
            '1': VARIABLE_ONE,
            'j0^1': { 'vars': j0, 'no_vars_used': basis_regular },
            'j0^2': { 'vars': j0^2, 'no_vars_used': basis_regular },
            'j0^3': { 'vars': j0^3, 'no_vars_used': basis_regular },
            'j1^1': { 'vars': j1, 'no_vars_used': basis_regular },
            'j1^2': { 'vars': j1^2, 'no_vars_used': basis_regular },
            'j1^3': { 'vars': j1^3, 'no_vars_used': basis_regular },
            'j0j1': { 'vars': j0*j1, 'no_vars_used': basis_regular }
        }
    else:
        # We don't count the 1
        variables = k + 1 \
            + k + 1 \
            + k + 1 \
            + k
        var_dict = {
            '1': {'vars': 1, 'no_vars_used': 1 },
            'j0^1': { 'vars': j0, 'no_vars_used': 1 },
            'j0^2': { 'vars': j0^2, 'no_vars_used': 1 },
            'j0^3': { 'vars': j0^3, 'no_vars_used': 1 },
            'j1^1': { 'vars': j1, 'no_vars_used': 1 },
            'j1^2': { 'vars': j1^2, 'no_vars_used': 1 },
            'j1^3': { 'vars': j1^3, 'no_vars_used': 1 },
            'j0j1': { 'vars': j0*j1, 'no_vars_used': 1 }
        }

    [new_constr, new_vars, new_nnz] = assert_square(in_fp, var_dict['j0^1'], var_dict['j0^2'])
    constraints += new_constr * (k + 1)
    non_zero_entries += new_nnz * (k + 1)
    variables += new_vars * (k + 1)

    [new_constr, new_vars, new_nnz] = assert_product(in_fp, var_dict['j0^2'], var_dict['j0^1'], var_dict['j0^3'])
    constraints += new_constr * (k + 1)
    non_zero_entries += new_nnz * (k + 1)
    variables += new_vars * (k + 1)

    [new_constr, new_vars, new_nnz] = assert_product(in_fp, var_dict['j0^1'], var_dict['j1^1'], var_dict['j0j1'])
    constraints += new_constr * k
    non_zero_entries += new_nnz * k
    variables += new_vars * k

    # Compute polynomial
    left_hand_side = scale_var(in_fp, -1488, var_dict['j0j1'])
    right_hand_side = sum_vars(in_fp, [
        [1, var_dict['j0^1']],
        [1, var_dict['j1^1']],
        [-1/1488, var_dict['j0j1']],
    ])
    output = sum_vars(in_fp, [
        [1, var_dict['j0^3']],
        [1, var_dict['j1^3']],
        [-162000, var_dict['j0^2']],
        [-162000, var_dict['j1^2']],
        [8748000000, var_dict['j0^1']],
        [8748000000, var_dict['j1^1']],
        [40773375, var_dict['j0j1']],
        [-157464000000000, var_dict['1']],
    ])
    [new_constr, new_vars, new_nnz] = assert_product(in_fp, left_hand_side, right_hand_side, output, -phi(j0=j0,j1=j1))
    constraints += new_constr * k
    non_zero_entries += new_nnz * k
    variables += new_vars * k

    return constraints, variables, non_zero_entries


# Compute the number of constraints, variables and NNZ
def simple_approach(ell, in_fp):
    phi, s, coefs, coefsp = phi_coefs_for_ell(ell)

    constraints = 0
    non_zero_entries = 0
    if in_fp:
        # We don't count the 1
        variables = 2 * k + 2 * (k + 1)
        # Variables for real part, variables for imaginary part, variables for sum
        var_dict = {
            '1': VARIABLE_ONE,
            'X^1': { 'vars': X, 'no_vars_used': basis_sum },
            'y0': { 'vars': j0 - coefs[1], 'no_vars_used': basis_sum },
            'y1': { 'vars': j1 - coefs[1], 'no_vars_used': basis_sum }
        }
    else:
        # We don't count the 1
        variables = k + k + 1
        var_dict = {
            '1': {'vars': 1, 'no_vars_used': 1 },
            'X^1': { 'vars': X, 'no_vars_used': 1 },
            'y0': { 'vars': j0 - coefs[1], 'no_vars_used': 1 },
            'y1': { 'vars': j1 - coefs[1], 'no_vars_used': 1 }
        }
    
    # Compute the powers up to ell, inclusive.
    # We first add the powers. We can always obtain these by squaring, except
    # for the biggest power, if its odd. We only do this over Fp, since
    # squaring is only cheaper here and the nnz even becomes more expensive
    # over Fp^2.
    #
    # We have an additional trick to minimize the number of non zero entries.
    # When using the squaring trick to obtain odd powers, we make sure to
    # obtain a linear combination that is immediately useful. This saves a
    # non-zero entry in the squaring, as well as in the computation of one of
    # the polynomials. We need to make sure that the power is high enough that
    # it is not needed in plain to compute another power.
    # I.e. we replace the middle power by the linear combination
    # X^{2i} + X^{2i + 1} + X^{2i + 2}
    # if we do not need 2i + 2i + 1. Since we allways compute ell directly, it suffices that
    # 4i + 1 is bigger than the next odd number, so
    # 4i + 1 > 2 * floor(ell / 2) - 1
    # or, we can replace power j = 2i + 1 by a linear combination with j + 1 if
    # j > floor(ell / 2)
    for i in range(2, ell + 1):
        if in_fp:
            variables += 2 * k
            if i == 2 and i == ell: # Only in this specific case the regular basis is better
                var_dict['X^' + str(i)] = {\
                    'vars': X^i,\
                    'no_vars_used': basis_regular,\
                }
            elif i != ell and i % 2 == 1 and i > floor(ell / 2):
                # We use the linear combination directly here, where the second
                # two terms are already correct for the first polynomial.
                coef_a = coefs[i + 1]
                coef_b = coefs[i + 2]
                var_dict['xX^' + str(i - 1) + ' + cX^' + str(i) + ' + cX^' + str(i + 1)] = {\
                    'vars': coef_a^2/coef_b/4*X^(i-1) + coef_a*X^(i) + coef_b*X^(i + 1),\
                    'no_vars_used': basis_sum,\
                }
            else:
                var_dict['X^' + str(i)] = {\
                    'vars': X^i,\
                    'no_vars_used': basis_sum,\
                }
        else:
            variables += 1 * k
            var_dict['X^' + str(i)] = {'vars': X^i, 'no_vars_used': 1 }

    # Now actually compute these powers. Again, we use squaring whereever
    # possible.
    for i in range(2, ell + 1):
        if i % 2 == 0:
        # A simple squaring!
            compute_from = var_dict['X^' + str(i / 2)]
            output = var_dict['X^' + str(i)]
            [new_constraints, new_no_variables, new_non_zero_entries] = assert_square(in_fp, compute_from, output)
        elif i == ell or not in_fp:
            # We do not have a higher even power to subtract, so multiply directly
            left_hand_side = var_dict['X^1']
            right_hand_side = var_dict['X^' + str(i - 1)]
            output = var_dict['X^' + str(i)]
            [new_constraints, new_no_variables, new_non_zero_entries] = assert_product(in_fp, left_hand_side, right_hand_side, output)
        else:
            # We use the following trick:
            # (x^i + x^{i + 1})^2 = x^{2i} + 2x^{2i + 1} + x^{2i + 2}
            #
            # Now, if the power is high enough, we just use the result to
            # express two powers correctly in the polynomial. Otherwise, we
            # subtract of the wrong powers. In practice, there will be only one
            # power, namely X^3 for ell=5 and X^5 for ell=7 or ell=13.
            if i > floor(ell / 2):
                coef_a = coefs[i + 1]
                coef_b = coefs[i + 2]
                inp = sum_vars(in_fp, [\
                    [coef_a / coef_b / 2, var_dict['X^' + str(floor(i / 2))]],\
                    [1, var_dict['X^' + str(ceil(i / 2))]],\
                ])
                outputs = [\
                    [1/coef_b, var_dict['xX^' + str(i - 1) + ' + cX^' + str(i) + ' + cX^' + str(i + 1)]],\
                ]
            else:
                inp = sum_vars(in_fp, [\
                    [1, var_dict['X^' + str(floor(i / 2))]],\
                    [1, var_dict['X^' + str(ceil(i / 2))]],\
                ])
                outputs = [\
                    [1, var_dict['X^' + str(i - 1)]],\
                    [2, var_dict['X^' + str(i)]],\
                    [1, var_dict['X^' + str(i + 1)]],\
                ]
            output = sum_vars(in_fp, outputs)
            [new_constraints, new_no_variables, new_non_zero_entries] = assert_square(in_fp, inp, output)

        # Update the variables
        constraints += new_constraints * k
        variables += new_no_variables * k
        non_zero_entries += new_non_zero_entries * k

    # Build the last two constraints

    # First constraint
    left_input = var_dict['X^1']
    right_input_list = [[-1, var_dict['y0']]]
    # Compute the polynomial on the right hand side.
    for i in range(1, ell + 1):
        if in_fp and i != ell and i % 2 == 1 and i > floor(ell / 2):
            # Here we can use this linear combination. We can omit the next
            # power, and the first power can be corrected for by an earlier power.
            right_input_list.append([1, var_dict['xX^' + str(i - 1) + ' + cX^' + str(i) + ' + cX^' + str(i + 1)]])
        elif in_fp and (i + 1) < ell and (i + 1) % 2 == 1 and (i + 1) > floor(ell / 2):
            # Compensate for the odd square trick
            right_input_list.append([\
                coefs[i + 1] - coefs[i + 2]^2 / coefs[i + 3] / 4,\
                var_dict['X^' + str(i)],\
            ])
        elif in_fp and (i - 1) < ell and (i - 1) % 2 == 1 and (i - 1) > floor(ell / 2):
            # Do nothing, this power has already been included in the previous linear combination!
            pass
        else:
            right_input_list.append([coefs[i + 1], var_dict['X^' + str(i)]])
    right_input = sum_vars(in_fp, right_input_list)
    output = scale_var(in_fp, -coefs[0], var_dict['1'])

    [new_constraints, new_no_variables, new_non_zero_entries] = assert_product(in_fp, left_input, right_input, output, phi(x=X,j=j0))
    constraints += new_constraints * k
    variables += new_no_variables * k
    non_zero_entries += new_non_zero_entries * k

    # Second constraint
    left_input = var_dict['X^' + str(ell)]
    right_input = sum_vars(in_fp, [\
        [coefsp[ell + 1], var_dict['X^1']],\
        [-1, var_dict['y1']],\
    ])
    output_list = [[-coefsp[0], var_dict['1']]]
    for i in range(1, ell - 1 + 1):
        if in_fp and i != ell and i % 2 == 1 and i > floor(ell / 2):
            output_list.append([-coefsp[i]/coefs[i + 1], var_dict['xX^' + str(i - 1) + ' + cX^' + str(i) + ' + cX^' + str(i + 1)]])
        # Compensate for the linear combination
        elif in_fp and i + 1 != ell and (i + 1) % 2 == 1 and (i + 1) > floor(ell / 2):
            output_list.append([-coefsp[i] + coefsp[i + 1] / coefs[i + 2]*coefs[i+2]^2/coefs[i + 3]/4, var_dict['X^' + str(i)]])
        elif in_fp and i - 1 != ell and (i - 1) % 2 == 1 and (i - 1) > floor(ell / 2):
            output_list.append([-coefsp[i] + coefsp[i - 1]/coefs[i]*coefs[i+1], var_dict['X^' + str(i)]])
        else:
            output_list.append([-coefsp[i], var_dict['X^' + str(i)]])
    output = sum_vars(in_fp, output_list)

    [new_constraints, new_no_variables, new_non_zero_entries] = assert_product(in_fp, left_input, right_input, output, X^(ell + 1)/ell^s * phi(x=ell^s/X,j=j1))
    constraints += new_constraints * k
    variables += new_no_variables * k
    non_zero_entries += new_non_zero_entries * k

    return constraints, variables, non_zero_entries

# A specialized method for ell=2 over Fp
def special_2_over_fp():
    ell = 2
    phi, s, coefs, coefsp = phi_coefs_for_ell(ell)

    y0 = coefs[1] - coefs[2]^2/4/coefs[3] - j0 - coefs[0]*coefs[2]/coefs[3]/2/ell^s
    y1 = coefs[1] - coefs[2]^2/4/coefs[3] - j1 - coefs[0]*coefs[2]/coefs[3]/2/ell^s


    constraints = 0
    non_zero_entries = 0
    variables = 2 * 2 * k + 2 * (k + 1)
    var_dict = {
        '1': VARIABLE_ONE,
        'Z+': { 'vars': coefs[2]/2/coefs[3] + X, 'no_vars_used': basis_regular },
        'Z-': { 'vars': coefs[2]/2/coefs[3]/ell^s + 1/X, 'no_vars_used': basis_regular },
        'y0': { 'vars': y0, 'no_vars_used': basis_regular },
        'y1': { 'vars': y1, 'no_vars_used': basis_regular },
    }

    # Compute product
    # X*X^{-1} = 1
    [n_constr, n_vars, n_nnz] = assert_product( \
        True, \
        sum_vars(True, [
            [1, var_dict['Z+']],\
            [-coefs[2]/2/coefs[3], var_dict['1']],\
        ]),
        sum_vars(True, [
            [1, var_dict['Z-']],\
            [-coefs[2]/2/coefs[3]/ell^s, var_dict['1']],\
        ]),
        var_dict['1']
    )
    constraints += n_constr * k
    variables += n_vars * k
    non_zero_entries += n_nnz * k

    # Compute first square
    [n_constr, n_vars, n_nnz] = assert_square(\
        True,\
        var_dict['Z+'],\
        sum_vars(\
            True,
            [# Constant term is included in y
                [-4096, var_dict['Z-']],\
                [-1, var_dict['y0']],\
            ]
        ),\
        phi(x=X,j=j0)/X
    )
    constraints += n_constr * k
    variables += n_vars * k
    non_zero_entries += n_nnz * k

    # Compute second square
    [n_constr, n_vars, n_nnz] = assert_square(\
        True,\
        var_dict['Z-'],\
        sum_vars(\
            True,
            [# Constant term is included in y
                [-1/16777216, var_dict['Z+']],\
                [-1/16777216, var_dict['y1']],\
            ]
        ),\
        X/coefs[3]/ell^(3*s)*phi(x=ell^s/X, j=j1)
    )
    constraints += n_constr * k
    variables += n_vars * k
    non_zero_entries += n_nnz * k

    return constraints, variables, non_zero_entries

# Compute the number of constraints, variables and NNZ
def special_3_over_fp():
    ell = 3
    phi, s, coefs, coefsp = phi_coefs_for_ell(ell)

    y0 = j0 - coefs[1]
    y1 = j1 - coefs[1]

    constraints = 0
    non_zero_entries = 0
    # We don't count the 1
    variables = 3 * (2 * k) + 2 * (k + 1)
    # Variables for real part, variables for imaginary part, variables for sum
    var_dict = {
        '1': VARIABLE_ONE,
        'y0': { 'vars': y0, 'no_vars_used': basis_sum },
        'y1': { 'vars': y1, 'no_vars_used': basis_sum },
        'X^1': { 'vars': X, 'no_vars_used': basis_sum },
        'X^2': { 'vars': X^2, 'no_vars_used': basis_sum },
        'X^3': { 'vars': X^3, 'no_vars_used': basis_sum },
    }

    [new_constraints, new_no_variables, new_non_zero_entries] = assert_square(True,\
            var_dict['X^1'],\
            var_dict['X^2'],\
    )
    # Update the variables
    constraints += new_constraints * k
    variables += new_no_variables * k
    non_zero_entries += new_non_zero_entries * k

    [new_constraints, new_no_variables, new_non_zero_entries] = assert_product(True,\
            var_dict['X^1'],\
            var_dict['X^2'],\
            var_dict['X^3'],\
    )
    # Update the variables
    constraints += new_constraints * k
    variables += new_no_variables * k
    non_zero_entries += new_non_zero_entries * k

    # Build the last two constraints

    # First constraint
    left_input = var_dict['X^1']
    right_input = sum_vars(True, [
        [-1, var_dict['y0']],\
        [1, var_dict['X^3']],\
    ])
    # Compute the polynomial on the right hand side.
    output = sum_vars(True, [\
        [-36, var_dict['X^3']],\
        [-270, var_dict['X^2']],\
        [-729, var_dict['1']],\
    ])
    [new_constraints, new_no_variables, new_non_zero_entries] = assert_product(True, left_input, right_input, output, phi(x=X,j=j0))
    constraints += new_constraints * k
    variables += new_no_variables * k
    non_zero_entries += new_non_zero_entries * k

    # Second constraint
    left_input = var_dict['X^3']
    right_input = sum_vars(True, [\
        [coefsp[ell + 1], var_dict['X^1']],\
        [-1, var_dict['y1']],\
    ])
    output = sum_vars(True, [\
        [-coefsp[0], var_dict['1']],\
        [-19131876, var_dict['X^1']],\
        [-196830, var_dict['X^2']],\
    ])
    [new_constraints, new_no_variables, new_non_zero_entries] = assert_product(True, left_input, right_input, output, X^(ell + 1)/ell^s * phi(x=ell^s/X,j=j1))
    constraints += new_constraints * k
    variables += new_no_variables * k
    non_zero_entries += new_non_zero_entries * k

    return constraints, variables, non_zero_entries


# Compute the number of constraints, variables and NNZ
def special_new_3_over_fp():
    ell = 3
    phi, s, coefs, coefsp = phi_coefs_for_ell(ell)

    y0 = j0 + 13946586480/729
    y1 = j1 + 13946586480/729

    constraints = 0
    non_zero_entries = 0
    # We don't count the 1
    variables = 3 * (2 * k) + 2 * (k + 1)
    # Variables for real part, variables for imaginary part, variables for sum
    var_dict = {
        '1': VARIABLE_ONE,
        'y0': { 'vars': y0, 'no_vars_used': basis_sum },
        'y1': { 'vars': y1, 'no_vars_used': basis_sum },
        'X^1': { 'vars': X, 'no_vars_used': basis_sum },
        'X^2': { 'vars': X^2, 'no_vars_used': basis_sum },
        'square': { 'vars': (X^2 + 18*X - 27)^2, 'no_vars_used': basis_regular },
    }

    [new_constraints, new_no_variables, new_non_zero_entries] = assert_square(True,\
            var_dict['X^1'],\
            var_dict['X^2'],\
    )
    # Update the variables
    constraints += new_constraints * k
    variables += new_no_variables * k
    non_zero_entries += new_non_zero_entries * k

    # First constraint, build the difference between polynomials to eliminate
    # the constant term and divide by X.
    constr1 = (phi(x=X,j=j0)*coefsp[0] - X^(ell+1)/ell^s*phi(x=ell^s/X, j=j1)*coefs[0])/X
    left_input = var_dict['X^2']
    right_input = sum_vars(True, [
        [729, var_dict['y1']],\
        [387419760, var_dict['X^1']],\
    ])
    # Compute the polynomial on the right hand side.
    output = sum_vars(True, [\
        [-104460042960, var_dict['X^1']],\
        [387420489, var_dict['y0']],\
        [-7412066808269760, var_dict['1']],\
    ])
    [new_constraints, new_no_variables, new_non_zero_entries] = assert_product(True, left_input, right_input, output, constr1)
    constraints += new_constraints * k
    variables += new_no_variables * k
    non_zero_entries += new_non_zero_entries * k

    [new_constraints, new_no_variables, new_non_zero_entries] = assert_square(True, sum_vars(True, [\
            [1, var_dict['X^2']],\
            [18, var_dict['X^1']],\
            [-27, var_dict['1']],\
        ]),
        var_dict['square'],
    )
    constraints += new_constraints * k
    variables += new_no_variables * k
    non_zero_entries += new_non_zero_entries * k

    [new_constraints, new_no_variables, new_non_zero_entries] = assert_product(True,
        var_dict['X^1'],
        sum_vars(True, [\
            [19132848, var_dict['1']],\
            [-1, var_dict['y0']],
        ]),
        sum_vars(True, [\
            [-1, var_dict['square']],
        ]),
        phi(x=X,j=j0)
    )
    constraints += new_constraints * k
    variables += new_no_variables * k
    non_zero_entries += new_non_zero_entries * k

    return constraints, variables, non_zero_entries

# Compute the number of constraints, variables and NNZ
def special_5_over_fp():
    ell = 5
    phi, s, coefs, coefsp = phi_coefs_for_ell(ell)

    y0 = j0 - coefs[1]
    y1 = j1 - coefs[1]

    constraints = 0
    non_zero_entries = 0
    # We don't count the 1
    variables = 5 * (2 * k) + 2 * (k + 1)
    # Variables for real part, variables for imaginary part, variables for sum
    var_dict = {
        '1': VARIABLE_ONE,
        'y0': { 'vars': y0, 'no_vars_used': basis_sum },
        'y1': { 'vars': y1, 'no_vars_used': basis_sum },
        'X^1': { 'vars': X, 'no_vars_used': basis_sum },
        'X^2': { 'vars': X^2, 'no_vars_used': basis_sum },
        'X^2 + X^3 + X^4': {
            'vars': coefs[3]^2 / coefs[4] / 4 * X^2 + coefs[3]*X^3 + coefs[4]*X^4,
            'no_vars_used': basis_sum
        },
        'X^4': { 'vars': X^4, 'no_vars_used': basis_sum },
        'X^5': { 'vars': X^5, 'no_vars_used': basis_sum },
    }

    # Now actually compute these powers. Again, we use squaring whereever
    # possible.
    for i in range(2, ell + 1):
        if i % 2 == 0:
            # A simple squaring!
            compute_from = var_dict['X^' + str(i / 2)]
            output = var_dict['X^' + str(i)]
            [new_constraints, new_no_variables, new_non_zero_entries] = assert_square(True, compute_from, output)
        elif i == 3:
            [new_constraints, new_no_variables, new_non_zero_entries] = assert_square(True,\
                sum_vars(True, [\
                    [coefs[3]/coefs[4]/2, var_dict['X^1']],\
                    [1, var_dict['X^2']],\
                ]),
                scale_var(True, 1/coefs[4], var_dict['X^2 + X^3 + X^4']),
            )
        else:
            assert i == 5
            [new_constraints, new_no_variables, new_non_zero_entries] = assert_product(True,\
                var_dict['X^1'],\
                var_dict['X^4'],\
                var_dict['X^5'],\
            )

        # Update the variables
        constraints += new_constraints * k
        variables += new_no_variables * k
        non_zero_entries += new_non_zero_entries * k

    # Build the last two constraints

    # First constraint
    left_input = var_dict['X^1']
    right_input = sum_vars(True, [
        [-1, var_dict['y0']],\
        [1, var_dict['X^5']],\
    ])
    # Compute the polynomial on the right hand side.
    output = sum_vars(True, [\
        [-30, var_dict['X^5']],\
        [-1, var_dict['X^2 + X^3 + X^4']],\
        [-14725/63, var_dict['X^2']],\
        [-125, var_dict['1']],\
    ])
    [new_constraints, new_no_variables, new_non_zero_entries] = assert_product(True, left_input, right_input, output, phi(x=X,j=j0))
    constraints += new_constraints * k
    variables += new_no_variables * k
    non_zero_entries += new_non_zero_entries * k

    # Second constraint
    left_input = var_dict['X^5']
    right_input = sum_vars(True, [\
        [coefsp[ell + 1], var_dict['X^1']],\
        [-1, var_dict['y1']],\
    ])
    output = sum_vars(True, [\
        [4725000, var_dict['X^4']],\
        [-coefsp[3]/coefs[3], var_dict['X^2 + X^3 + X^4']],\
        [-37439453125/63, var_dict['X^2']],\
        [-7324218750, var_dict['X^1']],\
        [-30517578125, var_dict['1']],\
    ])
    [new_constraints, new_no_variables, new_non_zero_entries] = assert_product(True, left_input, right_input, output, X^(ell + 1)/ell^s * phi(x=ell^s/X,j=j1))
    constraints += new_constraints * k
    variables += new_no_variables * k
    non_zero_entries += new_non_zero_entries * k

    return constraints, variables, non_zero_entries

# Compute the number of constraints, variables and NNZ, over Fp^2
def better_approach(ell):
    phi, s, coefs, coefsp = phi_coefs_for_ell(ell)

    constraints = 0
    non_zero_entries = 0
    t = (ell + 1) / 2

    y0 = j0 - coefs[1]
    y1 = j1 - coefs[1]

    # We don't count the 1
    variables = k + k + 1
    var_dict = {
        '1': {'vars': 1, 'no_vars_used': 1 },
        'X^1': { 'vars': X, 'no_vars_used': 1 },
        'y0': { 'vars': y0, 'no_vars_used': 1 },
        'y1': { 'vars': y1, 'no_vars_used': 1 }
    }

    # Store the powers of X
    # Unlike in the simple approach, we also need X^{t - 1} in plain.
    for i in range(2, t + 1):
        variables += 1 * k
        var_dict['X^' + str(i)] = {
            'vars': X^i,
            'no_vars_used': 1
        }

    # Compute the powers of X, inclusive
    for i in range(2, t + 1):
        if i % 2 == 0:
            # We can square
            compute_from = var_dict['X^' + str(i / 2)]
            output = var_dict['X^' + str(i)]
            [new_constraints, new_no_variables, new_non_zero_entries] = assert_square(false, compute_from, output)
        else:
            # We do not have a higher even power to subtract, so multiply directly
            left_hand_side = var_dict['X^1']
            right_hand_side = var_dict['X^' + str(i - 1)]
            output = var_dict['X^' + str(i)]
            [new_constraints, new_no_variables, new_non_zero_entries] = assert_product(false, left_hand_side, right_hand_side, output)

        # Update the variables
        constraints += new_constraints * k
        variables += new_no_variables * k
        non_zero_entries += new_non_zero_entries * k

    # Compute the crossterms with j
    left_hand_side = var_dict['X^1']
    right_hand_side = var_dict['y0']
    variables += 1 * k
    var_dict['X^1y0'] = { 'vars': X*y0, 'no_vars_used': 1 }
    output = var_dict['X^1y0']
    [n_constr, n_vars, nnnz] = assert_product(false, left_hand_side, right_hand_side, output)
    constraints += n_constr * k
    variables += n_vars * k
    non_zero_entries += nnnz * k

    left_hand_side = var_dict['X^' + str(t - 1)]
    right_hand_side = var_dict['y1']
    variables += 1 * k
    var_dict['X^' + str(t - 1) + 'y1'] = { 'vars': X^(t-1)*y1, 'no_vars_used': 1 }
    output = var_dict['X^' + str(t - 1) + 'y1']
    [n_constr, n_vars, nnnz] = assert_product(false, left_hand_side, right_hand_side, output)
    constraints += n_constr * k
    variables += n_vars * k
    non_zero_entries += nnnz * k

    # Build the last two constraints

    # First constraint
    left_input = var_dict['X^' + str(t)]
    right_input_list = [[coefs[t], var_dict['1']]]
    for i in range(1, t + 1):
        assert i - 1 != t
        right_input_list.append([coefs[i + t], var_dict['X^' + str(i)]])
    right_input = sum_vars(false, right_input_list)
    output_list = [\
        [-coefs[0], var_dict['1']],\
        [1, var_dict['X^1y0']],\
    ]
    for i in range(2, t - 1 + 1):
        output_list.append([-coefs[i], var_dict['X^' + str(i)]])
    output = sum_vars(false, output_list)
    [new_constraints, new_no_variables, new_non_zero_entries] = assert_product(false, left_input, right_input, output, phi(x=X,j=j0))
    constraints += new_constraints * k
    variables += new_no_variables * k
    non_zero_entries += new_non_zero_entries * k

    # Second constraint
    left_input = var_dict['X^' + str(t)]
    right_input_list = [\
        [coefsp[t], var_dict['1']],\
        [-1, var_dict['X^' + str(t - 1) + 'y1']],\
    ]
    for i in range(1, t + 1):
        if i == t - 1:
            # This is included in y'.
            pass
        else:
            right_input_list.append([coefsp[t + i], var_dict['X^' + str(i)]])
    right_input = sum_vars(false, right_input_list)
    output_list = [[-coefsp[0], var_dict['1']]]
    for i in range(1, t - 1 + 1):
        output_list.append([-coefsp[i], var_dict['X^' + str(i)]])
    output = sum_vars(false, output_list)
    [new_constraints, new_no_variables, new_non_zero_entries] = assert_product(false, left_input, right_input, output, X^(ell+1)/ell^s*phi(x=ell^s/X,j=j1))
    constraints += new_constraints * k
    variables += new_no_variables * k
    non_zero_entries += new_non_zero_entries * k

    return constraints, variables, non_zero_entries

def special_7_over_fp():
    ell = 7
    phi, s, coefs, coefsp = phi_coefs_for_ell(ell)

    constraints = 0
    non_zero_entries = 0
    t = (ell + 1) / 2
    
    y0 = 1728 - j0
    y1 = 1728 - j1

    # We don't count the 1
    variables = 2 * k + 2 * (k + 1)
    # Variables for real part, variables for imaginary part, variables for sum
    var_dict = {
        '1': VARIABLE_ONE,
        'X^1': { 'vars': X, 'no_vars_used': basis_sum },
        'y0': { 'vars': y0, 'no_vars_used': basis_sum },
        'y1': { 'vars': y1, 'no_vars_used': basis_sum }
    }

    # Store the powers of X
    # Unlike in the simple approach, we also need X^{t - 1} in plain.
    for i in range(2, t + 1):
        variables += 2 * k
        if i == 3:
            var_dict['xX^2 + cX^3 + cX^4'] = {
                'vars': X^4 + 14*X^3 + (14*14/4 + 63 - 14*14/4)*X^2,
                'no_vars_used': basis_sum
            }
        else:
            var_dict['X^' + str(i)] = {
                'vars': X^i,
                'no_vars_used': basis_sum
            }

    # Compute the powers of X, inclusive
    for i in range(2, t + 1):
        if i % 2 == 0:
            # We can square
            compute_from = var_dict['X^' + str(i / 2)]
            output = var_dict['X^' + str(i)]
            [new_constraints, new_no_variables, new_non_zero_entries] = assert_square(true, compute_from, output)
        else:
            assert i == 3
            # Use the squaring trick for odd powers!
            inputs = [\
                [7, var_dict['X^1']],\
                [1, var_dict['X^2']],\
            ]
            inp = sum_vars(true, inputs)
            output = sum_vars(true, [\
                [1, var_dict['xX^2 + cX^3 + cX^4']],\
                [-14, var_dict['X^2']],\
            ])
            [new_constraints, new_no_variables, new_non_zero_entries] = assert_square(true, inp, output)

        # Update the variables
        constraints += new_constraints * k
        variables += new_no_variables * k
        non_zero_entries += new_non_zero_entries * k

    # Compute the crossterms with j
    variables += 2 * k
    var_dict['X^1y0'] = { 'vars': X*y0, 'no_vars_used': basis_sum }

    left_hand_side = var_dict['X^1']
    right_hand_side = var_dict['y0']
    output = var_dict['X^1y0']
    [n_constr, n_vars, nnnz] = assert_product(true, left_hand_side, right_hand_side, output)
    constraints += n_constr * k
    variables += n_vars * k
    non_zero_entries += nnnz * k

    variables += 2 * k
    var_dict['1 + X^1 + X^2 + X^3j + X^3 + X^4'] = {
        'vars': coefsp[t] - coefsp[3]/14 + coefsp[t + 1] * X^1 + coefsp[t + 2] * X^2 + coefsp[t + 3] * X^3 + coefsp[t + 4] * X^4 - j1*X^3,
        'no_vars_used': basis_sum
    }
    left_hand_side = sum_vars(true, [\
        [-9/2, var_dict['X^2']],\
        [1/14, var_dict['xX^2 + cX^3 + cX^4']],\
        [-1/14, var_dict['X^4']],
    ])
    # Here we subtract the factor for X^7 away
    right_hand_side = sum_vars(true, [\
        [1, var_dict['y1']],\
        [coefsp[t + 3] - 1728, var_dict['1']],\
    ])
    output = sum_vars(true, [\
        [-coefsp[t] + coefsp[3]/14, var_dict['1']],\
        [-coefsp[t + 1], var_dict['X^1']],\
        [-coefsp[t + 2], var_dict['X^2']],\
        [-coefsp[t + 4], var_dict['X^4']],\
        [1, var_dict['1 + X^1 + X^2 + X^3j + X^3 + X^4']],\
    ])
    [n_constr, n_vars, nnnz] = assert_product(true, left_hand_side, right_hand_side, output)
    constraints += n_constr * k
    variables += n_vars * k
    non_zero_entries += nnnz * k

    # Build the last two constraints

    # First constraint
    inpt = sum_vars(true, [\
        [-7, var_dict['1']],\
        [70, var_dict['X^1']],\
        [1, var_dict['xX^2 + cX^3 + cX^4']],\
    ])
    output = scale_var(true, -1, var_dict['X^1y0'])
    [new_constraints, new_no_variables, new_non_zero_entries] = assert_square(true, inpt, output, phi(x=X, j=j0))
    constraints += new_constraints * k
    variables += new_no_variables * k
    non_zero_entries += new_non_zero_entries * k

    # Second constraint
    left_input = var_dict['X^4']
    right_input = var_dict['1 + X^1 + X^2 + X^3j + X^3 + X^4']
    output_list = [
        [-coefsp[0], var_dict['1']],\
        [-coefsp[1], var_dict['X^1']],\
        [-41564215210, var_dict['X^2']],\
        [-coefsp[3]/14, var_dict['xX^2 + cX^3 + cX^4']],\
    ]
    output = sum_vars(true, output_list)
    [new_constraints, new_no_variables, new_non_zero_entries] = assert_product(true, left_input, right_input, output, X^(ell + 1)/ell^2 * phi(x=ell^s/X, j=j1))
    constraints += new_constraints * k
    variables += new_no_variables * k
    non_zero_entries += new_non_zero_entries * k

    return constraints, variables, non_zero_entries

def special_13_over_fp():
    ell = 13
    phi, s, coefs, coefsp = phi_coefs_for_ell(ell)
    constraints = 0
    non_zero_entries = 0
    t = (ell + 1) / 2
    y0 = j0
    y1 = j1

    # We don't count the 1
    variables = 2 * k + 2 * (k + 1)
    # Variables for real part, variables for imaginary part, variables for sum
    var_dict = {
        '1': VARIABLE_ONE,
        'X^1': { 'vars': X, 'no_vars_used': basis_sum },
        'y0': { 'vars': y0 - 4464, 'no_vars_used': basis_sum },
        'y1': { 'vars': y1 - 4464, 'no_vars_used': basis_sum }
    }

    # Store the powers of X
    # Unlike in the simple approach, we also need X^{t - 1} in plain.
    for i in range(2, t + 1):
        variables += 2 * k
        if i == 5:
            # 260*x^4 + 78*x^5 + 13*x^6
            var_dict['xX^4 + cX^5 + cX^6'] = {
                'vars': -13 + 143*X + 468*X^2 + 494*X^3 + 260*X^4 + 78*X^5 + 13*X^6 + X^7,
                'no_vars_used': basis_sum
            }
        else:
            var_dict['X^' + str(i)] = {
                'vars': X^i,
                'no_vars_used': basis_sum
            }

    # Compute the powers of X, inclusive
    for i in range(2, t + 1):
        if i % 2 == 0:
            # We can square
            compute_from = var_dict['X^' + str(i / 2)]
            output = var_dict['X^' + str(i)]
            [new_constraints, new_no_variables, new_non_zero_entries] = assert_square(true, compute_from, output)
        elif i == 7:
            # We do not have a higher even power to subtract, so multiply directly
            left_hand_side = var_dict['X^1']
            right_hand_side = var_dict['X^6']
            output = var_dict['X^7']
            [new_constraints, new_no_variables, new_non_zero_entries] = assert_product(true, left_hand_side, right_hand_side, output)
        elif i == 3:
            inp = sum_vars(true, [\
                [1, var_dict['X^1']],\
                [1, var_dict['X^2']],\
            ])
            output = sum_vars(true, [\
                [1, var_dict['X^2']],\
                [2, var_dict['X^3']],\
                [1, var_dict['X^4']],\
            ])
            [new_constraints, new_no_variables, new_non_zero_entries] = assert_square(true, inp, output)
        elif i == 5:
            inp = sum_vars(true, [\
                [78/13/2, var_dict['X^2']],\
                [1, var_dict['X^3']],\
            ])
            output = sum_vars(true, [
                [1, var_dict['1']],\
                [-11, var_dict['X^1']],\
                [-36, var_dict['X^2']],\
                [-38, var_dict['X^3']],\
                [-11, var_dict['X^4']],\
                [1/13, var_dict['xX^4 + cX^5 + cX^6']],\
                [-1/13, var_dict['X^7']],\
            ])
            [new_constraints, new_no_variables, new_non_zero_entries] = assert_square(true, inp, output)

        # Update the variables
        constraints += new_constraints * k
        variables += new_no_variables * k
        non_zero_entries += new_non_zero_entries * k

    # Compute the crossterms with j
    variables += 2 * k
    var_dict['X(j + X + X^2 + X^3 + X^4)'] = {
        'vars': 104*X^5 + 988*X^4 + 3848*X^3 - j0*X + 6864*X^2 + 4464*X - 156,
        'no_vars_used': basis_regular
    }
    left_hand_side = var_dict['X^1']
    right_hand_side = sum_vars(true, [\
        [-1, var_dict['y0']],\
        [104, var_dict['X^4']],\
    ])
    output = sum_vars(true, [\
        [1, var_dict['X(j + X + X^2 + X^3 + X^4)']],\
        [156, var_dict['1']],\
        [-6864, var_dict['X^2']],\
        [-3848, var_dict['X^3']],\
        [-988, var_dict['X^4']],\
    ])

    [n_constr, n_vars, nnnz] = assert_product(true, left_hand_side, right_hand_side, output)
    constraints += n_constr * k
    variables += n_vars * k
    non_zero_entries += nnnz * k

    # Build the right hand side of the second constraint
    variables += 2 * k
    var_dict['X^6(j + 1)'] = {
        'vars': 3125503631770/3 + coefsp[8] * X + coefsp[9] * X^2 + coefsp[10] * X^3 + coefsp[11] * X^4 + coefsp[12] * X^5 + coefsp[13] * X^6 + coefsp[14]*X^7 -j1*X^6,
        'no_vars_used': basis_sum
    }

    left_hand_side = var_dict['X^6']
    right_hand_side = sum_vars(true, [\
        [-1, var_dict['y1']],\
        [-219193/6, var_dict['1']],\
    ])
    output = sum_vars(true, [\
        [-2083669153475/2, var_dict['1']],\
        [-1134742015693/6, var_dict['X^1']],\
        [-15273812710, var_dict['X^2']],\
        [-2333005961/3, var_dict['X^3']],\
        [-61331114/3, var_dict['X^4']],\
        [-196885/78, var_dict['xX^4 + cX^5 + cX^6']],\
        [1, var_dict['X^6(j + 1)']],\
        [15139/6, var_dict['X^7']],\
    ])

    [n_constr, n_vars, nnnz] = assert_product(true, left_hand_side, right_hand_side, output)
    constraints += n_constr * k
    variables += n_vars * k
    non_zero_entries += nnnz * k

    # Build the last two constraints

    # First constraint
    # (-13 + 143*x + 468*x^2 + 494*x^3 + 260*x^4 + 78*x^5 + 13*x^6 + x^7)^2
    inpt = var_dict['xX^4 + cX^5 + cX^6']
    output = scale_var(true, -1, var_dict['X(j + X + X^2 + X^3 + X^4)'])
    [new_constraints, new_no_variables, new_non_zero_entries] = assert_square(true, inpt, output, phi(x=X, j=j0))
    constraints += new_constraints * k
    variables += new_no_variables * k
    non_zero_entries += new_non_zero_entries * k

    # Second constraint
    left_input = var_dict['X^7']
    right_input = var_dict['X^6(j + 1)']
    output = sum_vars(true, [
        [-930788723466329/3, var_dict['1']],\
        [-1573453198968248/3, var_dict['X^1']],\
        [-316491283787185, var_dict['X^2']],\
        [-211580971490096/3, var_dict['X^3']],\
        [3223767809392/3, var_dict['X^4']],\
        [-44326807379140/78, var_dict['xX^4 + cX^5 + cX^6']],\
        [-7413360792448/3, var_dict['X^6']]
    ])
    [new_constraints, new_no_variables, new_non_zero_entries] = assert_product(true, left_input, right_input, output, X^(ell + 1)/ell^s*phi(x=ell^s/X, j=j1))
    constraints += new_constraints * k
    variables += new_no_variables * k
    non_zero_entries += new_non_zero_entries * k

    return constraints, variables, non_zero_entries

print('Original paper')
print('Over Fp2')
constraints, variables, non_zero_entries = original_paper(false)
print('- constraints: ', constraints)
print('- variables: ', variables)
print('- non_zero_entries: ', non_zero_entries)
print('Over Fp')
constraints, variables, non_zero_entries = original_paper(true)
print('- constraints: ', constraints)
print('- variables: ', variables)
print('- non_zero_entries: ', non_zero_entries)


print('\n\nSimple method:')

var('lamb')

def compute_relative_absolute(do_relative, ell):
    if do_relative:
        return lamb / log(ell, 2).n()
    else:
        return k

# True if the numbers should be computed for some security parameter lambda
# instead of the number of steps k directly.
do_relative = True

for ell in [2, 3, 5]:
    # Comment for concrete results
    k = compute_relative_absolute(do_relative, ell)
    print('ell = ', ell, ' Fp2')
    constraints, variables, non_zero_entries = simple_approach(ell, false)
    print('- constraints: ', constraints(k = k))
    print('- variables: ', variables(k = k))
    print('- non_zero_entries: ', non_zero_entries(k = k))

print('\n\nBetter method:')

for ell in [7, 13]:
    # Comment for concrete results
    k = compute_relative_absolute(do_relative, ell)
    print('ell = ', ell, ' Fp2')
    constraints, variables, non_zero_entries = better_approach(ell)
    print('- constraints: ', constraints(k = k))
    print('- variables: ', variables(k = k))
    print('- non_zero_entries: ', non_zero_entries(k = k))

# Relative and concrete results are equal
print('\n\nSpecial case for ell = 2 over Fp')
# Comment for concrete results
k = compute_relative_absolute(do_relative, 2)
constraints, variables, non_zero_entries = special_2_over_fp()
print('- constraints: ', constraints(k = k))
print('- variables: ', variables(k = k))
print('- non_zero_entries: ', non_zero_entries(k = k))

# Relative and concrete results are equal
print('\n\nSpecial case for ell = 3 over Fp')
# Comment for concrete results
k = compute_relative_absolute(do_relative, 3)
constraints, variables, non_zero_entries = special_3_over_fp()
print('- constraints: ', constraints(k = k))
print('- variables: ', variables(k = k))
print('- non_zero_entries: ', non_zero_entries(k = k))

# Relative and concrete results are equal
print('\n\nSpecial case for ell = 5 over Fp')
# Comment for concrete results
k = compute_relative_absolute(do_relative, 5)
constraints, variables, non_zero_entries = special_5_over_fp()
print('- constraints: ', constraints(k = k))
print('- variables: ', variables(k = k))
print('- non_zero_entries: ', non_zero_entries(k = k))

# Relative and concrete results are equal
print('\n\nSpecial case for ell = 7 over Fp')
# Comment for concrete results
k = compute_relative_absolute(do_relative, 7)
constraints, variables, non_zero_entries = special_7_over_fp()
print('- constraints: ', constraints(k = k))
print('- variables: ', variables(k = k))
print('- non_zero_entries: ', non_zero_entries(k = k))

# Relative and concrete results are equal
print('\n\nSpecial case for ell = 13 over Fp')
# Comment for concrete results
k = compute_relative_absolute(do_relative, 13)
constraints, variables, non_zero_entries = special_13_over_fp()
print('- constraints: ', constraints(k = k))
print('- variables: ', variables(k = k))
print('- non_zero_entries: ', non_zero_entries(k = k))

# Relative and concrete results are equal
print('\n\nAlternative method for ell = 3 over Fp')
# Comment for concrete results
k = compute_relative_absolute(do_relative, 3)
constraints, variables, non_zero_entries = special_new_3_over_fp()
print('- constraints: ', constraints(k = k))
print('- variables: ', variables(k = k))
print('- non_zero_entries: ', non_zero_entries(k = k))
