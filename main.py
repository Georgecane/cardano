import sympy as sp
import re
from decimal import Decimal, getcontext
import math
from itertools import permutations
import os
import pycfd.symbolicly.symbolic_regression as sre

getcontext().prec = 30

# Global storage for all mathematical objects
algebraic_structures = {}
tensors = {}
variables = {}
functions = {}
matrices = {}
user_sets = {}
number_sets = {
    'N': 'Natural numbers',
    'Z': 'Integers',
    'Q': 'Rational numbers',
    'R': 'Real numbers',
    'C': 'Complex numbers'
}

def clear_screen():
    """Clear the terminal screen"""
    try:
        # For Windows
        if os.name == 'nt':
            _ = os.system('cls')
        # For Unix/Linux/MacOS
        else:
            _ = os.system('clear')
        return "Screen cleared."
    except Exception as e:
        return f"Error clearing screen: {e}"

def preprocess_expression(expression):

    processed_expression = expression.replace('**', '^')
    processed_expression = re.sub(r'(\d)([a-zA-Z])', r'\1*\2', processed_expression)
    processed_expression = re.sub(r'([a-zA-Z])(\d)', r'\1*\2', processed_expression)
    processed_expression = re.sub(r'(\))(\d)', r'\1*\2', processed_expression)
    processed_expression = re.sub(r'(\d)(\()', r'\1*\2', processed_expression)

    trig_funcs = ['sin', 'cos', 'tan', 'cot', 'sec', 'csc', 'asin', 'acos', 'atan', 'acot', 'asec', 'acsc']
    for func in trig_funcs:
        processed_expression = processed_expression.replace(f'{func}(', f'sp.{func}(')

    processed_expression = processed_expression.replace('exp', 'sp.exp')
    processed_expression = processed_expression.replace('log', 'sp.log')
    processed_expression = processed_expression.replace('log10', 'sp.log10')
    processed_expression = processed_expression.replace('log2', 'sp.log2')

    processed_expression = re.sub(r'nroot\(([^,]+),\s*([^)]+)\)', r'compute_nroot(\1, \2)', processed_expression)


    processed_expression = processed_expression.replace('nPr', 'math.perm')
    processed_expression = processed_expression.replace('nCr', 'sp.binomial')
    processed_expression = processed_expression.replace('oo', 'sp.oo')
    processed_expression = processed_expression.replace('set', 'define_set')

    processed_expression = processed_expression.replace('pi', 'sp.pi')
    processed_expression = processed_expression.replace('mathe', '2.718281828459045235360287471352')
    processed_expression = processed_expression.replace("mp", "1.672621637 * 10 ^ -27")
    processed_expression = processed_expression.replace("mn", "1.674927211 * 10 ^ -27")
    processed_expression = processed_expression.replace("me", "9.10938215 * 10 ^ -31")
    processed_expression = processed_expression.replace("mμ", "1.8835313 * 10 ^ -28")
    processed_expression = processed_expression.replace('a_0', '5.291772086 * 10 ^ -11')
    processed_expression = processed_expression.replace('μN', '5.05078324 * 10 ^ -27')
    processed_expression = processed_expression.replace('μB', '9.27400915 * 10 -24')
    processed_expression = processed_expression.replace('ħ', '1.054571628 * 10 ^ -34')
    processed_expression = processed_expression.replace('Na', '6.02214076 * 10 ^ 23')
    processed_expression = processed_expression.replace('C_0', '300000000')
    processed_expression = processed_expression.replace('r_number', '0.6172839455060657')
    processed_expression = processed_expression.replace('I', 'sp.I')
    return processed_expression

def compute_derivative(expression, var):
    try:
        processed_expression = preprocess_expression(expression)
        expr = sp.sympify(processed_expression)
        var_symbol = sp.Symbol(var)
        derivative = sp.diff(expr, var_symbol)
        return clean_multiplication(str(preprocess_expression(str(derivative))))
    except Exception as e:
        return f"Error in compute_derivative: {e}"

def compute_partial_derivative(expression, var):
    try:
        processed_expression = preprocess_expression(expression)
        expr = sp.sympify(processed_expression)
        var_symbol = sp.Symbol(var)
        partial_derivative = sp.diff(expr, var_symbol)
        return clean_multiplication(str(preprocess_expression(str(partial_derivative))))
    except Exception as e:
        return f"Error in compute_partial_derivative: {e}"

def convert_to_input_format(expression):
    expression = str(expression)
    expression = re.sub(r'(\d)\*([a-zA-Z])', r'\1\2', expression)
    expression = re.sub(r'([a-zA-Z])\*([a-zA-Z])', r'\1\2', expression)
    expression = expression.replace('**', '^')
    return expression

def clean_multiplication(expression):
    expression = re.sub(r'(\d)\*([a-zA-Z])', r'\1\2', expression)
    expression = re.sub(r'([a-zA-Z])\*([a-zA-Z])', r'\1\2', expression)
    expression = expression.replace("**", "^")

    return expression

def evaluate_expression(expression):
    try:
        processed_expression = preprocess_expression(expression)

        if 'compute_nroot' in processed_expression:
            match = re.match(r'compute_nroot\(([^,]+),\s*([^)]+)\)', processed_expression)
            if match:
                n, k = match.groups()
                return compute_nroot(float(n), float(k))

        eval_dict = {
            'sp': sp, 'Decimal': Decimal, 'compute_nroot': compute_nroot,
            'compute_determinant': compute_determinant, 'compute_inverse': compute_inverse,
            'compute_eigenvalues': compute_eigenvalues, 'compute_eigenvectors': compute_eigenvectors,
            'multiply_matrices': multiply_matrices, 'add_matrices': add_matrices,
            'transpose_matrix': transpose_matrix, **variables, **functions, **number_sets, **user_sets, **algebraic_structures,
            'math': math
        }

        # Evaluate general expressions
        symbols = re.findall(r'\b[a-zA-Z_]\w*\b', processed_expression)
        for symbol in symbols:
            if symbol not in eval_dict:
                eval_dict[symbol] = sp.Symbol(symbol)

        evaluated_expression = sp.sympify(processed_expression, locals=eval_dict)

        if isinstance(evaluated_expression, sp.Basic):
            expanded_expression = sp.expand(evaluated_expression)
        else:
            expanded_expression = evaluated_expression

        if expanded_expression.is_number:
            if isinstance(evaluated_expression, sp.Rational):
                if evaluated_expression.q == 1 or evaluated_expression.q % 2 == 0 or evaluated_expression.q % 5 == 0:
                    return float(evaluated_expression)
            else:
                return str(expanded_expression.evalf())

        return clean_multiplication(str(expanded_expression))

    except sp.SympifyError as e:
        return f"Error in evaluating expression: {e}"
    except Exception as e:
        return f"Error in evaluate_expression: {e}"

def solve_equation(equation):
    """Solve mathematical equations"""
    try:
        # Split equation into left and right sides
        left_side, right_side = equation.split('=')

        # Process both sides
        left_expr = preprocess_expression(left_side)
        right_expr = preprocess_expression(right_side)

        # Convert to SymPy expressions
        left = sp.sympify(left_expr)
        right = sp.sympify(right_expr)

        # Move everything to left side
        equation = left - right

        # Find all symbols in equation
        symbols = list(equation.free_symbols)

        if not symbols:
            return "No variables to solve for"

        # Solve equation
        solution = sp.solve(equation, symbols[0])
        return f"{symbols[0]} = {solution}"
    except Exception as e:
        return f"Error solving equation: {e}"

def solve_system_of_equations(equations):
    """Solve a system of equations"""
    try:
        # Clean and process each equation
        processed_equations = []
        variables = set()

        for eq in equations:
            # Strip whitespace and ensure equation format
            eq = eq.strip()
            if '=' not in eq:
                eq += '= 0'

            # Split into LHS and RHS
            lhs, rhs = eq.split('=')
            lhs = lhs.strip()
            rhs = rhs.strip()

            # Convert to standard form (all terms on LHS)
            expr = f"{lhs}-({rhs})"
            processed_equations.append(expr)

            # Collect variables
            for var in ['x', 'y', 'z']:  # Add more variables if needed
                if var in eq:
                    variables.add(var)

        # Convert variables to SymPy symbols
        symbols = [sp.Symbol(var) for var in sorted(variables)]

        # Convert equations to SymPy expressions
        sympy_equations = [sp.sympify(eq) for eq in processed_equations]

        # Solve the system
        solution = sp.solve(sympy_equations, symbols)

        # Format the solution
        if isinstance(solution, dict):
            return ", ".join([f"{var} = {val}" for var, val in solution.items()])
        elif isinstance(solution, tuple):
            return ", ".join([f"{symbols[i]} = {val}" for i, val in enumerate(solution)])
        else:
            return str(solution)

    except Exception as e:
        return f"Error solving system: {e}"

def define_variable(expression):
    try:
        match = re.match(r'([a-zA-Z_]\w*)\s*:=\s*(.*)', expression)
        if match:
            var_name = match.group(1)
            var_value = match.group(2)
            variables[var_name] = sp.sympify(var_value)
            return f"Variable {var_name} defined."
        else:
            return "Invalid variable definition."
    except Exception as e:
        return f"Error in define_variable: {e}"

def clear_variable(var_name):
    try:
        if var_name in variables:
            del variables[var_name]
            return f"Variable {var_name} cleared."
        else:
            return f"Variable {var_name} does not exist."
    except Exception as e:
        return f"Error in clear_variable: {e}"

def clear_matrix(name):
    try:
        if name in matrices:
            del matrices[name]
            return f"Matrix {name} cleared."
        else:
            return f"Matrix {name} does not exist."
    except Exception as e:
        return f"Error in clear_matrix: {e}"

def define_matrix(name, matrix_str):
    try:
        matrix_str = matrix_str.strip()[1:-1]  # Remove the curly braces
        rows = matrix_str.split('}, {')
        matrix = []
        for row in rows:
            elements = row.split(',')
            matrix.append([sp.sympify(e.strip()) for e in elements])
        matrices[name] = sp.Matrix(matrix)
        return f"Matrix {name} defined."
    except Exception as e:
        return f"Error in define_matrix: {e}"

def compute_determinant(matrix_name):
    try:
        if matrix_name in matrices:
            return matrices[matrix_name].det()
        else:
            return "Matrix not defined."
    except Exception as e:
        return f"Error in compute_determinant: {e}"

def compute_inverse(matrix_name):
    try:
        if matrix_name in matrices:
            return matrices[matrix_name].inv()
        else:
            return "Matrix not defined."
    except Exception as e:
        return f"Error in compute_inverse: {e}"

def compute_eigenvalues(matrix_name):
    try:
        if matrix_name in matrices:
            return matrices[matrix_name].eigenvals()
        else:
            return "Matrix not defined."
    except Exception as e:
        return f"Error in compute_eigenvalues: {e}"

def compute_eigenvectors(matrix_name):
    try:
        if matrix_name in matrices:
            eigenvectors = matrices[matrix_name].eigenvects()
            eigenvectors_dict = {str(key): str(value) for key, value in eigenvectors}
            return eigenvectors_dict
        else:
            return "Matrix not defined."
    except Exception as e:
        return f"Error in compute_eigenvectors: {e}"

def multiply_matrices(matrix_name1, matrix_name2):
    try:
        if matrix_name1 in matrices and matrix_name2 in matrices:
            return matrices[matrix_name1] * matrices[matrix_name2]
        else:
            return "One or both matrices not defined."
    except Exception as e:
        return f"Error in multiply_matrices: {e}"

def add_matrices(matrix_name1, matrix_name2):
    try:
        if matrix_name1 in matrices and matrix_name2 in matrices:
            return matrices[matrix_name1] + matrices[matrix_name2]
        else:
            return "One or both matrices not defined."
    except Exception as e:
        return f"Error in add_matrices: {e}"

def transpose_matrix(matrix_name):
    try:
        if matrix_name in matrices:
            return matrices[matrix_name].transpose()
        else:
            return "Matrix not defined."
    except Exception as e:
        return f"Error in transpose_matrix: {e}"


def define_function(func_definition):
    try:
        # Regular expression to match multi-variable function definition like `f :: (x) = expr`
        match = re.match(r'(\w+)\s*::\s*\(([\w\s,]+)\)\s*=\s*(.*)', func_definition)

        if match:
            func_name = match.group(1)
            args = match.group(2).replace(" ", "").split(",")  # Split multiple arguments
            body = match.group(3)

            # Convert argument names to SymPy symbols
            arg_symbols = [sp.Symbol(arg) for arg in args]

            # Process and convert the function body, including handling of custom functions like isum
            body_expr = preprocess_expression(body)
            func_expr = sp.sympify(body_expr, locals={
                **variables, **functions, **number_sets,
                **{arg: sp.Symbol(arg) for arg in args}
            })

            # Create a Lambda function with arguments as a tuple
            functions[func_name] = sp.Lambda(tuple(arg_symbols), func_expr)
            return f"Function {func_name} defined with arguments ({', '.join(args)})."
        else:
            return "Invalid function definition. Use the format: func_name :: (arg1, arg2, ...) = body"
    except Exception as e:
        return f"Error in define_function: {e}"


def compute_nroot(n, k):
    try:
        return Decimal(n) ** (Decimal(1) / Decimal(k))
    except Exception as e:
        return f"Error in compute_nroot: {e}"


def clear_all():
    try:
        variables.clear()
        matrices.clear()
        functions.clear()
        tensors.clear()
        algebraic_structures.clear()
        user_sets.clear()
        return "All mathematical objects cleared."
    except Exception as e:
        return f"Error in clear_all: {e}"

def compute_indefinite_integral(expression, var):
    try:
        processed_expression = preprocess_expression(expression)
        expr = sp.sympify(processed_expression)
        var_symbol = sp.Symbol(var)
        integral = sp.integrate(expr, var_symbol)
        return clean_multiplication(str(preprocess_expression(str(integral))))
    except Exception as e:
        return f"Error in compute_indefinite_integral: {e}"

def compute_definite_integral(expression, var, a, b):
    try:
        processed_expression = preprocess_expression(expression)
        expr = sp.sympify(processed_expression)
        var_symbol = sp.Symbol(var)
        integral = sp.integrate(expr, (var_symbol, a, b))
        return clean_multiplication(str(preprocess_expression(str(integral))))
    except Exception as e:
        return f"Error in compute_definite_integral: {e}"

def compute_higher_order_derivative(expression, var, order):
    try:
        processed_expression = preprocess_expression(expression)
        expr = sp.sympify(processed_expression)
        var_symbol = sp.Symbol(var)
        derivative = sp.diff(expr, var_symbol, order)
        return clean_multiplication(str(preprocess_expression(str(derivative))))

    except Exception as e:
        return f"Error in compute_higher_order_derivative: {e}"

def get_help():
    return (
        "Available commands:\n"
        "1. `variable := expression` - Define a variable with a given expression.\n"
        "2. `matrix_name := {matrix_str}` - Define a matrix with given values.\n"
        "3. `clear variable_name` - Clear a specific variable.\n"
        "4. `clear matrix_name` - Clear a specific matrix.\n"
        "5. `det(matrix_name)` - Compute the determinant of a matrix.\n"
        "6. `inv(matrix_name)` - Compute the inverse of a matrix.\n"
        "7. `eigenvals(matrix_name)` - Compute the eigenvalues of a matrix.\n"
        "8. `eigenvects(matrix_name)` - Compute the eigenvectors of a matrix.\n"
        "9. `transpose(matrix_name)` - Transpose a matrix.\n"
        "10. `mul(matrix_name1, matrix_name2)` - Multiply two matrices.\n"
        "11. `add(matrix_name1, matrix_name2)` - Add two matrices.\n"
        "12. `d/dx(expression)` - Compute the derivative of an expression with respect to x.\n"
        "13. `dd/ddx(expression, variable)` - Compute the partial derivative of an expression.\n"
        "14. `indefinite integral(expression, variable)` - Compute the indefinite integral of an expression.\n"
        "15. `definite integral(expression, variable, a, b)` - Compute the definite integral of an expression from a to b.\n"
        "16. `derivative(expression, variable, order)` - Compute the nth order derivative of an expression.\n"
        "17. `exit` - Exit the CLI.\n"
        "18. `clear all` - Clear all variables, matrices, and functions.\n"
        "19. `cls` - Clear the terminal screen.\n"
        "20. `help` - Display this help message\n"
        "21. Algebraic Structure Commands:\n"
        "    - `group name {elements} {operation} {identity} {inverse}` - Define a group\n"
        "    - `ring name {elements} {operations} {properties}` - Define a ring\n"
        "    - `center structure_name` - Compute center of structure\n"
        "    - `generators structure_name` - Find generators\n"
        "    - `is_cyclic structure_name` - Check if structure is cyclic\n"
        "    - `quotient structure1 structure2` - Compute quotient structure\n"
        "    - `direct_product structure1 structure2` - Compute direct product\n"
    )

def compute_infinite_series(expression, var, n):
    try:
        processed_expression = preprocess_expression(expression)
        expr = sp.sympify(processed_expression)
        var_symbol = sp.Symbol(var)
        series = sp.Sum(expr, (var_symbol, 0, n)).doit()
        return clean_multiplication(str(preprocess_expression(str(series))))
    except Exception as e:
        return f"Error in compute_infinite_series: {e}"

def compute_large_sum(expression, var, start, end):
    try:
        processed_expression = preprocess_expression(expression)
        expr = sp.sympify(processed_expression)
        var_symbol = sp.Symbol(var)
        result = sp.Sum(expr, (var_symbol, start, end)).doit()
        return clean_multiplication(str(preprocess_expression(str(result))))
    except Exception as e:
        return f"Error in compute_large_sum: {e}"

def compute_large_product(expression, var, start, end):
    try:
        processed_expression = preprocess_expression(expression)
        expr = sp.sympify(processed_expression)
        var_symbol = sp.Symbol(var)
        result = sp.Product(expr, (var_symbol, start, end)).doit()
        return clean_multiplication(str(preprocess_expression(str(result))))
    except Exception as e:
        return f"Error in compute_large_product: {e}"


def compute_sum(expression, var, start):
    try:
        processed_expression = preprocess_expression(expression)
        expr = sp.sympify(processed_expression)
        var_symbol = sp.Symbol(var)

        # Check for convergence using the ratio test
        ratio_test = sp.limit(expr.subs(var_symbol, var_symbol + 1) / expr, var_symbol, sp.oo)

        if ratio_test < 1:
            result = sp.sympify(sp.summation(expr, (var_symbol, start, sp.oo)))
            return f"The series is convergent and approaches {clean_multiplication(str(result))}."
        else:
            return "The series is divergent."
    except Exception as e:
        return f"Error in compute_sum: {e}"


def compute_product(expression, var, start):
    try:
        processed_expression = preprocess_expression(expression)
        expr = sp.sympify(processed_expression)
        var_symbol = sp.Symbol(var)

        # Check for convergence using the limit comparison test
        limit_result = sp.limit(expr, var_symbol, sp.oo)

        if limit_result != 1:
            result = sp.sympify(sp.Product(expr, (var_symbol, start, sp.oo)).doit())
            return f"The product is convergent and approaches {clean_multiplication(str(result))}."
        else:
            return "The product is divergent."
    except Exception as e:
        return f"Error in compute_product: {e}"


def define_set(name, set_definition):
    try:
        set_definition = set_definition.strip()

        if not set_definition:
            return "Set definition cannot be empty."

        # Define finite sets like {1, 2, 3}
        if set_definition.startswith('{') and set_definition.endswith('}') and "|" not in set_definition:
            elements = set_definition[1:-1].split(',')
            elements = [sp.sympify(e.strip()) for e in elements]
            user_sets[name] = sp.FiniteSet(*elements)
            return f"Set {name} defined as {user_sets[name]}."

        # Define sets with conditions like {x | x ∈ R}
        elif '|' in set_definition:
            set_parts = set_definition.split('|')
            if len(set_parts) != 2:
                return "Invalid set definition format. Use '{x | x ∈ R}'."

            variable = set_parts[0].strip()
            condition = set_parts[1].strip()

            # Replace '∈' with Python-compatible 'in'
            condition = condition.replace('∈', 'in').replace('∉', 'not in')
            variable_symbol = sp.Symbol(variable)

            # Determine the universal set based on condition
            if 'R' in condition:
                universal_set = sp.S.Reals
            elif 'Z' in condition:
                universal_set = sp.S.Integers
            elif 'Q' in condition:
                universal_set = sp.S.Rationals
            elif 'N' in condition:
                universal_set = sp.S.Naturals
            else:
                return "Unknown universal set in condition."

            # Parse the condition expression
            try:
                condition_expr = sp.sympify(condition, locals={**number_sets, **user_sets})
            except Exception as e:
                return f"Error in parsing condition expression: {e}"

            user_sets[name] = sp.ConditionSet(variable_symbol, condition_expr, universal_set)
            return f"Set {name} defined with condition: {user_sets[name]}."

        else:
            return "Invalid set definition format. Use '{1, 2, 3}' or '{x | x ∈ R}'."
    except Exception as e:
        return f"Error in define_set: {e}"

def clear_set(name):
    """Clear a set from memory"""
    try:
        if name in user_sets:
            del user_sets[name]
            return f"Set {name} cleared."
        return f"Set {name} not found."
    except Exception as e:
        return f"Error clearing set: {e}"

def define_algebraic_structure(structure_type, name, elements_str):
    """Define algebraic structures (groups, rings, fields, monoids)"""
    try:
        # Parse elements and operations from string
        params = elements_str.split()
        elements = set()
        operations = {}
        properties = {}

        for param in params:
            if param.startswith('elements='):
                elements_part = param[9:].strip('{}')
                elements = {eval(x.strip()) for x in elements_part.split(',')}
            elif param.startswith('operation='):
                op_part = param[10:].strip('{}')
                operations['*'] = eval('{' + op_part + '}')
            elif param.startswith('addition='):
                add_part = param[9:].strip('{}')
                operations['+'] = eval('{' + add_part + '}')
            elif param.startswith('identity='):
                properties['identity'] = eval(param[9:])
            elif param.startswith('inverse='):
                inv_part = param[8:].strip('{}')
                properties['inverse'] = eval('{' + inv_part + '}')

        # Validate structure based on type
        if structure_type == 'group':
            if not validate_group(elements, operations.get('*', {}), properties):
                return "Invalid group structure"
        elif structure_type == 'ring':
            if not validate_ring(elements, operations, properties):
                return "Invalid ring structure"

        algebraic_structures[name] = {
            'type': structure_type,
            'elements': elements,
            'operations': operations,
            'properties': properties
        }
        return f"{structure_type.capitalize()} {name} defined."
    except Exception as e:
        return f"Error in define_algebraic_structure: {e}"

def algebraic_operation(name, op, name2=None):
    """Handle various algebraic operations"""
    try:
        if name not in algebraic_structures:
            return "Structure not found"

        structure = algebraic_structures[name]

        if op == 'order':
            return len(structure['elements'])
        elif op == 'center':
            return compute_center(structure)
        elif op == 'is_cyclic':
            return is_cyclic_structure(structure)
        elif op == 'generators':
            return find_generators(structure)
        elif op == 'quotient' and name2:
            if name2 not in algebraic_structures:
                return "Second structure not found"
            return compute_quotient(structure, algebraic_structures[name2])
        elif op == 'direct_product' and name2:
            if name2 not in algebraic_structures:
                return "Second structure not found"
            return compute_direct_product(structure, algebraic_structures[name2])

        return f"Unknown operation {op}"
    except Exception as e:
        return f"Error in algebraic operation: {e}"

def compute_center(structure):
    """Compute the center of an algebraic structure"""
    try:
        elements = structure['elements']
        operation = structure['operations'].get('*', {})
        center = set()

        for a in elements:
            is_central = True
            for b in elements:
                if operation[a][b] != operation[b][a]:
                    is_central = False
                    break
            if is_central:
                center.add(a)

        return f"Center: {center}"
    except Exception as e:
        return f"Error computing center: {e}"

def is_cyclic_structure(structure):
    """Check if structure is cyclic and find generators"""
    try:
        elements = structure['elements']
        operation = structure['operations'].get('*', {})

        for a in elements:
            generated = {a}
            current = a
            while True:
                current = operation[current][a]
                if current in generated:
                    break
                generated.add(current)
            if generated == elements:
                return f"Structure is cyclic. Generator: {a}"

        return "Structure is not cyclic"
    except Exception as e:
        return f"Error checking cyclicity: {e}"

def find_generators(structure):
    """Find all generators of the structure"""
    try:
        elements = structure['elements']
        operation = structure['operations'].get('*', {})
        generators = set()

        for a in elements:
            generated = {a}
            current = a
            power = 1
            while True:
                current = operation[current][a]
                if current in generated:
                    break
                generated.add(current)
                power += 1
            if generated == elements:
                generators.add(a)

        return f"Generators: {generators}"
    except Exception as e:
        return f"Error finding generators: {e}"

def compute_quotient(structure1, structure2):
    """Compute quotient structure"""
    try:
        cosets = compute_cosets(structure1, structure2)
        operation = structure1['operations'].get('*', {})

        # Define operation on cosets
        quotient_operation = {}
        for coset1 in cosets:
            quotient_operation[coset1] = {}
            for coset2 in cosets:
                quotient_operation[coset1][coset2] = operate_cosets(coset1, coset2, structure1)

        return quotient_operation
    except Exception as e:
        return f"Error computing quotient: {e}"

def compute_direct_product(structure1, structure2):
    """Compute direct product of structures"""
    try:
        elements1 = structure1['elements']
        elements2 = structure2['elements']
        op1 = structure1['operations'].get('*', {})
        op2 = structure2['operations'].get('*', {})

        # Create product elements
        product_elements = {(a, b) for a in elements1 for b in elements2}

        # Define operation on product
        product_operation = {}
        for (a1, b1) in product_elements:
            product_operation[(a1, b1)] = {}
            for (a2, b2) in product_elements:
                product_operation[(a1, b1)][(a2, b2)] = (op1[a1][a2], op2[b1][b2])

        return {
            'elements': product_elements,
            'operations': {'*': product_operation}
        }
    except Exception as e:
        return f"Error computing direct product: {e}"

def is_normal_substructure(sub, structure):
    """Check if substructure is normal"""
    try:
        elements = structure['elements']
        operation = structure['operations'].get('*', {})

        for n in sub['elements']:
            for g in elements:
                conjugate = operation[operation[g][n]][inverse(g, structure)]
                if conjugate not in sub['elements']:
                    return False
        return True
    except Exception:
        return False

def compute_cosets(structure1, structure2):
    """Compute cosets for quotient structures"""
    try:
        elements = structure1['elements']
        subgroup_elements = structure2['elements']
        operation = structure1['operations'].get('*', {})
        cosets = set()

        for g in elements:
            coset = frozenset(operation[g][h] for h in subgroup_elements)
            cosets.add(coset)

        return cosets
    except Exception as e:
        return f"Error computing cosets: {e}"

def operate_cosets(coset1, coset2, structure):
    """Operate on cosets"""
    try:
        operation = structure['operations'].get('*', {})
        result = set()

        for a in coset1:
            for b in coset2:
                result.add(operation[a][b])

        return frozenset(result)
    except Exception as e:
        return f"Error operating on cosets: {e}"

def clear_algebraic_structure(name):
    """Clear an algebraic structure from memory"""
    try:
        if name in algebraic_structures:
            del algebraic_structures[name]
            return f"Algebraic structure {name} cleared."
        return f"Algebraic structure {name} not found."
    except Exception as e:
        return f"Error clearing algebraic structure: {e}"

def handle_algebraic_error(operation, error):
    """Handle errors in algebraic operations"""
    error_messages = {
        'invalid_structure': "Invalid algebraic structure",
        'not_closed': "Operation not closed",
        'not_associative': "Operation not associative",
        'no_identity': "Identity element not found",
        'no_inverse': "Inverse elements not found",
        'not_commutative': "Operation not commutative",
        'not_distributive': "Distributive property fails"
    }
    return f"Error in {operation}: {error_messages.get(str(error), str(error))}"

def validate_ring(elements, operations, properties):
    """Validate ring axioms"""
    try:
        # Check if addition forms an abelian group
        add_op = operations.get('+', {})
        if not validate_group(elements, add_op, properties):
            return False

        # Check multiplication closure and associativity
        mul_op = operations.get('*', {})
        for a in elements:
            for b in elements:
                if mul_op[a][b] not in elements:
                    return False
                for c in elements:
                    if mul_op[mul_op[a][b]][c] != mul_op[a][mul_op[b][c]]:
                        return False

        # Check distributivity
        for a in elements:
            for b in elements:
                for c in elements:
                    # Left distributivity
                    left_dist = mul_op[a][add_op[b][c]]
                    right_dist = add_op[mul_op[a][b]][mul_op[a][c]]
                    if left_dist != right_dist:
                        return False
                    # Right distributivity
                    left_dist = mul_op[add_op[b][c]][a]
                    right_dist = add_op[mul_op[b][a]][mul_op[c][a]]
                    if left_dist != right_dist:
                        return False

        return True
    except Exception:
        return False

def convert_algebraic_to_matrix(structure_name):
    """Convert algebraic structure to matrix representation"""
    try:
        if structure_name not in algebraic_structures:
            return "Structure not found"

        structure = algebraic_structures[structure_name]
        elements = list(structure['elements'])
        operation = structure['operations'].get('*', {})

        matrix_name = f"M_{structure_name}"
        matrix = []
        for i in range(len(elements)):
            row = []
            for j in range(len(elements)):
                element = operation[elements[i]][elements[j]]
                idx = elements.index(element)
                row.append(idx)
            matrix.append(row)

        matrices[matrix_name] = sp.Matrix(matrix)
        return f"Structure {structure_name} converted to matrix {matrix_name}"
    except Exception as e:
        return f"Error in conversion: {e}"

def tensor_structure_product(tensor_name, structure_name):
    """Compute tensor product with algebraic structure"""
    try:
        if tensor_name not in tensors or structure_name not in algebraic_structures:
            return "Tensor or structure not found"

        tensor = tensors[tensor_name]
        structure = algebraic_structures[structure_name]

        # Implementation of tensor-structure product
        result = compute_tensor_structure_product(tensor, structure)

        return result
    except Exception as e:
        return f"Error in tensor-structure product: {e}"

def define_tensor(name, components, rank):
    """Define a tensor with given components and rank"""
    try:
        # Parse components string into a nested list/array
        components = eval(f"[{components}]")

        # Validate tensor structure based on rank
        if not validate_tensor_structure(components, rank):
            return "Invalid tensor structure for given rank"

        tensors[name] = {
            'components': components,
            'rank': rank
        }
        return f"Tensor {name} of rank {rank} defined."
    except Exception as e:
        return f"Error in define_tensor: {e}"

def validate_tensor_structure(components, rank):
    """Validate if components match the expected tensor rank"""
    try:
        if rank == 0:  # Scalar
            return isinstance(components, (int, float))
        elif rank == 1:  # Vector
            return isinstance(components, list)
        elif rank == 2:  # Matrix
            return isinstance(components, list) and all(isinstance(row, list) for row in components)
        # Add more rank validations as needed
        return True
    except Exception:
        return False

def compute_tensor_structure_product(tensor, structure):
    """Compute the tensor product with an algebraic structure"""
    try:
        tensor_components = tensor['components']
        structure_elements = structure['elements']
        operation = structure['operations'].get('*', {})

        # Result will have rank = tensor rank + 1
        result = []

        # For each tensor component
        for component in tensor_components:
            # For each structure element
            row = []
            for element in structure_elements:
                if isinstance(component, (int, float)):
                    # Scalar case
                    row.append(operation[component][element])
                else:
                    # Higher rank case - recursive computation needed
                    row.append([operation[c][element] for c in component])
            result.append(row)

        return result
    except Exception as e:
        return f"Error in tensor-structure product: {e}"

def structure_homomorphism(structure_name, func_name):
    """Check if a function is a homomorphism for the given structure"""
    try:
        if structure_name not in algebraic_structures or func_name not in functions:
            return "Structure or function not found"

        structure = algebraic_structures[structure_name]
        func = functions[func_name]

        # Check homomorphism property for each operation
        for op_name, operation in structure['operations'].items():
            for a in structure['elements']:
                for b in structure['elements']:
                    # f(a * b) = f(a) * f(b)
                    left_side = func(operation[a][b])
                    right_side = operation[func(a)][func(b)]
                    if left_side != right_side:
                        return f"Not a homomorphism for operation {op_name}"

        return "Function is a homomorphism"
    except Exception as e:
        return f"Error checking homomorphism: {e}"

def validate_group(elements, operation, properties):
    """Validate group axioms"""
    try:
        # Check closure
        for a in elements:
            for b in elements:
                if operation[a][b] not in elements:
                    return False

        # Check associativity
        for a in elements:
            for b in elements:
                for c in elements:
                    if operation[operation[a][b]][c] != operation[a][operation[b][c]]:
                        return False

        # Check identity
        identity = properties.get('identity')
        if not identity or identity not in elements:
            return False

        # Check inverses
        inverses = properties.get('inverse', {})
        for a in elements:
            if a not in inverses or operation[a][inverses[a]] != identity:
                return False

        return True
    except Exception:
        return False

def compute_automorphisms(structure):
    """Compute automorphisms of the structure"""
    try:
        elements = structure['elements']
        operation = structure['operations'].get('*', {})
        automorphisms = []

        # Generate all possible bijective mappings
        for perm in permutations(elements):
            mapping = dict(zip(elements, perm))
            if is_homomorphism(mapping, structure, structure):
                automorphisms.append(mapping)

        return f"Automorphisms: {automorphisms}"
    except Exception as e:
        return f"Error computing automorphisms: {e}"

def is_homomorphism(mapping, structure1, structure2):
    """Check if mapping is a homomorphism between structures"""
    try:
        # Check if mapping preserves operations
        for op_name, operation1 in structure1['operations'].items():
            operation2 = structure2['operations'].get(op_name, {})

            for a in structure1['elements']:
                for b in structure1['elements']:
                    if a not in mapping or b not in mapping:
                        return False
                    # Check if f(a * b) = f(a) * f(b)
                    left_side = mapping[operation1[a][b]]
                    right_side = operation2[mapping[a]][mapping[b]]
                    if left_side != right_side:
                        return False
        return True
    except Exception:
        return False

def clear_tensor(name):
    """Clear a tensor from memory"""
    try:
        if name in tensors:
            del tensors[name]
            return f"Tensor {name} cleared."
        return f"Tensor {name} not found."
    except Exception as e:
        return f"Error clearing tensor: {e}"

def inverse(element, structure):
    """Find the inverse of an element in a structure"""
    try:
        # Get the inverse from structure properties if available
        if 'properties' in structure and 'inverse' in structure['properties']:
            inverses = structure['properties']['inverse']
            if element in inverses:
                return inverses[element]

        # If not in properties, compute it
        elements = structure['elements']
        operation = structure['operations'].get('*', {})
        identity = structure['properties'].get('identity')

        # If no identity element, can't find inverse
        if not identity:
            return None

        # Search for inverse by checking operation result
        for e in elements:
            if (operation[element][e] == identity and
                operation[e][element] == identity):
                return e

        return None  # No inverse found
    except Exception as e:
        return f"Error finding inverse: {e}"

def detect_and_solve_equation(equation, conditions=None):
    """Automatically detect and solve differential equations with conditions"""
    try:
        # Parse equation
        if '=' in equation:
            lhs, rhs = equation.split('=')
            lhs, rhs = lhs.strip(), rhs.strip()
        else:
            lhs, rhs = equation.strip(), '0'

        # Analyze equation
        eq_type, variables = analyze_equation(lhs + '=' + rhs)

        # Parse conditions if provided
        parsed_conditions = []
        if conditions:
            for condition in conditions:
                if '=' in condition:
                    left, right = condition.split('=')
                    left = left.strip()
                    right = right.strip()

                    if '(' in left:
                        # Handle multi-variable conditions like y(x,0) = -x
                        func_name = left[:left.index('(')]
                        point_str = left[left.index('(')+1:left.index(')')]
                        points = [p.strip() for p in point_str.split(',')]

                        parsed_conditions.append({
                            'type': 'value',
                            'points': points,
                            'value': right
                        })
                    elif "'" in left:
                        order = left.count("'")
                        point = 0
                        parsed_conditions.append({
                            'type': 'derivative',
                            'order': order,
                            'point': point,
                            'value': float(right)
                        })

        # Force PDE if multiple variables detected in derivatives or conditions
        if any(',' in str(cond.get('points', '')) for cond in parsed_conditions):
            eq_type = 'pde'
        elif variables['derivatives'] and len(set(d[1] for d in variables['derivatives'])) > 1:
            eq_type = 'pde'

        # Solve based on type
        if eq_type == 'pde':
            return solve_pde(lhs, rhs, variables, parsed_conditions)
        else:
            return solve_ode(lhs, rhs, variables, parsed_conditions)

    except Exception as e:
        return f"Error in equation detection: {str(e)}"

def analyze_equation(equation):
    """Analyze equation to determine type and variables"""
    try:
        # Remove spaces
        eq = equation.replace(' ', '')

        # Initialize variables dictionary
        variables = {
            'dependent': None,
            'independent': set(),
            'derivatives': [],
            'functions': set()
        }

        # Find all derivatives
        derivatives = re.findall(r'd([a-zA-Z][a-zA-Z0-9]*)/d([a-zA-Z][a-zA-Z0-9]*)', eq)
        primed = re.findall(r'([a-zA-Z][a-zA-Z0-9]*)[\'"]+', eq)

        if derivatives:
            # Get variables from derivatives
            dep_vars = {d[0] for d in derivatives}
            indep_vars = {d[1] for d in derivatives}

            variables['dependent'] = list(dep_vars)[0]
            variables['independent'] = indep_vars
            variables['derivatives'] = derivatives

            # It's a PDE only if we have multiple DIFFERENT independent variables
            return ('pde' if len(indep_vars) > 1 else 'ode', variables)

        elif primed:
            # Handle primed notation (always ODE)
            variables['dependent'] = primed[0]
            all_vars = set(re.findall(r'[a-zA-Z][a-zA-Z0-9]*', eq)) - {'d'} - {variables['dependent']}
            variables['independent'] = {next(iter(all_vars))} if all_vars else {'x'}
            return ('ode', variables)

        else:
            # Default case
            all_vars = list(set(re.findall(r'[a-zA-Z][a-zA-Z0-9]*', eq)) - {'d'})
            if all_vars:
                variables['dependent'] = all_vars[0]
                variables['independent'] = {all_vars[1]} if len(all_vars) > 1 else {'x'}
            else:
                variables['dependent'] = 'y'
                variables['independent'] = {'x'}
            return ('ode', variables)

    except Exception as e:
        raise Exception(f"Error in analysis: {e}")

def solve_ode(lhs, rhs, variables, conditions=None):
    """Solve ordinary differential equations with conditions"""
    try:
        # Create symbols
        x = sp.Symbol(list(variables['independent'])[0])
        y = sp.Function(variables['dependent'])

        # Process equation
        expr = lhs
        dep_var = variables['dependent']
        indep_var = list(variables['independent'])[0]

        # Create local dictionary
        local_dict = {
            indep_var: x,
            dep_var: y
        }

        # Replace derivatives
        if variables['derivatives']:
            expr = expr.replace(f'd{dep_var}/d{indep_var}',
                              str(sp.Derivative(y(x), x)))

        # Replace function calls
        pattern = rf'\b{dep_var}\b(?!\()'
        expr = re.sub(pattern, f'{dep_var}({indep_var})', expr)

        # Create equation
        lhs_expr = sp.sympify(expr, locals=local_dict)
        rhs_expr = sp.sympify(rhs, locals=local_dict)
        eq = sp.Eq(lhs_expr, rhs_expr)

        # Solve ODE
        if conditions:
            ics = {}
            for cond in conditions:
                if cond['type'] == 'value':
                    ics[y(cond['point'])] = cond['value']
                elif cond['type'] == 'derivative':
                    ics[sp.Derivative(y(x), x, cond['order']).subs(x, cond['point'])] = cond['value']
            solution = sp.dsolve(eq, y(x), ics=ics)
        else:
            solution = sp.dsolve(eq, y(x))

        return f"ODE Solution: {solution}"

    except Exception as e:
        return f"Error in solve_ode: {str(e)}\nEquation parsed as: {expr}"

def solve_pde(lhs, rhs, variables, conditions=None):
    """Solve partial differential equations with conditions"""
    try:
        # Create symbols for all independent variables
        var_symbols = {var: sp.Symbol(var) for var in variables['independent']}
        dep_var = variables['dependent']
        u = sp.Function(dep_var)  # Use the actual dependent variable name
        F = sp.Function('F')  # Arbitrary function

        # Create local dictionary with all symbols
        local_dict = {
            **var_symbols,  # All independent variables
            dep_var: u,
            'F': F
        }

        # Process equation
        expr = lhs

        # Replace partial derivatives
        for d_var, i_var in variables['derivatives']:
            expr = expr.replace(
                f'd{d_var}/d{i_var}',
                str(sp.Derivative(u(*var_symbols.values()), var_symbols[i_var]))
            )

        # Replace function calls
        pattern = rf'\b{dep_var}\b(?!\()'
        expr = re.sub(pattern, f'{dep_var}({",".join(var_symbols.keys())})', expr)

        # Create equation
        lhs_expr = sp.sympify(expr, locals=local_dict)
        rhs_expr = sp.sympify(rhs, locals=local_dict)
        eq = sp.Eq(lhs_expr, rhs_expr)

        # Get general solution
        general = sp.pdsolve(eq, u(*var_symbols.values()))

        # If we have conditions, solve for F
        if conditions:
            # Process all conditions
            condition_eqs = []
            for cond in conditions:
                if cond['type'] == 'value':
                    try:
                        # Parse condition points and value
                        if isinstance(cond['points'], list):
                            # Full point specification y(x,t) = value
                            points = {var: sp.sympify(val, locals=local_dict)
                                    for var, val in zip(var_symbols.keys(), cond['points'])}
                        else:
                            # Single variable specification y(x) = f(t) or y(t) = f(x)
                            var = cond['points']
                            other_vars = set(var_symbols.keys()) - {var}
                            points = {v: var_symbols[v] for v in other_vars}
                            points[var] = 0  # Set specified variable to 0

                        # Parse value expression
                        value_expr = sp.sympify(cond['value'], locals=local_dict)

                        # Create condition equation
                        cond_eq = sp.Eq(general.rhs.subs(points), value_expr)
                        condition_eqs.append((cond_eq, points, value_expr))

                    except Exception as e:
                        print(f"Warning: Could not process condition: {e}")

            # Try to solve the system of conditions
            if condition_eqs:
                try:
                    # Get characteristic variable
                    char_var = next(arg for arg in general.rhs.find(F)
                                  if isinstance(arg, sp.Function)).args[0]

                    # Use all conditions to determine F
                    particular = general.rhs
                    for cond_eq, points, value_expr in condition_eqs:
                        char_val = char_var.subs(points)

                        # Update particular solution
                        particular = particular.subs(
                            F(char_var),
                            value_expr + F(char_var) - F(char_val)
                        )

                    return f"PDE Solution: Eq({dep_var}({','.join(var_symbols.keys())}), {particular})"

                except Exception as e:
                    print(f"Warning: Could not solve for particular solution: {e}")

        return f"PDE Solution: {general}"

    except Exception as e:
        return f"Error in solve_pde: {str(e)}\nEquation parsed as: {expr}"

def compute_roots(expression, variable):
    try:
        processed_expression = preprocess_expression(expression)
        var_symbol = sp.Symbol(variable)
        expr = sp.sympify(processed_expression, locals={**variables, **functions})

        solutions = sp.solve(expr, var_symbol)

        if solutions:
            cleaned_solutions = [clean_multiplication(str(sol)) for sol in solutions]
            return f"Roots: {', '.join(cleaned_solutions)}"
        else:
            return "No roots found."

    except Exception as e:
        return f"Error in compute_roots: {e}"

def process_commands():
    while True:
        try:
            user_input = input("Cardano CLI > ")

            # Basic commands
            if user_input.lower() == 'exit':
                break
            elif user_input.lower() == 'cls':
                print(clear_screen())
            elif user_input.lower() == 'help':
                print(get_help())
            elif user_input.lower() == 'clear all':
                print(clear_all())

            # Clear specific objects
            elif user_input.lower().startswith('clear '):
                name = user_input.split(' ')[1]
                if name in variables:
                    print(clear_variable(name))
                elif name in matrices:
                    print(clear_matrix(name))
                elif name in user_sets:
                    print(clear_set(name))
                elif name in algebraic_structures:
                    print(clear_algebraic_structure(name))
                elif name in tensors:
                    print(clear_tensor(name))
                else:
                    print(f"{name} not found.")

            # Matrix operations
            elif user_input.startswith('det('):
                matrix_name = user_input[4:-1]
                print(compute_determinant(matrix_name))
            elif user_input.startswith('inv('):
                matrix_name = user_input[4:-1]
                print(compute_inverse(matrix_name))
            elif user_input.startswith('eigenvals('):
                matrix_name = user_input[10:-1]
                print(compute_eigenvalues(matrix_name))
            elif user_input.startswith('eigenvects('):
                matrix_name = user_input[11:-1]
                print(compute_eigenvectors(matrix_name))
            elif user_input.startswith('transpose('):
                matrix_name = user_input[10:-1]
                print(transpose_matrix(matrix_name))
            elif user_input.startswith('mul('):
                matrices = user_input[4:-1].split(',')
                print(multiply_matrices(matrices[0].strip(), matrices[1].strip()))
            elif user_input.startswith('add('):
                matrices = user_input[4:-1].split(',')
                print(add_matrices(matrices[0].strip(), matrices[1].strip()))

            # Differential equations
            elif user_input.startswith('solve ode ') or user_input.startswith('solve pde '):
               parts = user_input.split(',')
               equation = parts[0].replace('solve ode ', '').replace('solve pde ', '').strip()
               conditions = [cond.strip() for cond in parts[1:]] if len(parts) > 1 else None
               print(detect_and_solve_equation(equation, conditions))

            # Calculus operations
            elif user_input.startswith('d/dx('):
                expr = user_input[5:-1]
                print(compute_derivative(expr, 'x'))
            elif user_input.startswith('dd/ddx('):
                parts = user_input[7:-1].split(',')
                print(compute_partial_derivative(parts[0], parts[1]))
            elif user_input.startswith('indefinite integral('):
                parts = user_input[19:-1].split(',')
                print(compute_indefinite_integral(parts[0], parts[1]))
            elif user_input.startswith('definite integral('):
                parts = user_input[17:-1].split(',')
                print(compute_definite_integral(parts[0], parts[1], parts[2], parts[3]))
            elif user_input.startswith('derivative('):
                parts = user_input[11:-1].split(',')
                print(compute_higher_order_derivative(parts[0], parts[1], int(parts[2])))

            # Series and sums
            elif user_input.startswith('series('):
                parts = user_input[7:-1].split(',')
                print(compute_infinite_series(parts[0], parts[1], int(parts[2])))
            elif user_input.startswith('sum('):
                parts = user_input[4:-1].split(',')
                print(compute_large_sum(parts[0], parts[1], int(parts[2]), int(parts[3])))
            elif user_input.startswith('product('):
                parts = user_input[8:-1].split(',')
                print(compute_large_product(parts[0], parts[1], int(parts[2]), int(parts[3])))
            elif user_input.startswith('isum('):
                parts = user_input[5:-1].split(',')
                print(compute_sum(parts[0].strip(), parts[1].strip(), int(parts[2].strip())))
            elif user_input.startswith('iproduct('):
                parts = user_input[9:-1].split(',')
                print(compute_product(parts[0].strip(), parts[1].strip(), int(parts[2].strip())))

            # Algebraic structures and operations
            elif user_input.startswith('center '):
                name = user_input.split()[1]
                print(compute_center(algebraic_structures.get(name)))
            elif user_input.startswith('generators '):
                name = user_input.split()[1]
                print(find_generators(algebraic_structures.get(name)))
            elif user_input.startswith('is_cyclic '):
                name = user_input.split()[1]
                print(is_cyclic_structure(algebraic_structures.get(name)))
            elif user_input.startswith('quotient '):
                name1, name2 = user_input.split()[1:]
                print(algebraic_operation(name1, 'quotient', name2))
            elif user_input.startswith('direct_product '):
                name1, name2 = user_input.split()[1:]
                print(algebraic_operation(name1, 'direct_product', name2))
            elif user_input.startswith('convert_matrix '):
                name = user_input.split()[1]
                print(convert_algebraic_to_matrix(name))
            elif user_input.startswith('tensor_product '):
                tensor_name, structure_name = user_input.split()[1:]
                print(tensor_structure_product(tensor_name, structure_name))
            elif user_input.startswith('homomorphism '):
                structure_name, func_name = user_input.split()[1:]
                print(structure_homomorphism(structure_name, func_name))

            # Set operations and definitions
            elif ':=' in user_input and '{' in user_input:
                name, set_def = user_input.split(':=')
                if 'tensor' in set_def.lower():
                    pattern = r'tensor\s*\{([^}]+)\}\s*rank\s*(\d+)'
                    match = re.match(pattern, set_def.strip())
                    if match:
                        components, rank = match.groups()
                        print(define_tensor(name.strip(), components, int(rank)))
                else:
                    print(define_set(name.strip(), set_def.strip()))

            # Function definitions
            elif '::' in user_input:
                print(define_function(user_input))

            # Algebraic structure definitions
            elif any(s in user_input for s in ['group', 'ring', 'field', 'monoid']):
                structure_type = next(s for s in ['group', 'ring', 'field', 'monoid']
                                   if s in user_input)
                match = re.match(
                    f'({structure_type})\s+([A-Za-z][A-Za-z0-9]*)\s+(.*)',
                    user_input,
                    re.IGNORECASE
                )
                if match:
                    _, name, params = match.groups()
                    print(define_algebraic_structure(structure_type, name, params))

            # Root computation
            elif user_input.startswith('nroot('):
                parts = user_input[6:-1].split(',')
                print(compute_nroot(parts[0].strip(), parts[1].strip()))

            # Equation solving
            elif '=' in user_input and 'solve' not in user_input:
                print(solve_equation(user_input))

            # System of equations
            elif user_input.startswith('solve '):
                equations = [eq.strip() for eq in user_input[6:].split(',')]
                if len(equations) > 1:
                    print(solve_system_of_equations(equations))
                else:
                    print(solve_equation(equations[0]))

            # Variable assignments
            elif ':=' in user_input and '::' not in user_input:
                name, expr = user_input.split(':=')
                if '{' in expr and '}' in expr:  # Matrix definition
                    print(define_matrix(name.strip(), expr.strip()))
                else:  # Variable definition
                    print(define_variable(name.strip(), expr.strip()))


            elif 'root' in user_input:
                match = re.match(r'root\((.*),\s*(\w+)\)', user_input)
                if match:
                    expression, variable = match.groups()
                    print(compute_roots(expression, variable))

            # Default expression evaluation
            else:
                print(evaluate_expression(user_input))

        except Exception as e:
            print(f"Error in command processing: {e}")

process_commands()
