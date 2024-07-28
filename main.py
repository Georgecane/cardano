import sympy as sp
import re
import os
import platform
from decimal import Decimal, getcontext

getcontext().prec = 30

variables = {}
matrices = {}
functions = {}

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

    processed_expression = re.sub(r'nroot\((\d+),\s*(\d+)\)', r'compute_nroot(\1, \2)', processed_expression)

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
    return expression

def evaluate_expression(expression):
    try:
        processed_expression = preprocess_expression(expression)
        print(f"Processed Expression: {processed_expression}")

        eval_dict = {'sp': sp, 'Decimal': Decimal, 'compute_nroot': compute_nroot,
                     'compute_determinant': compute_determinant,
                     'compute_inverse': compute_inverse, 'compute_eigenvalues': compute_eigenvalues,
                     'compute_eigenvectors': compute_eigenvectors, 'multiply_matrices': multiply_matrices,
                     'add_matrices': add_matrices, 'transpose_matrix': transpose_matrix,
                     **variables, **functions}

        # Define symbols
        symbols = re.findall(r'\b[a-zA-Z_]\w*\b', processed_expression)
        for symbol in symbols:
            if symbol not in variables and symbol not in functions:
                variables[symbol] = sp.Symbol(symbol)

        # Check for function calls
        match = re.match(r'(\w+)\((.*)\)', processed_expression)
        if match:
            func_name = match.group(1)
            func_arg = match.group(2)
            if func_name in functions:
                func = functions[func_name]
                evaluated_arg = sp.sympify(func_arg, locals=eval_dict)
                evaluated_expression = func(evaluated_arg)
            else:
                return f"Function {func_name} is not defined."
        else:
            evaluated_expression = sp.sympify(processed_expression, locals=eval_dict)

        expanded_expression = sp.expand(evaluated_expression)
        evaluated_numeric_expression = expanded_expression.evalf() if isinstance(expanded_expression, sp.Basic) else expanded_expression
        cleaned_numeric_result = clean_multiplication(str(evaluated_numeric_expression))

        print(f"Expanded Expression: {expanded_expression}")
        print(f"Numeric Evaluation: {cleaned_numeric_result}")
        expanded_expression = preprocess_expression(str(expanded_expression))
        return clean_multiplication(str(expanded_expression))

    except Exception as e:
        return f"Error in evaluate_expression: {e}"

def solve_equation(equation):
    try:
        # تنها معادلات ریاضی را پردازش کنید
        if '::' in equation or ':=' in equation:
            return "Invalid equation format. Use '=' for equations."

        lhs, rhs = equation.split('=')
        lhs_sym = preprocess_expression(lhs.strip())
        rhs_sym = preprocess_expression(rhs.strip())

        lhs_expr = sp.sympify(lhs_sym, locals={**variables, **functions})
        rhs_expr = sp.sympify(rhs_sym, locals={**variables, **functions})

        solutions = sp.solve(lhs_expr - rhs_expr)
        return f"Solutions: {solutions}"
    except Exception as e:
        return f"Error in solve_equation: {e}"

def solve_system_of_equations(equations):
    try:
        eq_list = [eq.strip() for eq in equations.strip('{}').split(',')]
        sympy_eqs = []
        variables_set = set()
        for eq in eq_list:
            if eq.count('=') != 1:
                return "Each equation must contain exactly one '=' sign."
            processed_eq = preprocess_expression(eq)
            lhs, rhs = processed_eq.split('=')
            lhs = sp.sympify(lhs, locals={**variables, **matrices})
            rhs = sp.sympify(rhs, locals={**variables, **matrices})
            sympy_eqs.append(sp.Eq(lhs, rhs))
            variables_set.update(lhs.free_symbols)
            variables_set.update(rhs.free_symbols)
        solutions = sp.solve(sympy_eqs, list(variables_set))
        if solutions:
            cleaned_solutions = {}
            for var, sol in solutions.items():
                cleaned_solutions[str(var)] = clean_multiplication(str(sol))
            return cleaned_solutions
        else:
            return "No solutions found."
    except Exception as e:
        return f"Error in solve_system_of_equations: {e}"

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
        matrix_str = matrix_str.strip()[1:-1]
        rows = matrix_str.split(', ')
        matrix = []
        for row in rows:
            elements = row.split()
            matrix.append([sp.sympify(e) for e in elements])
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
        match = re.match(r'(\w+)\s*::\s*(\w+)\s*=\s*(.*)', func_definition)
        if match:
            func_name = match.group(1)
            arg = match.group(2)
            body = match.group(3)

            arg_symbol = sp.Symbol(arg.strip())

            body_expr = preprocess_expression(body)
            func_expr = sp.sympify(body_expr, locals={**variables, **functions, arg: arg_symbol})

            functions[func_name] = sp.Lambda(arg_symbol, func_expr)
            return f"Function {func_name} defined."
        else:
            return "Invalid function definition. The correct format is: func_name :: arg = body"
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
        return "All variables, matrices, and functions cleared."
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
        "19. `help` - Display this help message\n"

    )

while True:
    try:
        user_input = input("Cardano CLI > ")

        if user_input.lower() == 'exit':
            break
        elif user_input.lower() == 'clear all':
            print(clear_all())
        elif user_input.lower().startswith('clear '):
            name = user_input.split(' ')[1]
            if name in variables:
                print(clear_variable(name))
            elif name in matrices:
                print(clear_matrix(name))
            else:
                print(f"{name} does not exist.")

        elif user_input.lower() == "help":
            print(get_help())

        elif '=' in user_input and '::' not in user_input and ':=' not in user_input:
            print(solve_equation(user_input))
        elif '{' in user_input and '}' in user_input:
            print(solve_system_of_equations(user_input))
        elif ':=' in user_input:
            print(define_variable(user_input))
        elif '::' in user_input:
            print(define_function(user_input))
        elif 'det' in user_input:
            match = re.match(r'det\((\w+)\)', user_input)
            if match:
                matrix_name = match.group(1)
                print(compute_determinant(matrix_name))
        elif 'inv' in user_input:
            match = re.match(r'inv\((\w+)\)', user_input)
            if match:
                matrix_name = match.group(1)
                print(compute_inverse(matrix_name))
        elif 'eigenvals' in user_input:
            match = re.match(r'eigenvals\((\w+)\)', user_input)
            if match:
                matrix_name = match.group(1)
                print(compute_eigenvalues(matrix_name))
        elif 'eigenvects' in user_input:
            match = re.match(r'eigenvects\((\w+)\)', user_input)
            if match:
                matrix_name = match.group(1)
                print(compute_eigenvectors(matrix_name))
        elif 'transpose' in user_input:
            match = re.match(r'transpose\((\w+)\)', user_input)
            if match:
                matrix_name = match.group(1)
                print(transpose_matrix(matrix_name))
        elif 'mul' in user_input:
            match = re.match(r'mul\((\w+),\s*(\w+)\)', user_input)
            if match:
                matrix_name1, matrix_name2 = match.groups()
                print(multiply_matrices(matrix_name1, matrix_name2))
        elif 'add' in user_input:
            match = re.match(r'add\((\w+),\s*(\w+)\)', user_input)
            if match:
                matrix_name1, matrix_name2 = match.groups()
                print(add_matrices(matrix_name1, matrix_name2))
        elif '{' in user_input and '}' in user_input:
            matrix_match = re.match(r'(\w+)\s*:=\s*{\s*([^}]+)\s*}', user_input)
            if matrix_match:
                matrix_name, matrix_str = matrix_match.groups()
                print(define_matrix(matrix_name, matrix_str))
        elif 'd/dx' in user_input:
            match = re.match(r'd/dx\((.*)\)', user_input)
            if match:
                expression = match.group(1)
                print(compute_derivative(expression, 'x'))
        elif 'dd/ddx' in user_input:
            match = re.match(r'dd/ddx\((.*),\s*(.*)\)', user_input)
            if match:
                expression, var = match.groups()
                print(compute_partial_derivative(expression, var))
        elif 'indefinite integral' in user_input:
            match = re.match(r'indefinite integral\((.*),\s*(.*)\)', user_input)
            if match:
                expression, var = match.groups()
                print(compute_indefinite_integral(expression, var))
        elif 'definite integral' in user_input:
            match = re.match(r'definite integral\((.*),\s*(.*),\s*(.*),\s*(.*)\)', user_input)
            if match:
                expression, var, a, b = match.groups()
                print(compute_definite_integral(expression, var, a, b))
        elif 'derivative' in user_input:
            match = re.match(r'derivative\((.*),\s*(.*),\s*(\d+)\)', user_input)
            if match:
                expression, var, order = match.groups()
                print(compute_higher_order_derivative(expression, var, int(order)))
        else:
            print(evaluate_expression(user_input))
    except Exception as e:
        print(f"Error in main loop: {e}")