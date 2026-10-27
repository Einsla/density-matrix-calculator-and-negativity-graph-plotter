import sympy as sp
import computation_mathematica as cm
from sympy.printing.str import sstr

# global variable
IJ_value = cm.loadIJ('MachinePrecision')
tau1, tau2 = sp.symbols('𝜏1 𝜏2', real=True)
O, l ,Ir, R = sp.symbols('Ω λ Λ R', real=True)
pi, i = sp.pi, sp.I
I, J = sp.symbols('I J')
A = cm.Detector('A', l, O)
B = cm.Detector('B', l, O)

class Term:
    def __init__(self, term_sign, integral_type, matrix, wightman, pairing):
        self.sign, self.coeff_euler, self.coeff_polar = encoder(matrix) # these variables are used to store magnitude of product of monopole, ui and p0
        self.coeff_euler *= term_sign
        self.unknown = IJ_encoder(integral_type, pairing, self.sign, wightman)
        self.wightman = wightman

    # packing term into form of {integral type}{detector1}{detector2}{monopole1}{monopole2}, eg JABpp
    def t_encoder(self, form="euler"):
        if form == "euler":
            expr = self.coeff_euler.multiply_elementwise(self.unknown)
        elif form == "polar":
            expr = self.coeff_polar.multiply_elementwise(self.unknown)
        else:
            raise ValueError("invalid form")
        return expr

    # evaluation of each term
    def evaluate(self):
        # create a substitution dictionary
        subs_dict = {
            sym: IJ_value[str(sym)]
            for sym in self.unknown.free_symbols
            if str(sym) in IJ_value
        }
        term_coupling = self.wightman.detector1.coupling_strength * self.wightman.detector2.coupling_strength
        term_coupling = term_coupling.subs(l, cm.lc)
        unknown = self.unknown.subs(subs_dict)
        expr = self.coeff_euler.multiply_elementwise(unknown)
        expr *= term_coupling
        return expr

# update IJ_value dictionary
def set_parameter(time_step, Rc, Precision='MachinePrecision'):
    global IJ_value
    cm.time_step = time_step
    cm.Rc = Rc
    IJ_value = cm.loadIJ(Precision)

# parameter adjust and retrieval
def get_parameter():
    return cm.time_step, cm.Rc

# detector monopole operator
def u(detector, time):
    operator = sp.Matrix.zeros(4, 4)
    if detector.label == 'A':
        operator[2, 0] = sp.exp(+1 * i * detector.Energy_gap * time)
        operator[3, 1] = sp.exp(+1 * i * detector.Energy_gap * time)
        operator[0, 2] = sp.exp(-1 * i * detector.Energy_gap * time)
        operator[1, 3] = sp.exp(-1 * i * detector.Energy_gap * time)
    elif detector.label == 'B':
        operator[1, 0] = sp.exp(+1 * i * detector.Energy_gap * time)
        operator[3, 2] = sp.exp(+1 * i * detector.Energy_gap * time)
        operator[0, 1] = sp.exp(-1 * i * detector.Energy_gap * time)
        operator[2, 3] = sp.exp(-1 * i * detector.Energy_gap * time)
    else:
        print("invalid detector")
        exit(0)
    return operator

# density matrix calculation
def partial_density(first, second, p0, form="euler"):
    pairing = first.label + second.label
    u1 = u(first, tau1)
    u2 = u(second, tau2)
    if form == "euler":
        p0 = euler_format(p0)

    # matrix calculation
    uup = (u1 * u2) * p0
    upu = (u1 * p0) * u2.H
    puu = (p0 * u2.H) * u1.H

    # Term assignment
    UUp = Term(-1, J, uup, cm.Wightman(first, second, False), pairing)
    UpU = Term(+1, I, upu, cm.Wightman(first, second, True), pairing)
    pUU = Term(-1, J, puu, cm.Wightman(first, second, True), pairing)

    # pack into container list
    pf = [UUp, UpU, pUU]
    return pf

# total 2nd order density matrix
def density_change(p0, evaluate=True, form='euler'):
    # couples all terms together in a single list
    pAA = partial_density(A, A, p0, form)
    pAB = partial_density(A, B, p0, form)
    pBA = partial_density(B, A, p0, form)
    pBB = partial_density(B, B, p0, form)
    pd2 = pAA + pAB + pBA + pBB

    # evaluate numerically or symbolically
    rhoT2 = sp.Matrix.zeros(4, 4)
    for term in pd2:
        if evaluate:
            rhoT2 += term.evaluate()
            rhoT2 = rhoT2.evalf()
        else:
            rhoT2 += term.t_encoder(form)
    return rhoT2

# encode matrix term(product of ui and p0)
def encoder(matrix):
    # Wildcards
    a = sp.Wild('a', exclude=[tau1, tau2, O, i], properties=[lambda x: True])
    b = sp.Wild('b', exclude=[tau1, tau2, O, i], properties=[lambda x: True])
    s1 = sp.Wild('s1', exclude=[tau2, tau1, O])
    s2 = sp.Wild('s2', exclude=[tau1, tau2, O])

    # initialization
    rows, cols = matrix.shape
    sign_matrix = [[['x','x'] for _ in range(cols)] for _ in range(rows)]
    euler_coeff_matrix = sp.zeros(rows, cols)
    polar_coeff_matrix = sp.zeros(rows, cols)

    for im in range(rows):
        for jm in range(cols):
            expr = matrix[im, jm]
            factors = expr.args if expr.is_Mul else [expr]

            coeff = 1  # everything not containing tau1/tau2
            exp_factors = []  # everything containing tau1/tau2

            # coefficient matching
            for f in factors:
                if f.has(tau1) or f.has(tau2):
                    exp_factors.append(f)
                else:
                    coeff *= f
            coeff_e = sp.expand_complex(coeff)
            coeff_p = coeff

            # sign matching
            sign1, sign2 = "=", "="
            for exp_term in exp_factors:
                arg = exp_term.args[0]  # exponent
                if arg.has(tau1):
                    if arg.has(-1):
                        sign1 = -1
                    else:
                        sign1 = 1
                if arg.has(tau2):
                    if arg.has(-1):
                        sign2 = -1
                    else:
                        sign2 = 1

            sign1 = cm.sign_symbol[sign1] if sign1 in [1, -1] else "="
            sign2 = cm.sign_symbol[sign2] if sign2 in [1, -1] else "="
            s12 = [sign1, sign2]

            sign_matrix[im][jm] = s12
            euler_coeff_matrix[im, jm] = coeff_e
            polar_coeff_matrix[im, jm] = coeff_p
    return sign_matrix, euler_coeff_matrix, polar_coeff_matrix

# reformat density matrix to polar form
def euler_format(matrix):
    matrix_euler = matrix.applyfunc(sp.expand_complex)
    matrix_euler = matrix_euler.applyfunc(lambda x: x.evalf())
    return matrix_euler

def ket_to_density(a, v):
    k0 = sp.Matrix([a[j] * sp.exp(i * v[j]) for j in range(4)])
    p0 = k0 * k0.H

    return p0

# encode the integral pairing and sign together
def IJ_encoder(integral_type, pairing, sign_matrix, wightman):
    row = len(sign_matrix)
    col = len(sign_matrix[0])
    expr = sp.zeros(row, col)
    sw = cm.m if wightman.flip else cm.p
    for im in range(row):
        for jm in range(col):
            sign = list(sign_matrix[im][jm])
            ref = sp.Symbol(f"{integral_type}{sw}{pairing}{sign[0]}{sign[1]}", commutative=False)
            expr[im, jm] = ref
    return expr

# close kernel for cm
def close():
    if cm == cm:
        cm.close()

# for display testing purpose
# format numerical value in expression to n sf
def round_sf(expr, sf=4):
    return sstr(
        expr.xreplace({
            n: sp.Float(n, sf)
            for n in expr.atoms(sp.Float)
        })
    )

# display equation in encoder form
def display(matrix_sum,interface=False):
    matrix_expr = sp.Matrix.zeros(4, 4)

    for matrix in matrix_sum:
        matrix_expr += matrix.t_encoder()
    # matrix_expr = sp.expand_complex(matrix_expr)

    cell_width = 80
    if interface:
        print("{")
    for im in range(0, matrix_expr.rows):
        row_str = ""
        for jm in range(0, matrix_expr.cols):
            cell = round_sf(matrix_expr[im, jm])
            mathematica_char = {"(":"[", ")":"]", "exp":"Exp", "e-":"*10^-", "e+":"*10^+"}
            if interface:
                for old, new in mathematica_char.items():
                    cell=cell.replace(old, new)
            cell = cell.replace("*", " * ")
            if interface and jm != matrix_expr.cols - 1:
                cell += ","
            row_str += f"{cell:^{cell_width}}"
        if interface:
            row_str =f"{{{row_str}}}"
        print(row_str)
    if interface:
        print("}")

# display matrix product
def matrix_display(matrix):
    cell_width = 20
    sign_matrix, coeff_matrix_euler, coeff_matrix_polar = encoder(matrix)
    for im in range(matrix.rows):
        row_str = ""
        for jm in range(matrix.cols):
            interaction_phase = sp.Integer(0)
            if "=" not in sign_matrix[im][jm]:
                t = [tau1, tau2]
                for n in range(2):
                    sn = cm.symbol_sign[sign_matrix[im][jm][n]]
                    ipn = sn * i * O * t[n]
                    interaction_phase += ipn
            expr = round_sf(coeff_matrix_euler[im, jm]) * sp.exp(interaction_phase)
            cell = f"{expr}"
            row_str += f"{cell:^{cell_width}}"
        print(row_str)