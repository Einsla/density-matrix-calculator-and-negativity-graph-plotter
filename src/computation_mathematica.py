import sympy as sp
import itertools
import mathematica_coding as mc
import time as t
import pickle
from pathlib import Path

# constant value
Oc = 1.0
lc = 0.01
Irc = 0.01
Rc = 4.0
time_step = 0.5
tau0c = 0

# sympy symbol declaration
p, m= sp.symbols('p m')
O, l ,Ir, R = sp.symbols('Ω λ Λ R', real=True)
E = sp.symbols('ε', positive=True)
tau, tau0, tau1, tau2 = sp.symbols('𝜏 𝜏0 𝜏1 𝜏2', real=True)
pi, i, infty = sp.pi, sp.I, sp.oo
I, J = sp.symbols('I J')

# dictionary definition
integral_boundary = {
    J: [[tau0, tau],[tau0, tau1]],
    I: [[tau0, tau],[tau0, tau]]}
sign_symbol = {+1: p, -1: m}
symbol_sign = {v: k for k, v in sign_symbol.items()}  # backward retrieval of sign_symbol

# define parent file directory
project_root = Path(__file__).resolve().parent.parent
data_folder = project_root / "data"
data_folder.mkdir(exist_ok=True)

# starting kernel
mc.warm_up()

class Detector:
    def __init__(self, label, coupling_strength=l, Energy_gap=O):
        if label == 'A':
            self.position = 0
        elif label == 'B':
            self.position = R
        else:
            print("invalid detector")
            exit(0)
        self.label = label
        self.coupling_strength = coupling_strength
        self.Energy_gap = Energy_gap


    # get position,x of detector at time, t for wightman calculation
    def x(self, time):
        return self.position

    # get time, t of detector for wightman calculation
    def t(self, time):
        return time

class Wightman:
    def __init__(self, first, second, flip=False):
        # flip means order flip and take adjoint of it
        self.flip = flip
        if flip:
            self.time1 = tau2
            self.time2 = tau1
            self.detector1 = second
            self.detector2 = first
        else:
            self.time1 = tau1
            self.time2 = tau2
            self.detector1 = first
            self.detector2 = second

    # for evaluation purpose
    def evaluate(self):
        t1 = self.detector1.t(self.time1)
        t2 = self.detector2.t(self.time2)
        x1 = self.detector1.x(self.time1)
        x2 = self.detector2.x(self.time2)
        phase = sp.Piecewise(
            (0, (x1 - x2) ** 2 - (t1 - t2) ** 2 > 0),
            (pi, sp.And(((x1 - x2) ** 2 - (t1 - t2) ** 2 < 0), (t1 - t2 > 0))),
            (-pi, sp.And(((x1 - x2) ** 2 - (t1 - t2) ** 2 < 0), (t1 - t2 < 0))),
        )
        expr = (-1 / (4 * pi)) * (sp.log(Ir ** 2) + sp.log(sp.Abs((t1 - t2) ** 2 - (x1 - x2) ** 2)) + i * phase)
        return expr

# for code IJ value
def codingIJ(integral_type, signW, pairing, sign1, sign2, Precision):
    # initialization
    d1 = Detector(pairing[0])
    d2 = Detector(pairing[1])
    wightman = Wightman(d1, d2, flip=(signW == m))
    boundary = integral_boundary[integral_type]
    integrand = sp.exp(i * O * tau1 * symbol_sign[sign1] + i * O * tau2 * symbol_sign[sign2]) * wightman.evaluate()

    # substitute constant
    const = {O: Oc, l: lc, Ir: Irc, R: Rc, tau0: tau0c, tau: time_step}
    T1, T2 = [[pos.subs(const) for pos in layer] for layer in boundary]
    integrand_num = integrand.subs(const).evalf()

    # translate to mathematica language
    coding = mc.mathematica_translator(integrand_num, T1, T2, Precision)
    return coding

# create data
def createIJ(Precision):
    # initialization
    t1 = t.time()
    integral = [I, J]
    detector = ['A', 'B']
    sign = [p, m]
    expr = []
    data = {}
    n = 0

    # iterating possible symbol combination
    for it, sw, d1, d2, s1, s2 in itertools.product(integral, sign, detector, detector, sign, sign):
        if it == I and sw == p:
            continue
        symbol = str(it) + str(sw) + d1 + d2 + str(s1) + str(s2)
        data[symbol] = 0
        expr.append(codingIJ(it, sw, d1 + d2, s1, s2, Precision))
    final_expr = "List[" + ",".join(expr) + "]"
    result = mc.numerical_integrate(final_expr)

    # built dictionary
    for key in data.keys():
        data[key] = result[n]
        n += 1

    # save data in pkl format
    subfolder = f'{Precision}'
    file_path = data_folder / subfolder / f'data{{R = {Rc:.2f}, Δt = {time_step:.2f}}}.pkl'
    with open(file_path, 'wb') as file:
        pickle.dump(data, file)
    t2 = t.time()
    print("new data created")
    print(f"data creation took {t2-t1:.2f}s")

# retrieve data from stored file
def loadIJ(Precision):
    subfolder = f'{Precision}'
    file_path = data_folder / subfolder / f'data{{R = {Rc:.2f}, Δt = {time_step:.2f}}}.pkl'
    # check if dataset exist
    if not file_path.exists():
        print("Creating data...")
        createIJ(Precision)

    # load the dataset
    with open(file_path, "rb") as f:
        data = pickle.load(f)
    return data

# close kernel
def close():
    mc.session.terminate()
