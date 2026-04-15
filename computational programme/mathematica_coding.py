from wolframclient.evaluation import WolframLanguageSession
from wolframclient.language import wlexpr
from sympy.printing.mathematica import mathematica_code, MCodePrinter
import sympy as sp
import time

# define mathematica file location
kernel_path = r"C:\Program Files\Wolfram Research\Wolfram\14.3\WolframKernel.exe"
session = WolframLanguageSession(kernel_path)

# adding piecewise function to python mathematica translator
def _print_Piecewise(self, expr):
    """Print Piecewise in Mathematica syntax."""
    pieces = []
    for (subexpr, cond) in expr.args:
        # Convert each condition and expression to Mathematica form
        expr_str = self._print(subexpr)
        cond_str = self._print(cond)
        pieces.append(f"{{{expr_str}, {cond_str}}}")
    pieces_str = ",".join(pieces)
    return f"Piecewise[{{{pieces_str}}}]"

# adding and function to python mathematica translator
def _print_And(self, expr):
    # recursively print all arguments
    args_str = [self._print(arg) for arg in expr.args]
    return " && ".join(args_str)

MCodePrinter._print_Piecewise = _print_Piecewise
MCodePrinter._print_And = _print_And

# warm up calculation
def warm_up():
    print("Starting Wolfram kernel...")
    session.evaluate(wlexpr("1+1"))  # trivial evaluation
    print("Kernel ready.")

# mathematica coding translator
def mathematica_translator(expr, T1, T2, Precision):
    tau1, tau2 = sp.symbols('𝜏1 𝜏2', real=True)
    t1l, t1u = T1
    t2l, t2u = T2

    # convert python input to mathematica input
    expr_str = mathematica_code(expr)
    coding = f"""
            NIntegrate[
                {expr_str},
                {{{tau1}, {t1l}, {t1u}}},
                {{{tau2}, {t2l}, {t2u}}},
                WorkingPrecision -> {Precision}
            ]"""
    return coding

# mathematica processing input and output back to python
def numerical_integrate(expr):
    result = session.evaluate(
        wlexpr(expr)
    )
    result = list(result)
    for i in range(len(result)):
        # convert mathematica output to python output
        if 'Complex' in str(result[i]):
            real = float(result[i][0])
            imag = float(result[i][1])
        else:
            real = float(result[i])
            imag = float(0.0)
        result[i] = complex(real, imag)
    return result

# for testing the parameter pass to mathematica in notebook file format
def export_to_mathematica_file(expr_str, filename="test.nb"):
    nb_content = f"""Notebook[{{Cell[TextData[{{"{expr_str}"}}], "Input"]}}]"""
    with open(filename, "w", encoding="utf-8") as f:
        f.write(nb_content)
    print(f"Mathematica notebook saved to: {filename}")
    exit(0)