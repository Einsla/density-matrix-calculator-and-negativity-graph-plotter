import density_calculation as dc
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
import pandas as pd
import time
from pathlib import Path
import contextlib
from scipy.interpolate import interp1d


# define parent file directory and filename
project_root = Path(__file__).resolve().parent.parent
data_folder = project_root / "result"
data_folder.mkdir(exist_ok=True)

# compute partial transpose of a density matrix
def partial_transpose(rho, dA, dB, subsystem="B"):
    """
    Partial transpose of a bipartite density matrix in SymPy.

    Parameters:
        rho : sp.Matrix of size (dA*dB, dA*dB)
        dA, dB : dimensions of subsystem A and B
        subsystem : "A" or "B"

    Returns:
        sp.Matrix : partially transposed density matrix
    """
    rho_pt = sp.zeros(dA*dB, dA*dB)

    for iA in range(dA):
        for iB in range(dB):
            for jA in range(dA):
                for jB in range(dB):
                    row = iA*dB + iB
                    col = jA*dB + jB

                    if subsystem == "A":
                        # swap iA and jA
                        new_row = jA*dB + iB
                        new_col = iA*dB + jB
                    elif subsystem == "B":
                        # swap iB and jB
                        new_row = iA*dB + jB
                        new_col = jA*dB + iB
                    else:
                        raise ValueError("subsystem must be 'A' or 'B'")

                    rho_pt[new_row, new_col] = rho[row, col]
    return rho_pt

# compute eigenvalue for dataset naming purpose
def compute_eigenvalue(rho, dims=(2,2)):
    dA, dB = dims


    # partial transpose
    rho_pt = partial_transpose(rho, 2, 2)

    # compute eigenvalue
    eigenvalues_dict = rho_pt.eigenvals()

    # convert to list
    eigenvalues = []
    for val, mult in eigenvalues_dict.items():
        eigenvalues.extend([val.evalf()] * mult)
    eigenvalues = np.array(eigenvalues)
    return eigenvalues

def compute_negativity(rho, dims=(2,2)):
    eigenvalues = compute_eigenvalue(rho, dims)
    trace_norm = np.sum(np.abs(eigenvalues))
    negativity = (trace_norm - 1.0) / 2
    return negativity

def iteration_data(n, p, time_step, Rc):
    dc.set_parameter(time_step=time_step, Rc=Rc)
    p0 = p
    data = []
    eigenvalues = compute_eigenvalue(p0)
    negativity0 = compute_negativity(p0)
    data.append([0, 0, negativity0])
    for i in range(n):
        t = (i + 1) * time_step
        dp = dc.density_change(p0)
        p0 += dp
        p0 = p0.evalf()
        purity = sp.Trace(p0 * p0).doit().evalf()

        negativity = compute_negativity(p0)
        progress = (i + 1) * 100.0 / n
        dataset = [progress, t, negativity]
        data.append(dataset)
        print(t)
        print("progress:", progress, "%")
        print("negativity:", negativity)
        print("purity:", purity)
        print()
    return data, eigenvalues

def continuous_data(n, p, time_step, Rc, Precision):
    p0 = p
    data = []
    eigenvalues = compute_eigenvalue(p0)
    negativity0 = compute_negativity(p0)
    data.append([0, 0, negativity0])
    for i in range(n):
        t = (i + 1) * time_step
        dc.set_parameter(time_step=t, Rc=Rc, Precision=Precision)
        dp = dc.density_change(p0)
        pf = p0 + dp
        pf = pf.evalf()
        purity = sp.Trace(pf * pf).doit().evalf()

        negativity = compute_negativity(pf)
        progress = (i + 1) * 100.0 / n
        dataset = [progress, t, negativity]
        data.append(dataset)
        print(t)
        print("progress:", progress, "%")
        print("negativity:", negativity)
        print("purity:", purity)
        print()
    return data, eigenvalues

def data_save(args, filename, p):
    filename += ".txt"
    filepath = data_folder / filename
    progress, t, negativity = args
    df = pd.DataFrame({
        "progress": progress,
        "time": t,
        "negativity": negativity,
    })
    cell_widths = [10, 10, 20]
    with open(filepath, "w", encoding="utf-8") as f:
        # Write initial matrix using your function
        f.write("Initial Matrix:\n")
        with contextlib.redirect_stdout(f):
            dc.matrix_display(p)  # prints go directly to file
        f.write("\n\n")

        # Write pandas DataFrame as table
        f.write("Data:\n")

        # header
        header = ["progress", "time", "negativity"]
        f.write("".join(f"{h:>{w}}" for h, w in zip(header, cell_widths)) + "\n")

        # data rows
        for row in df.itertuples(index=False):
            f.write("".join(f"{val:>{w}.6g}" for val, w in zip(row, cell_widths)) + "\n")
        f.write("\n")

def plot_negativity(datatype, n, p, time_step, Precision, Rc=4.0):
    # data processing
    if datatype == "iteration":
        data, eigenvalues = iteration_data(n, p, time_step, Rc)
    elif datatype == "continuous":
        data, eigenvalues = continuous_data(n, p, time_step, Rc, Precision)
    else:
        raise ValueError("datatype must be either 'iteration' or 'continuous'")
    data = np.array(data, dtype=object)
    progress = np.array(data[:, 0], dtype=float)
    t = np.array(data[:, 1], dtype=float)
    negativity = np.array(data[:, 2], dtype=float)

    # create smooth curve using cubic interpolation
    f = interp1d(t, negativity, kind='cubic')  # cubic spline
    t_smooth = np.linspace(t.min(), t.max(), 500)  # fine t grid
    neg_smooth = f(t_smooth)

    # generate filename
    filename = "result(" + ",".join(f"{x:.3f}" for x in eigenvalues)
    filename += f"){{R = {Rc:.2f}, Δt = {time_step:.2f}}}"
    data_save((progress, t, negativity), filename, p)

    # plotting
    fig, ax = plt.subplots()
    ax.plot(t_smooth, neg_smooth, linestyle='-', color='blue', label='Negativity')
    ax.set_xlabel('t')
    ax.set_ylabel('Negativity')
    ax.set_title('Negativity vs t (smooth)')
    ax.legend()
    ax.grid(True)

    # save plot
    filename += ".png"
    filepath = data_folder / filename
    fig.savefig(filepath)
    fig.show()
    return fig
