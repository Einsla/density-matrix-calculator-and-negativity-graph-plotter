import streamlit as st
import sympy as sp
import negativity_plotter as neg
import density_calculation as dc
import itertools
import pypandoc
import tempfile
import os


st.set_page_config(layout="wide")
st.title("🧮 2-Qubit Quantum Entanglement Calculator")

amp = sp.symbols("A1 A2 A3 A4", real=True)
phi = sp.symbols("θ1 θ2 θ3 θ4", real=True)
integral = sp.symbols("I J")
sign = sp.symbols("p m")
flip_sign = {sign[0]:sign[1], sign[1]:sign[0]}
detector = sp.symbols("A B")

IJ_conv = {}
for it, ws, d1, d2, s1, s2 in itertools.product(integral, sign, detector, detector, sign, sign):
    key = f"{it}{ws}{d1}{d2}{s1}{s2}"
    if (it == integral[0]) and (ws == sign[0]):
        continue
    symbol = sp.symbols(f"{it}^{{{d1}{d2}}}_{{{s1}{s2}}}", commutative=False)
    if (it == integral[1]) and (ws == sign[1]):
        symbol = sp.symbols(f"{it}^{{{d1}{d2}}}_{{{flip_sign[s1]}{flip_sign[s2]}}}", commutative=False)
        symbol = sp.conjugate(symbol)
    IJ_conv[key] = symbol

allowed_symbols = {
    "sqrt": sp.sqrt,
    "pi": sp.pi,
    "I": sp.I,
}

# -----------------------------
# Mode Selection
# -----------------------------
mode = st.radio(
    "Choose input mode:",
    ["Ket vector", "Density matrix"]
)

# -----------------------------
# Option 1: From Ket
# -----------------------------
if mode == "Ket vector":

    st.latex(r"""
    |\psi\rangle = A_1e^{i\theta_1}|00\rangle + A_2e^{i\theta_2}|01\rangle + A_3e^{i\theta_3}|10\rangle + A_4e^{i\theta_4}|11\rangle
    """)

    col1, col2 = st.columns(2)
    with col1:
        st.subheader("Amplitude")
        a = []
        for i in range(4):
            ai = st.text_input(rf"$A_{i+1}$", amp[i])
            a.append(ai)
    with col2:
        p = []
        st.subheader("Phase")
        for i in range(4):
            pi = st.text_input(rf"$\theta_{i+1}$", phi[i])
            p.append(pi)

    try:
        for i in range(4):
            a[i] = amp[i] if a[i] == str(amp[i]) else sp.sympify(a[i])
            p[i] = phi[i] if p[i] == str(phi[i]) else sp.sympify(p[i])
        psi = sp.Matrix([
            [a[0] * sp.exp(sp.I * p[0])],
            [a[1] * sp.exp(sp.I * p[1])],
            [a[2] * sp.exp(sp.I * p[2])],
            [a[3] * sp.exp(sp.I * p[3])]
        ])

    except(TypeError, ValueError):
        st.error("Invalid input.")
        st.stop()
    st.session_state['psi'] = psi

# -----------------------------
# Option 2: Manual Matrix Input
# -----------------------------
else:

    rho_0 = dc.ket_to_density(amp, phi)
    rho = sp.zeros(4,4)

    for i in range(4):
        cols = st.columns(4)
        for j in range(4):
            value = cols[j].text_input(f"$\\rho_{{{i+1}{j+1}}}$", rho_0[i, j])
            try:
                rho[i, j] = rho_0[i, j] if value == str(rho_0[i, j]) else sp.sympify(value)

            except(TypeError, ValueError):
                st.error(f"Invalid input at ({i+1},{j+1})")
                st.stop()
    st.session_state['rho_in'] = rho

if st.button("compute"):
    if mode == "Ket vector":
        psi = st.session_state['psi']
        norm = sp.sqrt(sum(sp.Abs(ai) ** 2 for ai in a)).evalf()
        psi = psi / norm
        if norm.is_Number and norm != 1:
            a = [ai / norm for ai in a]
            a_norm = [f"{ai:.2f}" for ai in a]
            st.warning(
                f"""⚠ Ket parameter provided is not normalized yet
            normalization output: [{", ".join(a_norm)}]
            """)

        rho = dc.ket_to_density(a, p)
        st.success("Density matrix computed successfully!")

    else:
        rho = st.session_state['rho_in']
        norm = sp.Trace(rho).doit().evalf()
        if norm != sum(ai ** 2 for ai in amp):
            rho = rho / norm
        if norm.is_Number and norm != 1:
            st.warning(
                f"""⚠ matrix entry provided is not normalized yet
                                normalization output as in density matrix shown below compute
                                """)

    st.markdown("---")
    # -----------------------------
    # Display Density Matrix
    # -----------------------------
    st.subheader("Density Matrix ρ")
    for i in range(4):
        col = st.columns(4)
        for j in range(4):
            col[j].latex(sp.latex(rho[i,j]))

    st.markdown("---")
    # -----------------------------
    # Display initial eigenvalue
    # -----------------------------
    st.subheader("Eigenvalues")

    # noinspection SpellCheckingInspection
    eigvals = neg.compute_eigenvalue(rho)

    for i in range(len(eigvals)):
        col = st.columns([1, 9])
        col[0].latex(f"\lambda_{i+1}")
        col[1].latex(sp.latex(eigvals[i]))

    st.markdown("---")
    try:
        # compute rho, eigenvalues
        st.session_state['rho'] = rho
        st.session_state['eigvals'] = eigvals
        st.success("Density matrix & eigenvalues computed")
    except:
        st.error("Computation failed")
    st.session_state['nested button'] = True

# -----------------------------
# button operation
# -----------------------------
if st.session_state.get('nested button', False):
    col1, col2 = st.columns(2)
    rho = st.session_state['rho']
    eigvals = st.session_state['eigvals']
    if col2.button("final density matrix expression"):
        rhodT = dc.density_change(rho, False, "polar")
        rhoT = rho + rhodT
        for i in range(4):
            col = st.columns(4)
            for j in range(4):
                subs_dict = {
                    sym: IJ_conv[str(sym)]
                    for sym in rhoT[i,j].free_symbols
                    if str(sym) in IJ_conv
                }
                rhoT[i,j] = rhoT[i,j].subs(subs_dict)
                with col[j]:
                    with st.container(border=True):
                        st.latex(sp.latex(rhoT[i,j]))

        # preparation for data download
        latex_rows = [" & ".join([sp.latex(rhoT[i, j]) for j in range(rhoT.cols)]) + r" \\"
                      for i in range(rhoT.rows)]
        latex_matrix = r"\begin{bmatrix}" + "\n".join(latex_rows) + r"\end{bmatrix}"

        # Add some document preamble
        latex_doc = r"""
        \documentclass{article}
        \usepackage{amsmath}
        \begin{document}
        Here is the matrix:
        \[
        """ + latex_matrix + r"""
        \]
        \end{document}
        """

        # Streamlit download
        st.download_button(
            label="Download LaTeX file",
            data=latex_doc,
            file_name="matrix.tex",
            mime="text/plain"
        )

    if col1.button("plot negativity"):
        plotable = True
        for val in eigvals:
            # Try evaluating to a numeric complex number
            try:
                v_numeric = complex(val.evalf())  # converts both real & imaginary parts to float
            except (TypeError, ValueError):
                plotable = False
                st.error("Ensure all input is numeric (no symbolic variables like A1, A2...)")
                st.stop()
        if plotable:
            st.session_state['show_plot_form'] = True
    if st.session_state.get('show_plot_form', False):
        with st.form("plot negativity"):
            graph_mode = st.radio(
                "Choose graph mode:",
                ["continuous", "iteration"]
            )
            col1, col2, col3 = st.columns(3)
            with col1:
                R = st.text_input(f"Separation, $R$", 4.0)
            with col2:
                t = st.text_input(f"time step, $\Delta t$", 0.1)
            with col3:
                n = st.text_input(f"number of nodes, $n$", 40)
            R = float(R)
            t = float(t)
            n = int(n)
            submitted = st.form_submit_button("Plot")
            if submitted:
                if (n * t) > R:
                    st.error(
                        """the plot include timelike separation
                        Please scale up R or scale down n or Δt""")
                    st.stop()
                plt = neg.plot_negativity(graph_mode, n, rho, t, "MachinePrecision", R)  # call your plotting function
                st.pyplot(plt)
dc.close()