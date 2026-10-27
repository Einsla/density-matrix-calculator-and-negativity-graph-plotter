# Density Matrix Calculator and Negativity Graph Plotter

A Python/Streamlit application for exploring two-qubit density matrices in a quantum entanglement harvesting/degradation model. It can construct an initial density matrix from a ket vector or matrix entries, calculate eigenvalues of the partial transpose, plot negativity as a function of the interaction switch-off time, and display a symbolic second-order density-matrix expression.

> **Research code notice:** This repository is a computational research project, not a validated scientific package. Check the derivation, assumptions, numerical settings, and output before using results in a report or publication. The instructions below describe the current code structure; they do not imply that every workflow has been tested on a clean machine.

## Features

- Enter a two-qubit pure state using four amplitudes and four phases, or enter a (4\times4) density matrix.
- Construct/display the density matrix and calculate eigenvalues of its partial transpose.
- Plot negativity using **continuous** or **iteration** modes.
- Display the symbolic final density-matrix expression in terms of integral symbols.
- Use Wolfram Mathematica for numerical integration.

## Requirements

- Python 3.10 or newer is recommended.
- Wolfram Mathematica, with a local Wolfram Kernel installation.
- Python packages listed in [requirements.txt](requirements.txt).

The Mathematica kernel path has a Windows default in `src/mathematica_coding.py`. You can override it without editing the source by setting `WOLFRAM_KERNEL_PATH` to your Wolfram Kernel executable path in the terminal before starting Streamlit. For example, in PowerShell:

```powershell
$env:WOLFRAM_KERNEL_PATH = "C:\Program Files\Wolfram Research\Wolfram\14.3\WolframKernel.exe"
```

Change the example path/version to match your installation.

## Setup (Windows)

1. Install Python and Wolfram Mathematica.
2. Download or clone this repository.
3. Open PowerShell or Command Prompt in the repository directory.
4. Create and activate a virtual environment:

   ```powershell
   python -m venv .venv
   .\.venv\Scripts\Activate.ps1
   ```

   In Command Prompt, activate it with `.venv\Scripts\activate.bat`.

5. Install Python dependencies:

   ```powershell
   python -m pip install --upgrade pip
   pip install -r requirements.txt
   ```

6. Start the Streamlit interface from the code directory:

   ```powershell
   streamlit run src/interface.py
   ```

7. Streamlit should print a local URL (usually `http://localhost:8501`). Open it in your browser.

You can also create a Windows `.bat` launcher containing:

```bat
@echo off
cd /d "C:\path\to\density-matrix-calculator-and-negativity-graph-plotter\src"
streamlit run interface.py
pause
```

Replace the example path with the actual location on your computer.

## How to use the interface

1. Select **Ket vector** or **Density matrix**.
2. Enter the state parameters or matrix entries.
3. Click **compute** to display the density matrix and partial-transpose eigenvalues.
4. Choose **final density matrix expression** to display the symbolic expression, or **plot negativity** to configure a plot.
5. For a negativity plot, choose **continuous** or **iteration**, then set the detector separation (R), time step, and number of nodes.

The two plot modes follow different update procedures in the current implementation. Review the code and the physical model before interpreting them as equivalent.

## Project layout

```text
.
├── README.md
├── requirements.txt
├── .gitignore
├── src/
│   ├── interface.py                 # Streamlit user interface
│   ├── density_calculation.py       # Density-matrix and perturbative calculations
│   ├── negativity_plotter.py        # Partial transpose, negativity, plots and output
│   ├── computation_mathematica.py   # Integral setup and cached integral data
│   └── mathematica_coding.py        # Wolfram Kernel communication and translation
├── data/                            # Generated integral cache (created at runtime)
└── result/                          # Generated plot/data output (created at runtime)
```

The application creates `data/` and `result/` at runtime. Keep generated data if you want to reuse cached integrals or retain previous outputs.

## Important notes

- **Mathematica is required** for the numerical-integration workflow. Installing the Python dependencies alone is not sufficient.
- The kernel executable path is machine-specific; configure it for your installation.
- The first plot may take a while because integral data may need to be calculated and cached.
- Numeric inputs are needed for plotting; symbolic parameters cannot be numerically evaluated by that workflow.
- The code uses local pickle files for cached integral values. Only load pickle files that you trust.
- Generated data and plots may contain research outputs. Back them up before deleting or regenerating them.

## Troubleshooting

- **Wolfram Kernel fails to start:** verify Mathematica is installed and set `WOLFRAM_KERNEL_PATH` to the correct executable path, or update the default in `src/mathematica_coding.py`.
- **Missing Python module:** activate the intended virtual environment and run `pip install -r requirements.txt`.
- **Streamlit command not found:** install dependencies in the active environment, then run `python -m streamlit run src/interface.py` from the repository root.
- **Plotting/interpolation errors:** check that the requested number of nodes is sufficient for the interpolation method and inspect the numerical results.

## License

No license file was present when this README was prepared. Add a license if you want to specify how others may use, modify, or redistribute this project.
