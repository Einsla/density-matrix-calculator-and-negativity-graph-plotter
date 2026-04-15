# density-matrix-calculator-and-negativity-graph-plotter
This code is used to calculate density matrix symbollically under quantum entanglement harvesting and degradation setup using dyson expansion operator up to 2nd order. To visualize how the negativity of the state change under this effect, this coding also allow user to sketch the negativity graph with switch off time

Startup guide:
To use the interface directly without python compiler:
1. create a new txt file
2. type as follow according to directory of folder:
   @echo off
   cd /d "(directory of the downloaded folder)"
   streamlit run interface.py
   pause
3. save as .bat file format

Operation guide:
- The interface provide 2 type of method to create initial density matrix:
  1. ket state (8 parameter)
  2. density matrix (16 parameter)
- The compute functionality will automatically normalized the ket state and display the calculated initial density matrix and eigenvalue
(This function does not support symbolic expression of parameter except for the default parameter set)
- The plot negativity functionality allows user to visualized the negativity trend as the state evolve in 2 different treatment:
  1. Continuous  ：all data point is calculated using initial density matrix by changing switch off time gradually by the time step defined by user
  2. Iteration   : data point are evealuated using final density matrix of previous result iteratively with defined time step as the switch off time of integral.
- The final density matrix compute symbolically the final density matrix expression in term of individual integral symbol 


