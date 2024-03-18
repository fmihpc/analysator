How to run local runs to calculate Diffusion coefficients and plots
ATTENTION: All paths are hardcoded and need to be changed

RUN IN ORDER
-- writeConfig.py : Creates config (.cfg) files and .sh files to run local runs based on Temperature Anisotropy and Parallel Beta. Time is based on cyclotron gyration.
- Will create the folders to individually run each run. the 'indx' variable makes sure to not run all runs at the same time (currently coded to run a batch of 100 runs out of 400).
- Will change all variables written as '$variable$' in LossCone.cfg and LossCone.sh.
- Will queue runs by itself.
- DO NOT RUN IN interactive node. Run as 'python writeConfig.py indx'.

When batch of runs are completed
-- writeVariablesSh.py : Creates the .sh files to get the variables (1dmuspace, Taniso, Beta, etc.) needed to compute Dmumu and run them.
- Computes the box averaged variables between the specified interval
- Requires indx similar to writeConfig.py
- Creates the path to store the data
- Based on 'variables.sh' which itself uses 'getVariables.py'
- DO NOT RUN IN interactive node. Run as 'python writeVariablesSh.py indx bulkStart bulkEnd interval'.

When variables batch are completed
-- writeDmumuSh.py : Creates the .sh files to compute Dmumu and run them.
- Based on 'DmumuFit.sh' which itself runs 'solverIvascenko.py'
- Requires the same input than the previous script. In addition, requires 'twindow'.
- DO NOT RUN IN interactive node. Run as 'python writeVariablesSh.py indx bulkStart bulkEnd interval twindow'.

When Dmumu batch are completed
-- plotTimeDep.py : plots the (betaPara,Taniso,nu0,t) map.
- Has 'rerun' keyword to avoid having to load all the data again (False/True)
- Creates the dat file needed for SubGrid model as well. 
- Run as 'python plotTimeDep.py run bulkStart bulkEnd interval rerun'
- example of output in 'ExampleTimeDep.png'

--------------------------------------------------------------------------------------------------------------------------
Other scripts
-- getvariables.py : Compute the box averaged variables between time intervals provided by writeVariablesSh.py

-- solverIvascenko.py : Compute Dmumu(nu0) values based on Ivascenko et al. (2016) and Dubart et al. (2022). Used by writeDmumuSh.py
- saves Dmumu, nu0Analytic (see Dubart et al. (2022) Annex), nu0Fit, and the mean between the 2.

-- plot_Dmumu_Ivascenko : Plots fmu, VDFs, Dmumu and the fits calculated by solverIvascenko.py for 1 run.
- example of output in 'ExamplePlotDmumu.png'
