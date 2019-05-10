# 
# This file is part of Analysator.
# Copyright 2013-2016 Finnish Meteorological Institute
# Copyright 2017-2018 University of Helsinki
# 
# For details of usage, see the COPYING file and read the "Rules of the Road"
# at http://www.physics.helsinki.fi/vlasiator/
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
# 

activepopulation='proton' # default

speciesdict ={
    'avgs': 'p',
    'proton': 'p',
    'helium': 'He',
    'oxygen': 'O',
    'electron': 'e',
}
speciesamu ={
    'avgs': 1,
    'proton': 1,
    'helium': 4,
    'oxygen': 16,
    'electron': 5.488e-4,
}

# Define some units for intrinsic values
unitsdict = {
    'Rhom': 'kg/m3',
    'Rhoq': 'C/m3',
    'rho': '1/m3',
    'RhoBackstream': '1/m3',
    'RhoNonBackstream': '1/m3',
    'rho_v': '1/m2s',
    'RhoVBackstream': '1/m2s',
    'RhoVNonBackstream': '1/m2s',
    'V': 'm/s',
    'VBackstream': 'm/s',
    'VNonBackstream': 'm/s',
    'B': 'T',
    'B_vol': 'T',
    'background_B': 'T',
    'perturbed_B': 'T',
    'BGB': 'T',
    'PERB': 'T',
    'PERB_vol': 'T',
    'E': 'V/m',
    'E_vol': 'V/m',
    'EXHALL_000_100': 'V/m',
    'EXHALL_001_101': 'V/m',
    'EXHALL_010_110': 'V/m',
    'EXHALL_011_111': 'V/m',
    'EYHALL_000_010': 'V/m',
    'EYHALL_001_011': 'V/m',
    'EYHALL_100_110': 'V/m',
    'EYHALL_101_111': 'V/m',
    'EZHALL_000_001': 'V/m',
    'EZHALL_010_011': 'V/m',
    'EZHALL_100_101': 'V/m',
    'EZHALL_110_111': 'V/m',
    'pressure': 'Pa',
    'pressure_dt2': 'Pa',
    'pressure_r': 'Pa',
    'pressure_v': 'Pa',
    'PTensorDiagonal': 'Pa',
    'PTensorOffDiagonal': 'Pa',
    'PTensorBackstreamDiagonal': 'Pa',
    'PTensorBackstreamOffDiagonal': 'Pa',
    'PTensorNonBackstreamDiagonal': 'Pa',
    'PTensorNonBackstreamOffDiagonal': 'Pa',
    'max_v_dt': 's',
    'max_r_dt': 's',
    'max_fields_dt': 's',
    'MinValue': 's3/m6',
    'EffectiveSparsityThreshold': 's3/m6',
    'rho_loss_adjust': '1/m3',
    'EnergyDensity': 'eV/cm3',
    'PrecipitationDiffFlux': '1/(cm2 sr s eV)'
    }

# Define some LaTeX markup names for intrinsic values
latexdict = {
    'Rhom': r'$\rho_m$',
    'Rhoq': r'$\rho_q$',
    'rho': r'$n_\mathrm{p}$',
    'RhoBackstream': r'$n_\mathrm{p,st}$',
    'RhoNonBackstream': r'$n_\mathrm{p,th}$',
    'rho_v': r'$\Gamma_\mathrm{p}$',
    'RhoVBackstream': r'$\Gamma_\mathrm{p,st}$',
    'RhoVNonBackstream': r'$\Gamma_\mathrm{p,th}$',
    'V': r'$V$',
    'VBackstream': r'$V_\mathrm{p,st}$',
    'VNonBackstream': r'$V_\mathrm{p,th}$',
    'B': r'$B$',
    'B_vol': r'$B_\mathrm{vol}$',
    'background_B': r'$B_\mathrm{bg}$',
    'perturbed_B': r'$B_\mathrm{pert}$',
    'BGB': r'$B_\mathrm{bg}$',
    'PERB': r'B_\mathrm{pert}$',
    'PERB_vol': r'B_\mathrm{vol,pert}$',
    'E': r'$E$',
    'E_vol': r'$E_\mathrm{vol}$',
    'EXHALL_000_100': r'$E_\mathrm{Hall,000,100}$',
    'EXHALL_001_101': r'$E_\mathrm{Hall,001,101}$',
    'EXHALL_010_110': r'$E_\mathrm{Hall,010,110}$',
    'EXHALL_011_111': r'$E_\mathrm{Hall,011,111}$',
    'EYHALL_000_010': r'$E_\mathrm{Hall,000,010}$',
    'EYHALL_001_011': r'$E_\mathrm{Hall,001,011}$',
    'EYHALL_100_110': r'$E_\mathrm{Hall,100,110}$',
    'EYHALL_101_111': r'$E_\mathrm{Hall,101,111}$',
    'EZHALL_000_001': r'$E_\mathrm{Hall,000,001}$',
    'EZHALL_010_011': r'$E_\mathrm{Hall,010,011}$',
    'EZHALL_100_101': r'$E_\mathrm{Hall,100,101}$',
    'EZHALL_110_111': r'$E_\mathrm{Hall,110,111}$',
    'pressure': r'$P$',
    'pressure_dt2': r'$P_{\mathrm{d}t/2}}$',
    'pressure_r': r'$P_r$',
    'pressure_v': r'$P_v$',
    'PTensorDiagonal': r'$\mathcal{P}_\mathrm{diag}$',
    'PTensorOffDiagonal': r'$\mathcal{P}_\mathrm{off-diag}$',
    'PTensorBackstreamDiagonal': r'$\mathcal{P}_\mathrm{st,diag}$',
    'PTensorBackstreamOffDiagonal': r'$\mathcal{P}_\mathrm{st,off-diag}$',
    'PTensorNonBackstreamDiagonal': r'$\mathcal{P}_\mathrm{th,diag}$',
    'PTensorNonBackstreamOffDiagonal': r'$\mathcal{P}_\mathrm{th,off-diag}$',
    'max_v_dt': r'$\Delta t_{\mathrm{max},v}$',
    'max_r_dt': r'$\Delta t_{\mathrm{max},r}$',
    'max_fields_dt': r'$\Delta t_\mathrm{max,FS}$',
    'MinValue': r'$f_\mathrm{Min}$',
    'EffectiveSparsityThreshold': r'$f_\mathrm{Min}$',
    'rho_loss_adjust': r'$\Delta_\mathrm{loss} n_\mathrm{p}$',
    }

# Define some LaTeX markup names for intrinsic values
latexdictmultipop = {
    'rho': r'$n_\mathrm{REPLACEPOP}$',
    'RhoBackstream': r'$n_\mathrm{REPLACEPOP,st}$',
    'RhoNonBackstream': r'$n_\mathrm{REPLACEPOP,th}$',
    'rho_v': r'$\Gamma_\mathrm{REPLACEPOP}$',
    'RhoVBackstream': r'$\Gamma_\mathrm{REPLACEPOP,st}$',
    'RhoVNonBackstream': r'$\Gamma_\mathrm{REPLACEPOP,th}$',
    'V': r'$V_\mathrm{REPLACEPOP}$',
    'VBackstream': r'$V_\mathrm{REPLACEPOP,st}$',
    'VNonBackstream': r'$V_\mathrm{REPLACEPOP,th}$',
    'pressure': r'$P_\mathrm{REPLACEPOP}$',
    'pressure_dt2': r'$P_{\mathrm{REPLACEPOP},\mathrm{d}t/2}}$',
    'pressure_r': r'$P_{\mathrm{REPLACEPOP},r}$',
    'pressure_v': r'$P_{\mathrm{REPLACEPOP},v}$',
    'PTensorDiagonal': r'$\mathcal{P}_\mathrm{REPLACEPOP,diag}$',
    'PTensorOffDiagonal': r'$\mathcal{P}_\mathrm{REPLACEPOP,off-diag}$',
    'PTensorBackstreamDiagonal': r'$\mathcal{P}_\mathrm{REPLACEPOP,st,diag}$',
    'PTensorBackstreamOffDiagonal': r'$\mathcal{P}_\mathrm{REPLACEPOP,st,off-diag}$',
    'PTensorNonBackstreamDiagonal': r'$\mathcal{P}_\mathrm{REPLACEPOP,th,diag}$',
    'PTensorNonBackstreamOffDiagonal': r'$\mathcal{P}_\mathrm{REPLACEPOP,th,off-diag}$',
    'MinValue': r'$f_\mathrm{REPLACEPOP,Min}}$',
    'EffectiveSparsityThreshold': r'$f_\mathrm{REPLACEPOP,Min}$',
    'EnergyDensity': r'$U_\mathrm{REPLACEPOP}$',
    'PrecipitationDiffFlux': r'$\mathcal{F}_\mathrm{REPLACEPOP}$'
    }

# Define some LaTeX markup units for intrinsic values
latexunitsdict = {
    'Rhom': r'$\mathrm{kg}\,\mathrm{m}^{-3}$',
    'Rhoq': r'$\mathrm{C}\,\mathrm{m}^{-3}$',
    'rho': r'$\mathrm{m}^{-3}$',
    'RhoBackstream': r'$\mathrm{m}^{-3}$',
    'RhoNonBackstream': r'$\mathrm{m}^{-3}$',
    'rho_v': r'$\mathrm{m}^{-2}$s',
    'RhoVBackstream': r'$\mathrm{m}^{-2}$s',
    'RhoVNonBackstream': r'$\mathrm{m}^{-2}$s',
    'V': r'$\mathrm{m}\,\mathrm{s}^{-1}$',
    'VBackstream': r'$\mathrm{m}\,\mathrm{s}^{-1}$',
    'VNonBackstream': r'$\mathrm{m}\,\mathrm{s}^{-1}$',
    'B': r'T',
    'B_vol': r'T',
    'background_B': r'T',
    'perturbed_B': r'T',
    'BGB': r'T',
    'PERB': r'T',
    'PERB_vol': r'T',
    'E': r'$\mathrm{V}\,\mathrm{m}^{-1}$',
    'E_vol': r'$\mathrm{V}\,\mathrm{m}^{-1}$',
    'EXHALL_000_100': r'$\mathrm{V}\,\mathrm{m}^{-1}$',
    'EXHALL_001_101': r'$\mathrm{V}\,\mathrm{m}^{-1}$',
    'EXHALL_010_110': r'$\mathrm{V}\,\mathrm{m}^{-1}$',
    'EXHALL_011_111': r'$\mathrm{V}\,\mathrm{m}^{-1}$',
    'EYHALL_000_010': r'$\mathrm{V}\,\mathrm{m}^{-1}$',
    'EYHALL_001_011': r'$\mathrm{V}\,\mathrm{m}^{-1}$',
    'EYHALL_100_110': r'$\mathrm{V}\,\mathrm{m}^{-1}$',
    'EYHALL_101_111': r'$\mathrm{V}\,\mathrm{m}^{-1}$',
    'EZHALL_000_001': r'$\mathrm{V}\,\mathrm{m}^{-1}$',
    'EZHALL_010_011': r'$\mathrm{V}\,\mathrm{m}^{-1}$',
    'EZHALL_100_101': r'$\mathrm{V}\,\mathrm{m}^{-1}$',
    'EZHALL_110_111': r'$\mathrm{V}\,\mathrm{m}^{-1}$',
    'pressure': r'Pa',
    'pressure_dt2': r'Pa',
    'pressure_r': r'Pa',
    'pressure_v': r'Pa',
    'PTensorDiagonal': r'Pa',
    'PTensorOffDiagonal': r'Pa',
    'PTensorBackstreamDiagonal': r'Pa',
    'PTensorBackstreamOffDiagonal': r'Pa',
    'PTensorNonBackstreamDiagonal': r'Pa',
    'PTensorNonBackstreamOffDiagonal': r'Pa',
    'max_v_dt': r's',
    'max_r_dt': r's',
    'max_fields_dt': r's',
    'MinValue': r'$\mathrm{m}^{-6}\,\mathrm{s}^{3}$',
    'EffectiveSparsityThreshold': r'$\mathrm{m}^{-6}\,\mathrm{s}^{3}$',
    'rho_loss_adjust': r'$\mathrm{m}^{-3}$',
    'EnergyDensity': r'$\mathrm{eV}\,\mathrm{cm}^{-3}$',
    'PrecipitationDiffFlux': r'$\mathrm{cm}^{-2} \,\mathrm{sr}^{-1}\,\mathrm{s}^{-1}\,\mathrm{eV}^{-1}$'
    }

