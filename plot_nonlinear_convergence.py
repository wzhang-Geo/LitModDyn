# -*- coding: utf-8 -*-
"""
Created on Thu Aug 20 21:42:05 2026

@author: wzhang
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Post-processing script for nonlinear Stokes convergence in LitModDyn.

Reads:
    N_ext_0.mat
    N_ext_1.mat
    N_ext_2.mat
    ...

Required variables in each .mat file:
    vx1, vy1, MEII, MEII0

Definitions:
    C_v^k =
        ||v^k - v^(k-1)||_2 / ||v^k||_2

    C_epsilon^k =
        ||MEII^k - MEII0^k||_2 / ||MEII^k||_2

where MEII0^k is the strain-rate invariant stored at the beginning
of nonlinear iteration k and MEII^k is the updated value after the
Stokes solve.

Author: Wentao Zhang / LitModDyn convergence post-processing
"""

import os
import re
import glob
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat


# ============================================================
# User settings
# ============================================================

# Prefix of saved MAT files:
# Northern profile: 'N_ext_'
# Southern profile: 'S_ext_'
PREFIX = 'N_ext_'

if PREFIX == 'N_ext_':
    Title_fig = 'a) Northern profile'
else:
    Title_fig = 'b) Southern profile'

# Directory containing the MAT files.
# Use '.' if this script is in the same folder as the MAT files.
DATA_DIR = '.'

# Convergence tolerance shown in the figure
TOL = 1.0e-2

# Save figure and numerical convergence history
SAVE_FIGURE = True
SAVE_DATA = True


# ============================================================
# Helper functions
# ============================================================

def iteration_number(filename, prefix):
    """
    Extract nonlinear iteration number from files such as:
        N_ext_0.mat
        N_ext_1.mat
        ...
    """
    basename = os.path.basename(filename)

    pattern = rf'^{re.escape(prefix)}(\d+)\.mat$'
    match = re.match(pattern, basename)

    if match is None:
        raise ValueError(
            f'Cannot extract iteration number from file: {basename}'
        )

    return int(match.group(1))


def relative_l2_change(current, previous):
    """
    Relative L2-norm change:
        ||current - previous||_2 / ||current||_2
    """
    current = np.asarray(current, dtype=float)
    previous = np.asarray(previous, dtype=float)

    numerator = np.linalg.norm(
        (current - previous).ravel(),
        ord=2
    )

    denominator = np.linalg.norm(
        current.ravel(),
        ord=2
    )

    if denominator == 0.0:
        return np.nan

    return numerator / denominator


def velocity_relative_change(vx_current, vy_current,
                             vx_previous, vy_previous):
    """
    Relative L2-norm change of the complete 2-D velocity solution:

        C_v =
        sqrt(||dvx||_2^2 + ||dvy||_2^2)
        --------------------------------
        sqrt(||vx ||_2^2 + ||vy ||_2^2)
    """

    dvx = np.asarray(vx_current, dtype=float) - \
          np.asarray(vx_previous, dtype=float)

    dvy = np.asarray(vy_current, dtype=float) - \
          np.asarray(vy_previous, dtype=float)

    numerator = np.sqrt(
        np.linalg.norm(dvx.ravel(), ord=2)**2
        +
        np.linalg.norm(dvy.ravel(), ord=2)**2
    )

    denominator = np.sqrt(
        np.linalg.norm(
            np.asarray(vx_current, dtype=float).ravel(),
            ord=2
        )**2
        +
        np.linalg.norm(
            np.asarray(vy_current, dtype=float).ravel(),
            ord=2
        )**2
    )

    if denominator == 0.0:
        return np.nan

    return numerator / denominator


# ============================================================
# Find and sort MAT files
# ============================================================

search_pattern = os.path.join(
    DATA_DIR,
    PREFIX + '*.mat'
)

mat_files = glob.glob(search_pattern)

# Keep only files with exact format PREFIX<number>.mat
valid_files = []

for filename in mat_files:
    basename = os.path.basename(filename)

    if re.match(
        rf'^{re.escape(PREFIX)}\d+\.mat$',
        basename
    ):
        valid_files.append(filename)

mat_files = sorted(
    valid_files,
    key=lambda f: iteration_number(f, PREFIX)
)

if len(mat_files) == 0:
    raise FileNotFoundError(
        f'No files found matching {PREFIX}<iteration>.mat '
        f'in directory: {os.path.abspath(DATA_DIR)}'
    )


print('=' * 68)
print(' LitModDyn nonlinear convergence post-processing')
print('=' * 68)

print(f'\nData directory : {os.path.abspath(DATA_DIR)}')
print(f'File prefix    : {PREFIX}')
print(f'Files detected : {len(mat_files)}\n')

for f in mat_files:
    print('  ', os.path.basename(f))


# ============================================================
# Calculate convergence metrics
# ============================================================

iterations = []
Cv = []
Cepsilon = []

previous_vx = None
previous_vy = None


for filename in mat_files:

    iteration0 = iteration_number(
        filename,
        PREFIX
    )

    # Display nonlinear iteration starting from 1
    iteration = iteration0 + 1

    data = loadmat(filename)

    required_vars = [
        'vx1',
        'vy1',
        'MEII',
        'MEII0'
    ]

    missing_vars = [
        var for var in required_vars
        if var not in data
    ]

    if missing_vars:
        raise KeyError(
            f'{os.path.basename(filename)} is missing variables: '
            + ', '.join(missing_vars)
        )

    vx = np.asarray(
        data['vx1'],
        dtype=float
    )

    vy = np.asarray(
        data['vy1'],
        dtype=float
    )

    MEII = np.asarray(
        data['MEII'],
        dtype=float
    )

    MEII0 = np.asarray(
        data['MEII0'],
        dtype=float
    )


    # --------------------------------------------------------
    # Strain-rate convergence
    #
    # MEII0 = value at beginning of the current iteration
    # MEII  = updated value after the current Stokes solution
    # --------------------------------------------------------

    ceps = relative_l2_change(
        MEII,
        MEII0
    )


    # --------------------------------------------------------
    # Velocity convergence
    #
    # First saved iteration has no previous velocity solution,
    # therefore Cv is undefined (NaN).
    # --------------------------------------------------------

    if previous_vx is None:
        cv = np.nan
    else:
        cv = velocity_relative_change(
            vx,
            vy,
            previous_vx,
            previous_vy
        )


    iterations.append(iteration)
    Cv.append(cv)
    Cepsilon.append(ceps)


    cv_text = (
        'N/A'
        if np.isnan(cv)
        else f'{cv:.6e}'
    )

    print(
        f'Iteration {iteration:2d}: '
        f'C_v = {cv_text:>12s}, '
        f'C_epsilon = {ceps:.6e}'
    )


    previous_vx = vx.copy()
    previous_vy = vy.copy()


iterations = np.asarray(
    iterations,
    dtype=int
)

Cv = np.asarray(
    Cv,
    dtype=float
)

Cepsilon = np.asarray(
    Cepsilon,
    dtype=float
)


# ============================================================
# Check convergence
# ============================================================

print('\n' + '-' * 68)
print(f'Convergence tolerance: {TOL:.1e}')

converged_iteration = None

for i in range(len(iterations)):

    if (
        np.isfinite(Cv[i])
        and np.isfinite(Cepsilon[i])
        and Cv[i] < TOL
        and Cepsilon[i] < TOL
    ):
        converged_iteration = iterations[i]
        break


if converged_iteration is not None:

    print(
        'Both convergence criteria are first satisfied at '
        f'iteration {converged_iteration}:'
    )

    idx = np.where(
        iterations == converged_iteration
    )[0][0]

    print(
        f'    C_v       = {Cv[idx]:.6e}'
    )

    print(
        f'    C_epsilon = {Cepsilon[idx]:.6e}'
    )

else:

    print(
        'The two convergence criteria are not simultaneously '
        'below the prescribed tolerance in the available files.'
    )

print('-' * 68)


# ============================================================
# Save convergence history
# ============================================================

if SAVE_DATA:

    output_txt = os.path.join(
        DATA_DIR,
        PREFIX + 'nonlinear_convergence.txt'
    )

    table = np.column_stack(
        (
            iterations,
            Cv,
            Cepsilon
        )
    )

    np.savetxt(
        output_txt,
        table,
        fmt=['%d', '%.8e', '%.8e'],
        header=(
            'Iteration    '
            'C_v    '
            'C_epsilon\n'
            'C_v = ||v(k)-v(k-1)||_2 / ||v(k)||_2\n'
            'C_epsilon = ||MEII(k)-MEII0(k)||_2 / ||MEII(k)||_2'
        )
    )

    print(
        '\nConvergence data saved to:\n  '
        + output_txt
    )


# ============================================================
# Plot convergence history
# ============================================================

fig, ax = plt.subplots(
    figsize=(6.6, 4.8)
)


# C_v starts from the second nonlinear iteration
mask_v = (
    np.isfinite(Cv)
    & (Cv > 0.0)
)

ax.semilogy(
    iterations[mask_v],
    Cv[mask_v]/(iterations[mask_v])**0,
    '-o',
    linewidth=1.8,
    markersize=5.5,
    label=r'$C_v$'
)


# C_epsilon can be calculated from the first iteration
mask_e = (
    np.isfinite(Cepsilon)
    & (Cepsilon > 0.0)
)

ax.semilogy(
    iterations[mask_e],
    Cepsilon[mask_e]/(iterations[mask_e])**0,
    '-s',
    linewidth=1.8,
    markersize=5.5,
    label=r'$C_{\dot{\varepsilon}}$'
)


# Convergence threshold
ax.axhline(
    TOL,
    linestyle='--',
    linewidth=1.3,
    # label=rf'Tolerance = ${TOL:.0e}$'
    label=rf'Tolerance = $10^{{{int(np.log10(TOL))}}}$'
)


ax.set_xlabel(
    'Nonlinear iteration',
    fontsize=12
)

ax.set_ylabel(
    r'Relative $L_2$-norm change',
    fontsize=12
)


ax.set_title(
    Title_fig
)


ax.set_xticks(
    iterations
)

ax.tick_params(
    axis='both',
    labelsize=10
)

ax.grid(
    True,
    which='both',
    linestyle=':',
    linewidth=0.7
)

ax.legend(
    frameon=False,
    fontsize=10
)

fig.tight_layout()


if SAVE_FIGURE:

    output_png = os.path.join(
        DATA_DIR,
        PREFIX + 'nonlinear_convergence.png'
    )

    output_pdf = os.path.join(
        DATA_DIR,
        PREFIX + 'nonlinear_convergence.pdf'
    )

    fig.savefig(
        output_png,
        dpi=300,
        bbox_inches='tight'
    )

    fig.savefig(
        output_pdf,
        bbox_inches='tight'
    )

    print(
        '\nFigures saved to:\n'
        f'  {output_png}\n'
        f'  {output_pdf}'
    )


plt.show()