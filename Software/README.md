# Baseline LM Solver (MATLAB)

## 1. What this folder provides

This folder provides the **baseline solver** used in the paper:
- Magnetic dipole model
- Levenberg-Marquardt optimization (`lsqnonlin`)
- Serial data parsing and pose estimation

## 2. One-click entrypoint

Run in MATLAB from the `Software` folder:

```matlab
run_baseline_solver
```

## 3. File structure

- `run_baseline_solver.m`: main entrypoint + all runtime helper functions in one file
- `compute_NT.mlx`: optional calibration/analysis notebook

## 4. Keyboard controls during runtime

- `q`: quit
- `t`: enable trajectory drawing
- `c`: clear trajectory
- `r`: append one row of current result to `magnet_data.xlsx`

## 5. Output data

When pressing `r`, the script appends:

1. `X (mm)`
2. `Y (mm)`
3. `Z (mm)`
4. `Alpha (deg)`
5. `Beta (deg)`
6. `Gamma (deg)`
7. `Solve Time (s)`

## 6. Environment notes

- MATLAB with Optimization Toolbox (`lsqnonlin`) is required.
- Default serial config is `COM3`, `115200`; edit `defaultConfig()` in `run_baseline_solver.m` if needed.
- `NT` is magnetic moment magnitude. Default is `5.8`; different magnets require different values.

## 7. Serial data format and parsing logic

The solver expects one full serial frame as **three consecutive lines**:

1. `U` line: starts with `U`, followed by `N` numeric values
2. `V` line: starts with `V`, followed by `N` numeric values
3. `W` line: starts with `W`, followed by `N` numeric values

For a `5x5` sensor array, `N = 25`.  
The parser reads `U/V/W` and stacks them into an `N x 3` matrix:

- column 1 = `U`
- column 2 = `V`
- column 3 = `W`

If any axis line does not contain exactly `N` values, the frame is discarded and retried (up to `serialMaxRetry` times).

