# SIMPLE-CFD
This C++ code implements the Semi-Implicit Method for Pressure Linked Equations.

# Input Formatting

The following table specifies all the input fields and the expected data type. The input file consists of the field name and value on each line, seperated by a space.

| VAR | TYPE |
| --- | --- |
| mu | double |
| rho | double |
| relax | double |
| relaxp | double |
| NX | double |
| NY | double |
| LX | double |
| LY | double |
| FIELD_T | int (bool) |
| FIELD_L | int (bool) |
| FIELD_R | int (bool) |
| FIELD_B | int (bool) |
| U_T | double |
| U_L | double |
| U_R | double |
| U_B | double |
| V_T | double |
| V_L | double |
| V_R | double |
| V_B | double |
| P_T | double |
| P_L | double |
| P_R | double |
| P_B | double |



# Compiling
## Requirements
This code depends on Eigen3 for linear algebra and multidimensional array support.

## Run
```
$ ./SIMPLE < INPUT_FILE.txt
```

## Plot
Plotting is done with ```plot.py```. 

# TODO LIST
- Add RANS option
- Add compressible flow option