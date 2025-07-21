import numpy as np

# Step-wise deformation gradients
F_A1_A2 = np.array([[0.19246862, 0.        ],
                    [0.00329900, 0.89423077]])

F_A2_A3 = np.array([[1.        , 0.        ],
                    [0.02197288, 1.01075269]])

F_A3_A4 = np.array([[ 1.08695652, -0.13829787],
                    [-0.01086957,  0.87234043]])

# Total deformation gradient from A1 to A4
F_total = F_A3_A4 @ F_A2_A3 @ F_A1_A2

# Green-Lagrange strain tensor calculation
I = np.identity(2)
C_total = F_total.T @ F_total
E_total = 0.5 * (C_total - I)

# Output results
print("📐 Total Deformation Gradient F_total:\n", F_total)
print("\n💥 Total Green-Lagrange Strain Tensor E_total:\n", E_total)
