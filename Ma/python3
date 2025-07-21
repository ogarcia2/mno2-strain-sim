import numpy as np

# Atom coordinates for each structure (3 atoms used to define the triangle)
structures = {
    "delta": [
        [15.09821, -2.29574, 8.78582],
        [7.31011, -2.29574, 8.78582],
        [6.5590607, -3.56558899, 16.54733386],
    ],
    "TS1": [
        [15.50684, -2.60366, 8.29335],
        [7.79384, -2.60366, 8.29335],
        [6.8021954, -3.68996651, 16.64557779],
    ],
    "MS1": [
        [15.94859, -2.33862, 7.85623],
        [7.79929, -2.33862, 7.85623],
        [6.297646, -3.11887829, 16.57725618],
    ],
    "TS2": [
        [18.44607, -0.46791, 6.9194],
        [10.23657, -0.46791, 6.9194],
        [8.4037451, -4.49515367, 14.74076286],
    ],
    "alpha": [
        [18.3474, -0.18058, 6.91377],
        [10.0347, -0.18058, 6.91377],
        [7.97224, -4.07039, 13.84072],
    ],
}

def build_2D_basis(points):
    """Drop y-axis and return x–z 2D basis from 3 atoms"""
    a = np.array(points[0])[[0, 2]]
    b = np.array(points[1])[[0, 2]]
    c = np.array(points[2])[[0, 2]]
    v1 = a - b
    v2 = c - b
    return np.column_stack((v1, v2))

def compute_F_E2D(A_pts, B_pts):
    A = build_2D_basis(A_pts)
    B = build_2D_basis(B_pts)
    F = B @ np.linalg.inv(A)
    E = 0.5 * (F.T @ F - np.eye(2))
    return F, E

# Perform calculations for each transition
labels = ["delta", "TS1", "MS1", "TS2", "alpha"]
F_list, E_list = [], []

for i in range(len(labels) - 1):
    F, E = compute_F_E2D(structures[labels[i]], structures[labels[i+1]])
    F_list.append(F)
    E_list.append(E)
    print(f"\n=== {labels[i]} → {labels[i+1]} ===")
    print("F_2D =\n", np.round(F, 6))
    print("E0_2D =\n", np.round(E, 6))

# Chain total F and E0
F_total = F_list[-1] @ F_list[-2] @ F_list[-3] @ F_list[-4]
E_total = 0.5 * (F_total.T @ F_total - np.eye(2))

print("\n=== Overall delta → alpha (2D) ===")
print("F_total_2D =\n", np.round(F_total, 6))
print("E0_total_2D =\n", np.round(E_total, 6))

