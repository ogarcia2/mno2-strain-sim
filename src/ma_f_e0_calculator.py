import numpy as np

# Atom coordinates: 3 per structure (used to define local triangle basis)
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

def build_basis(pts):
    origin = np.array(pts[1])
    v1 = np.array(pts[0]) - origin
    v2 = np.array(pts[2]) - origin
    v3 = np.cross(v1, v2)
    return np.column_stack((v1, v2, v3))

def compute_F_E(A_pts, B_pts):
    A = build_basis(A_pts)
    B = build_basis(B_pts)
    F = B @ np.linalg.inv(A)
    E = 0.5 * (F.T @ F - np.eye(3))
    return F, E

# Compute all stepwise transforms
labels = ["delta", "TS1", "MS1", "TS2", "alpha"]
F_list, E_list = [], []

for i in range(len(labels) - 1):
    F, E = compute_F_E(structures[labels[i]], structures[labels[i + 1]])
    F_list.append(F)
    E_list.append(E)
    print(f"\n=== {labels[i]} → {labels[i+1]} ===")
    print("F =\n", np.round(F, 6))
    print("E₀ =\n", np.round(E, 6))

# Compute overall transform
F_total = F_list[-1] @ F_list[-2] @ F_list[-3] @ F_list[-4]
E_total = 0.5 * (F_total.T @ F_total - np.eye(3))

print("\n=== Overall delta → alpha ===")
print("F_total =\n", np.round(F_total, 6))
print("E_total =\n", np.round(E_total, 6))

