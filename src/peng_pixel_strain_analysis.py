import numpy as np

def average_box_dimensions(sides):
    """Takes pixel lengths in order: top, left, right, bottom → returns width, height"""
    width = np.mean([sides[0], sides[3]])
    height = np.mean([sides[1], sides[2]])
    return np.array([width, height])

# Raw pixel counts per structure (top, left, right, bottom)
pixels = {
    "A1_raw": [398, 409, 391, 399],
    "A2": [363, 353, 364, 360],
    "A3": [368, 381, 380, 370],
    "A4": [365, 440, 440, 365],
}

# Shrink A1 by the known 10.9% scaling factor
scale_factor = 1 / 1.109
pixels["A1"] = [p * scale_factor for p in pixels["A1_raw"]]

# Convert to 2D lattice vectors (width, height)
shapes = {k: average_box_dimensions(v) for k, v in pixels.items()}

# Stepwise transitions
transitions = ["A1", "A2", "A3", "A4"]
F_list = []
E_list = []

for i in range(len(transitions) - 1):
    A = shapes[transitions[i]]
    B = shapes[transitions[i + 1]]
    F = np.diag(B / A)
    E = 0.5 * (F.T @ F - np.eye(2))
    F_list.append(F)
    E_list.append(E)

# Total transformation
F_total = F_list[2] @ F_list[1] @ F_list[0]
E_total = 0.5 * (F_total.T @ F_total - np.eye(2))

# Print results
np.set_printoptions(precision=6, suppress=True)
print("📐 Total Deformation Gradient F_total:\n", F_total)
print("\n💥 Total Green-Lagrange Strain Tensor E_total:\n", E_total)

