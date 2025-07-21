import numpy as np
import sys
import argparse

def compute_deformation_gradient(A_coords, B_coords):
    """
    Computes deformation gradient F using vectors from bottom-left to bottom-right and top-left.
    Parameters:
        A_coords: list of 3 tuples from original lattice: [top-left, bottom-left, bottom-right]
        B_coords: list of 3 tuples from deformed lattice: same order
    Returns:
        F: 2x2 numpy array
    """
    # Reorder to origin at bottom-left
    A_top, A_botL, A_botR = [np.array(p) for p in A_coords]
    B_top, B_botL, B_botR = [np.array(p) for p in B_coords]

    A_vec1 = A_botR - A_botL  # x-direction
    A_vec2 = A_top - A_botL   # y-direction

    B_vec1 = B_botR - B_botL
    B_vec2 = B_top - B_botL

    A_matrix = np.column_stack((A_vec1, A_vec2))
    B_matrix = np.column_stack((B_vec1, B_vec2))

    F = B_matrix @ np.linalg.inv(A_matrix)
    return F

def compute_e0(F):
    """Computes strain energy tensor E₀ = 0.5 * (FᵀF - I)"""
    I = np.identity(2)
    C = F.T @ F
    E = 0.5 * (C - I)
    return E

def parse_args():
    parser = argparse.ArgumentParser(description="Calculate deformation gradient F and strain tensor E₀ from 2D lattice corner coordinates.")
    parser.add_argument("--coords", nargs=12, type=float, metavar=('Ax1','Ay1','Ax2','Ay2','Ax3','Ay3','Bx1','By1','Bx2','By2','Bx3','By3'),
                        help="Coordinates for A and B lattices: top-left, bottom-left, bottom-right of A, then same for B.")
    return parser.parse_args()

def prompt_for_coords():
    print("\n🔧 Enter pixel coordinates of lattice corners in this order:")
    print("  1. Top-left\n  2. Bottom-left\n  3. Bottom-right\n")
    A = []
    B = []
    for label in ['Original (A)', 'Deformed (B)']:
        print(f"\n→ {label} lattice:")
        for point in ['Top-left', 'Bottom-left', 'Bottom-right']:
            x = float(input(f"  {point} X: "))
            y = float(input(f"  {point} Y: "))
            (A if label == 'Original (A)' else B).append((x, y))
    return A, B

def main():
    args = parse_args()

    if args.coords:
        coords = args.coords
        if len(coords) != 12:
            print("❌ Error: Must provide exactly 12 numbers (6 points).")
            sys.exit(1)
        A_coords = [(coords[0], coords[1]), (coords[2], coords[3]), (coords[4], coords[5])]
        B_coords = [(coords[6], coords[7]), (coords[8], coords[9]), (coords[10], coords[11])]
    else:
        A_coords, B_coords = prompt_for_coords()

    F = compute_deformation_gradient(A_coords, B_coords)
    E = compute_e0(F)

    print("\n📐 Deformation Gradient F:\n", F)
    print("\n💥 Strain Tensor E₀:\n", E)

if __name__ == "__main__":
    main()
