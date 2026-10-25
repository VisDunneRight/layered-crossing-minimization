"""
Adaptive breakpoints for a piecewise-linear (PWL) approximation of
phi = arctan(s), for use with Gurobi's addGenConstrPWL.
 
Each piece is made as long as possible while keeping the exact worst-case
angle error at or below eps. Pieces therefore come out long where arctan is
nearly linear (around s = 0) and short where it bends (the steep ends), so
you hit a given accuracy with far fewer breakpoints than uniform spacing.
"""
 
import numpy as np
 
 
def chord_max_error(s1, s2):
    """Exact max |chord - arctan| on [s1, s2], with the chord through the
    two endpoints on the curve.
 
    The error e(s) = chord(s) - arctan(s) is zero at both endpoints, and its
    derivative m - 1/(1+s^2) vanishes only at s = +/- sqrt(1/m - 1), where m
    is the chord slope. So the max is at one of those points if it lies
    inside the interval.
    """
    if s2 - s1 <= 0.0:
        return 0.0
    p1, p2 = np.arctan(s1), np.arctan(s2)
    m = (p2 - p1) / (s2 - s1)          # always 0 < m <= 1 for arctan
    if m >= 1.0:
        return 0.0
    r = np.sqrt(1.0 / m - 1.0)
    err = 0.0
    for c in (-r, r):
        if s1 < c < s2:
            err = max(err, abs(p1 + m * (c - s1) - np.arctan(c)))
    return err
 
 
def adaptive_breakpoints(s_min, s_max, eps):
    """Greedy breakpoints on [s_min, s_max] with max angle error <= eps (rad).
 
    Starting from the left end, each new breakpoint is pushed as far right
    as bisection allows while the piece still meets the tolerance.
    """
    pts = [s_min]
    a = s_min
    while True:
        if chord_max_error(a, s_max) <= eps:
            pts.append(s_max)
            break
        lo, hi = a, s_max                  # lo always meets the tolerance
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            if chord_max_error(a, mid) <= eps:
                lo = mid
            else:
                hi = mid
        if lo <= a:                        # safety guard; should not happen
            raise RuntimeError("Tolerance too tight to make progress.")
        pts.append(lo)
        a = lo
    return np.array(pts)
 
 
def uniform_angle_breakpoints(s_min, s_max, K):
    """K breakpoints evenly spaced in angle (the previous approach)."""
    return np.tan(np.linspace(np.arctan(s_min), np.arctan(s_max), K))
 
 
def max_error(s_pts):
    """Worst-case angle error (rad) of a breakpoint set."""
    return max(chord_max_error(a, b) for a, b in zip(s_pts[:-1], s_pts[1:]))


def add_arctan_pwl(model, s_var, phi_var, s_pts):
    """Add phi_var = PWL(arctan(s_var)) to a Gurobi model.
    Gurobi extrapolates a PWL beyond its first and last points, so make sure
    s_var's bounds lie inside [s_pts[0], s_pts[-1]].
    """
    s_pts = np.asarray(s_pts, dtype=float)
    model.addGenConstrPWL(s_var, phi_var,
                          s_pts.tolist(), np.arctan(s_pts).tolist())

 
if __name__ == "__main__":
    # Example: y-values in [0, 6] give slopes in [-6, 6] (about +/- 80.5 deg).
    s_min, s_max = -30, 30
 
    for tol_deg in [0.5, 0.1, 0.05, 0.01]:
        eps = np.radians(tol_deg)
        ad = adaptive_breakpoints(s_min, s_max, eps)
 
        # Smallest uniform-in-angle K that meets the same tolerance.
        K = 2
        while max_error(uniform_angle_breakpoints(s_min, s_max, K)) > eps:
            K += 1
 
        print(f"tol {tol_deg:>5} deg per line: adaptive K = {len(ad):>3}, "
              f"uniform-angle K = {K:>3}  "
              f"(adaptive max err {np.degrees(max_error(ad)):.4f} deg)")
 
    # Show where the breakpoints land, in degrees.
    ad = adaptive_breakpoints(s_min, s_max, np.radians(0.1))
    print("\nBreakpoint angles for 0.1 deg tolerance:")
    print(np.round(np.degrees(np.arctan(ad)), 1))
 
    # Gurobi usage sketch:
    # import gurobipy as gp
    # m = gp.Model()
    # s   = m.addVar(lb=s_min, ub=s_max, name="s")
    # phi = m.addVar(lb=np.arctan(s_min), ub=np.arctan(s_max), name="phi")
    # m.addConstr(s == y2 - y1)
    # add_arctan_pwl(m, s, phi, adaptive_breakpoints(s_min, s_max, eps))