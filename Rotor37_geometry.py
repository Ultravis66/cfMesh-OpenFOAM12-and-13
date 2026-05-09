"""
================================================================================
NASA Rotor 37 — Blade Geometry & CFD Domain Builder
================================================================================
Author: Mitchell R. Stolk
License: MIT
Reads blade profile coordinates from GrabCAD .curve files and builds a complete
single-passage CFD domain geometry for use with cfMesh + OpenFOAM.

Source geometry files required (place in same directory):
    profile.curve   — 6 spanwise blade profiles (0%, 20%, 40%, 60%, 80%, 100%)
    hub.curve       — hub endwall meridional curve
    shroud.curve    — shroud endwall meridional curve

Units in source files: centimeters — converted to meters throughout.

Output STL files:
    rotor37_blade_capped.stl    — blade surface with hub/tip caps
    rotor37_hub_surface.stl     — hub endwall (10 deg revolution)
    rotor37_shroud_surface.stl  — shroud endwall (10 deg revolution)
    rotor37_inlet.stl           — inlet annular face
    rotor37_outlet.stl          — outlet annular face
    rotor37_periodic1.stl       — periodic face at theta_start
    rotor37_periodic2.stl       — periodic face at theta_start + 10 deg
    rotor37_cfmesh.stl          — named multi-region STL for cfMesh

Passage geometry:
    36 blades → 10 deg passage arc
    One blade passage with periodic boundary conditions on circumferential faces
    Zero tip clearance (standard for initial validation runs)

References:
    Reid, L. and Moore, R.D. (1978). Design and Overall Performance of Four
    Highly Loaded, High-Speed Inlet Stages for an Advanced High-Pressure-Ratio
    Core Compressor. NASA Technical Paper 1337.

================================================================================
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from stl import mesh as stl_mesh


# ============================================================
# SECTION 1 — READ CURVE FILES
# ============================================================

def read_curve_file(filepath):
    """
    Read a .curve file and return a dict of named profiles.

    Handles files with # section headers and files with raw numbers only.
    Each profile is a numpy array of shape (N, 3) — X, Y, Z in meters.
    Source files are in centimeters; division by 100 converts on read.
    """
    profiles = {}
    current_label = None
    current_points = []

    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('#'):
                if current_label is not None and current_points:
                    profiles[current_label] = np.array(current_points) / 100.0
                current_label = line[1:].strip()
                current_points = []
            else:
                try:
                    vals = [float(v) for v in line.split()]
                    if len(vals) == 3:
                        current_points.append(vals)
                except ValueError:
                    continue

    # Handle files with no section headers — treat entire file as one profile
    if current_label is None and current_points:
        label = filepath.replace('\\', '/').split('/')[-1]
        profiles[label] = np.array(current_points) / 100.0
    elif current_label is not None and current_points:
        profiles[current_label] = np.array(current_points) / 100.0

    return profiles


# Read all three curve files
blade_profiles = read_curve_file('profile.curve')
hub_profile    = read_curve_file('hub.curve')
shroud_profile = read_curve_file('shroud.curve')

print("=" * 60)
print("  Geometry Files Loaded")
print("=" * 60)
print(f"  Blade profiles: {len(blade_profiles)}")
for label, pts in blade_profiles.items():
    print(f"    {label}: {len(pts)} points")
print(f"  Hub sections:    {len(hub_profile)}")
print(f"  Shroud sections: {len(shroud_profile)}")


# ============================================================
# SECTION 2 — VISUALISE BLADE SECTIONS
# ============================================================

fig = plt.figure(figsize=(14, 8))
ax  = fig.add_subplot(111, projection='3d')

colors = plt.cm.viridis(np.linspace(0, 1, len(blade_profiles)))
for (label, pts), color in zip(blade_profiles.items(), colors):
    ax.plot(pts[:, 2], pts[:, 1], pts[:, 0], color=color, linewidth=1.5, label=label)

for label, pts in hub_profile.items():
    ax.plot(pts[:, 2], pts[:, 1], pts[:, 0], 'k-', linewidth=2.0, label='Hub')

for label, pts in shroud_profile.items():
    ax.plot(pts[:, 2], pts[:, 1], pts[:, 0], 'r-', linewidth=2.0, label='Shroud')

ax.set_xlabel('Z — Axial [m]')
ax.set_ylabel('Y — Tangential [m]')
ax.set_zlabel('X — Radial [m]')
ax.set_title('NASA Rotor 37 — Blade Sections with Hub and Shroud')
ax.view_init(elev=25, azim=-50)
ax.legend(fontsize=7, loc='upper left')
plt.tight_layout()
plt.savefig('rotor37_blade_sections.png', dpi=150)
plt.show()
print("\n  Plot saved: rotor37_blade_sections.png")


# ============================================================
# SECTION 3 — STL UTILITY FUNCTIONS
# ============================================================

def revolve_curve(pts, angle):
    """Rotate a set of points by angle (radians) around the Z axis."""
    cos_a   = np.cos(angle)
    sin_a   = np.sin(angle)
    rotated = pts.copy()
    rotated[:, 0] = pts[:, 0] * cos_a - pts[:, 1] * sin_a
    rotated[:, 1] = pts[:, 0] * sin_a + pts[:, 1] * cos_a
    return rotated


def profiles_to_stl(profiles_dict, filename):
    """
    Triangulate a blade surface from ordered spanwise profiles and save as STL.
    Adjacent profiles are connected with a strip of triangles.
    """
    profile_list = list(profiles_dict.values())
    n_profiles   = len(profile_list)
    n_points     = len(profile_list[0])
    n_triangles  = (n_profiles - 1) * (n_points - 1) * 2
    blade_mesh   = stl_mesh.Mesh(np.zeros(n_triangles, dtype=stl_mesh.Mesh.dtype))

    tri_idx = 0
    for i in range(n_profiles - 1):
        p0 = profile_list[i]
        p1 = profile_list[i + 1]
        for j in range(n_points - 1):
            blade_mesh.vectors[tri_idx][0] = p0[j]
            blade_mesh.vectors[tri_idx][1] = p1[j]
            blade_mesh.vectors[tri_idx][2] = p0[j + 1]
            tri_idx += 1
            blade_mesh.vectors[tri_idx][0] = p1[j]
            blade_mesh.vectors[tri_idx][1] = p1[j + 1]
            blade_mesh.vectors[tri_idx][2] = p0[j + 1]
            tri_idx += 1

    blade_mesh.update_normals()
    blade_mesh.save(filename)
    print(f"  Saved: {filename} ({n_triangles} triangles)")


def add_caps_to_blade(profiles_dict, filename):
    """
    Build a watertight blade STL by adding hub and tip end caps.
    Each cap is a fan of triangles from the profile centroid.
    """
    profile_list = list(profiles_dict.values())
    hub_section  = profile_list[0]    # 0% span
    tip_section  = profile_list[-1]   # 100% span
    n            = len(hub_section)

    hub_centre = hub_section.mean(axis=0)
    tip_centre = tip_section.mean(axis=0)

    n_cap_tris   = (n - 1) * 2
    n_blade_tris = (len(profile_list) - 1) * (n - 1) * 2
    n_total      = n_blade_tris + n_cap_tris * 2

    full_mesh = stl_mesh.Mesh(np.zeros(n_total, dtype=stl_mesh.Mesh.dtype))
    tri_idx   = 0

    # Blade surface
    for i in range(len(profile_list) - 1):
        p0 = profile_list[i]
        p1 = profile_list[i + 1]
        for j in range(n - 1):
            full_mesh.vectors[tri_idx][0] = p0[j]
            full_mesh.vectors[tri_idx][1] = p1[j]
            full_mesh.vectors[tri_idx][2] = p0[j + 1]
            tri_idx += 1
            full_mesh.vectors[tri_idx][0] = p1[j]
            full_mesh.vectors[tri_idx][1] = p1[j + 1]
            full_mesh.vectors[tri_idx][2] = p0[j + 1]
            tri_idx += 1

    # Hub cap
    for j in range(n - 1):
        full_mesh.vectors[tri_idx][0] = hub_centre
        full_mesh.vectors[tri_idx][1] = hub_section[j + 1]
        full_mesh.vectors[tri_idx][2] = hub_section[j]
        tri_idx += 1

    # Tip cap
    for j in range(n - 1):
        full_mesh.vectors[tri_idx][0] = tip_centre
        full_mesh.vectors[tri_idx][1] = tip_section[j]
        full_mesh.vectors[tri_idx][2] = tip_section[j + 1]
        tri_idx += 1

    full_mesh.update_normals()
    full_mesh.save(filename)
    print(f"  Saved: {filename} ({n_total} triangles, capped)")


def curve_to_surface_stl(pts, angle, filename):
    """
    Create a ruled surface by revolving a meridional curve by angle (radians).
    Used to build hub and shroud endwall surfaces.
    """
    pts2    = revolve_curve(pts, angle)
    n       = len(pts)
    n_tris  = (n - 1) * 2
    s       = stl_mesh.Mesh(np.zeros(n_tris, dtype=stl_mesh.Mesh.dtype))
    tri_idx = 0
    for j in range(n - 1):
        s.vectors[tri_idx][0] = pts[j]
        s.vectors[tri_idx][1] = pts[j + 1]
        s.vectors[tri_idx][2] = pts2[j]
        tri_idx += 1
        s.vectors[tri_idx][0] = pts2[j]
        s.vectors[tri_idx][1] = pts[j + 1]
        s.vectors[tri_idx][2] = pts2[j + 1]
        tri_idx += 1
    s.save(filename)
    print(f"  Saved: {filename} ({n_tris} triangles)")


def make_annular_face_stl(r_hub_val, r_tip_val, z_val, angle, filename,
                          n_r=20, n_theta=20, theta_start=0.0):
    """
    Create an annular face at constant Z for inlet or outlet patches.
    Spans from r_hub_val to r_tip_val over the passage angle.
    """
    r_vals     = np.linspace(r_hub_val, r_tip_val, n_r)
    theta_vals = np.linspace(theta_start, theta_start + angle, n_theta)
    n_tris     = (n_r - 1) * (n_theta - 1) * 2
    s          = stl_mesh.Mesh(np.zeros(n_tris, dtype=stl_mesh.Mesh.dtype))
    tri_idx    = 0
    for i in range(n_r - 1):
        for j in range(n_theta - 1):
            p00 = np.array([r_vals[i]   * np.cos(theta_vals[j]),   r_vals[i]   * np.sin(theta_vals[j]),   z_val])
            p10 = np.array([r_vals[i+1] * np.cos(theta_vals[j]),   r_vals[i+1] * np.sin(theta_vals[j]),   z_val])
            p01 = np.array([r_vals[i]   * np.cos(theta_vals[j+1]), r_vals[i]   * np.sin(theta_vals[j+1]), z_val])
            p11 = np.array([r_vals[i+1] * np.cos(theta_vals[j+1]), r_vals[i+1] * np.sin(theta_vals[j+1]), z_val])
            s.vectors[tri_idx][0] = p00
            s.vectors[tri_idx][1] = p10
            s.vectors[tri_idx][2] = p01
            tri_idx += 1
            s.vectors[tri_idx][0] = p10
            s.vectors[tri_idx][1] = p11
            s.vectors[tri_idx][2] = p01
            tri_idx += 1
    s.save(filename)
    print(f"  Saved: {filename} ({n_tris} triangles)")


def make_periodic_face_stl(hub_pts, shr_pts, filename, n_pts=200):
    """
    Create a ruled surface connecting the hub edge to the shroud edge.
    Used for periodic boundary faces on the circumferential sides.
    Resamples both curves to n_pts to handle mismatched point counts.
    """
    def resample(pts, n):
        t_old = np.linspace(0, 1, len(pts))
        t_new = np.linspace(0, 1, n)
        return np.column_stack([
            np.interp(t_new, t_old, pts[:, 0]),
            np.interp(t_new, t_old, pts[:, 1]),
            np.interp(t_new, t_old, pts[:, 2]),
        ])

    hub_r   = resample(hub_pts, n_pts)
    shr_r   = resample(shr_pts, n_pts)
    n_tris  = (n_pts - 1) * 2
    s       = stl_mesh.Mesh(np.zeros(n_tris, dtype=stl_mesh.Mesh.dtype))
    tri_idx = 0
    for j in range(n_pts - 1):
        s.vectors[tri_idx][0] = hub_r[j]
        s.vectors[tri_idx][1] = shr_r[j]
        s.vectors[tri_idx][2] = hub_r[j + 1]
        tri_idx += 1
        s.vectors[tri_idx][0] = shr_r[j]
        s.vectors[tri_idx][1] = shr_r[j + 1]
        s.vectors[tri_idx][2] = hub_r[j + 1]
        tri_idx += 1
    s.save(filename)
    print(f"  Saved: {filename} ({n_tris} triangles)")


def write_named_stl(meshes_dict, filename):
    """
    Write a single ASCII STL with named solid regions for cfMesh.

    cfMesh requires one STL file with named 'solid' blocks — one per patch.
    Triangle winding is flipped so normals point into the fluid domain.
    """
    with open(filename, 'w') as f:
        for name, m in meshes_dict.items():
            f.write(f"solid {name}\n")
            for tri in m.vectors:
                v0, v1, v2 = tri[0], tri[2], tri[1]   # flipped winding
                n    = np.cross(v1 - v0, v2 - v0)
                norm = np.linalg.norm(n)
                if norm > 0:
                    n = n / norm
                f.write(f"  facet normal {n[0]:.6e} {n[1]:.6e} {n[2]:.6e}\n")
                f.write(f"    outer loop\n")
                f.write(f"      vertex {v0[0]:.6e} {v0[1]:.6e} {v0[2]:.6e}\n")
                f.write(f"      vertex {v1[0]:.6e} {v1[1]:.6e} {v1[2]:.6e}\n")
                f.write(f"      vertex {v2[0]:.6e} {v2[1]:.6e} {v2[2]:.6e}\n")
                f.write(f"    endloop\n")
                f.write(f"  endfacet\n")
            f.write(f"endsolid {name}\n")
    print(f"  Saved: {filename}")


# ============================================================
# SECTION 4 — BLADE STL
# ============================================================

print("\n" + "=" * 60)
print("  Exporting Blade STL")
print("=" * 60)

profiles_to_stl(blade_profiles, 'rotor37_blade.stl')
add_caps_to_blade(blade_profiles, 'rotor37_blade_capped.stl')


# ============================================================
# SECTION 5 — PASSAGE DOMAIN
# ============================================================

# Rotor 37: 36 blades → 10 degree passage arc
N_BLADES    = 36
THETA       = 2 * np.pi / N_BLADES   # 10 degrees in radians
Z_INLET     = -0.0419                 # axial position of inlet face [m]
Z_OUTLET    =  0.1067                 # axial position of outlet face [m]

print("\n" + "=" * 60)
print("  Building Passage Domain")
print("=" * 60)

hub_pts = list(hub_profile.values())[0]
shr_pts = list(shroud_profile.values())[0]

# ----------------------------------------------------------
# Angular alignment
# ----------------------------------------------------------
# The hub, shroud, and blade hub section must all start at
# the same circumferential angle. The blade defines the
# reference; hub and shroud are rotated to match.
blade_pts_ref     = list(blade_profiles.values())[0]
blade_theta_start = np.arctan2(blade_pts_ref[0, 1], blade_pts_ref[0, 0])
hub_correction    = blade_theta_start - np.arctan2(hub_pts[0, 1], hub_pts[0, 0])
shr_correction    = blade_theta_start - np.arctan2(shr_pts[0, 1], shr_pts[0, 0])

print(f"  Blade reference angle: {np.degrees(blade_theta_start):.3f} deg")
print(f"  Hub correction:        {np.degrees(hub_correction):.3f} deg")
print(f"  Shroud correction:     {np.degrees(shr_correction):.3f} deg")

hub_pts_aligned = revolve_curve(hub_pts, hub_correction)
shr_pts_aligned = revolve_curve(shr_pts, shr_correction)

# ----------------------------------------------------------
# Endwall surfaces
# ----------------------------------------------------------
curve_to_surface_stl(hub_pts_aligned, THETA, 'rotor37_hub_surface.stl')
curve_to_surface_stl(shr_pts_aligned, THETA, 'rotor37_shroud_surface.stl')

# ----------------------------------------------------------
# Inlet and outlet annular faces
# ----------------------------------------------------------
r_hub_inlet  = np.sqrt(hub_pts_aligned[0,  0]**2 + hub_pts_aligned[0,  1]**2)
r_tip_inlet  = np.sqrt(shr_pts_aligned[0,  0]**2 + shr_pts_aligned[0,  1]**2)
r_hub_outlet = np.sqrt(hub_pts_aligned[-1, 0]**2 + hub_pts_aligned[-1, 1]**2)
r_tip_outlet = np.sqrt(shr_pts_aligned[-1, 0]**2 + shr_pts_aligned[-1, 1]**2)

theta_inlet  = np.arctan2(hub_pts_aligned[0,  1], hub_pts_aligned[0,  0])
theta_outlet = np.arctan2(hub_pts_aligned[-1, 1], hub_pts_aligned[-1, 0])

print(f"  Inlet:  r_hub={r_hub_inlet:.4f}  r_tip={r_tip_inlet:.4f}  theta={np.degrees(theta_inlet):.3f} deg")
print(f"  Outlet: r_hub={r_hub_outlet:.4f}  r_tip={r_tip_outlet:.4f}  theta={np.degrees(theta_outlet):.3f} deg")

make_annular_face_stl(r_hub_inlet,  r_tip_inlet,  Z_INLET,  THETA, 'rotor37_inlet.stl',  theta_start=theta_inlet)
make_annular_face_stl(r_hub_outlet, r_tip_outlet, Z_OUTLET, THETA, 'rotor37_outlet.stl', theta_start=theta_outlet)

# ----------------------------------------------------------
# Periodic faces
# ----------------------------------------------------------
# Periodic 1 — at passage start (theta_inlet)
make_periodic_face_stl(hub_pts_aligned, shr_pts_aligned, 'rotor37_periodic1.stl')

# Periodic 2 — at passage end (theta_inlet + 10 deg)
hub_pts_rotated = revolve_curve(hub_pts_aligned, THETA)
shr_pts_rotated = revolve_curve(shr_pts_aligned, THETA)
make_periodic_face_stl(hub_pts_rotated, shr_pts_rotated, 'rotor37_periodic2.stl')

print("\n  Passage domain complete.")


# ============================================================
# SECTION 6 — NAMED MULTI-REGION STL FOR cfMesh
# ============================================================

print("\n" + "=" * 60)
print("  Exporting named multi-region STL for cfMesh")
print("=" * 60)

named_meshes = {
    'blade':     stl_mesh.Mesh.from_file('rotor37_blade_capped.stl'),
    'hub':       stl_mesh.Mesh.from_file('rotor37_hub_surface.stl'),
    'shroud':    stl_mesh.Mesh.from_file('rotor37_shroud_surface.stl'),
    'inlet':     stl_mesh.Mesh.from_file('rotor37_inlet.stl'),
    'outlet':    stl_mesh.Mesh.from_file('rotor37_outlet.stl'),
    'periodic1': stl_mesh.Mesh.from_file('rotor37_periodic1.stl'),
    'periodic2': stl_mesh.Mesh.from_file('rotor37_periodic2.stl'),
}
write_named_stl(named_meshes, 'rotor37_cfmesh.stl')


# ============================================================
# SECTION 7 — DIAGNOSTICS
# ============================================================

print("\n" + "=" * 60)
print("  Angular Alignment Diagnostic")
print("=" * 60)
for fname in ['rotor37_blade_capped.stl', 'rotor37_hub_surface.stl', 'rotor37_shroud_surface.stl']:
    m      = stl_mesh.Mesh.from_file(fname)
    pts    = m.vectors.reshape(-1, 3)
    angles = np.degrees(np.arctan2(pts[:, 1], pts[:, 0]))
    print(f"  {fname}")
    print(f"    Z:     {pts[:,2].min():.4f} to {pts[:,2].max():.4f} m")
    print(f"    Theta: {angles.min():.2f} to {angles.max():.2f} deg")

print("\n" + "=" * 60)
print("  Coordinate Summary")
print("=" * 60)
hub_pts_raw = list(hub_profile.values())[0]
shr_pts_raw = list(shroud_profile.values())[0]
bld_pts_raw = list(blade_profiles.values())[0]
print(f"  Hub    X: {hub_pts_raw[:,0].min():.4f} to {hub_pts_raw[:,0].max():.4f} m  |  Z: {hub_pts_raw[:,2].min():.4f} to {hub_pts_raw[:,2].max():.4f} m")
print(f"  Shroud X: {shr_pts_raw[:,0].min():.4f} to {shr_pts_raw[:,0].max():.4f} m  |  Z: {shr_pts_raw[:,2].min():.4f} to {shr_pts_raw[:,2].max():.4f} m")
print(f"  Blade  Z: {bld_pts_raw[:,2].min():.4f} to {bld_pts_raw[:,2].max():.4f} m")
print(f"  Hub outlet theta: {np.degrees(theta_outlet):.3f} deg")