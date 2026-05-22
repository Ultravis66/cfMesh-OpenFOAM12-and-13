# cfMesh for OpenFOAM 13

> **USABLE BETA** — active development on OpenFOAM Foundation v13. Quality system significantly improved over original cfMesh. Run-to-run variance exists due to residual OMP nondeterminism in base cfMesh pipeline — being addressed incrementally.

A working port of cfMesh to OpenFOAM Foundation version 13, with significant boundary layer robustness improvements beyond the original cfMesh.

Original cfMesh developed by Franjo Juretic / Creative Fields, licensed under GPL v3.

Tested on **NASA Rotor 37** two-passage turbomachinery sector geometry with boundary layers.

---

## Branch status

| Branch | OpenFOAM Version | Status |
|--------|-----------------|--------|
| `master` | OpenFOAM 13 | **Active development** |
| `openfoam12-archive` | OpenFOAM 12 | Archived — no further updates |

---

## Quick start

**Tested with:** OpenFOAM Foundation v13, Ubuntu 22.04 / WSL2

```bash
source /opt/openfoam13/etc/bashrc
git clone https://github.com/Ultravis66/cfMesh-OpenFOAM13.git cfMesh
cd cfMesh
./Allwclean
./Allwmake
```

Verify the build:
```bash
ls $FOAM_USER_APPBIN | grep -E "cartesianMesh|tetMesh|pMesh"
```

Run the validation case:
```bash
cd tutorials/cartesianMesh/unitCube
cartesianMesh
checkMesh
```

---

## What works
- `cartesianMesh` — Cartesian hex-dominant mesher, primary supported executable, tested on NASA Rotor 37
- `tetMesh` — Builds successfully, experimental, not fully validated
- `pMesh` — Builds successfully, experimental, not fully validated
- `libmeshLibrary.so` — Core mesh library

---

## API changes from original cfMesh (compile fixes)

| Original | Updated |
|----------|---------|
| `Pstream::blocking/scheduled/nonBlocking` | `Pstream::commsTypes::blocking/scheduled/nonBlocking` |
| `eigenVector(T, scalar)` | `eigenVectors(T)` with direct component access |
| `circumCentre()/circumRadius()` | `circumSphere().first()/.second()` |
| `isFile()` | `std::filesystem::exists()` — requires filesystem include and c++17 |
| `rm/isDir/rmDir` | requires `OSspecific.H` include in affected files |
| edgeMesh library | include paths updated to use meshTools |
| `unallocLabelList` | `labelList` |

---

## Runtime fixes

| Bug | Fix |
|-----|-----|
| `face::centre()` garbage for cfMesh face ordering — phantom inversions | Replace with vertex average everywhere |
| `face::normal()` garbage — wrong BL extrusion direction | Replace with triangle fan normal everywhere |
| Inversion check on unmapped octree mesh — infinite loop | Map vertices onto surface BEFORE `preMapVertices()` |
| `preMapVertices()` internal untangleSurface hangs | Use `preMapVertices(0)` |
| `findLeafContainingVertex` throws `char const*` | Replace throw with return -1 |
| `newPositionSurfaceOptimizer` garbage vertex positions | Displacement bound + vecY near-zero guard |
| `pointFieldPMG::setSize` uninitialized over-allocation | Zero-initialize over-allocated region |
| Octree points uninitialized | Zero-initialize `octreePointsPtr_` |
| `checkFacePatchesGeometry` unbounded loop | Iteration limit (`nGeomIter < 5`) |
| `findNearestSurfacePoint` 100-iteration hang | Reduce to 10 iterations |
| `optimizePoint` FPE from zero division | Guard zero divisions in surfaceOptimizer |
| `untangleSurface` too many global iterations | Limit to 2 global iterations |
| Garbage denormalized points in output | Sanitize before write |
| `face::areaAndCentre` hangs in OF13 | Replace with vertex average/triangle fan normal |
| `omp_get_num_procs()` over-threads OMP sections | Replace with `omp_get_max_threads()` |
| `meshOctreeCube::leavesInBox` infinite recursion | Add size guard and self-check |
| `meshSurfaceMapperCornersAndEdges` 20-iteration hang | Reduce to 5 iterations |
| `surfaceOptimizer` 100-iteration hang | Reduce to 10 iterations |
| `addLayerForAllPatches` O(n²) BL quality rollback hangs | Removed |
| `mapEdgeNodes` OMP nondeterminism — intermittent hangs and bad face explosions at full thread count | Deterministic `parallelBndNodes` build: compute in parallel, move and sync in serial |

---

## Boundary Layer Generation Improvements

### Patch-based BL control
Patches with `nLayers 0` (inlet, outlet, periodic) are fully excluded from BL generation. Previously cfMesh required hardcoded patch names — now reads directly from `meshDict boundaryLayers/patchBoundaryLayers`. Three code paths fixed: `findPatchesToBeTreatedTogether`, `addLayerForPatch`, `addLayerForAllPatches`.

### BL/BL junction topology infrastructure (C1)

Persistent metadata for sharp BL/BL junction points is now carried through the boundary-layer generation pipeline:

- `blblJunctionPoints_` — `labelHashSet` on `boundaryLayers` capturing boundary points on sharp BL/BL feature curves.
- Junction points use `layerScale_ = 0.0` to suppress spurious extrusion at the singular feature line.
- Sharp BL/BL junctions are detected using patch-normal angle comparison rather than hardcoded patch names.

This provides groundwork for a future transition-topology implementation at sharp BL/BL intersections. The remaining planned work is explicit wedge/pyramid/O-grid style topology emission for junction regions where ordinary prism extrusion becomes ill-conditioned.

### layerScale ramp field
Smooth suppression zone at BL/no-BL patch boundaries:

| Zone | Scale |
|------|-------|
| Transition edge | 0.0 |
| Ring 1 | 0.25 |
| Ring 2 | 0.50 |
| Ring 3 | 0.75 |
| Ring 4 | 1.00 |

Previously ring1 was incorrectly 0.0, creating a flat zone and steep gradient that caused warts at patch junctions.

## Surface quality system (v0.1–v0.5)

A commercial-mesher-style surface projection and repair pipeline, developed incrementally on top of the original cfMesh surface mapping code.

### Architecture fix — patch-identity-first projection
The original cfMesh pipeline ran `mapMeshToSurface` (global surface projection) before `extractPatches` (patch identity assignment). This meant every boundary point was projected using a global nearest-triangle search with no patch constraint — causing points near concave multi-patch corners to snap to the wrong patch surface, producing visible protrusions.

Fix: `extractPatches` now runs before `mapMeshToSurface`, establishing patch identity on the raw Cartesian hex mesh before any projection occurs. A second `extractPatches` call after projection corrects any misassignments from the unprojected hex geometry. This is the primary fix for surface protrusions.

### Boundary point topology classifier
Every boundary point is classified into one of five topology classes before projection:

| Class | Description |
|-------|-------------|
| `SINGLE_PATCH` | Interior surface point — standard projection |
| `TWO_PATCH_EDGE` | Feature curve point — constrained to edge |
| `MULTI_PATCH_CORNER` | 3+ patch junction — constrained to corner endpoint |
| `BL_NOBL` | BL/no-BL transition — locked |
| `NONMANIFOLD` | Non-manifold geometry — locked and reported |

VTK diagnostic output available for each class for ParaView inspection.

### Feature curve extraction and corner snap
Dual-mode feature curve extraction (face-patch data when available, pointPatches fallback for early pipeline stages). Extracts ~10,260 feature edge segments. 24/24 multi-patch corner points snap correctly to nearest incident feature curve endpoint.

### Pure quality evaluator — `proposedMoveIsValid`
Reusable per-point quality gate used by all movement logic. Evaluates a proposed new position purely geometrically — no mesh mutation, no `const_cast`, thread-safe by construction.

Checks per incident face:
1. Movement distance cap (relative to local cell size)
2. Face area (reject near-zero collapse)
3. Pyramid height (relative to face area — scale-invariant)
4. Skewness proxy (rejects moves that worsen local face skewness beyond threshold)

### Constrained surface smoothing
Patch-aware Laplacian smoothing pass after corner/edge snap. Two-phase update (compute from frozen geometry, apply in fixed order):
- `SINGLE_PATCH` points: Laplacian average of same-patch neighbors → re-project to own patch → `proposedMoveIsValid` gate
- `TWO_PATCH_EDGE` and `MULTI_PATCH_CORNER`: locked
- `BL_NOBL`: locked

### Targeted repair of validity-rejected points
After global projection, points that failed the validity check are collected. Their 2-ring neighborhood is expanded, classified, and repaired using constrained patch-aware projection. Typically repairs ~1,250 of ~1,500 neighborhood points per run, targeting the blade/periodic junction regions where protrusions concentrate.

### Bisection backtracking (experimental, gated off)
Binary search along the snap vector for points that fail full projection. Gated behind `useBisectionBacktracking = false` pending stronger quality bundle (skewness + non-ortho comparison against old position required before enabling).

### Version history

| Tag | Description |
|-----|-------------|
| `v0.1-beta-stable` | Pipeline reorder + relative validity thresholds |
| `v0.2` | Skewness proxy added to quality evaluator |
| `v0.35` | Validity-rejected point VTK diagnostics |
| `v0.4` | Targeted repair of validity-rejected points |
| `v0.5` | Pure quality evaluator — no const_cast, thread-safe, tight variance |

### BL/no-BL transition edge constraint system

Generic protection system for points where BL patches meet no-BL patches (wall meets inlet/outlet/periodic). Fully geometry-agnostic — driven entirely by `nLayersForPatch_` from meshDict, no patch names hardcoded.

| Component | Description |
|-----------|-------------|
| `detectBLNoBlTransitionEdges()` | Lightweight pre-detection before edge extraction. Runs before any surface projection. No mesh modification. |
| `blNoBlEdgePoints_` / `blNoBlEdges_` | Generic layer-termination edge and point sets |
| `blNoBlPointPatch_` | BL-side patch ownership per interface point |
| `protectedPoints_` + `protectedPointPatches_` | Per-mapper protection: interface points redirected to `findNearestSurfacePointInRegion` on their own patch only, preventing projection across feature curves onto the wrong patch |
| Guards | `mapEdgeNodes`, `mapCorners`, `mapVerticesOntoSurface`, `mapVerticesOntoSurfacePatches` |
| `meshSurfaceEdgeExtractorNonTopo` | Overloaded constructor accepts and passes through protected point sets |


### meshDict geometry controls

| Parameter | Default | Description |
|-----------|---------|-------------|
| `blblFeatureAngleDeg` | 40 | Angle threshold for detecting sharp BL/BL junction curves |
| `layerScaleRing1` | 0.25 | BL scale at ring 1 (nearest to transition edge) |
| `layerScaleRing2` | 0.50 | BL scale at ring 2 |
| `layerScaleRing3` | 0.75 | BL scale at ring 3 |
| `layerScaleRing4` | 1.00 | BL scale at ring 4 (full BL thickness) |
| `cornerSnapRelaxation` | 0.25 | Corner vertex snap damping (0=locked, 1=full snap) |
| `postBLSnap` | false | Post-BL surface re-projection (experimental) |
| `gapDetection` | false | Enable gap detection pre-pass (diagnostic only) |
| `gapSuppression` | false | Suppress BL at detected gap-risk points |
| `gapDetectionDistance` | 1e-4 | Absolute gap width threshold for detection |
| `gapProbeDistance` | 1e-6 | Inward normal probe distance for gap measurement |

### BL/BL sharp-junction suppression
Points where two BL patches meet at a sharp angle (blade+hub, blade+shroud) are detected using patch normal angle comparison and suppressed. Prevents degenerate layer cells at blade-endwall junctions in turbomachinery geometry.

### refineFace graceful handling
`refineBoundaryLayers::refineFace` previously called `FatalError` when no split direction could be found for transition-zone faces. Now returns the original face unrefined instead of crashing.

### Two-pass candidate quality rollback (experimental infrastructure)
Initial infrastructure exists for detecting degenerate layer-cell candidates and reducing local layer complexity, but the detection/rollback policy is still under development:
- `mixedDisplacementFaces_` member on `boundaryLayers`
- `forceSingleLayerFaces_` + `forceSingleLayerAtFaces()` on `refineBoundaryLayers`
- Layer count smoothing on boundary-face adjacency graph (max jump of 1 between neighbors)
- Phase 1 non-mutating diagnostic in progress

### Known remaining issues

**Residual OMP nondeterminism** — Run-to-run variance in mesh quality metrics (non-ortho errors, severely non-ortho face count) due to dynamic OMP scheduling in base cfMesh projection passes. Skewness is now stable (9.3–9.5). Non-ortho errors vary 0–1 typical. Being addressed incrementally by making projection passes fully deterministic.

**Residual spikes at blade_4/periodic junction** — Small cluster of cells at the blade trailing edge / periodic boundary intersection still show elevated non-orthogonality. Root cause: points near this junction are reverted by the validity checker during projection and cannot be fully repaired in one pass. Targeted repair pass reduces but does not eliminate this.

**BL/no-BL interface protrusions** — BL prism cells extrude through adjacent no-BL patch surfaces at corner junctions. Fix in progress: clamp extrusion normals at `blNoBlEdgePoints_` to BL-side faces only.

**2 high-aspect-ratio cells at blade-endwall junctions** — Machine epsilon volume at geometric singularity. Planned fix: wedge/pyramid topology emission at sharp BL/BL junctions.

> **Diagnostic note:** SnappyHexMesh produces clean geometry at the same patch junctions using the same STL, confirming the protrusions are cfMesh-specific and not a geometry defect.

---
## Mesh quality on NASA Rotor 37 (v0.5)

Tested with boundary layers active (3-layer blade, 5-layer hub/shroud).
Results shown as typical range across runs — residual variance from base cfMesh OMP scheduling.

| Metric | v0.5 typical | Original cfMesh |
|--------|-------------|-----------------|
| Cells | ~2.1M | ~2.1M |
| Average non-orthogonality | 10.6° | 10.5° |
| Max non-orthogonality | 90–100° | 147° |
| Non-orthogonality errors | 0–1 | 1+ |
| Negative volume cells | 0–2 | 5–7 |
| Bad pyramid faces | 1–9 | 13–32 |
| Max skewness | 9.3–9.5 | 15.6+ |
| Severely non-ortho faces (>70°) | 1,800–2,100 | 2,058 |
| High aspect ratio cells | 3–5 | 13 |
| Runtime | ~75–100s | ~75s |
| Visible surface protrusions | Significantly reduced | Severe at patch junctions |

---
## Running the regression test

A automated validation script is included that runs `cartesianMesh` on the Rotor 37 case and checks all key quality metrics against known-good baselines.

```bash
~/cfMesh/tutorials/rotor37/regression_test.sh ~/rotor37_mesh
```

Expected output:
```
======================================================
cfMesh Rotor 37 Regression Test
...
[ 1/3 ] Running cartesianMesh...
PASS: cartesianMesh completed in ~75s
PASS: Bad faces after untangling: 0 (limit: 10)
[ 2/3 ] Running checkMesh...
PASS: Cell count: ~2.0M (baseline: 2100000 +/-10%)
PASS: Average non-orthogonality: ~10deg (limit: 15deg)
PASS: High aspect ratio cells: 3-5 (known limit: 10)
PASS: Negative volume cells: 0-2 (known limit: 5)
[ 3/3 ] Summary
PASSED: 6
FAILED: 0
RESULT: PASS
```

Run this after any code change to catch regressions before committing.

## Troubleshooting

**Build fails with missing header:** Make sure OpenFOAM 13 is sourced before building. Run `echo $WM_PROJECT_VERSION` — should return `13`.

**`cartesianMesh` not found after build:** Check `$FOAM_USER_APPBIN` is in your PATH. Run `source /opt/openfoam13/etc/bashrc` again.

**Mesh hangs:** Should not occur in current version — OMP nondeterminism in `mapEdgeNodes` was fixed in this port. If you encounter a hang, please open a GitHub issue with your STL and meshDict.

---

## Tested on
- Ubuntu 22.04 WSL2, OpenFOAM Foundation v13
- NASA Rotor 37 two-passage turbomachinery sector with boundary layers

## NASA Rotor 37 boundary-layer mesh preview

<p align="center">
  <img src="NASA_Rotor_37_1.png" alt="NASA Rotor 37 cfMesh OpenFOAM 13 boundary-layer mesh preview 1" width="48%">
  <img src="NASA_Rotor_37_2.png" alt="NASA Rotor 37 cfMesh OpenFOAM 13 boundary-layer mesh preview 2" width="48%">
</p>

<p align="center">
  <img src="NASA_Rotor_37_3.png" alt="NASA Rotor 37 cfMesh OpenFOAM 13 boundary-layer mesh preview 3" width="48%">
  <img src="NASA_Rotor_37_4.png" alt="NASA Rotor 37 cfMesh OpenFOAM 13 boundary-layer mesh preview 4" width="48%">
</p>

<p align="center">
  <em>NASA Rotor 37 two-passage turbomachinery sector generated with the OpenFOAM 13 cfMesh port.</em>
</p>

---

## Original cfMesh
- Source: https://sourceforge.net/projects/cfmesh/
- Authors: Franjo Juretic, Philippose Rajan, Ivor Clifford

---
## Tutorial: meshing NASA Rotor 37

A complete worked example using the included Rotor 37 two-passage turbomachinery sector geometry.

### Prerequisites
- cfMesh built and sourced (see Quick start above)
- ParaView installed (for visualization)

### 1. Set up the case

```bash
mkdir -p ~/rotor37_mesh/constant/triSurface
mkdir -p ~/rotor37_mesh/system
cp /path/to/Rotor37_cfmesh.stl ~/rotor37_mesh/constant/triSurface/
```

### 2. Copy the included meshDict
```bash
cp /path/to/cfMesh/tutorials/rotor37/system/meshDict ~/rotor37_mesh/system/
```

The included `meshDict` is pre-configured for the Rotor 37 geometry with 3-layer blade BL and 5-layer hub/shroud BL.

### 3. Run the mesher

```bash
cd ~/rotor37_mesh
cartesianMesh 2>&1 | tee mesh.log
```

Expected runtime: ~60 seconds on a modern CPU. Watch for: (7950X3d in my case)
```
Finished renaming boundary patches
```

### 4. Check mesh quality

```bash
checkMesh -case ~/rotor37_mesh 2>&1 | grep -E "cells:|average non-orth|Failed|aspect|negative"
```

Expected output:
cells:            ~2.1M
Mesh non-orthogonality Max: ~100 degrees  average: 10.5 degrees
***High aspect ratio cells found, number of cells 3-5
***Zero or negative cell volume detected, Number of negative volume cells: 0-2
Failed 1-3 mesh checks.

### 5. Visualize in ParaView

```bash
paraFoam -case ~/rotor37_mesh &
```

Or create a `.foam` file and open manually:
```bash
touch ~/rotor37_mesh/rotor37.foam
# Open rotor37.foam in ParaView
# Apply → select patches: blade_1-4, hub, shroud
# Color by: p or U after running a solver
```
To inspect boundary layers: in ParaView, use **Filters → Extract Block** to isolate wall patches, then apply **Surface with Edges** representation. The 3–5 layer prism stacks on blade and endwall surfaces should be clearly visible.

---

## License and Attribution
This project is licensed under the GNU General Public License v3.0, consistent with the original cfMesh license.
Original cfMesh authorship and copyright remain with the original cfMesh authors, including Franjo Juretić / Creative Fields.
This fork includes OpenFOAM 13 compatibility updates, build fixes, and meshing robustness modifications by Mitch Stolk.
Modifications in this fork Copyright (C) 2026 Mitch Stolk.
If you use this fork in academic, research, commercial evaluation, or public engineering work, please cite or acknowledge this repository:
Mitchell Stolk, cfMesh-OpenFOAM13: OpenFOAM 13 port and robustness-focused fork of cfMesh, 2026.

## Disclaimer
Unofficial community port. Not affiliated with Creative Fields or the OpenFOAM Foundation.


