# cfMesh for OpenFOAM 13

v0.8.2 BETA: active development on OpenFOAM Foundation v13. Core OMP apply-loop races fixed in v0.8; major Rotor 37 quality metrics are now repeatable run-to-run, with minor residual variance still under investigation. 

### NOTE: NASA Rotor 37 is being used as a stress-test case because it is a sufficiently complex turbomachinery geometry, with tight blade/hub/shroud junctions, periodic sector boundaries, and boundary-layer growth challenges. The intent is not to make the mesher Rotor-37-specific, but to use this geometry to expose robustness issues that should improve performance on other external and internal CFD geometries. Limited testing on CFD, a few tests have been run just to see if the mesh would work. 


### Known remaining issues
**Residual minor variance** Core OMP apply-loop nondeterminism fixed in v0.8. Small run-to-run variance in severely non-ortho face count (~1,650–1,760) remains from dynamic scheduling in the untangling termination logic. Under investigation.

Sharp corners where blade, Hub/Shroud/Periodic boundries meet create jank protrusions, currently working on virtual tropology fix for this.
---
Original cfMesh developed by Franjo Juretic / Creative Fields, licensed under GPL v3.

Tested on **NASA Rotor 37** two-passage turbomachinery sector geometry with boundary layers.

---

## Branch status

| Branch | OpenFOAM Version | Status |
|--------|-----------------|--------|
| `master` | OpenFOAM 13 | **Active development** |
| `openfoam12-archive` | OpenFOAM 12 | Archived, no further updates |

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
- `cartesianMesh` Cartesian hex-dominant mesher, primary supported executable, tested on NASA Rotor 37
- `tetMesh` Builds successfully, experimental, not fully validated
- `pMesh` Builds successfully, experimental, not fully validated
- `libmeshLibrary.so` Core mesh library

---

## API changes from original cfMesh (compile fixes)

| Original | Updated |
|----------|---------|
| `Pstream::blocking/scheduled/nonBlocking` | `Pstream::commsTypes::blocking/scheduled/nonBlocking` |
| `eigenVector(T, scalar)` | `eigenVectors(T)` with direct component access |
| `circumCentre()/circumRadius()` | `circumSphere().first()/.second()` |
| `isFile()` | `std::filesystem::exists()` requires filesystem include and c++17 |
| `rm/isDir/rmDir` | requires `OSspecific.H` include in affected files |
| edgeMesh library | include paths updated to use meshTools |
| `unallocLabelList` | `labelList` |

---

## Runtime fixes

| Bug | Fix |
|-----|-----|
| `face::centre()` garbage for cfMesh face ordering phantom inversions | Replace with vertex average everywhere |
| `face::normal()` garbage wrong BL extrusion direction | Replace with triangle fan normal everywhere |
| Inversion check on unmapped octree mesh infinite loop | Map vertices onto surface BEFORE `preMapVertices()` |
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
| `mapEdgeNodes` OMP nondeterminism intermittent hangs and bad face explosions at full thread count | Deterministic `parallelBndNodes` build: compute in parallel, move and sync in serial |
| `smoothLaplacianFC` apply loop OMP race concurrent writes to shared mesh point storage | Serialize apply loop; compute remains parallel |
| `smoothSurfaceOptimizer` apply loop OMP race | Serialize apply loop; compute remains parallel |
| `nodeDisplacementLaplacianParallel` constructs `meshSurfaceEngineModifier` per point inside loop | Hoist modifier outside loop; serialize apply |
| `nodeDisplacementLaplacianFCParallel` apply loop OMP race | Serialize apply loop |
| `edgeNodeDisplacementParallel` apply loop OMP race + `std::map::operator[]` mutating map during parallel read | Serialize apply loop; replace `operator[]` with `find()` |
| `lp.coordinates() / lp.pointLabel()` divide-by-zero in parallel Laplacian smoothers | Guard `pointLabel() <= 0`; fallback to original position where needed |
| `globalToLocal[key]` unguarded on received parallel data possible silent map insertion / invalid lookup | Add `globalToLocal.found()` guards with correct packet consumption in receive loops |
| `calculateBoundaryOwners` uses patch 0 start offset regardless of `activePatch_` | Use `boundaries[activePatch_].patchStart()` when `activePatch_ >= 0` |
| `calculateBoundaryFacePatches` iterates all patches when `activePatch_ >= 0` out-of-bounds write risk | Add `activePatch_` branch; fill selected patch label only |
| `moveBoundaryVertex` dereferences `faceNormalsPtr_` assuming non-null when `pointNormalsPtr_` is set | Guard with `pointNormalsPtr_ && faceNormalsPtr_` |
| Face normal/centre computation has no size guard on degenerate faces | `bf.size() < 3` gives zero normal; `bf.empty()` gives zero centre |
| `bp[e.start()]` / `bp[point]` used without `-1` guard throughout edge addressing | Add `bps < 0` / `bpI < 0` continue guards |
| `checkAndFixCellGroupsAtBndVertices` performs topology surgery under OMP dynamic scheduling nondeterministic face connectivity | Serialize repair loop when modifying topology |
| `checkEdgeFaceConnections` / `checkFaceGroupsAtBndVertices` missing `clearAddressingData()` after `removeCells()` | Add `mesh_.clearAddressingData()` after each `removeCells()` |
| `calcGlobalBoundaryPointLabels` receive loops use remote face/point indices without bounds checks | Guard face index, point index, and `bp[]` result before use |
| `findLeafContainingVertex` throws raw `char const*` for points near octree boundary — process terminates on fine meshes | Replace `FatalErrorIn` + `throw` with `WarningIn` + `return -1`; all callers already handle `-1` gracefully |
| Near-coincident STL vertices below mesh resolution cause protrusions at sharp trailing-edge corners | meshDict-controlled geometry preprocessing pass: octree-based vertex weld before octree generation (geometryPreprocessing { active true; weldTolerance 1e-4; }) |
| BL+BL+neutral 3-way corners (blade/shroud/periodic, blade/hub/periodic) extruded with averaged bisector normal causing protrusions | Acute corner classifier measures patch-normal dot product; acute corners (dot < threshold) get 0.02 layerScale + projection guard; mild corners get 0.15 + ring taper |


---

## Boundary Layer Generation Improvements

### Patch-based BL control
Patches with `nLayers 0` (inlet, outlet, periodic) are fully excluded from BL generation. Previously cfMesh required hardcoded patch names now reads directly from `meshDict boundaryLayers/patchBoundaryLayers`. Three code paths fixed: `findPatchesToBeTreatedTogether`, `addLayerForPatch`, `addLayerForAllPatches`.

### BL/BL junction topology infrastructure (C1)

Persistent metadata for sharp BL/BL junction points is now carried through the boundary-layer generation pipeline:

- `blblJunctionPoints_` `labelHashSet` on `boundaryLayers` capturing boundary points on sharp BL/BL feature curves.
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

### Architecture fix patch-identity-first projection
The original cfMesh pipeline ran `mapMeshToSurface` (global surface projection) before `extractPatches` (patch identity assignment). This meant every boundary point was projected using a global nearest-triangle search with no patch constraint causing points near concave multi-patch corners to snap to the wrong patch surface, producing visible protrusions.

Fix: `extractPatches` now runs before `mapMeshToSurface`, establishing patch identity on the raw Cartesian hex mesh before any projection occurs. A second `extractPatches` call after projection corrects any misassignments from the unprojected hex geometry. This is the primary fix for surface protrusions.

### Boundary point topology classifier
Every boundary point is classified into one of five topology classes before projection:

| Class | Description |
|-------|-------------|
| `SINGLE_PATCH` | Interior surface point standard projection |
| `TWO_PATCH_EDGE` | Feature curve point constrained to edge |
| `MULTI_PATCH_CORNER` | 3+ patch junction constrained to corner endpoint |
| `BL_NOBL` | BL/no-BL transition locked |
| `NONMANIFOLD` | Non-manifold geometry locked and reported |

VTK diagnostic output available for each class for ParaView inspection.

### Feature curve extraction and corner snap
Dual-mode feature curve extraction (face-patch data when available, pointPatches fallback for early pipeline stages). Extracts ~10,260 feature edge segments. 24/24 multi-patch corner points snap correctly to nearest incident feature curve endpoint.

### Pure quality evaluator `proposedMoveIsValid`
Reusable per-point quality gate used by all movement logic. Evaluates a proposed new position purely geometrically no mesh mutation, no `const_cast`, thread-safe by construction.

Checks per incident face:
1. Movement distance cap (relative to local cell size)
2. Face area (reject near-zero collapse)
3. Pyramid height (relative to face area scale-invariant)
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
| `v0.5` | Pure quality evaluator no const_cast, thread-safe, tight variance |
| `v0.6` | Partial deterministic optimizer — 8 unsafe OMP loops serialized |
| `v0.7` | BL/neutral edge detection + BL-side projection constraint + VTK diagnostics |
| `v0.7.1` | BL/neutral extrusion constraint + crossing clamp + layer rollback scaffold |
| `v0.8` | Full codebase audit — 60+ bug fixes: OMP serialization, globalToLocal guards, activePatch_ bugs, face size guards, bounds guards throughout |
| `v0.8.1` | `feature/bl-transition-topology` | BL/BL junction extrusion fix: per-patch conservative normal selection in createNewVertex. Max skewness drops from 11.1 to 8.74. |
| `v0.8.2` | `feature/bl-transition-topology` | Geometry preprocessing: meshDict-controlled near-coincident vertex welding. BL+BL+neutral corner classification. Bad pyramids 6-9→1, non-ortho errors 1→0. |

### BL/no-BL transition edge constraint system

Generic protection system for points where BL patches meet no-BL patches (wall meets inlet/outlet/periodic). Fully geometry-agnostic driven entirely by `nLayersForPatch_` from meshDict, no patch names hardcoded.

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
| `geometryPreprocessing { active true; weldTolerance 1e-4; }` | off | Pre-mesh STL near-coincident vertex welding. Fixes knife-edge TE defects below mesh resolution. |
| `blblCornerAcuteThreshold` | 0.3 | Dot product threshold for acute BL+BL+neutral corner classification. 0.3=73°, 0.0=90°, -0.3=107°. |
| `lockAcuteCornerPoints` | true | Lock acute BL+BL+neutral corners in final surface optimizer. Prevents post-BL smoothing from pushing corners across patch boundaries. |
| `layerScaleRing1` | 0.10 | Updated default — softer taper start |
| `layerScaleRing2` | 0.30 | Updated default |
| `layerScaleRing3` | 0.60 | Updated default |

### BL/BL sharp-junction suppression
Points where two BL patches meet at a sharp angle (blade+hub, blade+shroud) are detected using patch normal angle comparison and suppressed. Prevents degenerate layer cells at blade-endwall junctions in turbomachinery geometry.

### refineFace graceful handling
`refineBoundaryLayers::refineFace` previously called `FatalError` when no split direction could be found for transition-zone faces. Now returns the original face unrefined instead of crashing.

### new 3-gate rollback in development
Initial infrastructure exists for detecting degenerate layer-cell candidates and reducing local layer complexity, but the detection/rollback policy is still under development:
- `mixedDisplacementFaces_` member on `boundaryLayers`
- `forceSingleLayerFaces_` + `forceSingleLayerAtFaces()` on `refineBoundaryLayers`
- Layer count smoothing on boundary-face adjacency graph (max jump of 1 between neighbors)
- Phase 1 non-mutating diagnostic in progress

**Blade TE singularity** One residual high-aspect cell and one near-zero volume cell at the blade trailing-edge geometric singularity (knife-edge meeting shroud+periodic). This is a CAD/STL limitation — the feature is below mesh resolution. The geometry preprocessing weld mitigates but cannot fully resolve a true zero-thickness trailing edge.

**BL/BL feature-edge rounding** Blade/hub and blade/shroud feature curves can be slightly rounded by the surface smoother. Fix in progress: feature-curve locking during smooth passes.

> **Diagnostic note:** SnappyHexMesh produces clean geometry at the same patch junctions using the same STL, confirming the protrusions are cfMesh-specific and not a geometry defect.

---
## What's new in v0.7.1

### Patch-role classification system
Patches are now classified into three roles driven by meshDict:
- `BL_PATCH` wall patches with active boundary layers
- `TERMINATION_PATCH` inlet/outlet (explicit via `terminationPatches`)  
- `NEUTRAL_PATCH` periodic/symmetry (no BL suppression triggered)

```
boundaryLayers
{
    terminationPatches (inlet outlet);
    ...
}
```

### BL/neutral edge constraint (v0.7)
Blade/periodic junction points where a BL wall patch meets a periodic/symmetry patch are now fully constrained:
- **Projection**: locked to BL-side patch only via `findNearestSurfacePointInRegion`
- **Extrusion normal**: recomputed from BL-side faces only, preventing tilt toward periodic plane
- **Crossing clamp**: extruded points that cross the periodic plane are clamped back to surface
- **VTK diagnostic**: `blNeutralEdgePoints_predetect.vtk` written per run for ParaView inspection

### Topology-aware BL suppression
Corner points (3+ patches) at BL/termination junctions are fully suppressed (`layerScale_ = 0.0`). Sharp two-patch BL/termination edges are suppressed based on patch normal angle (75° default threshold). Periodic/neutral patches never trigger suppression regardless of angle.

### Full deterministic optimizer (v0.8)
All known race-sensitive OMP apply loops in the surface and mesh optimizer have been serialized. Major run-to-run quality spikes have been eliminated on Rotor 37; minor residual variance in severe non-orthogonality counts remains under investigation.
---
## Mesh quality on NASA Rotor 37 (v0.8)

Tested with boundary layers active (3-layer blade, 5-layer hub/shroud).

| Metric | v0.8.2 current | v0.8.1 | v0.8 | v0.7.1 best case |
|--------|----------------|--------|------|------------------|
| Cells | ~2.28M | ~2.28M | ~2.28M | ~2.1M |
| Average non-orthogonality | 9.88° | 9.86° | 9.86° | 10.6° |
| Max non-orthogonality | ~90° | ~90° | ~90° | 90–100° |
| Non-ortho errors | 0 | 1 | 1 | 0–1 |
| Negative volume cells | 1 | 1 | 1 | 0–2 |
| Bad pyramid faces | **0** | 6 | 8–11 | 1–9 |
| Max skewness | **8.74** | 8.74 | 11.1 | 9.3–9.5 |
| Severely non-ortho faces (>70°) | ~1,696 | ~1,960 | ~1,700 | 1,800–2,100 |
| Runtime | ~83s | ~83s | ~85–91s | ~75–100s |
| Deterministic / repeatable major metrics | Yes | Yes | Yes | No |

---

## Version history

| Tag | Branch | Description |
|-----|--------|-------------|
| `v0.8` | `feature/bl-transition-topology` | Full codebase audit 60+ bug fixes: OMP apply-loop serialization, globalToLocal guards, activePatch_ owner/patch bugs, face size guards, bounds guards throughout. Deterministic results on Rotor 37. |
| `v0.7.1` | `feature/bl-transition-topology` | BL/neutral extrusion constraint + crossing clamp + layer rollback scaffold |
| `v0.7` | `feature/bl-transition-topology` | BL/neutral edge detection + BL-side projection constraint + VTK diagnostics |
| `v0.6` | `feature/deterministic-optimizer` | Partial deterministic optimizer 8 unsafe OMP loops serialized |
| `v0.5` | `openfoam13` | Pure quality evaluator, thread-safe, tight skewness variance |
| `v0.4` | `openfoam13` | Targeted repair of validity-rejected points |
| `v0.35` | `openfoam13` | VTK diagnostics for validity-rejected points |
| `v0.2` | `openfoam13` | Skewness proxy in quality evaluator |
| `v0.1-beta-stable` | `openfoam13` | Pipeline reorder + relative validity thresholds |
| `v0.8.2` | `feature/bl-transition-topology` | Geometry preprocessing (STL vertex welding), BL+BL+neutral corner classifier (acute/mild, angle-based), acute corner projection guard (`lockAcuteCornerPoints`), softer 4-ring taper. Zero bad pyramids, zero non-ortho errors. |

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
PASS: cartesianMesh completed in ~85–95s
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

**Build fails with missing header:** Make sure OpenFOAM 13 is sourced before building. Run `echo $WM_PROJECT_VERSION` should return `13`.

**`cartesianMesh` not found after build:** Check `$FOAM_USER_APPBIN` is in your PATH. Run `source /opt/openfoam13/etc/bashrc` again.

**Mesh hangs:** Should not occur in current version OMP nondeterminism in `mapEdgeNodes` was fixed in this port. If you encounter a hang, please open a GitHub issue with your STL and meshDict.

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

Expected runtime: ~85–95 seconds on a Ryzen 7950X3D / WSL2 test setup.
```
Finished renaming boundary patches
```

### 4. Check mesh quality

```bash
checkMesh -case ~/rotor37_mesh 2>&1 | grep -E "cells:|average non-orth|Failed|aspect|negative"
```

Expected output:
cells:            ~2.1M
Mesh non-orthogonality Max: ~90 degrees  average: ~9.9 degrees
***High aspect ratio cells found, number of cells 3
***Zero or negative cell volume detected, Number of negative volume cells: 1
Failed 1 mesh checks.

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
If this fork is useful in academic, research, commercial evaluation, or public engineering work, attribution to this repository is appreciated.
Mitchell Stolk, cfMesh-OpenFOAM13: OpenFOAM 13 port and robustness-focused fork of cfMesh, 2026.

## Disclaimer
Unofficial community port. Not affiliated with Creative Fields or the OpenFOAM Foundation.


