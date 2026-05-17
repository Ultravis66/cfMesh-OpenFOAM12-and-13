# cfMesh for OpenFOAM 13

> **WORK IN PROGRESS** — active development on OpenFOAM 13 only. To keep development focused and maintainable, active development is now limited to OpenFOAM Foundation v13. If you would like to patch OF12 with the changes made here, feel free — the `openfoam12-archive` branch is preserved as-is.

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

### BL/BL sharp-junction suppression
Points where two BL patches meet at a sharp angle (blade+hub, blade+shroud) are detected using patch normal angle comparison and suppressed. Prevents degenerate layer cells at blade-endwall junctions in turbomachinery geometry.

### Geometry preservation pipeline
Commercial-mesher-style surface projection with validity checking:

- **Targeted centroid point projection** — `correctEdgesBetweenPatches::patchCorrection()` creates new centroid points when decomposing boundary faces at feature edges. These 89 points are now tracked and projected onto the STL surface immediately after creation.
- **Validity-checked surface projection** — `mapVerticesOntoSurface()` now stores old vertex positions before each projection pass and runs a serial validity check afterward. Any move that inverts an adjacent face pyramid is reverted. ~8000 invalid moves reverted per run.

- **Validity-checked edge snapping** — `mapEdgeNodes()` stores old positions and runs a serial validity check after projection. Invalid moves (130+ per run) are reverted before geometry update.
- **Corner-aware projection** — `mapCorners()` snaps 3+ patch intersection vertices toward the nearest STL corner using a damped relaxation factor (default 0.25, tunable via `cornerSnapRelaxation` in meshDict). A 1e-10 pyramid height threshold catches near-degenerate faces that SMALL misses. 20 invalid corner moves reverted per run.
- **Ordered snap sequence** — corners → edges → surface, each with validity rollback. Corners are processed first as the most constrained vertices, followed by two-patch edge points, then surface vertices.
- **`blblJunctionPoints_` wiring** — Junction point metadata flows end-to-end from `boundaryLayers` → `cartesianMeshGenerator` → `refineBoundaryLayers`, enabling future wedge topology emission.
- **Post-BL snap infrastructure** — Local octree rebuild after BL generation for targeted outer boundary projection. Currently disabled pending tuning; controlled via `postBLSnap` in meshDict.

### refineFace graceful handling
`refineBoundaryLayers::refineFace` previously called `FatalError` when no split direction could be found for transition-zone faces. Now returns the original face unrefined instead of crashing.

### Two-pass candidate quality rollback (experimental infrastructure)
Initial infrastructure exists for detecting degenerate layer-cell candidates and reducing local layer complexity, but the detection/rollback policy is still under development:
- `mixedDisplacementFaces_` member on `boundaryLayers`
- `forceSingleLayerFaces_` + `forceSingleLayerAtFaces()` on `refineBoundaryLayers`
- Layer count smoothing on boundary-face adjacency graph (max jump of 1 between neighbors)
- Phase 1 non-mutating diagnostic in progress

### Known remaining issue
2 high-aspect-ratio cells at blade-endwall junctions (reduced from 13 via geometry preservation pipeline). (`1.68119e+147` aspect ratio, `-9.44996e-16` negative volume — machine epsilon). These occur at blade+hub and blade+shroud geometric singularities where the BL/BL junction feature curve is physically singular. The planned fix is proper wedge/pyramid topology emission in `generateNewCellsPrism` at sharp BL/BL junctions — architectural groundwork (C1 detection, layerScale suppression, C3 degenerate edge handling) is already in place on the `openfoam13` branch.

---

## Mesh quality on NASA Rotor 37

| Metric | Value |
|--------|-------|
| Cells | ~2.1M |
| Average non-orthogonality | 10.5° |
| High aspect ratio cells | 2 (down from 13) |
| Negative volume cells | 2 (down from 7) |
| Failed mesh checks | 5 |
| Runtime | ~75s |
| Bad faces after untangling | 0 |
| Regression test | 6/6 PASS |
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
PASS: High aspect ratio cells: 2 (known limit: 13)
PASS: Negative volume cells: 2 (known limit: 10)
[ 3/3 ] Summary
PASSED: 6
FAILED: 0
RESULT: PASS
```

Run this after any code change to catch regressions before committing.

## Troubleshooting

**Build fails with missing header:** Make sure OpenFOAM 13 is sourced before building. Run `echo $WM_PROJECT_VERSION` — should return `13`.

**`cartesianMesh` not found after build:** Check `$FOAM_USER_APPBIN` is in your PATH. Run `source /opt/openfoam13/etc/bashrc` again.

**Mesh hangs:** cfMesh uses OpenMP. Set `OMP_NUM_THREADS=4` (or your core count) before running.

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
***High aspect ratio cells found, number of cells 2
***Zero or negative cell volume detected, Number of negative volume cells: 2
Failed 5 mesh checks.

The 2 high-aspect and 2 negative-volume cells are a known limitation

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

## License

This project is licensed under the GNU General Public License v3.0, consistent with the original cfMesh license.

## Disclaimer
Unofficial community port. Not affiliated with Creative Fields or the OpenFOAM Foundation.


