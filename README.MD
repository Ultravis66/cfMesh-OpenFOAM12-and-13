# cfMesh for OpenFOAM 13

> **WORK IN PROGRESS** — active development on OpenFOAM 13 only. Since this got complicated extremely fast, I will only be working on 13 going forward. If you would like to patch 12 with the changes I made to 13, feel free to do it.

A working port of cfMesh to OpenFOAM Foundation version 13, with significant boundary layer robustness improvements beyond the original cfMesh.

Original cfMesh developed by Franjo Juretic / Creative Fields, licensed under GPL v3.

Tested on **NASA Rotor 37** two-passage turbomachinery sector geometry with boundary layers.

---

## Branch status

| Branch | OpenFOAM Version | Status |
|--------|-----------------|--------|
| `master` | OpenFOAM 13 | **Active development** |
| `openfoam12-archive` | OpenFOAM 12 | Archived — no further updates |

**OpenFOAM 12 support has been dropped.** All development is now on OpenFOAM 13. The OF12 branch is preserved as-is for reference.

---

## What works
- `cartesianMesh` — Cartesian hex-dominant mesher, tested on NASA Rotor 37
- `tetMesh` — Tetrahedral mesher (experimental)
- `pMesh` — Polyhedral mesher (experimental)
- `libmeshLibrary.so` — Core mesh library

---

## Build instructions

```bash
source /opt/openfoam13/etc/bashrc
git clone https://github.com/Ultravis66/cfMesh-OpenFOAM13.git cfMesh
cd cfMesh/meshLibrary && wmake libso
cd ../executables/cartesianMesh && wmake
cd ../executables/tetMesh && wmake
cd ../executables/pMesh && wmake
```

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

### BL/BL sharp-junction suppression
Points where two BL patches meet at a sharp angle (blade+hub, blade+shroud) are detected using patch normal angle comparison and suppressed. Prevents degenerate layer cells at blade-endwall junctions in turbomachinery geometry.

### refineFace graceful handling
`refineBoundaryLayers::refineFace` previously called `FatalError` when no split direction could be found for transition-zone faces. Now returns the original face unrefined instead of crashing.

### Two-pass candidate quality rollback (infrastructure complete, detection pending)
Architecture for detecting degenerate layer cells and forcing single-layer treatment at source boundary faces:
- `mixedDisplacementFaces_` member on `boundaryLayers`
- `forceSingleLayerFaces_` + `forceSingleLayerAtFaces()` on `refineBoundaryLayers`
- Layer count smoothing on boundary-face adjacency graph (max jump of 1 between neighbors)
- Phase 1 non-mutating diagnostic in progress

### Known remaining issue
13 high-aspect-ratio cells at blade-endwall junctions (`1.68119e+147` aspect ratio, `-9.44996e-16` negative volume — machine epsilon). These are at blade+hub and blade+shroud geometric singularities where simple prism extrusion is topologically unsolvable without butterfly/O-grid topology. Two-pass quality rollback is the planned fix.

---

## Mesh quality on NASA Rotor 37

| Metric | Value |
|--------|-------|
| Cells | ~2.17M |
| Average non-orthogonality | 10.5° |
| Bad faces after untangling | 0 |
| Runtime | ~60s |

---

## Tested on
- Ubuntu 22.04 WSL2, OpenFOAM 13
- NASA Rotor 37 two-passage turbomachinery sector with boundary layers

---

## Original cfMesh
- Source: https://sourceforge.net/projects/cfmesh/
- Authors: Franjo Juretic, Philippose Rajan, Ivor Clifford

---

## Disclaimer
Unofficial community port. Not affiliated with Creative Fields or the OpenFOAM Foundation.
