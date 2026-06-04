# cfMesh for OpenFOAM 13

v0.10 beta: active development on OpenFOAM Foundation v13. This fork focuses on OpenFOAM 13 compatibility, deterministic parallel behavior, boundary-layer robustness, and surface-conformance fixes for complex STL-driven Cartesian/polyhedral meshing. Rotor 37 is currently used as the primary stress-test case. Serial runs are deterministic; parallel repeatability is substantially improved, with minor residual variance still under investigation.

### NOTE: NASA Rotor 37 is used as a stress-test case because it combines tight blade/hub/shroud junctions, periodic sector boundaries, acute trailing-edge features, refinement transitions, and boundary-layer growth challenges. The intent is not to make the mesher Rotor-37-specific; Rotor 37 is used to expose general robustness issues that should benefit other external, internal, and turbomachinery CFD geometries. CFD validation is still limited and should be treated as preliminary.

---
### Known remaining issues
**Residual minor variance** Major OMP races fixed in v0.8–v0.9. Serial (OMP=1) runs fully deterministic. OMP=32 variance substantially reduced — primary races (updateType, otherVrts_, frontPoints) fixed. Minor residual variance remains under investigation.

**optimiseBoundaryVolumeOptimizer call-site bug (fixed v0.9)** Original code passed `true`/`false` as iteration count to a function expecting `(nIterations, nonShrinking)`. Early boundary motion was unconstrained instead of non-shrinking; late boundary fallback was doing zero iterations. Fixed to `(1, true)` and `(1, false)`.

**Visible surface protrusions at patch boundaries** — Partially resolved in v0.10. Root cause: volume optimizer moves boundary points off STL surface with no constraint. Surface-constrained optimizer (Layer 3) now reduces visible-scale drift (>100 micron points) by ~46%. Residual drift on periodic patches remains under active development.

**Runtime increase with surface-constrained optimizer** — Serial octree pre-pass in `partTetMesh::updateOrigMesh` adds ~2x wall-clock time vs unconstrained. Can be disabled via `constrainOptimizerBoundaryMotion false`. Full parallel implementation planned.

### Current limitations

- Rotor 37 is a stress test, not a general validation suite.
- CFD validation is preliminary.
- Boundary-layer generation is much improved but still under active development near acute blade/endwall/periodic singularities.
- Surface-constrained optimization reduces boundary drift but currently increases runtime.
- A small number of negative-volume cells, bad pyramid faces, and non-orthogonality errors may remain in current development runs.
- The mesher is not yet production-certified; users should run `checkMesh`, inspect boundary surfaces visually, and validate solver behavior for their own cases.

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
  - Current Rotor 37 results are promising but still under active development; this should not yet be treated as a production-certified mesher.
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
| `partTetMesh::updateVerticesSMP` updateType |= CELLCENTRE/FACECENTRE race on shared tet centre nodes | Serialize updateType marking with omp single — primary cause of OMP=32 vs serial quality divergence |
| `partTriMesh::updateVerticesSMP` / `updateVertices` updateType |= FACECENTRE/SMOOTH race | Serialize updateType marking |
| `boundaryLayersCreateVertices` concurrent std::map find/insert/[] on otherVrts_ during OMP parallel BL vertex creation | Serialize otherVrts_ access; add FatalError guard for missing map entries |
| `extrudeLayer` frontPoints[f[pI]] |= FRONTVERTEX race on shared point indices | Serialize front point marking — was corrupting BL front layer detection |
| `meshSurfaceEngineModifier` updateFaces/updateBndPoint shared boolList marking races | Serialize marking loops |
| `optimiseBoundaryVolumeOptimizer` call-site bug — true/false passed as nIterations not nonShrinking | Fix to (1,true) and (1,false) |
| `tetMeshOptimisation::optimiseUsingMeshUntangler` missing helper[tetI]=true — false convergence declaration | Add missing assignment |
| Untangle stagnation deadlock at small residual bad-face count | Detect stuck iterations, switch to unconstrained fallback after 4 unchanged iterations with ≤8 bad faces |
| `mapVerticesOntoSurface` unconstrained global nearest-triangle projection causes visible protrusions at patch boundaries | Patch-constrained projection: validate returned patch against `pointPatches`; re-project to correct patch if mismatch |
| Volume optimizer moves boundary points up to 1.3mm off STL surface with no constraint | Surface-constrained `partTetMesh::updateOrigMesh`: limited-displacement single-patch boundary re-projection during optimizer write-back; current implementation uses a serial octree query phase with parallel write-back. |
| `volumeOptimizerEvaluateGradients` `FatalError` on degenerate zero-volume tets | Replace fatal with counted skip + end-of-function report |
| `finalUntangleRejected_` flag incorrectly set when untangle skipped due to pre-existing negVol | Skip rollback check when untangle did not run; BL proceeds normally |


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

A topology-aware surface projection and repair pipeline, developed incrementally on top of the original cfMesh surface mapping code.

### Architecture fix patch-identity-first projection
The original cfMesh pipeline ran `mapMeshToSurface` (global surface projection) before `extractPatches` (patch identity assignment). This meant every boundary point was projected using a global nearest-triangle search with no patch constraint causing points near concave multi-patch corners to snap to the wrong patch surface, producing visible protrusions.

Fix: `extractPatches` now runs before `mapMeshToSurface`, establishing patch identity on the raw Cartesian hex mesh before projection. A second `extractPatches` call after projection corrects misassignments from the unprojected hex geometry. This addressed one class of patch-boundary protrusions. Later diagnostics showed that a larger remaining source was boundary-point drift during volume optimization, addressed in v0.10 by the surface-constrained optimizer.

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

### Version

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


### meshDict controls

| Parameter | Default | Description |
|-----------|---------|-------------|
| `blblFeatureAngleDeg` | 40 | Junction angle threshold for BL/BL sharp curve detection (degrees) |
| `blblCornerAcuteThreshold` | 0.3 | Dot product threshold for acute corner classification. 0.5=60°, 0.3=73°, 0.0=90° |
| `lockAcuteCornerPoints` | true | Lock acute corners in final surface optimizer |
| `acuteCornerRing0-6` | 0.0/0.0/0.05/0.15/0.35/0.60/1.00 | layerScale taper seeded from acute corner points |
| `acuteCornerCapLayers` | false | Cap layer count near acute corners (experimental) |
| `virtualTopologyExclusion` | false | Enable commercial-style BL exclusion face-rings at acute junctions |
| `virtualTopoRing0-2` | 0.05/0.15/0.40 | layerScale for VT face-ring suppression zones |
| `detectGaps` | false | Enable gap/contact proximity detection |
| `gapThreshold` | 5e-4 | Proximity distance threshold in metres |
| `gapPatchPairs` | none | Explicit patch pairs to check for proximity |
| `gapSuppressPatches` | none | Patches to suppress BL on (others are protected) |
| `gapFaceRingExclusion` | true | Enable face-ring BFS around gap points (suppresses prism topology) |
| `gapFaceRing0-2Scale` | 0.02/0.05/0.20 | layerScale for gap face-ring suppression zones |
| `layerScaleRing1-6` | 0.10/0.25/0.45/0.65/0.80/1.00 | BL/no-BL transition ramp scales |
| `terminationPatches` | none | Patches with no BL (inlet/outlet) |
| `postBLSnap` | false | Post-BL surface re-projection (experimental) |
| `postRefineBLOptimisation` | true | Run surface optimizer after BL refinement |


### BL/BL sharp-junction suppression
Points where two BL patches meet at a sharp angle (blade+hub, blade+shroud) are detected using patch normal angle comparison and suppressed. Prevents degenerate layer cells at blade-endwall junctions in turbomachinery geometry.

## Surface-constrained optimizer (v0.10)

Layer 3 surface-constrained optimizer architecture: boundary points are constrained to their STL patch surface during volume optimization, preventing drift accumulation.

### Architecture
The volume optimizer (`optimizeMeshFV`, `untangleMeshFV`, `optimizeLowQualityFaces`) previously moved boundary points freely in 3D with no surface awareness. Points could drift up to 1.3mm from the STL surface, causing visible protrusions that the downstream snap could not correct without inverting cells.

Fix: `partTetMesh::updateOrigMesh` now applies a limited-displacement re-projection for single-patch boundary points before committing optimizer results:

1. **Serial pre-pass**: query octree for each drifted boundary point (drift > `1e-5m`), compute limited target position (1% step toward STL)
2. **Parallel write-back**: apply precomputed targets in OMP — no octree access in parallel path

### Results on Rotor 37
- Visible-scale drift (>100 micron points): 2863 → 1540 (46% reduction)
- Max drift: 1.19mm → 0.52mm
- BL coverage: 95–99% across all patches

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
---
## What's new

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

### Parallel hardening and deterministic-optimizer work (v0.8–v0.9)
Many known race-sensitive OMP apply loops have been serialized or guarded across the meshLibrary. v0.9 adds: negativeNode tet-marking loops (primary source of OMP=32 nondeterminism), isExitFace/isBndPoint/removeBndPoint shared-write races, and 100+ globalToLocal guards across 30+ files. Serial (OMP=1) runs are now fully deterministic. Minor residual variance at OMP=32 under continued investigation.
---
## Mesh quality on NASA Rotor 37 (v0.9)

Tested with boundary layers active. Current development testing commonly uses 3-layer blade boundary layers and 3–5 layer hub/shroud configurations depending on the experiment.

| Metric | v0.10 (current) | v0.9 |
|--------|-----------------|------|
| Cells | ~1.77M | ~1.77M |
| Average non-orthogonality | ~12.5° | TBD |
| Non-ortho errors | 3–9 | TBD |
| Negative volume cells | 2–5 | TBD |
| Bad pyramid faces | 11–57 | TBD |
| Max skewness | ~8.5–11.6 | TBD |
| Visible-scale drift (>100µm pts) | ~1540 | ~2863 |
| Runtime (OMP=4) | ~580s wall | ~90s |
| Runtime (OMP=32) | ~204s wall | ~90s |
| Deterministic | OMP=4 stable | Improving |
---

## Version history

| Tag | Branch | Description |
|-----|--------|-------------|
| `v0.10` | `feature/bl-transition-topology` | Surface-constrained optimizer (Layer 3): patch-constrained projection in `mapVerticesOntoSurface`, limited-displacement `partTetMesh::updateOrigMesh` re-projection, two-phase serial/parallel write-back, post-optimizer incremental snap, pre-BL snap with displacement filter and pyramid-tolerant rollback. 46% reduction in visible-scale boundary drift. |
| `v0.9` | `feature/bl-transition-topology` | Comprehensive parallel hardening: 100+ globalToLocal guards, partTetMesh/partTriMesh updateType races, otherVrts_ std::map race in BL vertex creation, frontPoints extrusion race, meshSurfaceEngineModifier marking races, optimiseBoundaryVolumeOptimizer call-site bug, MeshUntangler convergence fix, untangle stagnation detection. Serial runs fully deterministic; OMP=32 quality substantially improved. |
| `v0.8.3` | `feature/bl-transition-topology` | Virtual topology face-ring BFS (commercial-style BL exclusion zones at acute junctions), VT/BL-ramp priority hierarchy, side-aware gap closure with face-level suppressLayerAtBndFace_ mask, configurable gapPatchPairs/gapSuppressPatches. 0 bad pyramids, skew 9.77. |
| `v0.8.2` | `feature/bl-transition-topology` | Geometry preprocessing (STL vertex welding), BL+BL+neutral corner classifier (acute/mild, angle-based), acute corner projection guard (`lockAcuteCornerPoints`), softer 4-ring taper. Zero bad pyramids, zero non-ortho errors. |
| `v0.8` | `feature/bl-transition-topology` | Full codebase audit 60+ bug fixes: OMP apply-loop serialization, globalToLocal guards, activePatch_ owner/patch bugs, face size guards, bounds guards throughout. Deterministic results on Rotor 37. |
| `v0.7.1` | `feature/bl-transition-topology` | BL/neutral extrusion constraint + crossing clamp + layer rollback scaffold |
| `v0.7` | `feature/bl-transition-topology` | BL/neutral edge detection + BL-side projection constraint + VTK diagnostics |
| `v0.6` | `feature/deterministic-optimizer` | Partial deterministic optimizer 8 unsafe OMP loops serialized |
| `v0.5` | `openfoam13` | Pure quality evaluator, thread-safe, tight skewness variance |
| `v0.4` | `openfoam13` | Targeted repair of validity-rejected points |
| `v0.35` | `openfoam13` | VTK diagnostics for validity-rejected points |
| `v0.2` | `openfoam13` | Skewness proxy in quality evaluator |
| `v0.1-beta-stable` | `openfoam13` | Pipeline reorder + relative validity thresholds |

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
[ 2/3 ] Running checkMesh...
PASS: Cell count: ~1.8 (baseline: 2100000 +/-10%)
PASS: Average non-orthogonality: ~10deg (limit: 15deg)
PASS: High aspect ratio cells: 3-5 (known limit: 10)
[ 3/3 ] Summary
PASSED: 6
FAILED: 0
RESULT: PASS
```

Run this after any code change to catch regressions before committing.

## Troubleshooting

**Build fails with missing header:** Make sure OpenFOAM 13 is sourced before building. Run `echo $WM_PROJECT_VERSION` should return `13`.

**`cartesianMesh` not found after build:** Check `$FOAM_USER_APPBIN` is in your PATH. Run `source /opt/openfoam13/etc/bashrc` again.

**Mesh hangs or stalls:** Several known OMP and topology-repair hangs have been fixed, but this is still an active-development mesher. If a case hangs, rerun with `OMP_NUM_THREADS=1` to separate serial logic from parallel nondeterminism, then open an issue with the STL, `meshDict`, thread count, and relevant log output.

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

Expected runtime depends strongly on whether `constrainOptimizerBoundaryMotion` is enabled. On a Ryzen 7950X3D / WSL2 test setup, unconstrained runs are significantly faster; v0.10 surface-constrained runs are currently slower because the octree projection pre-pass is still being optimized.

```
Finished renaming boundary patches
```

### 4. Check mesh quality

```bash
checkMesh -case ~/rotor37_mesh 2>&1 | grep -E "cells:|average non-orth|Failed|aspect|negative"
```

Expected output:
cells:            ~1.77M
Mesh non-orthogonality Max: ~89 degrees  average: ~8.6 degrees
*Number of severely non-orthogonal (> 70 degrees) faces: ~600
***Max skewness = ~8.2-8.4

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


