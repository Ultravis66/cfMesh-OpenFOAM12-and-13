---

## API changes from original cfMesh (compile fixes)

| Original | Updated |
|----------|---------|
| `Pstream::blocking/scheduled/nonBlocking` | `Pstream::commsTypes::blocking/scheduled/nonBlocking` |
| `eigenVector(T, scalar)` | `eigenVectors(T)` with direct component access |
| `circumCentre()/circumRadius()` | `circumSphere().first()/.second()` |
| `isFile()` | `std::filesystem::exists()` requires `#include <filesystem>` and `-std=c++17` |
| `rm/isDir/rmDir` | requires `#include "OSspecific.H"` in affected files |
| `edgeMesh` library | include paths updated to use `meshTools` |
| `unallocLabelList` | `labelList` |

---

## Runtime fixes

Previous ports compile but crash or hang at runtime. This port fixes:

| Bug | File(s) | Fix |
|-----|---------|-----|
| `face::centre()` garbage for cfMesh face ordering → phantom inversions | `meshSurfaceEngineCalculateBoundaryNodesAndFaces`, `meshSurfaceEngineModifier`, `edgeExtractor`, `decomposeFaces`, many others | Replace with vertex average everywhere |
| `face::normal()` garbage → wrong BL extrusion direction | `extrudeLayer`, `boundaryLayersCreateVertices`, `boundaryLayerOptimisationNormals`, `symmetryPlaneOptimisation`, many others | Replace with triangle fan normal everywhere |
| Inversion check on unmapped octree mesh → infinite loop | `cartesianMeshGenerator` | Map vertices onto surface BEFORE `preMapVertices()` |
| `preMapVertices()` internal `untangleSurface` hangs | `cartesianMeshGenerator`, `meshSurfaceMapperPremapVertices` | Use `preMapVertices(0)` |
| `findLeafContainingVertex` throws `char const*` | `meshOctreeNeighbourSearches` | Replace throw with `return -1` |
| `newPositionSurfaceOptimizer` garbage vertex positions | `meshSurfaceOptimizerI` | Displacement bound + vecY near-zero guard |
| `pointFieldPMG::setSize` uninitialized over-allocation | `pointFieldPMGI` | Zero-initialize over-allocated region |
| Octree points uninitialized | `meshOctreeAddressingCreation` | Zero-initialize `octreePointsPtr_` |
| `checkFacePatchesGeometry` unbounded loop | `edgeExtractor` | Iteration limit (`nGeomIter < 5`) |
| `findNearestSurfacePoint` 100-iteration hang | `meshOctreeFindNearestSurfacePoint` | Reduce to 10 iterations |
| `optimizePoint` FPE from zero division | `surfaceOptimizer` | Guard zero divisions |
| `untangleSurface` too many global iterations | `meshSurfaceOptimizerOptimizeSurface` | Limit to 2 global iterations |
| Garbage denormalized points in output | `cartesianMeshGenerator` | Sanitize before write |

---

## Boundary layer improvements (beyond original cfMesh)

These are new capabilities not present in the original cfMesh:

### Feature-aware vertex normals
Replaced simple face-normal averaging with angle-weighted normals that respect sharp features. At blade trailing edges, leading edges, and blade-endwall junctions, the new code detects the feature angle and uses only the dominant-patch normals instead of averaging across the sharp edge. Result: **96% reduction in bad BL cells** on NASA Rotor 37 (229,000 → 9,646).

### Gap distance detection
For each wall vertex being extruded, ray-casts along the surface normal to detect opposing surfaces in the 2-ring face neighborhood. Limits BL thickness to 50% of the detected gap distance. Prevents cell overlap in narrow passages such as blade trailing edges and blade-endwall junctions.

### Concave edge layer termination
Implemented the previously stub-only `terminateLayersAtConcaveEdges` flag using `meshSurfaceCheckEdgeTypes`. Identifies concave edge vertices and skips BL extrusion there. Enabled by default in `cartesianMeshGenerator`.

---

## Tested on
- Ubuntu 22.04 WSL2, OpenFOAM 12 and 13
- NASA Rotor 37 two-passage turbomachinery sector with boundary layers
- Unit cube validation case

---

## Original cfMesh
- Source: https://sourceforge.net/projects/cfmesh/
- Authors: Franjo Juretic, Philippose Rajan, Ivor Clifford

---

## Disclaimer
Unofficial community port. Not affiliated with Creative Fields or the OpenFOAM Foundation.
EOF

cd ~/cfmesh-repo
git add README.MD
git commit -m "Update README with BL improvements and full fix table"
git push origin master
git checkout openfoam13
git cherry-pick HEAD~0 2>/dev/null || git merge master --no-edit
git push origin openfoam13
git checkout master
