cat > ~/cfmesh-repo/README.MD << 'EOF'
# cfMesh OpenFOAM 12 & 13 Port

A working port of cfMesh to OpenFOAM 12 and 13 (openfoam.org).
Original cfMesh by Franjo Juretic, Creative Fields. Licensed under GPL v3.

Successfully tested on **NASA Rotor 37** turbomachinery passage geometry.

---

## What works
- `cartesianMesh`  Cartesian hex-dominant mesher
- `tetMesh`  Tetrahedral mesher
- `pMesh`  Polyhedral mesher
- `libmeshLibrary.so`  Core mesh library

---

## Branches

| Branch | OpenFOAM Version | Status |
|--------|-----------------|--------|
| `master` | OpenFOAM 12 | ✅ Working |
| `openfoam13` | OpenFOAM 13 | ✅ Working |

---

## Build instructions

**OpenFOAM 12:**
source /opt/openfoam12/etc/bashrc
git clone -b master https://github.com/Ultravis66/cfMesh-OpenFOAM12-and-13.git cfMesh
cd cfMesh/meshLibrary && wmake libso
cd ../executables/cartesianMesh && wmake
cd ../executables/tetMesh && wmake
cd ../executables/pMesh && wmake

**OpenFOAM 13:**
source /opt/openfoam13/etc/bashrc
git clone -b openfoam13 https://github.com/Ultravis66/cfMesh-OpenFOAM12-and-13.git cfMesh
cd cfMesh/meshLibrary && wmake libso
cd ../executables/cartesianMesh && wmake
cd ../executables/tetMesh && wmake
cd ../executables/pMesh && wmake

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
| `unallocLabelList` | `labelList` (fixed in FMSToVTK and mergeSurfacePatches) |

---

## Runtime fixes (this port goes beyond just compiling)

Previous ports of cfMesh to modern OpenFOAM compile successfully but crash or
hang at runtime. This port fixes the following runtime bugs:

| Bug | File | Fix |
|-----|------|-----|
| `face::centre()` returns garbage for cfMesh face ordering → phantom 1946 inverted faces | `meshSurfaceEngineCalculateBoundaryNodesAndFaces.C`, `meshSurfaceEngineModifier.C`, `edgeExtractor.C` | Replace with simple vertex average |
| Inversion check runs on unmapped octree mesh → infinite untangler loop | `cartesianMeshGenerator.C` | Map vertices onto surface BEFORE `preMapVertices()` |
| `preMapVertices()` internal `untangleSurface` hangs on complex geometry | `cartesianMeshGenerator.C`, `meshSurfaceMapperPremapVertices.C` | Use `preMapVertices(0)` to skip internal untangle |
| `findLeafContainingVertex` throws `char const*` → uncaught crash | `meshOctreeNeighbourSearches.C` | Replace throw with `return -1` |
| `newPositionSurfaceOptimizer` produces out-of-domain vertex positions | `meshSurfaceOptimizerI.H` | Add displacement bound check and vecY near-zero guard |
| `pointFieldPMG::setSize` over-allocates without zero-init → garbage points in output | `pointFieldPMGI.H` | Zero-initialize over-allocated region |
| Octree points uninitialized for unreferenced nodes | `meshOctreeAddressingCreation.C` | Zero-initialize `octreePointsPtr_` |
| `checkFacePatchesGeometry` unbounded `do-while` loop → infinite loop | `edgeExtractor.C` | Add iteration limit (`nGeomIter < 5`) |
| `findNearestSurfacePoint` 100-iteration expansion → hang for out-of-STL vertices | `meshOctreeFindNearestSurfacePoint.C` | Reduce limit to 10 |
| `optimizePoint` FPE from zero `funcBefore` division | `surfaceOptimizer.C` | Guard zero divisions |
| `untangleSurface` global iteration limit too high for complex geometry | `meshSurfaceOptimizerOptimizeSurface.C` | Limit to 2 global iterations |

---

## Tested On
- Ubuntu 22.04 WSL2 with OpenFOAM 12 and 13
- NASA Rotor 37 two-passage turbomachinery sector (232,000+ cells)
- Unit cube validation case

---

## Original cfMesh
- Source: https://sourceforge.net/projects/cfmesh/
- Authors: Franjo Juretic, Philippose Rajan, Ivor Clifford

---

## Disclaimer
This is an unofficial community port and is not affiliated with Creative Fields
or the OpenFOAM Foundation.
EOF
