# cfMesh OpenFOAM 12 & 13 Port

A working port of cfMesh to OpenFOAM Foundation versions 12 and 13.

Original cfMesh was developed by Franjo Juretic / Creative Fields and is licensed under GPL v3.

This port has been successfully built and tested on:
- OpenFOAM 12
- OpenFOAM 13
- Ubuntu 22.04 / WSL2
- Unit cube validation case
- NASA Rotor 37 two-passage turbomachinery sector geometry

---

## What works
- `cartesianMesh` Cartesian hex-dominant mesher
- `tetMesh` Tetrahedral mesher
- `pMesh` Polyhedral mesher
- `libmeshLibrary.so` Core mesh library

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

Previous ports of cfMesh to modern OpenFOAM compile but crash or hang at runtime.
This port fixes the following runtime bugs discovered through debugging on a
NASA Rotor 37 turbomachinery geometry:

| Bug | Fix |
|-----|-----|
| `face::centre()` returns garbage for cfMesh face ordering → phantom 1946 inverted faces | Replace with simple vertex average in `calculateFaceCentres()`, `meshSurfaceEngineModifier`, `edgeExtractor` |
| Inversion check on unmapped octree mesh → infinite untangler loop | Map vertices onto surface BEFORE `preMapVertices()` in `mapMeshToSurface()` |
| `preMapVertices()` internal `untangleSurface` hangs on complex geometry | Use `preMapVertices(0)` to skip internal untangle |
| `findLeafContainingVertex` throws `char const*` → uncaught crash | Replace throw with `return -1` |
| `newPositionSurfaceOptimizer` produces out-of-domain vertex positions | Add displacement bound and vecY near-zero guard |
| `pointFieldPMG::setSize` uninitialized over-allocation → garbage points | Zero-initialize over-allocated region |
| Octree points uninitialized for unreferenced nodes | Zero-initialize `octreePointsPtr_` |
| `checkFacePatchesGeometry` unbounded loop → infinite loop | Add iteration limit (`nGeomIter < 5`) |
| `findNearestSurfacePoint` 100-iteration hang for out-of-STL vertices | Reduce limit to 10 |
| `optimizePoint` FPE from zero division | Guard zero divisions in `surfaceOptimizer` |
| `untangleSurface` too many global iterations for complex geometry | Limit to 2 global iterations |
| Point sanitizer | Zero out near-zero denormalized values before writing mesh |

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
or the OpenFOAM Foundation
