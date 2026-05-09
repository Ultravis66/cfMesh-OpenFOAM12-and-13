# cfMesh — OpenFOAM 12 Port

A working port of cfMesh to OpenFOAM 12 (openfoam.org).
Original cfMesh by Franjo Juretic, Creative Fields. Licensed under GPL v3.

---

## What works

- `cartesianMesh` — Cartesian hex-dominant mesher
- `tetMesh` — Tetrahedral mesher
- `pMesh` — Polyhedral mesher
- `libmeshLibrary.so` — Core mesh library

---

## Build instructions

```bash
source /opt/openfoam12/etc/bashrc
git clone https://github.com/Ultravis66/cfMesh-OpenFOAM12.git cfMesh
cd cfMesh/meshLibrary && wmake libso
cd ../executables/cartesianMesh && wmake
cd ../executables/tetMesh && wmake
cd ../executables/pMesh && wmake
```

---

## API changes from original cfMesh

| Original | Updated |
|----------|---------|
| `Pstream::blocking/scheduled/nonBlocking` | `Pstream::commsTypes::blocking/scheduled/nonBlocking` |
| `eigenVector(T, scalar)` | `eigenVectors(T)` with direct component access |
| `circumCentre()/circumRadius()` | `circumSphere().first()/.second()` |
| `isFile()` | `std::filesystem::exists()` — requires `#include <filesystem>` and `-std=c++17` |
| `rm/isDir/rmDir` | requires `#include "OSspecific.H"` in affected files |
| `edgeMesh` library | include paths updated to use `meshTools` |
| `unallocLabelList` | `labelList` (fixed in FMSToVTK and mergeSurfacePatches) |

---

## Known issues
- Tested on Ubuntu 22.04 WSL2 with OpenFOAM 12 only

---

## Original cfMesh

- Source: https://sourceforge.net/projects/cfmesh/
- Authors: Franjo Juretic, Philippose Rajan, Ivor Clifford
