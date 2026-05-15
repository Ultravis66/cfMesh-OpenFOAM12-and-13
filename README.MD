python3 << 'PYEOF'
content = """# cfMesh OpenFOAM 12 & 13 Port

> **WORK IN PROGRESS** — active development, more patches coming.

A working port of cfMesh to OpenFOAM Foundation versions 12 and 13, with significant boundary layer robustness improvements beyond the original cfMesh.

Original cfMesh developed by Franjo Juretic / Creative Fields, licensed under GPL v3.

Tested on **NASA Rotor 37** two-passage turbomachinery sector geometry with boundary layers.

---

## What works
- `cartesianMesh` — Cartesian hex-dominant mesher tested
- `tetMesh` — Tetrahedral mesher (experimental)
- `pMesh` — Polyhedral mesher (experimental)
- `libmeshLibrary.so` — Core mesh library

---

## Branches

| Branch | OpenFOAM Version | Status |
|--------|-----------------|--------|
| `master` | OpenFOAM 12 | Working |
| `openfoam13` | OpenFOAM 13 | Working |

---

## Build instructions

OpenFOAM 12:

    source /opt/openfoam12/etc/bashrc
    git clone -b master https://github.com/Ultravis66/cfMesh-OpenFOAM12-and-13.git cfMesh
    cd cfMesh/meshLibrary && wmake libso
    cd ../executables/cartesianMesh && wmake
    cd ../executables/tetMesh && wmake
    cd ../executables/pMesh && wmake

OpenFOAM 13:

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
| Pstream::blocking/scheduled/nonBlocking | Pstream::commsTypes::blocking/scheduled/nonBlocking |
| eigenVector(T, scalar) | eigenVectors(T) with direct component access |
| circumCentre()/circumRadius() | circumSphere().first()/.second() |
| isFile() | std::filesystem::exists() - requires filesystem include and c++17 |
| rm/isDir/rmDir | requires OSspecific.H include in affected files |
| edgeMesh library | include paths updated to use meshTools |
| unallocLabelList | labelList |

---

## Runtime fixes

Previous ports compile but crash or hang at runtime. This port fixes:

| Bug | Fix |
|-----|-----|
| face::centre() garbage for cfMesh face ordering - phantom inversions | Replace with vertex average everywhere |
| face::normal() garbage - wrong BL extrusion direction | Replace with triangle fan normal everywhere |
| Inversion check on unmapped octree mesh - infinite loop | Map vertices onto surface BEFORE preMapVertices() |
| preMapVertices() internal untangleSurface hangs | Use preMapVertices(0) |
| findLeafContainingVertex throws char const* | Replace throw with return -1 |
| newPositionSurfaceOptimizer garbage vertex positions | Displacement bound + vecY near-zero guard |
| pointFieldPMG::setSize uninitialized over-allocation | Zero-initialize over-allocated region |
| Octree points uninitialized | Zero-initialize octreePointsPtr_ |
| checkFacePatchesGeometry unbounded loop | Iteration limit (nGeomIter < 5) |
| findNearestSurfacePoint 100-iteration hang | Reduce to 10 iterations |
| optimizePoint FPE from zero division | Guard zero divisions in surfaceOptimizer |
| untangleSurface too many global iterations | Limit to 2 global iterations |
| Garbage denormalized points in output | Sanitize before write |

---

## Boundary layer improvements (beyond original cfMesh)

### Feature-aware vertex normals
Replaced simple face-normal averaging with angle-weighted normals that respect sharp features. At blade trailing edges, leading edges, and blade-endwall junctions, the new code detects the feature angle and uses only the dominant-patch normals instead of averaging across the sharp edge. Result: 96% reduction in bad BL cells on NASA Rotor 37 (229,000 to 9,646).

### Gap distance detection
For each wall vertex being extruded, ray-casts along the surface normal to detect opposing surfaces in the 2-ring face neighborhood. Limits BL thickness to 50% of the detected gap distance. Prevents cell overlap in narrow passages such as blade trailing edges and blade-endwall junctions.

### Concave edge layer termination
Implemented the previously stub-only terminateLayersAtConcaveEdges flag using meshSurfaceCheckEdgeTypes. Identifies concave edge vertices and skips BL extrusion there. Enabled by default in cartesianMeshGenerator.

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
"""

with open('/root/cfmesh-repo/README.MD', 'w') as f:
    f.write(content)
print("README written")
PYEOF

cd ~/cfmesh-repo
git add README.MD
git commit -m "Update README with BL improvements and full fix table"
git push origin master
git checkout openfoam13
git merge master --no-edit
git push origin openfoam13
git checkout master
