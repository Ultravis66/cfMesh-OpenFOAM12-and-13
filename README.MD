# cfMesh for OpenFOAM Foundation 13

An actively developed fork of cfMesh for OpenFOAM Foundation v13.

**NOTE:** This release candidate replaces the earlier rc1 development state. The previous build contained boundary-layer quality-recovery behavior that could generate meshes which completed meshing checks but were not sufficiently robust for practical CFD execution. The current rc1 corrects those issues and has now been validated through actual OpenFOAM CFD execution.

This fork focuses on:

* OpenFOAM 13 compatibility
* Robust Cartesian and polyhedral mesh generation
* Boundary-layer generation near difficult topology
* Deterministic and race-hardened OpenMP behavior
* Transactional mesh-quality repair
* Topology and addressing validation
* Reliable, lossless mesh serialization
* Surface conformance for STL-driven CFD workflows
> **Current status:** `v1.0.0-rc1` release candidate  
> **Current branch:** `release/v1.0.0-rc1`  
> This is a research and development mesher, not a production-certified industrial tool.

The original cfMesh was developed by Franjo Juretić and Creative Fields. This fork preserves the original GPL licensing and authorship while adding extensive OpenFOAM 13 compatibility and robustness work.

## Boundary-layer architecture selection

This fork contains three selectable boundary-layer generation paths. The architecture is selected in `system/meshDict`, allowing users to retain the original cfMesh behavior, use the earlier enhanced implementation, or enable the new constraint-planner pipeline.

| Architecture | Description | Recommended use |
|---|---|---|
| `classicCfMesh` | Preserved classic cfMesh boundary-layer implementation, adapted for OpenFOAM Foundation 13. | Compatibility testing, regression comparisons, and geometries that worked reliably with the original cfMesh approach. |
| `legacyEnhanced` | Heavily modified continuation of the original boundary-layer implementation, including robustness, projection, OpenMP, and transition-handling improvements. | Existing cases developed with earlier versions of this fork. |
| `constraintPlanner` | New CFMitch constraint-planning architecture. It computes explicit per-face layer plans, enforces neighboring-layer compatibility, tracks repair provenance, and supports transactional refinement and rollback. | Current development path and recommended architecture for difficult boundary-layer topology. |

> **Validation status:** The Rotor37 and propeller results documented below were generated using `constraintPlanner`. In a four-layer `singleOrifice` regression test, `legacyEnhanced` produced a stock-`checkMesh`-clean mesh. The preserved `classicCfMesh` path completed successfully with valid topology and positive cell volumes, but stock `checkMesh` reported one incorrectly oriented face pyramid.

### CFMitch constraint-planner architecture

To use the current CFMitch architecture:

```cpp
boundaryLayers
{
    boundaryLayerArchitecture constraintPlanner;

    constraintPlanner
    {
        maxLayerStep           2;
        neutralProjectionFloor 0.0;

        heightPlanner
        {
            enabled          false;
            maxEdgeRatio     1.50;
            maxMoveFraction  0.05;
            iterations       1;
        }
    }

    // Patch roles and layer definitions follow...
}
```

### Classic cfMesh architecture

To use the preserved classic implementation:

```cpp
boundaryLayers
{
    boundaryLayerArchitecture classicCfMesh;

    // Classic boundary-layer settings follow...
}
```

### Legacy-enhanced architecture

To use the earlier enhanced implementation:

```cpp
boundaryLayers
{
    boundaryLayerArchitecture legacyEnhanced;

    // Enhanced boundary-layer settings follow...
}
```

> **Compatibility note:** These architectures are not expected to produce identical meshes. They represent distinct boundary-layer planning and topology-generation paths. Always run the stock OpenFOAM `checkMesh` utility and inspect boundary-layer junctions before using a generated mesh for CFD.

## CFMitch infrastructure

CFMitch identifies the new robustness, repair, and diagnostic infrastructure developed specifically for this fork. It is not a separate executable. It operates through `cartesianMesh` when the corresponding architecture and controls are selected in `meshDict`.

The CFMitch pipeline includes:

- Explicit boundary-layer constraint planning
- Neighboring-face layer-count compatibility
- Topology and face-incidence validation
- OpenFOAM-equivalent candidate-quality gates for signed volume, face pyramids, non-orthogonality, and skewness
- Transactional repair with acceptance or rollback
- Boundary-layer provenance tracking
- Provenance-guided adaptive local layer retreat, including post-Q3 hard-quality recovery
- Addressing validation before and after renumbering
- Safe handling of cells with more than 255 faces
- Lossless mesh-coordinate serialization

The original cfMesh implementation remains available both as a compatibility reference and as a practical fallback for cases where its simpler behavior is preferable.

## Development scope

This is a substantial engineering fork rather than a minimal OpenFOAM compatibility patch. The repository contains approximately 180,000 lines of source code, with heavily modified or added codebase during the OpenFOAM 13 port and robustness-development effort.

The modifications span boundary-layer planning and refinement, mesh topology generation, surface projection, geometry optimization, OpenMP race hardening, mesh addressing, diagnostic infrastructure, transactional repair, and mesh serialization.

These figures describe the scale of the development effort; mesh quality is evaluated independently using stock OpenFOAM utilities and reproducible validation cases.

## What this fork can generate

### Large 16-layer propeller validation mesh

<p align="center">
  <img src="Propeller_1.png" alt="Large cfMesh propeller boundary-layer mesh overview" width="32%">
  <img src="Propeller_2.png" alt="Detailed propeller boundary-layer mesh" width="32%">
  <img src="Propeller_3.png" alt="Propeller surface and boundary-layer mesh detail" width="32%">
</p>

<p align="center">
  <em>Large propeller stress-test mesh: 13.14 million cells, 16 boundary layers, zero negative-volume cells, zero incorrect face pyramids, and a clean stock OpenFOAM 13 checkMesh result.</em>
</p>

### NASA Rotor 37 turbomachinery stress test

<p align="center">
  <img src="NASA_Rotor_37_1.png" alt="NASA Rotor 37 cfMesh boundary-layer mesh overview" width="48%">
  <img src="NASA_Rotor_37_2.png" alt="NASA Rotor 37 blade and endwall mesh" width="48%">
</p>

<p align="center">
  <img src="NASA_Rotor_37_3.png" alt="NASA Rotor 37 boundary-layer mesh detail" width="48%">
  <img src="NASA_Rotor_37_4.png" alt="NASA Rotor 37 periodic-sector mesh detail" width="48%">
</p>

<p align="center">
  <em>NASA Rotor 37 two-passage sector used to expose difficult blade/endwall junctions, periodic boundaries, acute trailing edges, and boundary-layer transition topology. The current result has zero negative-volume cells, correctly oriented face pyramids, and no non-orthogonality errors; localized skewness remains the primary mesh-quality target.</em>

</p>

---

## Current release candidate

The `v1.0.0-rc1` development candidate introduces a new boundary-layer planning and repair architecture together with several important fixes to core mesh handling.

### Large propeller validation case

The current release candidate produced a clean OpenFOAM mesh for a large propeller geometry using a 16-layer boundary-layer configuration with a nominal first-layer height of approximately 40 μm.

| Metric                          |      Result |
| ------------------------------- | ----------: |
| Points                          |  13,408,984 |
| Faces                           |  39,692,917 |
| Internal faces                  |  39,498,813 |
| Cells                           |  13,141,866 |
| Hexahedra                       |  13,018,764 |
| Negative-volume cells           |           0 |
| Incorrect face pyramids         |           0 |
| Maximum non-orthogonality       |      82.35° |
| Severely non-orthogonal faces   |          13 |
| Maximum skewness                |        3.72 |
| Propeller BL effective coverage |      96.66% |
| Stock OpenFOAM `checkMesh`      | **Mesh OK** |

The written mesh was independently read and validated by the stock OpenFOAM Foundation v13 `checkMesh` utility.

### NASA Rotor 37 stress test

NASA Rotor 37 remains the primary difficult-topology development case. It combines:

* Blade/hub/shroud junctions
* Periodic sector boundaries
* Acute trailing-edge geometry
* BL/periodic and BL/no-BL interfaces
* Strong refinement transitions
* Multiple interacting boundary-layer patches

The current high-resolution Rotor37 run produced a valid, connected OpenFOAM mesh with zero negative-volume cells, correctly oriented face pyramids, and no non-orthogonality errors. Stock OpenFOAM `checkMesh` reports one remaining failed quality check for localized skewness.

| Metric                             |                   Current result |
| ---------------------------------- | -------------------------------: |
| Points                             |                        5,552,053 |
| Faces                              |                       15,557,523 |
| Internal faces                     |                       14,957,239 |
| Cells                              |                        5,002,615 |
| Number of regions                  |                           1 (OK) |
| Negative-volume cells              |                                0 |
| Interface-audit failures           |                                0 |
| Invalid owner/neighbour addressing |                                0 |
| Incorrect face pyramids            |                                0 |
| Non-orthogonality errors           |                                0 |
| Maximum non-orthogonality          |                           89.46° |
| Severely non-orthogonal faces      |                            1,306 |
| Maximum skewness                   |                            16.85 |
| Highly skew faces                  |                              636 |
| Stock OpenFOAM `checkMesh`         | 1 failed quality check: skewness |

Rotor37 is therefore a valid mesh-generation result and is runnable in OpenFOAM. It is not yet a CFD-performance validation claim: 1,306 faces remain severely non-orthogonal above 70°—although none exceed 90°—and 636 faces exceed the stock skewness threshold. Pressure ratio, efficiency, mass flow, and convergence validation remain pending.


---

## Major improvements in `v1.0.0-rc1`

### Constraint-planned boundary layers

The new boundary-layer constraint planner computes a compatible refinement plan before topology generation.

It supports:

* Patch-specific layer counts and heights
* BL, termination, and neutral patch roles
* Layer-count compatibility across adjacent faces
* Maximum layer-step constraints
* Height-first repair planning
* Local layer compression
* Local hard retreat where valid layered topology cannot be generated
* Preservation of classic cfMesh behavior as a selectable path

The planner is geometry-driven and does not depend on Rotor37-specific patch names.

### Transactional boundary-layer repair

Boundary-layer repair is evaluated through mesh snapshots rather than applying irreversible mutations.

The repair pipeline can:

1. Preserve the valid pre-refinement mesh
2. Generate and evaluate an initial BL result
3. Build repair plans from failed-face provenance
4. Restore the original topology
5. Retry with locally modified layer constraints
6. Accept only a safer result
7. Roll back rejected attempts

Quality comparisons include:

* Negative-volume cells
* Incorrect face pyramids
* Signed-volume magnitude
* Minimum cell volume
* Non-orthogonality
* Skewness
* Open-cell detection
* Unused-point checks

### Adaptive local BL retreat

When compression cannot produce a valid layered configuration, the mesher can locally suppress boundary-layer refinement on selected faces.

The current implementation includes:

* Failed-face provenance
* Cumulative retreat seeds
* Co-parent closure of related boundary faces
* Local maximum-layer overrides
* Compression context retained during retreat
* Acceptance against the best valid baseline
* Complete rollback of rejected candidates

On the propeller validation case, co-parent closure expanded 16 retreat seeds to 18 faces and produced a clean accepted mesh.
The final Rotor37 repair eliminated all incorrectly oriented face pyramids and all non-orthogonality errors while preserving zero negative-volume cells. It also reduced maximum skewness from 41.10 to 16.85, leaving skewness as the only remaining stock `checkMesh` failure.
### Generated-child topology audits

New diagnostic gates validate refined cell shells before and after face relabeling.

Checks include:

* Every child edge used exactly twice
* No malformed faces
* No repeated consecutive face vertices
* Valid cell-to-face references
* Closed generated child shells
* Closed cells before final relabeling
* Closed cells after final relabeling
* Correct final face incidence
* No unused, singly misclassified, or multiply owned internal faces

These checks cover refinement types 0–3 and appended child cells.

### Final face-incidence validation

Before committing the generated topology, the mesher verifies that:

* Every internal face is referenced by exactly two cells
* Every boundary face is referenced by exactly one cell
* No face is unused
* No face is referenced by three or more cells
* Every cell face reference is in range
* The derived internal-face boundary agrees with patch metadata

### Owner/neighbour addressing coherency

Snapshot and retry operations can replace points, faces, and cells wholesale. Earlier cache invalidation paths did not always discard cached owner/neighbour arrays associated with the previous topology.

The current implementation:

* Clears owner/neighbour addressing after topology restoration
* Rebuilds addressing before renumbering
* Audits addressing before renumbering
* Audits addressing after renumbering
* Audits addressing after boundary replacement
* Aborts instead of writing a corrupted mesh if final addressing is invalid

### `replaceBoundary` wide-counter fix

The original boundary replacement implementation used:

```cpp
List<direction> nFacesInCell
```

A `direction` is too narrow for cells containing more than 255 faces. Counter overflow silently truncated cell face lists during boundary replacement.

This became visible in the propeller stress test, which contains two 291-face polyhedral cells. Each affected counter wrapped at 256, corrupting the written mesh even though the pre-replacement topology was valid.

The counter now uses:

```cpp
List<label> nFacesInCell
```

This fix restored all cell face references and eliminated millions of apparently illegal cells reported by stock `checkMesh`.

### Lossless mesh serialization

Low-precision ASCII output previously rounded sensitive coordinates enough to introduce face-pyramid errors that were absent from the in-memory mesh.

The mesh writer now enforces at least:

```cpp
std::numeric_limits<scalar>::max_digits10
```

For a double-precision OpenFOAM build, this results in 17 significant digits. The requested output format—ASCII or binary—is preserved.

This provides numerical round-trip consistency between:

* The final in-memory mesh
* The written `constant/polyMesh`
* The mesh reconstructed by stock OpenFOAM utilities

### OpenFOAM primitive-mesh parity diagnostics

The generator includes a parity audit modeled on stock OpenFOAM calculations for:

* Cell centres
* Signed cell volumes
* Face pyramids
* Internal-face non-orthogonality
* Face skewness

This distinguishes actual generator defects from differences caused by serialization or by cfMesh/OpenFOAM geometric conventions.

---

## Additional robustness work

### Surface-constrained optimization

Boundary points can be constrained to their appropriate STL regions during optimization, reducing surface drift near patch interfaces.

The implementation includes:

* Patch-aware surface projection
* Feature-edge and corner classification
* BL/no-BL protection
* BL/neutral protection
* Limited-displacement reprojection
* Transactional acceptance and rollback
* Surface-drift diagnostics

### Parallel hardening

Numerous race-sensitive shared writes have been serialized or restructured, including:

* Boundary-node construction
* Laplacian smoother apply loops
* Surface-optimizer apply loops
* Boundary-layer vertex maps
* Front-point marking
* Tet and triangle update-type marking
* Boundary addressing flags
* Topology-changing repair loops
* Unguarded parallel map access

`OMP_NUM_THREADS=1` remains the preferred deterministic debugging and regression mode. Parallel repeatability has improved substantially, but high-thread-count validation remains ongoing.

### Geometry and numerical guards

The fork includes guards and recovery behavior for:

* Degenerate faces
* Zero-length vectors
* Zero-volume tetrahedra
* Invalid boundary-point addressing
* Missing map entries
* Octree search failures
* Excessive recursion
* Uninitialized over-allocation
* Invalid split-edge metadata
* Stagnating mesh untangling
* Unsafe optimizer results

---

## Supported executables

| Component           | Status                                    |
| ------------------- | ----------------------------------------- |
| `cartesianMesh`     | Primary supported mesher; actively tested |
| `tetMesh`           | Builds successfully; experimental         |
| `pMesh`             | Builds successfully; experimental         |
| `libmeshLibrary.so` | Core cfMesh library                       |

---

## Branch status

| Branch               | OpenFOAM version | Status                    |
| -------------------- | ---------------- | ------------------------- |
| `release/v1.0.0-rc1` | Foundation v13   | Current release candidate |
| `master`             | Foundation v13   | Previous development line |
| `openfoam12-archive` | OpenFOAM 12      | Archived                  |

---

## Quick start

### Requirements

* OpenFOAM Foundation v13
* Ubuntu 22.04 or compatible Linux environment
* GCC toolchain compatible with the OpenFOAM installation
* C++17 filesystem support
* WSL2 is supported for development and testing

### Clone and build the release candidate

```bash
source /opt/openfoam13/etc/bashrc

git clone https://github.com/Ultravis66/cfMesh-OpenFOAM13.git cfMesh
cd cfMesh

git switch release/v1.0.0-rc1

./Allwclean
./Allwmake
```

Verify the build:

```bash
ls "$FOAM_USER_APPBIN" |
grep -E 'cartesianMesh|tetMesh|pMesh'

ls "$FOAM_USER_LIBBIN/libmeshLibrary.so"
```

### Run the unit-cube tutorial

```bash
cd tutorials/cartesianMesh/unitCube

export OMP_NUM_THREADS=1

cartesianMesh 2>&1 |
tee log.cartesianMesh

checkMesh -constant 2>&1 |
tee log.checkMesh
```

---

## Recommended validation workflow

Always inspect the text produced by `checkMesh`. On some OpenFOAM versions, `checkMesh` can return process exit code zero even when geometric checks fail.

```bash
checkMesh -constant 2>&1 |
tee log.checkMesh

grep -E \
'Mesh OK|Failed|negative|non-orthogonality|face pyramids|skewness' \
log.checkMesh
```

For development comparisons, use:

```bash
export OMP_NUM_THREADS=1
```

before running `cartesianMesh`.

A generated mesh should not be accepted solely because the mesher exits successfully. At minimum, verify:

* Boundary definition
* Cell-to-face addressing
* Point usage
* Upper-triangular ordering
* Region count
* Cell volumes
* Face pyramids
* Non-orthogonality
* Skewness
* Surface conformity in ParaView

---

## Boundary-layer controls

Boundary-layer behavior is configured in `system/meshDict`.

The current architecture supports concepts such as:

```text
boundaryLayers
{
    boundaryLayerArchitecture constraintPlanner;

    terminationPatches
    (
        inlet
        outlet
    );

    patchBoundaryLayers
    {
        // Patch-specific layer settings
    }
}
```

Important controls include:

| Control                            | Purpose                                                |
| ---------------------------------- | ------------------------------------------------------ |
| `boundaryLayerArchitecture`        | Select constraint-planned or classic behavior          |
| `patchBoundaryLayers`              | Patch-specific BL layer counts and heights             |
| `terminationPatches`               | Identify inlet/outlet-style BL termination patches     |
| `preRefBLAutoRepair`               | Enable pre-refinement snapshotting and transactional BL repair |
| `constrainOptimizerBoundaryMotion` | Constrain boundary motion to surface geometry          |
| `postBLSnap`                       | Enable the experimental post-optimization surface snap |
| `geometryPreprocessing`            | Configure optional near-coincident STL vertex handling |

Not every control is appropriate for every geometry. Validate each configuration using stock OpenFOAM utilities and visual inspection.

---

## Current limitations

* Rotor37 now passes the topology, cell-volume, face-pyramid, and non-orthogonality-error checks, remains one connected region, and contains zero negative-volume cells.
* The current high-resolution Rotor37 mesh still contains 1,306 severely non-orthogonal faces above 70° and 636 highly skew faces. Skewness is its only remaining failed stock `checkMesh` criterion.
* CFD performance validation—including pressure ratio, efficiency, mass flow, and convergence behavior—is still pending.
* Very small first-layer heights and large layer counts require case-specific validation.
* Surface-constrained optimization can significantly increase runtime.
* OMP=1 remains the reference mode for deterministic regression testing.
* High-thread-count OpenMP repeatability remains under investigation.
* `tetMesh` and `pMesh` are not validated to the same level as `cartesianMesh`.
* This project is not production-certified.
* Users must validate every mesh independently.

---

## NASA Rotor37 stress-test workflow

A Rotor37 setup may be prepared using:

```bash
mkdir -p ~/rotor37_mesh/constant/triSurface
mkdir -p ~/rotor37_mesh/system

cp /path/to/Rotor37_cfmesh.stl \
    ~/rotor37_mesh/constant/triSurface/

cp tutorials/rotor37/system/meshDict \
    ~/rotor37_mesh/system/
```

Run deterministically:

```bash
cd ~/rotor37_mesh

export OMP_NUM_THREADS=1

cartesianMesh 2>&1 |
tee log.Rotor37

checkMesh -constant 2>&1 |
tee log.checkMesh.Rotor37
```

Rotor37 results are sensitive to geometry preparation, refinement settings, first-layer height, layer count, and active repair controls. Historical metric ranges should not be treated as universal acceptance thresholds.

---

## Selected OpenFOAM 13 API updates

| Original usage                      | OpenFOAM 13 update                      |
| ----------------------------------- | --------------------------------------- |
| `Pstream::blocking`                 | `Pstream::commsTypes::blocking`         |
| `Pstream::scheduled`                | `Pstream::commsTypes::scheduled`        |
| `Pstream::nonBlocking`              | `Pstream::commsTypes::nonBlocking`      |
| `eigenVector(T, scalar)`            | `eigenVectors(T)` with component access |
| `circumCentre()` / `circumRadius()` | `circumSphere().first()` / `.second()`  |
| `isFile()`                          | `std::filesystem::exists()`             |
| `unallocLabelList`                  | `labelList`                             |
| Legacy edge-mesh includes           | Updated `meshTools` include paths       |

---

## Release history

| Version      | Summary                                                                                                                                                                                                                            |
| ------------ | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `v1.0.0-rc1` | Constraint-planned BL refinement, transactional repair, OpenFOAM-equivalent quality gating, post-Q3 adaptive retreat, generated-child closure audits, face-incidence and addressing validation, wide cell-face counters, and lossless mesh serialization.
                                            |
| `v0.10`      | Surface-constrained optimizer and limited-displacement boundary reprojection                                                                                                                                                       |
| `v0.9`       | Expanded parallel hardening and optimizer/untangler fixes                                                                                                                                                                          |
| `v0.8.x`     | BL/BL and BL/neutral transition infrastructure, geometry preprocessing, and broad codebase safety audit                                                                                                                            |
| `v0.7.x`     | BL/no-BL and BL/neutral edge constraint system                                                                                                                                                                                     |
| `v0.6`       | Initial deterministic optimizer hardening                                                                                                                                                                                          |
| `v0.1–v0.5`  | OpenFOAM 13 port, surface mapping fixes, quality evaluation, and targeted repair                                                                                                                                                   |

`v1.0.0` will remain reserved for a final release after the release candidate has completed broader regression and documentation review.

---

## Original cfMesh

* Original project: https://sourceforge.net/projects/cfmesh/
* Original authors include Franjo Juretić, Philippose Rajan, and Ivor Clifford
* Original copyright and authorship remain with their respective owners

---

## License and attribution

This project is licensed under the GNU General Public License v3.0, consistent with the original cfMesh license.

This fork includes OpenFOAM Foundation v13 compatibility updates and meshing-robustness modifications by Mitchell Stolk.

Modifications in this fork:

```text
Copyright (C) 2026 Mitchell Stolk
```

Suggested attribution:

```text
Mitchell Stolk, cfMesh-OpenFOAM13:
OpenFOAM Foundation 13 port and robustness-focused cfMesh fork, 2026.
```

If this fork is useful in academic research, engineering evaluation, commercial assessment, or public technical work, attribution to this repository is appreciated.

---

## Disclaimer

This is an unofficial community fork. It is not affiliated with Creative Fields or the OpenFOAM Foundation.

The software is provided without warranty. Users are responsible for validating mesh topology, geometry, numerical quality, and solver suitability for every application.
