# Linear Network Flow (LNF) Model Implementation

This document describes the implementation of the Linear Network Flow (LNF) model in MODFLOW 6. The LNF model is designed for simulating flow in linear network systems (such as pipe networks or channel networks) using a simplified approach with static conductance values.

## Table of Contents

1. [Overview](#overview)
2. [Architecture](#architecture)
3. [File Structure](#file-structure)
4. [Core Components](#core-components)
5. [Package Descriptions](#package-descriptions)
6. [IDM Integration](#idm-integration)
7. [Build System Integration](#build-system-integration)
8. [Testing](#testing)
9. [Key Implementation Details](#key-implementation-details)

---

## Overview

The LNF model extends MODFLOW 6's `NumericalModelType` base class to provide a lightweight model for linear network flow simulation. Unlike the GWF (Groundwater Flow) model which computes conductance from hydraulic properties, LNF uses pre-specified static conductance values for network connections.

### Key Features

- **Static Conductance**: Flow between nodes is computed using user-specified conductance values
- **Unstructured Grid**: Uses DISU (unstructured discretization) to define arbitrary network topology
- **Standard Boundary Conditions**: Supports CHD (Constant Head) boundaries
- **Full Output Control**: HEAD and BUDGET output to binary files

### Flow Equation

The flow between connected nodes follows Darcy's law in discrete form:

```
Q_ij = C_ij * (h_i - h_j)
```

Where:
- `Q_ij` = flow from node i to node j
- `C_ij` = conductance of connection between nodes i and j
- `h_i`, `h_j` = heads at nodes i and j

---

## Architecture

The LNF model follows the standard MODFLOW 6 model architecture:

```
NumericalModelType (base class)
    └── LnfModelType
            ├── DisBaseType (DISU - discretization)
            ├── LnfCndType (CND - conductance package)
            ├── LnfOcType (OC - output control)
            └── BndType list (boundary packages like CHD)
```

### Model Lifecycle

The LNF model implements all required lifecycle procedures:

| Procedure | Purpose |
|-----------|---------|
| `lnf_df` | Define - create packages from input |
| `lnf_ac` | Add connections to sparse matrix structure |
| `lnf_mc` | Map connections after matrix structure is set |
| `lnf_ar` | Allocate and read - setup arrays and read data |
| `lnf_rp` | Read and prepare - process stress period data |
| `lnf_ad` | Advance - advance to next time step |
| `lnf_fc` | Formulate coefficients - fill matrix and RHS |
| `lnf_cq` | Calculate flows between cells |
| `lnf_bd` | Budget - accumulate flow terms |
| `lnf_ot` | Output - write results to files |
| `lnf_da` | Deallocate - cleanup memory |

---

## File Structure

### Source Files

```
src/Model/LinearNetworkFlow/
├── lnf.f90          # Main LNF model module
├── lnf-cnd.f90      # CND (Conductance) package
└── lnf-oc.f90       # OC (Output Control) package
```

### Definition Files (DFN)

```
doc/mf6io/mf6ivar/dfn/
├── lnf-nam.dfn      # NAM (Name file) definitions
├── lnf-disu.dfn     # DISU (Discretization) definitions
├── lnf-cnd.dfn      # CND (Conductance) definitions
├── lnf-chd.dfn      # CHD (Constant Head) definitions
└── lnf-oc.dfn       # OC (Output Control) definitions
```

### Generated IDM Files

```
src/Idm/
├── lnf-namidm.f90   # Generated from lnf-nam.dfn
├── lnf-disuidm.f90  # Generated from lnf-disu.dfn
├── lnf-cndidm.f90   # Generated from lnf-cnd.dfn
├── lnf-chdidm.f90   # Generated from lnf-chd.dfn
└── lnf-ocidm.f90    # Generated from lnf-oc.dfn
```

### Test Files

```
autotest/
└── test_lnf_01.py   # 3-node linear network test
```

---

## Core Components

### LnfModelType (lnf.f90)

The main model type that orchestrates all LNF functionality.

**Key Members:**
```fortran
type, extends(NumericalModelType) :: LnfModelType
  type(LnfCndType), pointer :: cnd => null()     ! Conductance package
  type(LnfOcType), pointer :: oc => null()       ! Output control
  integer(I4B), pointer :: incnd => null()       ! CND unit number
  integer(I4B), pointer :: inoc => null()        ! OC unit number
  real(DP), dimension(:), pointer :: flowja      ! Intercell flows
end type
```

**Key Procedures:**

1. **lnf_fc** - Formulates the coefficient matrix:
   - Calls `cnd%cnd_fc()` to fill AMAT with conductance terms
   - Calls boundary package `bnd_fc()` for each boundary

2. **lnf_cq** - Calculates intercell flows:
   - Calls `cnd%cnd_cq()` to compute `flowja` array
   - Calls boundary package `bnd_cq()` for boundary flows

3. **lnf_ot** - Handles output:
   - Saves `flowja` to CBC file via `cnd%cnd_save_model_flows()`
   - Saves boundary flows via `bnd_ot_model_flows()`
   - Saves HEAD via `oc%oc_ot()`
   - Prints budget summary

---

## Package Descriptions

### DISU Package (Discretization - Unstructured)

The LNF model uses the standard DISU package from MODFLOW 6 to define network topology.

**Key Input:**
- `NODES` - Number of nodes in network
- `NJA` - Number of connections (including diagonal)
- `IAC` - Number of connections per node
- `JA` - Connection array (CSR format)
- `TOP`, `BOT` - Node elevations
- `AREA` - Node areas
- `CL12`, `HWVA` - Connection geometry (optional for LNF)

### CND Package (Conductance)

**Purpose:** Defines static conductance values for all network connections.

**Key Input:**
- `NJA` - Number of entries in conductance array
- `CONDUCTANCE` - Array of conductance values (NJA entries, CSR format)

**Key Procedures:**

```fortran
subroutine cnd_fc(this, amat, idxglo, rhs, x)
  ! Fill coefficient matrix with conductance terms
  ! For each connection: amat(idiag) += cond, amat(ioff) -= cond
end subroutine

subroutine cnd_cq(this, x, flowja)
  ! Calculate intercell flows: flowja(ipos) = cond * (hn - hm)
end subroutine

subroutine cnd_save_model_flows(this, flowja, icbcfl, icbcun)
  ! Write FLOW-JA-FACE to CBC file
end subroutine
```

### CHD Package (Constant Head)

Uses the standard MODFLOW 6 CHD boundary package framework. Specifies fixed head values at selected nodes.

**Input Format:**
```
BEGIN PERIOD 1
  node_number  head_value
END PERIOD
```

### OC Package (Output Control)

Controls saving of HEAD and BUDGET to binary output files.

**Key Input:**
```
BEGIN OPTIONS
  HEAD FILEOUT model.hds
  BUDGET FILEOUT model.cbc
END OPTIONS

BEGIN PERIOD 1
  SAVE HEAD ALL
  SAVE BUDGET ALL
END PERIOD
```

---

## IDM Integration

The Input Data Model (IDM) system in MODFLOW 6 handles reading input files and storing data in memory. The LNF model integrates with IDM through:

### 1. Definition Files (DFN)

Each package has a `.dfn` file that defines input blocks, parameters, and data types. Example from `lnf-cnd.dfn`:

```
# --------------------- lnf cnd options ---------------------
block options
name print_input
type keyword
...

# --------------------- lnf cnd griddata ---------------------
block griddata
name conductance
type double precision
shape (nja)
reader readarray
```

### 2. Memory Path Convention

LNF uses context names ≤16 characters (LENCONTEXTNAME limit):
- Model context: `LNF01` (model name)
- Package mempath: `__INPUT__/LNF01/CND`

### 3. IDM Registration

The model type is registered in:
- `src/Utilities/Idm/DefinitionSelect.f90` - Package definition selection
- `src/Utilities/Idm/IdmLoad.f90` - Input loading
- `src/Utilities/Idm/SourceLoad.F90` - Source loading context

---

## Build System Integration

### Meson Build

LNF source files are automatically included via glob patterns in `src/meson.build`:

```meson
modflow6_sources = files(
  ...
)

# LinearNetworkFlow files included via directory glob
```

### Visual Studio (mf6core.vfproj)

For Windows builds, add to the project file:

```xml
<Filter Name="LinearNetworkFlow">
  <File RelativePath="..\..\src\Model\LinearNetworkFlow\lnf.f90"/>
  <File RelativePath="..\..\src\Model\LinearNetworkFlow\lnf-cnd.f90"/>
  <File RelativePath="..\..\src\Model\LinearNetworkFlow\lnf-oc.f90"/>
</Filter>
```

### Regenerating IDM Files

After modifying DFN files, regenerate IDM modules:

```bash
cd utils/idmloader/scripts
python dfn2f90.py
```

---

## Testing

### Test Case: test_lnf_01.py

A 3-node linear network test that validates basic functionality.

**Network Topology:**
```
Node 1 (CHD h=1.0) ---- Node 2 (active) ---- Node 3 (CHD h=0.0)
        C=1.0                    C=1.0
```

**Expected Results:**
- Head at Node 2: 0.5 (midpoint between boundary heads)
- CHD flow at Node 1: -0.5 (outflow)
- CHD flow at Node 3: +0.5 (inflow)
- Budget balance: IN = OUT = 0.5

**Running the Test:**
```bash
cd autotest
pytest test_lnf_01.py -v
```

---

## Key Implementation Details

### 1. Matrix Assembly (CSR Format)

The LNF model uses Compressed Sparse Row (CSR) format for the coefficient matrix:

- `IA(n)` - Index to first entry of row n in JA/AMAT
- `JA(k)` - Column index for entry k
- `AMAT(k)` - Matrix value for entry k

For conductance:
```fortran
! Diagonal: sum of all conductances connected to node
amat(idiag) = amat(idiag) + cond

! Off-diagonal: negative conductance
amat(ipos) = amat(ipos) - cond
```

### 2. Flow Calculation

Intercell flows stored in `flowja` array (same structure as JA):

```fortran
do n = 1, nodes
  do ipos = ia(n) + 1, ia(n + 1) - 1
    m = ja(ipos)
    cond = this%conductance(ipos)
    flowja(ipos) = cond * (x(n) - x(m))
  end do
end do
```

### 3. Memory Management

Arrays are allocated via the Memory Manager for proper tracking:

```fortran
call mem_allocate(this%flowja, this%dis%nja, 'FLOWJA', this%memoryPath)
```

### 4. Pointer Setup

The solution sets pointers for X (head), RHS, and IBOUND:

```fortran
call this%set_xptr(this%x, this%x_old, 'X', 'XOLD', this%name)
call this%set_rhsptr(this%rhs, this%name)
call this%set_iboundptr(this%ibound, this%name)
```

### 5. Output Control Integration

The `oc_rp` call in `lnf_rp` is critical for loading PERIOD block settings:

```fortran
subroutine lnf_rp(this)
  if (this%inoc > 0) call this%oc%oc_rp()  ! Load save/print settings
  ! ... boundary package RP calls
end subroutine
```

---

## Future Enhancements

Potential future additions to the LNF model:

1. **Additional Boundary Packages**: WEL (wells), DRN (drains), RIV (rivers)
2. **Time-Varying Conductance**: Support for transient conductance values
3. **Observations Package**: Support for head and flow observations
4. **Newton-Raphson**: Nonlinear formulation for advanced applications
5. **Parallel Support**: PRT-style parallelization for large networks

---

## References

- MODFLOW 6 Developer Documentation
- `src/Model/GroundWaterFlow/gwf.f90` - Reference implementation
- `src/Model/ChannelFlow/chf.f90` - Similar network-style model
- `doc/mf6io/mf6ivar/` - Input variable documentation system

---

*Document created: February 5, 2026*
*Last updated: February 5, 2026*
