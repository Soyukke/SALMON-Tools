# SALMON-Tools

Tools for SALMON (<https://salmon-tddft.jp/>)

# 1. poscar2salmon.jl

## Description

This script `poscar2salmon.jl` that converts POSCAR, which is a structural file for VASP, to an input file for SALMON.

## Requirement

- Julia 1.1.0

## Usage

Prepare POSCAR file and pseudopotential file.

Download the pseudopotential file from here.

<https://www.abinit.org/sites/default/files/PrevAtomicData/psp-links/lda_fhi>

Example of operation. file organization

```
.
├── examples
│   ├── POSCAR
│   └── psps
│       ├── 01-H.LDA.fhi
│       └── 08-O.LDA.fhi
└── poscar2salmon.jl
```

Execute the following.

```
cd examples
julia ../poscar2salmon.jl POSCAR
```

Load POSCAR and `psps/*.fhi` to generate the SALMON INPUT file ` salmon.inp`.