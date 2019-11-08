Lecture\_12
================

Prepare protein structure for docking
=====================================

We want to download the 1HSG PDB sturcture and then produce a "protein-only" and "ligand-only" new seperate PDB files

``` r
library(bio3d)

get.pdb("1HSG")
```

    ## Warning in get.pdb("1HSG"): ./1HSG.pdb exists. Skipping download

    ## [1] "./1HSG.pdb"

Produce a "1hsg\_protein.pdb" and "1hsg\_ligand.pdb" file
=========================================================

``` r
pdb <- read.pdb("1HSG.pdb")
pdb
```

    ## 
    ##  Call:  read.pdb(file = "1HSG.pdb")
    ## 
    ##    Total Models#: 1
    ##      Total Atoms#: 1686,  XYZs#: 5058  Chains#: 2  (values: A B)
    ## 
    ##      Protein Atoms#: 1514  (residues/Calpha atoms#: 198)
    ##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
    ## 
    ##      Non-protein/nucleic Atoms#: 172  (residues: 128)
    ##      Non-protein/nucleic resid values: [ HOH (127), MK1 (1) ]
    ## 
    ##    Protein sequence:
    ##       PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYD
    ##       QILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNFPQITLWQRPLVTIKIGGQLKE
    ##       ALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTP
    ##       VNIIGRNLLTQIGCTLNF
    ## 
    ## + attr: atom, xyz, seqres, helix, sheet,
    ##         calpha, remark, call

``` r
ligand <- atom.select(pdb,"ligand", value=T)
write.pdb(ligand, file="1hsg_ligand.pdb")
```

``` r
protein <- atom.select(pdb,"protein", value=T)
write.pdb(protein, file="1hsg_protein.pdb")
```
