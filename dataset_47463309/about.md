
## Links
- [link to item](https://github.com/orgs/open-reaction-database/projects/3/views/1?pane=issue&itemId=47463309)
- [link to paper](https://doi.org/10.1038/s41467-023-42446-5)
- [link to paper SI](https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-023-42446-5/MediaObjects/41467_2023_42446_MOESM1_ESM.pdf)

## Workflow
- I downloaded SI `22826312.zip` from data source and extracted [Fig. 3.xlsx](Fig.%203.xlsx) from the downloaded file.
- I saved the first sheet as `sheet0.csv`. This sheet contains all reaction info.
- I found chemical structures of P1-P13 in SI Fig 7, and the names in Fig 3a
  - both `molscribe` and `mathpix` do not work, use vendor search for smiles
```json
[
  {
    "P1": "[Ir(dFCF3ppy)2(dtbbpy)]PF6",
    "P2": "Ir(ppy)3",
    "P3": "[Ir(ppy)2(bpy)]PF6",
    "P4": "Ir(dF-ppy)3",
    "P5": "[Ru(bpy)3]Cl2•6H2O",
    "P6": "[Ru(bpy)3](PF6)2",
    "P7": "Mes-Acr-ClO4",
    "P8": "Na2-Eosin Y",
    "P9": "Rhodamine B",
    "P10": "Methylene blue",
    "P11": "p-Me TPT",
    "P12": "p-OMe TPT",
    "P13": "Riboflavin tetrabutyrate"
  },
  {
    "P1": "F[P](F)(F)(F)(F)F.CC(C)(C)C1=CC=[N@H]2C(=C1)C3=CC(=CC=[N@@H]3[Ir]2456c7cc(F)cc(F)c7C8=CC=C(C=[N]48)C(F)(F)F)C(C)(C)C.Fc9cc(F)c(C%10=[N]5C=C(C=C%10)C(F)(F)F)c6c9",
    "P2": "C1=CC=C([C-]=C1)C2=CC=CC=N2.C1=CC=C([C-]=C1)C2=CC=CC=N2.C1=CC=C([C-]=C1)C2=CC=CC=N2.[Ir+3]",
    "P3": "C1=CC=C([C-]=C1)C2=CC=CC=N2.C1=CC=C([C-]=C1)C2=CC=CC=N2.C1=CC=NC(=C1)C2=CC=CC=N2.F[P-](F)(F)(F)(F)F.[Ir+3]",
    "P4": "FC1=CC([Ir]C2=C(C3=NC=CC=C3)C(F)=CC(F)=C2)=C(C4=NC=CC=C4)C(F)=C1.FC5=CC(F)=CC=C5C6=NC=CC=C6",
    "P5": "[H]O[H].[H]O[H].[H]O[H].[H]O[H].[H]O[H].[H]O[H].Cl[Ru]Cl.c1ccc(nc1)-c2ccccn2.c3ccc(nc3)-c4ccccn4.c5ccc(nc5)-c6ccccn6",
    "P6": "[Ru++].F[P-](F)(F)(F)(F)F.F[P-](F)(F)(F)(F)F.c1ccc(nc1)-c2ccccn2.c3ccc(nc3)-c4ccccn4.c5ccc(nc5)-c6ccccn6",
    "P7": "Cc1cc(C)c(-c2c3ccccc3[n+](C)c3ccccc23)c(C)c1.[O-]Cl(=O)(=O)=O",
    "P8": "C1=CC=C(C(=C1)C2=C3C=C(C(=O)C(=C3OC4=C(C(=C(C=C24)Br)O)Br)Br)Br)C(=O)O.[Na].[Na]",
    "P9": "CCN(CC)C1=CC2=C(C=C1)C(=C3C=CC(=[N+](CC)CC)C=C3O2)C4=CC=CC=C4C(=O)O.[Cl-]",
    "P10": "CN(C)C1=CC2=C(C=C1)N=C3C=CC(=[N+](C)C)C=C3S2.[Cl-]",
    "P11": "Cc1ccc(-c2cc(-c3ccc(C)cc3)[o+]c(-c3ccc(C)cc3)c2)cc1.[B-](F)(F)(F)F",
    "P12": "COc1ccc(-c2cc(-c3ccc(OC)cc3)[o+]c(-c3ccc(OC)cc3)c2)cc1.[B-](F)(F)(F)F",
    "P13": "CCCC(=O)OCC(C(C(CN1C2=C(C=C(C(=C2)C)C)N=C3C1=NC(=O)NC3=O)OC(=O)CCC)OC(=O)CCC)OC(=O)CCC"
  }
]
```
- I found chemical structures of L1-L11 in Fig 3a
  - use `mathpix` on [L1_to_L11.PNG](L1_to_L11.PNG)
    - it got L3 wrong, use `molscribe` for it instead
```json
{
  "L1": "C1COCCN1",
  "L2": "C[Si](C)(C)OC(c1ccccc1)(c1ccccc1)[C@@H]1CCCN1",
  "L3": "C[C@@H]1[NH2+][C@@H](C(C)(C)C)N(C)C1=O.[O-]S(=O)(=O)C(F)(F)F",
  "L4": "C[Si](C)(C)OC(c1cc(C(F)(F)F)cc(C(F)(F)F)c1)(c1cc(C(F)(F)F)cc(C(F)(F)F)c1)[C@H]1CCCN1",
  "L5": "O=C(O)[C@H]1CCCN1",
  "L6": "CN1C(=O)[C@H](Cc2ccccc2)N[C@@H]1C(C)(C)C",
  "L7": "O=C(O)[C@@H]1C[C@@H]2CCCC[C@@H]2N1",
  "L8": "O=C(O)[C@H]1CCCCN1",
  "L9": "NC(=O)C1CCCN1",
  "L10": "O=C(O)[C@@H]1Cc2ccccc2N1",
  "L11": "O=C(O)C1CSCN1"
}
```
- I found chemical structures of S1-S10 in Fig 3a
```json
{
  "S1": "c1ccc(C(=O)CBr)cc1",
  "S2": "Cc1ccc(C(=O)CBr)cc1",
  "S3": "Oc1ccc(C(=O)CBr)cc1",
  "S4": "COc1ccc(C(=O)CBr)cc1",
  "S5": "c2ccccc2c1ccc(C(=O)CBr)cc1",
  "S6": "FC(F)(F)c1ccc(C(=O)CBr)cc1",
  "S7": "Fc1ccc(C(=O)CBr)cc1",
  "S8": "Clc1ccc(C(=O)CBr)cc1",
  "S9": "Brc1ccc(C(=O)CBr)cc1",
  "S10": "FC(F)(F)Oc1ccc(C(=O)CBr)cc1"
}
```
- I found compound 4 is used as the internal standard
```json
{
  "name": "1-phenylprop-2-yn-1-ol",
  "smiles": "OC(C#C)c1ccccc1"
}
```
- I identified the reaction template from Fig. 3a
```python
lhs = "[C:1]#[C:2][C:3][C:4][C:5][C:6][C:7]=[O:8].[c:9]1[c:10][c:11][c:12]([C:13](=[O:14])[C:15][Br:16])[c:17][c:18]1"
rhs = "[c:9]1[c:10][c:11][c:12]([C:13](=[O:14])[C:15][C:6]([C:5][C:4][C:3][C:2]#[C:1])([C:7]=[O:8]))[c:17][c:18]1"
```
- The reported *ee* value was calculated s.t. positive - R config, negative - S config
  ```text
  Among these catalysts, the chiral products obtained with L2 and L3 are mainly R-configuration, 
  while those obtained with L4 and L10 are mainly S-configuration.
  ```