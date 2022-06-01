# bcfw

Supplementary Sage code that accompanies the paper:

"The BCFW Amplituhedron Triangulation" 
Chaim Even-Zohar, Tsviqa Lakrec, and Ran J. Tessler.
Preprint available from: https://arxiv.org/abs/2112.02703

The code in "bcfw.sage" implements constructions that are used in the proof that the BCFW cells triangulate the amplituhedron. Namely, the construction of a domino matrix given a chord diagram, and finding a separating functionary between the BCFW cells that correspond to two given chord diagrams. 

## Example
```
sage
> load("bcfw.sage")
> list_separators(7)
```
