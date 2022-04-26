# ABACUS-kit
A toolkit for generating ABACUS input files

Author: Inabahu

ABACUS is a impressively ab-initio program developed by Chinese group. It is a wheel designed for fastly generating input files for [ABACUS program](https://github.com/deepmodeling/abacus-develop). It's still under programming and I hope it can become a useful toolkit as vaspkit for vasp.
If anyone have questions or comments, you can ask in the issue page or contact me through [email](Inabahu@tju.edu.cn) or you can design more I/O interface to help expand the program funcitons.

## Prerequisites
Python3
numpy

User Customized Variables
- Pesudopotential database foldpath: "PP_path_custom"
- Basis sets foldpath: "BS_path_custom"

> Attention: It is just a wheel and user should check your input files to obtain reasonable computaiton results. And for some special parameters (like `nspin` tag for spin polarization), you should customize these by yourself.

## Tutorial
```
POSCAR_filepath="./biphen-P.vasp"
[atom_spec,atom_nu,a,b,c]=STRU_gen_from_POSCAR(POSCAR_filepath,"SG15_ONCV","pw",PP_path_custom,BS_path_custom)
a=[float(a_i) for a_i in a.strip().split("         ")]
Lattice_a=np.sqrt(a[0]**2+a[1]**2+a[2]**2)
b=[float(b_i) for b_i in b.strip().split("        ")]
Lattice_b=np.sqrt(b[0]**2+b[1]**2+b[2]**2)
c=[float(c_i) for c_i in c.strip().split("        ")]
Lattice_c=np.sqrt(c[0]**2+c[1]**2+c[2]**2)
KPT_gen('auto','gamma',Lattice_a,Lattice_b,Lattice_c,0.04)
total_atom_num=sum([float(nu) for nu in atom_nu])
INPUT_gen("relax_pw",PP_path_custom,len(atom_spec),int(8*total_atom_num))
```

- Specify the "PP_path_custom" and "BS_path_custom" at the beginning if you have special Pesudopotential database foldpath
- Firstly, you should give the POSCAR filepath for `POSCAR_filepath`
- then you can choose the calculation type ("pw" or "lcao") (now this tool only support "SG15_ONCV" Pesudopotential. Other type is updating)
- Change `0.04` in the KPT_gen function to your customized value
  - 0.1  Gamma only Coarse calculation
  - 0.03~0.04   Middle (suitable for geometry relaxation)
  - 0.01~0.02  High (suitable for self-consistent calcualtion)
- Choose the computation type
  - "relax_pw": geometry relaxation using pw method
  - "relax_lcao": geometry relaxation using lcao method
  - "scf_pw": self-consistent calcualtion using pw method
  - "scf_lcao": self-consistent calcualtion using lcao method
  - "dos_pw": density of state (DOS) calcualtion using pw method (two-step calculation: INPUT-1 and INPUT-2)
  - "dos_lcao": density of state (DOS) calcualtion using lcao method (two-step calculation: INPUT-1 and INPUT-2)
  - "band_pw": band structure (DOS) calcualtion using pw method (two-step calculation: INPUT-1 and INPUT-2)
  - "band_lcao": band structure (DOS) calcualtion using lcao method (two-step calculation: INPUT-1 and INPUT-2)

