# Molecular-Docking-pipeline-Deleting-toxic-compounds.-

## Introduction

Azo Reductases, which degrade Azo dyes into colorless amines through reductive cleavage. This process requires a low molecular weight-reducing molecule such as FADH or NADH as an electron donor. 
Finally, the enzyme cleaves the azo-type bond (-N = N-) and transfers four electrons as a reducing equivalent. The azo dye works as an electron acceptor causing its discoloration and the formation of an intermediate (a highly toxic aromatic amine) that is subsequently degraded by an aerobic process (Singh RL, or microaerophilically, when conditions are anaerobic, membrane-bound azo reductases use a redox mediator as an electron shuttle.

The objective of the project is to degrade the colorimetric molecules indigo blue and acid orange 10, whose origin is the waste of the textile industry. It is currently known from the literature that azoreductases are a family of enzymes with a high potential to degrade a wide variety of polluting compounds. Inside the group of molecules that this enzyme can degrade are the azo dyes and indigoids. However, they have not been characterized specifically for the degradation of indigo blue and acid orange 10. This means that it would be necessary to study the interaction that laccase and azoreductase enzymes may have with both molecules.

We applied computational methods to study the binding and stability properties of indigo blue and acid orange 10. The azoreductase test was from Pseudomonas putida.

## Methodology 

Here I developed a bioinformatics pipeline to make molecular docking with the azoreductase protein.

1)	The sequences of the protein selected were searched in NCBI and its accession number of the crystallographic structures was retrieved.
   
2)	For the chemical structure of the two compounds to study, both were searched and retrieved from chemical databases.
   
3)	Both ligands and proteins were optimized before their preparation for molecular docking.
   
4)	The active site for the enzyme was searched in the literature. The box the proteins were placed where its active site is known to be located.
    
5)	The results of the molecular docking were evaluated considering the coupling energy and the type of interaction that the complex enzyme-substrate had.
    
6)	The best result of the docking was evaluated through a molecular dynamics simulation to determine the stability and bond of the complex

## Step 1.  LIGANDS STRUCTURE, OPTIMIZATION, AND PREPARATION

The structures of the colorimetric compounds were searched in the public database PubChem (PubChem blue indigo CID 10215 & acid orange CID 16015)
 
The molecular structures of the two dyes were downloaded and optimized using electronic structure calculations with Python software (Rdkit library) with the Merck Molecular Force Field (MMFF).


``` bash

#Install RDkit with anaconda environment 
conda install -c conda-forge rdkit
```

``` Python 
#import libraries from my browser

from rdkit import Chem
from rdkit.Chem import AllChem

#Molecule to minimize (Blue indigo)

m= Chem.MolFromSmiles('C1=CC=C2C(=C1)C(=C(N2)C3=NC4=CC=CC=C4C3=O)O')
m

#Add Hidrogens
m2=Chem.AddHs(m)
AllChem.EmbedMolecule(m2)

#Optimize molecule -> MMFF (Merck Molecular Force Field)

AllChem.MMFFOptimizeMolecule(m2)

#Import molecule optimized

from rdkit.Chem import Draw
img = Draw.MolToImage(m2)
img.show()

```
MMF is a force field to calculi energy and geometry in chemistry molecules. These approaches allow us to describe energy boundaries, flexion energy, solvation energy, etc....  (you can repeat the same approach for the Acid orange 10 ligand) 

How can we ensure that the geometries obtained correspond to a minimum on the potential energy surface, thus indicating a molecular geometry of the structures at their minimum energy?
 In addition to the optimization calculations, in both cases, electronic frequency calculations were performed using the Psi4 program. If you want to calculate this change, you will need to use the Psi4 program and not the Rdkit library (because does not make this step.) 

``` Python
import psi4

# molecule in Psi4
psi4_molecule = psi4.geometry("""
0 1
C
H 1 1.08
H 1 1.08 2 104.5
""")

# Set options for the calculation
psi4.set_options ({
    'basis': 'cc-pvdz',
    'scf_type': 'df',
    'e_convergence': 1e-8,
    'd_convergence': 1e-8
})

# Perform a frequency calculation
psi4_freqs = psi4.frequencies()

print ("Vibrational Frequencies:", psi4_freqs)
```
After the optimization, the ligands were prepared for docking using the Autodock Tool. The polar hydrogens and the partial charges were added, as the tutorial of Autodock Vina says. The ligands were prepared flexibly in both structures.







