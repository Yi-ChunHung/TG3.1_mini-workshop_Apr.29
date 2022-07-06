# 🧭 ftn58ToolBox

The tool box contains codes that improve quality of life while dealing with ftn58 or ftn58sparse. The codes are functions with their function be indicated either in their name or as their comments in second line.

## 🚩 [RSPCD](./RSPCD)

The code plots the distribution of a wave function resolving to atoms.

## 🚩 [RS_Hamiltonian](./RS_Hamiltonian.m)

The code tidies up the matrix elements of ftn58sparse for each translation.

## 🚩 [Stn58sparse2ftn58](./Sftn58sparse2ftn58.m)

The code tranforms a Sftn58sparse into a ftn58.

## 🚩 [check_onsite](./check_onsite.m)

The code displays the onsite energies of a ftn58sparse.

## 🚩 [ftn58sparse2ftn58](./ftn58sparse2ftn58.m)

The code transforms a ftn58sparse into a ftn58.

## 🚩 [ftn58sparse_transform_representation](./ftn58sparse_transform_representation.m)

The code transforms the representation of a ftn58sparse from Bloch representation to Wannier representation and vice versa.

## 🚩 [non_repeat_ftn58sparse](./non_repeat_ftn58sparse.m)

The code makes a ftn58sparse to have a unique [ij,dd]. Using such a ftn58sparse can accelerate most codes in this repository.

## 🚩 [readAinfo](./readAinfo.m)

The code displays the .Ainfo of a ftn58sparse as a table, which enables users to more efficiently access the information in .Ainfo.

## 🚩 [reconstruct_ftn58spasre](./reconstruct_ftn58sparse.m)

The code removes orbitals and their related hoppings on given atoms in a ftn58spasre.
