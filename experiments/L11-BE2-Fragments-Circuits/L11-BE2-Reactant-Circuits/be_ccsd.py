from functools import reduce
import numpy
from pyscf import gto, scf, ao2mo, cc
from pyscf import symm
from pyscf.tools import fcidump

ccsd_energies = []
# ccsd_t_energies = []
mf_energies = []
for i in range(15):
    mf = fcidump.to_scf(f'L11_reactant_BE2_f{i}/fcidump.txt', molpro_orbsym=True)
    # mf.mol.nelectron = mf.mol.nelectron*2
    mf.mol.verbose = 0
    mf.run()
    mf_energies.append(mf.e_tot)
    assert mf.converged
    mycc = cc.CCSD(mf)
    mycc.kernel()
    assert mycc.converged
    ccsd_energies.append(mycc.e_tot)
    ccsd_energy = mycc.e_tot
    # et = mycc.ccsd_t()
    # ccsd_t_energies.append(ccsd_energy + et)
    print(f'L11_reactant_BE2_f{i} CCSD energy: {ccsd_energy}')
    # print(f'L11_reactant_BE2_f{i} CCSD(T) energy: {ccsd_energy + et}')

print(ccsd_energies)

# Write a table of the energies
# Print a nicely formatted table of energies
header = f"{'i':^3} | {'HF energy':^18} | {'CCSD energy':^18}"
sep = f"{'-'*3}-+-{'-'*18}-+-{'-'*18}-+"
print(header)
print(sep)
for i in range(15):
    print(f"{i:^3} | {mf_energies[i]:>18.10f} | {ccsd_energies[i]:>18.10f}")