from ase.io import read
from ase.io.utils import PlottingVariables

atoms = read('mg2nih4-28.pwi', format='espresso-in')
pv = PlottingVariables(atoms, rotation='30x,-70y,15z')
print(pv.T)
