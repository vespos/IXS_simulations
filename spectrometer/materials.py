
#Materials and crystal planes

import sys
# sys.path.append(r'C:\Users\espov\Documents\Python\xrt')
sys.path.append('/mnt/c/Users/espov/Documents/Python/xrt')
import xrt.backends.raycing.materials as rmats

matRh = rmats.Material(
    elements=r"Rh",
    kind=r"mirror",
    rho=12.41,
    table=r"Chantler",
    name=None)

Si111_300K = rmats.CrystalSi(
    tK=300.0,
    name=None)

Si111_125K = rmats.CrystalSi(
    tK=125.0,
    name=None)

Si555_300K = rmats.CrystalSi(
    tK=300.0,
    hkl=[5, 5, 5],
    name=None)


Quartz102_300K = rmats.CrystalFromCell(
    name=r"alphaQuartz",
    hkl=[1, 0, 2],
    a=4.91304,
    c=5.40463,
    gamma=120,
    atoms=[14, 14, 14, 8, 8, 8, 8, 8, 8],
    atomsXYZ=[[0.4697, 0.0, 0.0], [-0.4697, -0.4697, 0.3333333333333333], [0.0, 0.4697, 0.6666666666666666], [0.4125, 0.2662, 0.1188], [-0.1463, -0.4125, 0.4521], [-0.2662, 0.1463, -0.2145], [0.1463, -0.2662, -0.1188], [-0.4125, -0.1463, 0.2145], [0.2662, 0.4125, 0.5479]],
    tK=300.0,
    table=r"Chantler")

C111_300K = rmats.CrystalDiamond(
    d=2.0595,
    elements=r"C",
    rho=3.51,
    name=r"Diamond")

C333_300K = rmats.CrystalDiamond(
    d=2.0595/3,
    elements=r"C",
    rho=3.51,
    name=r"Diamond")

Si844_125K = rmats.CrystalSi(
    tK=125.0,
    hkl=[8, 4, 4],
    name=None)

Si753_125K = rmats.CrystalSi(
    tK=125.0,
    hkl=[7, 5, 3],
    name=None)

Si220_300K = rmats.CrystalSi(
    tK=300,
    hkl=[2, 2, 0],
    name=None)

Si220_125K = rmats.CrystalSi(
    tK=125,
    hkl=[2, 2, 0],
    name=None)

Si333_125K = rmats.CrystalSi(
    tK=125,
    hkl=[3, 3, 3],
    name=None)

Si333_300K = rmats.CrystalSi(
    tK=125,
    hkl=[3, 3, 3],
    name=None)

Si311_300K = rmats.CrystalSi(
    tK=300,
    hkl=[3, 1, 1],
    name=None)

Be = rmats.Material(
    elements=r"Be",
    kind=r"lens",
    rho=1.85,
    name=None)

Si400_300K = rmats.CrystalSi(
    tK=300,
    hkl=[4, 0, 0],
    name=None)

Si400_125K = rmats.CrystalSi(
    tK=125,
    hkl=[4, 0, 0],
    name=None)

Si444_300K = rmats.CrystalSi(
    tK=300,
    hkl=[4, 4, 4],
    name=None)

Si800_300K = rmats.CrystalSi(
    tK=300,
    hkl=[8, 0, 0],
    name=None)

Si844_300K = rmats.CrystalSi(
    tK=300,
    hkl=[8, 4, 4],
    name=None)

Si777_300K = rmats.CrystalSi(
    tK=300,
    hkl=[7, 7, 7],
    name=None)

Si884_300K = rmats.CrystalSi(
    tK=300,
    hkl=[8, 8, 4],
    name=None)

Si999_300K = rmats.CrystalSi(
    tK=300,
    hkl=[9, 9, 9],
    name=None)

Si555_300K = rmats.CrystalSi(
    tK=300,
    hkl=[5, 5, 5],
    name=None)

Si553_300K = rmats.CrystalSi(
    tK=300,
    hkl=[5, 5, 3],
    name=None)

Si880_300K = rmats.CrystalSi(
    tK=300,
    hkl=[8, 8, 0],
    name=None)

Si440_300K = rmats.CrystalSi(
    tK=300,
    hkl=[4, 4, 0],
    name=None)

Si664_300K = rmats.CrystalSi(
    tK=300,
    hkl=[6, 6, 4],
    name=None)

Si660_300K = rmats.CrystalSi(
    tK=300,
    hkl=[6, 6, 0],
    name=None)

Ge400_300K = rmats.CrystalDiamond(
    hkl=[4, 0, 0],
    d=1.414456,
    elements=r"Ge",
    name=None)

Ge800_300K = rmats.CrystalDiamond(
    hkl=[8, 0, 0],
    d=0.707228,
    elements=r"Ge",
    name=None)

Ge111_300K = rmats.CrystalDiamond(
    hkl=[1, 1, 1],
    d=3.2665,
    elements=r"Ge",
    name=None)

Ge220_300K = rmats.CrystalDiamond(
    hkl=[2, 2, 0],
    d=2.0003428585799985,
    elements=r"Ge",
    name=None)

Ge311_300K = rmats.CrystalDiamond(
    hkl=[3, 1, 1],
    d=1.7058981216243407,
    elements=r"Ge",
    name=None)

Ge933_300K = rmats.CrystalDiamond(
    hkl=[9, 3, 3],
    d=0.5686327072081135,
    elements=r"Ge",
    name=None)


Nickel = rmats.Material(
    elements=r"Ni",
    rho=8.9,
    name=r"Nickel")

Carbon = rmats.Material(
    elements=r"C",
    rho=2.2,
    name=r"Carbon")

SiliconSubstrate = rmats.Material(
    elements=r"Si",
    rho=2.33,
    name=r"SiliconSubsrate")
