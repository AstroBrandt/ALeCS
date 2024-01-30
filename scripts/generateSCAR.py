"""
Usage: generateSCAR.py [options]

Options:
    -h --help               Show this screen
    --geometry=<file>       Geometry file for the molecule. The script is written to use .xyz or .pdb files
    --ALECS=<dir>           ALECS directory
    --kmax=<kmax>           Maximum k value [default: 6]
    --lEmin=<lEmin>         Minimum log10(E) [default: 1]
    --lEmax=<lEmax>         Maximum log10(E) [default: 6]
"""
import numpy as np
import scipy.special as sis
import scipy.interpolate as si
import matplotlib.pyplot as plt

from string import digits
from glob import glob
from docopt import docopt 

kmax = 6
atomBasisAllowed = ['H', 'He', 'C', 'N', 'O', 'P', 'S', 'Ar']
a0 = 5.29E-9
a0A = 5.29E-1
Ipot = {"H":13.59844, "He":24.58738, "B":8.29803, "C":11.26030, "N":14.53414, "O":13.61806, "P":10.48669, "S":10.36001, 
        "Ar":15.75962, "Ne":21.5646,}

remove_digits = str.maketrans('', '', digits)


def distance_matrix(mol):
    Natoms = mol['Natoms']
    dmatrix = np.zeros((Natoms, Natoms))
    i = 0
    while i < Natoms:
        ci = mol['coords'][i]
        j = 0
        while j < Natoms:
            cj = mol['coords'][j]
            dc = ci - cj
            dmatrix[i,j] = np.sqrt(dc[0]*dc[0] + dc[1]*dc[1] + dc[2]*dc[2])*1E-8
            j += 1
        i += 1
    return dmatrix

cacheDict = {}
def calcSCARCoeff(atomFits, N, dmat, atoms, E, i, k):
    if (i,k) in cacheDict.keys():
        return cacheDict[(i,k)]
    else:
        if k == 1:
            return 1
        fac = (N - k + 1.)/(N - 1.)
        sumFac = 0
        j = 1
        sigi = a0*a0*atomFits[atoms[i-1]](E)
        while j <= N:
            if j != i:
                sigj = a0*a0*atomFits[atoms[j-1]](E)
                dij = dmat[i-1,j-1]
                alphaij = max(4*np.pi*dij*dij, sigi, sigj)
                scri = calcSCARCoeff(atomFits, N, dmat, atoms, E, j, k-1)
                sumFac += sigj * scri / alphaij
            j += 1
        cacheDict[(i,k)] = fac*sumFac
        return fac*sumFac
    
def calcSCAR(atomFits, mol, E, kmax=None, benchmark=False):
    cacheDict.clear()
    Natoms = mol['Natoms']
    dmat = mol['dmatrix']
    atoms = mol['atoms']
    xsec = 0.0
    iat = 1
    if kmax == None:
        kmax = Natoms
    if kmax > Natoms:
        kmax = Natoms
    while iat <= Natoms:
        k = 1
        si = 0.
        kmaxIter = Natoms
        if benchmark:
            kmaxIter = kmax
        while (k <= kmaxIter):
            si += (-1)**(k+1) * calcSCARCoeff(atomFits, Natoms, dmat, atoms, E, iat, k)/sis.factorial(k)
            k += 1
        xsec += si * a0*a0*atomFits[atoms[iat-1]](E)
        iat += 1
    return xsec

if __name__ == '__main__':
    args = docopt(__doc__)
    geometry = args['--geometry']
    ALECS = args['--ALECS']
    kmax = int(args['--kmax'])
    lEmin = float(args['--lEmin'])
    lEmax = float(args['--lEmax'])

    print("Reading geometry file: ", geometry)
    print("ALECS directory: ", ALECS)
    print("Maximum k value: ", kmax)

    atomData = {}
    atomFits = {}
    for atom in atomBasisAllowed:
        atomData[atom] = np.loadtxt(ALECS + "/ion_xs/atoms/" + atom + ".xs", skiprows=2, comments="#")
        atomFits[atom] = si.interp1d(atomData[atom][:,0], atomData[atom][:,1], kind='linear')
    
    indFold = geometry.rfind('/')
    mol = geometry[indFold+1:geometry.rfind('.')]
    extension = geometry[geometry.rfind('.')+1:]

    molData = {}
    molData['Natoms'] = 0
    molData['atoms'] = []
    molData['coords'] = []
    if extension == 'xyz':
        with open(geometry, 'r') as f:
            for line in f.readlines():
                if line[0] == '#':
                    continue
                s = line.split()
                if len(s) < 4:
                    continue
                ai = s[0].translate(remove_digits)
                ci = np.array([float(si) for si in s[1:]])
                if ai not in atomBasisAllowed:
                    print("Molecule has atoms not in the database. Either add the cross sections or request them on the ALeCS github.")
                    print("Exiting...")
                    import sys
                    sys.exit()
                molData['atoms'].append(ai)
                molData['coords'].append(ci)
                molData['Natoms'] += 1
        molData['dmatrix'] = distance_matrix(molData)
    if extension == 'pdb':
        with open(geometry, 'r') as f:
            for line in f.readlines():
                if 'HETATM' in line:
                    s = line.split()
                    if s[2] not in atomBasisAllowed:
                        print("Molecule has atoms not in the database. Either add the cross sections or request them on the ALeCS github.")
                        print("Exiting...")
                        import sys
                        sys.exit()
                    molData['Natoms'] += 1
                    molData['atoms'].append(s[2])
                    molData['coords'].append(np.array([float(s[4]), float(s[5]), float(s[6])]))
        molData['dmatrix'] = distance_matrix(molData)
    
    Nene = 256
    eArr = np.logspace(lEmin, lEmax, Nene)
    sigArr = np.zeros(Nene)
    for (j, ei) in enumerate(eArr):
        sigArr[j] = calcSCAR(atomFits, molData, ei, kmax=kmax)/(a0*a0)

    with open(mol+".xs", 'w+') as f:
        f.write("# %s cross section\n"%mol)
        f.write("%d\n"%Nene)
        for (ei, ji) in zip(eArr, sigArr):
            f.write("%.6e %.6e\n" % (ei, ji))
        
    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111)
    ax.plot(eArr, sigArr)

    ax.set_xscale('linear')
    ax.set_xlim(0, 300)
    ax.set_xlabel("Energy (eV)", fontsize=16)
    ax.set_ylabel("Cross section (a$_0^2$)", fontsize=16)
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.grid(visible=True, which='major', axis='both', linestyle='--', 
            linewidth=0.75, color='k', alpha=0.5)
    fig.savefig(mol+".pdf", bbox_inches="tight")
