#!/usr/bin/env python

import numpy, ctypes
from pyscf import lib
from pyscf import gto, scf, dft
from pyscf.geomopt import berny_solver
from berny import Berny, geomlib, Logger, optimize as optimize_berny

_loaderpath = '/home/jluis/src/pyscf/pruebas/dftd3/lib'
libdftd3 = numpy.ctypeslib.load_library('libdftd3.so', _loaderpath)

name = 'c2h4_c2h4'

mol = gto.Mole()
mol.atom = '''
C  -0.471925  -0.471925  -1.859111
C   0.471925   0.471925  -1.859111
H  -0.872422  -0.872422  -0.936125
H   0.872422   0.872422  -0.936125
H  -0.870464  -0.870464  -2.783308
H   0.870464   0.870464  -2.783308
C  -0.471925   0.471925   1.859111
C   0.471925  -0.471925   1.859111
H  -0.872422   0.872422   0.936125
H   0.872422  -0.872422   0.936125
H  -0.870464   0.870464   2.783308
H   0.870464  -0.870464   2.783308
'''
mol.basis = 'aug-cc-pvdz'
mol.charge = 0
mol.spin = 0
mol.symmetry = 0
mol.verbose = 4
mol.build()

mf_grad_scan = dft.RKS(mol).nuc_grad_method().as_scanner()
mf_grad_scan.base = scf.addons.remove_linear_dep_(mf_grad_scan.base)
mf_grad_scan.base.level_shift = 0.2
mf_grad_scan.base.verbose = 4
mf_grad_scan.base.xc = 'pbe0'
mf_grad_scan.base.grids.level = 3
mf_grad_scan.base.grids.prune = dft.gen_grid.nwchem_prune
mf_grad_scan.grid_response = False

def mf_grad_with_dftd3(geom):
    e_tot, g_rhf = mf_grad_scan(geom)
    mol = mf_grad_scan.mol
    func = 'pbe0'
    version = 4
    tz = 0
    coords = mol.atom_coords()
    itype = numpy.zeros(mol.natm, dtype=numpy.int32)
    for ia in range(mol.natm):
        symb = mol.atom_pure_symbol(ia)
        itype[ia] = lib.parameters.NUC[symb]
    edisp = numpy.zeros(1)
    grad = numpy.zeros((mol.natm,3)) 
    libdftd3.dftd3(ctypes.c_int(mol.natm),
             coords.ctypes.data_as(ctypes.c_void_p),
             itype.ctypes.data_as(ctypes.c_void_p),
             ctypes.c_char_p(func),
             ctypes.c_int(version),
             ctypes.c_int(tz),
             edisp.ctypes.data_as(ctypes.c_void_p),
             grad.ctypes.data_as(ctypes.c_void_p))
    lib.logger.info(mf_grad_scan,"* Disp Energy [au]: %12.8f" % edisp[0])
    lib.logger.info(mf_grad_scan,"* Disp Gradients [au]:")
    atmlst = range(mol.natm)
    for k, ia in enumerate(atmlst):
        symb = mol.atom_pure_symbol(ia)
        lib.logger.info(mf_grad_scan,"* %d %s %12.8f %12.8f %12.8f" \
        % (ia, symb, grad[k,0], grad[k,1], grad[k,2]))
    e_tot += edisp[0]
    g_rhf += grad
    lib.logger.info(mf_grad_scan,"* Total Energy [au]: %12.8f" % e_tot)
    lib.logger.info(mf_grad_scan,"* Total Gradients [au]:")
    atmlst = range(mol.natm)
    for k, ia in enumerate(atmlst):
        symb = mol.atom_pure_symbol(ia)
        lib.logger.info(mf_grad_scan,"* %d %s %12.8f %12.8f %12.8f" \
        % (ia, symb, g_rhf[k,0], g_rhf[k,1], g_rhf[k,2]))
    return e_tot, g_rhf

mf = berny_solver.as_pyscf_method(mol, mf_grad_with_dftd3)
mol = berny_solver.kernel(mf, maxsteps=50, assert_convergence=True)
xyzfile = name + '_opt.xyz'
fspt = open(xyzfile,'w')
coords = mol.atom_coords()*lib.param.BOHR
fspt.write('%d \n' % mol.natm)
fspt.write('%d %d\n' % (mol.charge, (mol.spin+1)))
for ia in range(mol.natm):
    symb = mol.atom_pure_symbol(ia)
    fspt.write('%s  %12.6f  %12.6f  %12.6f\n' % (symb, \
    coords[ia][0],coords[ia][1], coords[ia][2]))

