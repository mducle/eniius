#!/usr/bin/env python3
import unittest
import numpy as np
import tempfile
import os
import nexusformat.nexus as nexus
import eniius

class EniiusTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.tmpdir = tempfile.TemporaryDirectory()
        cls.rootdir = os.path.dirname(os.path.realpath(eniius.__file__))
        cls.detdat = os.path.join(cls.rootdir, 'instruments', 'detector.dat')

    @classmethod
    def tearDownClass(cls):
        with open('success', 'w') as f:
            f.write('success')
        cls.tmpdir.cleanup()

    def test_save_nxs_from_mcstas(self):
        nxsfile = os.path.join(self.tmpdir.name, 'mcstas.nxs')
        instrfile = os.path.join(self.rootdir, 'instruments', 'isis_merlin.instr')
        wrapper = eniius.Eniius.from_mcstas(instrfile, self.detdat)
        wrapper.to_icp(nxsfile)
        with nexus.nxload(nxsfile) as nxs:
            self.assertTrue('mantid_workspace_1' in nxs)
            self.assertTrue('instrument' in nxs['mantid_workspace_1'])
            nxinst = nxs['mantid_workspace_1/instrument']
            self.assertTrue(isinstance(nxinst, nexus.NXinstrument))
            self.assertTrue('physical_detectors' in nxinst)
            self.assertTrue('physical_monitors' in nxinst)
            comp1 = nxinst[dir(nxinst)[0]]
            self.assertTrue('mcstas' in comp1)
            self.assertTrue(isinstance(comp1['mcstas'].nxvalue, str))
            self.assertTrue(isinstance(comp1['mcstas']._value, np.ndarray))
            self.assertTrue(isinstance(comp1['mcstas']._value, np.ndarray))
            self.assertEqual(comp1['mcstas']._value.dtype.kind, 'S')
            self.assertTrue('transforms' in comp1)
            self.assertTrue(isinstance(comp1['transforms'], nexus.NXtransformations))

    def test_save_nxspe_from_merlin(self):
        Ei = 180.
        nxspefile = os.path.join(self.tmpdir.name, 'horace.nxspe')
        wrapper = eniius.Eniius(eniius.horace.merlin_instrument(Ei), self.detdat)
        wrapper.to_nxspe(nxspefile)
        with nexus.nxload(nxspefile) as nxspe:
            root = nxspe[dir(nxspe)[0]]
            self.assertTrue('NXSPE_info' in root)
            self.assertEqual(root['NXSPE_info/fixed_energy'].nxvalue, Ei)
            self.assertTrue('instrument' in root)
            self.assertTrue('fermi' in root['instrument'])
            self.assertEqual(root['instrument/fermi/energy'].nxvalue, Ei)
            self.assertTrue('sample' in root)


    def test_save_nxspe_from_let(self):
        Ei = 3.7
        nxspefile = os.path.join(self.tmpdir.name, 'horace.nxspe')
        wrapper = eniius.Eniius(eniius.horace.let_instrument(Ei), self.detdat)
        wrapper.to_nxspe(nxspefile)
        with nexus.nxload(nxspefile) as nxspe:
            root = nxspe[dir(nxspe)[0]]
            self.assertTrue('NXSPE_info' in root)
            self.assertEqual(root['NXSPE_info/fixed_energy'].nxvalue, Ei)
            self.assertTrue('instrument' in root)
            self.assertTrue('fermi' in root['instrument'])
            self.assertEqual(root['instrument/fermi/energy'].nxvalue, Ei)
            self.assertTrue('sample' in root)


if __name__ == '__main__':
    unittest.main()

