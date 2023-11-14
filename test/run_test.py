#!/usr/bin/env python3
import unittest
import numpy as np
import tempfile
import os
import nexusformat.nexus as nexus
import eniius
import json

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
            self.assertTrue('transformations' in comp1)
            self.assertTrue(isinstance(comp1['transformations'], nexus.NXtransformations))

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

    def test_save_json_from_mcstas(self):
        jsonfile = os.path.join(self.tmpdir.name, 'mcstas.json')
        instrfile = os.path.join(self.rootdir, 'instruments', 'isis_merlin.instr')
        wrapper = eniius.Eniius.from_mcstas(instrfile, self.detdat)
        wrapper.to_json(jsonfile)
        with open(jsonfile) as file:
            data = json.load(file)

        nx_class_fields = ('name', 'type', 'children', 'attributes')

        self.assertTrue(len(data) == 1)
        self.assertTrue('children' in data)
        children = data['children']
        self.assertTrue(len(children) == 1)

        root = children[0]
        for field in nx_class_fields:
            self.assertTrue(field in root)
        self.assertTrue('entry' == root['name'])
        self.assertTrue(len(root['children']) == 1)
        self.assertTrue(len(root['attributes']) == 1)
        self.assertTrue(root['attributes'][0]['values'] == 'NXentry')

        instrument = root['children'][0]
        for field in nx_class_fields:
            self.assertTrue(field in instrument)
        self.assertTrue(len(instrument['children']) == 14)
        self.assertTrue(len(instrument['attributes']) == 1)
        self.assertTrue(instrument['attributes'][0]['values'] == 'NXinstrument')

        named_children = [x for x in instrument['children'] if 'name' in x]
        self.assertTrue(len(named_children) == 13)
        self.assertTrue('mcstas' in [x['name'] for x in named_children])
        mcstas_children = [x for x in named_children if x['name'] == 'mcstas']
        self.assertTrue(len(mcstas_children) == 1)

        for field in nx_class_fields:
            self.assertTrue(field in mcstas_children[0])
        named_mcstas_children = [x for x in mcstas_children[0]['children'] if 'name' in x]
        self.assertTrue('declare' in [x['name'] for x in named_mcstas_children])
        declare = [x for x in named_mcstas_children if x['name'] == 'declare']
        self.assertTrue(len(declare) == 1)

        for field in nx_class_fields:
            self.assertTrue(field in declare[0])
        parameters = declare[0]['children']
        for name in ('slit_curv', 'num_slits', 'width', 'len', 'phase_time', 'E_min', 'E_max'):
            self.assertTrue(len([x for x in parameters if x['config']['name'] == name]) == 1)
        # Check that extraneous empty parameters do not get saved:
        self.assertTrue(len([x for x in parameters if x['config']['name'] == '']) == 0)

    def test_save_slit_json_from_mcstas(self):
        jsonfile = os.path.join(self.tmpdir.name, 'mcstas.json')
        instrfile = os.path.join(os.path.dirname(__file__), 'one_slit_explicit.instr')
        wrapper = eniius.Eniius.from_mcstas(instrfile)
        wrapper.to_json(jsonfile)
        with open(jsonfile) as file:
            data = json.load(file)

        nx_class_fields = ('name', 'type', 'children', 'attributes')

        self.assertTrue(len(data) == 1)
        self.assertTrue('children' in data)
        children = data['children']
        self.assertTrue(len(children) == 1)

        root = children[0]
        for field in nx_class_fields:
            self.assertTrue(field in root)
        self.assertTrue('entry' == root['name'])
        self.assertTrue(len(root['children']) == 1)
        self.assertTrue(len(root['attributes']) == 1)
        self.assertTrue(root['attributes'][0]['values'] == 'NXentry')

        instrument = root['children'][0]
        for field in nx_class_fields:
            self.assertTrue(field in instrument)
        self.assertTrue(len(instrument['children']) == 5)
        self.assertTrue(len(instrument['attributes']) == 1)
        self.assertTrue(instrument['attributes'][0]['values'] == 'NXinstrument')

        named_children = [x for x in instrument['children'] if 'name' in x]
        self.assertTrue(len(named_children) == 4)
        self.assertTrue('mcstas' in [x['name'] for x in named_children])
        mcstas_children = [x for x in named_children if x['name'] == 'mcstas']
        self.assertTrue(len(mcstas_children) == 1)

        for field in nx_class_fields:
            self.assertTrue(field in mcstas_children[0])
        named_mcstas_children = [x for x in mcstas_children[0]['children'] if 'name' in x]
        self.assertFalse('declare' in [x['name'] for x in named_mcstas_children])

    def test_one_arm_json_from_mcstas(self):
        # Tests that coordinate transformation involing an Arm() component is correct
        jsonfile = os.path.join(self.tmpdir.name, 'mcstas.json')
        instrfile = os.path.join(os.path.dirname(__file__), 'simple_one_arm.instr')
        wrapper = eniius.Eniius.from_mcstas(instrfile)
        wrapper.to_json(jsonfile, absolute_depends_on=True)
        with open(jsonfile) as file:
            data = json.load(file)
        # Address of component is 'entry/instrument/monitor'
        instrument = data['children'][0]['children'][0]
        self.assertEqual(instrument['name'], 'instrument')
        idx = [ii for ii,v in enumerate(instrument['children']) if 'name' in v and v['name'] == 'monitor']
        self.assertEqual(len(idx), 1)
        monitor = instrument['children'][idx[0]]
        idx = [ii for ii,v in enumerate(monitor['children']) if 'name' in v and v['name'] == 'transformations']
        self.assertEqual(len(idx), 1)
        transforms = monitor['children'][idx[0]]
        rot = transforms['children'][0]['config']['values']
        attrs = {v['name']:v['values'] for v in transforms['children'][0]['attributes']}
        self.assertTrue(np.abs(rot - 125.0) < 1.e-5)
        self.assertTrue(np.sum(np.abs(np.array(attrs['vector']) - [0., 1., 0.])) < 1.e-5)
        self.assertTrue(np.sum(np.abs(np.array(attrs['offset']) - [0., 0., 0.])) < 1.e-5)
        self.assertEqual(attrs['transformation_type'], 'rotation')


if __name__ == '__main__':
    unittest.main()

