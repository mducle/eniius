import numpy as np
import scipy.io
import copy
import sys
import os

import nexusformat.nexus

# Monkey patch the nexus write function to use fixed-width ASCII NX_class labels
# because IBEX [ISISICP] only supports this format (as it uses the old napi.c)
old_wg = nexusformat.nexus.tree.NXFile._writegroup
def new_wg(self, group):
    links = old_wg(self, group)
    if group.nxpath != '' and group.nxpath != '/':
        path = self.nxpath + '/' + group.nxname
        if group.nxclass and group.nxclass != 'NXgroup':
            self[path].attrs['NX_class'] = np.array(group.nxclass, dtype='S')
    return links
nexusformat.nexus.tree.NXFile._writegroup = new_wg

from nexusformat.nexus import NXFile, NXfield, NXdata, nxload, NXtransformations, \
                              NXslit, NXsource, NXmoderator, NXnote, NXbeam, \
                              NXinstrument, NXfermi_chopper, NXdisk_chopper, \
                              nxopen, NXentry, NXcollection, NXsample, NXdetector, NXroot

def create_inst_nxs(outfile, inst_fun, ei, root_name='w1', det_file=None):
    n = 5
    with nxopen(outfile, 'w') as root:
        root[root_name] = NXentry()
        if (root_name == 'w1'):
            root[root_name + '/NXSPE_info'] = NXcollection(fixed_energy=NXfield(ei, units='meV'),
                                                           ki_over_kf_scaling=True,
                                                           psi=NXfield(np.nan, units='degrees'))
            root[root_name + '/definition'] = NXfield('NXSPE', version='1.3')
        root[root_name + '/program_name'] = NXfield('eniius', version='0.1')
        root[root_name + '/sample'] = NXsample()
        root[root_name + '/data'] = create_data(n, ei)
        root[root_name + '/instrument'] = inst_fun(ei)
        if (root_name == 'mantid_workspace_1'):
            inst_name = root['mantid_workspace_1/instrument/name'].nxvalue + ' '
            root['mantid_workspace_1/instrument/name'].replace(
                    NXfield(np.array(inst_name, dtype='S'), shape=(1,)))
        if det_file:
            # Get a list of instrument components to save
            cmps = ','.join([k for k in root[root_name + '/instrument'] if k != 'name']) + ' '
            root[root_name + '/instrument/components_to_save'] = NXnote( \
                    data=NXfield(np.array(cmps, dtype='S'), shape=(1,)),
                    description='Instrument components to be saved by IBEX/ISISICP',
                    type='text/plain')
            # Assumes ISIS format; cols=(det#, delta, L2, code, theta, phi, W_xyz, a_xyz, det_123)
            detdat = np.loadtxt(det_file, skiprows=3, usecols=(0, 2, 3, 4, 5))
            for fnm, idn in zip(['detectors', 'monitors'], [2, 1]):
                idx = np.where(detdat[:,2] == idn)[0]
                fd = {f'number_of_{fnm}': NXfield(np.array(len(idx), dtype='uint64'))}
                fd['azimuthal_angle'] = NXfield(detdat[idx,4], units='degree')
                fd['polar_angle'] = NXfield(detdat[idx,3], units='degree')
                fd['distance'] = NXfield(detdat[idx,1], units='metre')
                fd['detector_number'] = NXfield(detdat[idx,0].astype(np.int32))
                root[root_name + '/instrument/physical_' + fnm] = NXdetector(**fd)

def create_data(n, ei):
    dmat = np.random.rand(n, n)
    emat = np.random.rand(n, n) / 10.
    th = np.linspace(-25, 120, n)
    phi = np.zeros(n)
    en = np.linspace(ei/20, ei*0.9, n)
    wd = np.ones(n) / 10.
    dd = np.ones(n) * 4.
    dat = NXdata(data=NXfield(dmat, axes='polar:energy', signal=1), error=emat,
                 errors=emat, energy=NXfield(en, units='meV'),
                 azimuthal=th, azimuthal_width=wd, polar=phi, polar_width=wd, distance=dd)
    return dat

def let_instrument(ei):
    INST_TABLES = scipy.io.loadmat('inst_tables.mat')

    inst = NXinstrument(fermi=NXfermi_chopper(energy=ei))
    inst['name'] = NXfield(value='LET', short_name='LET')
    m_ = {'units':'metre'}
    # Defines the Source
    inst['source'] = NXsource(Name='ISIS', type='Spallation Neutron Source',
                              frequency=NXfield(10, units='hertz'), target_material='W')
    # Defines the moderator
    d_mod = NXtransformations(MOD_T_AXIS=NXfield(-25., transformation_type='translation',
                                                 vector=[0.,0.,1.], depends_on='.', **m_),
                              MOD_R_AXIS=NXfield(0., transformation_type='rotation',
                                                 vector=[0.,1.,0.], depends_on='MOD_T_AXIS', units='degree'))
    pulse = NXnote(type='ikcarp', data=[42.1304,0,0],
                   description='Empirical Ikeda-Carpenter type moderator pulse model')
    inst['moderator'] = NXmoderator(type='Liquid H2', temperature=NXfield(17., units='kelvin'),
                                    empirical_pulse_shape=pulse, transforms=d_mod)
    # Defines the divergences
    hdiv = NXdata(signal=NXfield(INST_TABLES['let_div'][:,1], unit='', name='Normalised Beam Profile'),
                  axes=NXfield(INST_TABLES['let_div'][:,0], unit='degree', name='Horizontal Divergence'))
    vdiv = NXdata(signal=NXfield(INST_TABLES['let_div'][:,3], unit='', name='Normalised Beam Profile'),
                  axes=NXfield(INST_TABLES['let_div'][:,2], unit='degree', name='Vertical Divergence'))
    inst['horiz_div'] = NXbeam(data=hdiv)
    inst['vert_div'] = NXbeam(data=vdiv)
    # Defines the choppers
    d_ch1 = NXtransformations(CH1_T_AXIS=NXfield(-17.17, transformation_type='translation',
                                                 vector=[0.,0.,1.], depends_on='.', **m_))
    sl = np.rad2deg(np.arctan2(0.04, 0.28) / 2.)
    inst['shaping_chopper'] = NXdisk_chopper(rotation_speed=NXfield(40., units='hertz'), radius=NXfield(0.28, **m_),
                                             type='contra_rotating_pair', slit_edges=[-sl, sl], transforms=d_ch1)
    d_ch5 = NXtransformations(CH5_T_AXIS=NXfield(-1.5, transformation_type='translation',
                                                 vector=[0.,0.,1.], depends_on='.', **m_))
    sl = np.rad2deg(np.arctan2(0.031, 0.28) / 2.)
    inst['mono_chopper'] = NXdisk_chopper(rotation_speed=NXfield(240., units='hertz'), radius=NXfield(0.28, **m_),
                                          type='contra_rotating_pair', slit_edges=[-sl, sl], transforms=d_ch5)
    return inst

def maps_instrument(ei):
    INST_TABLES = scipy.io.loadmat('inst_tables.mat')

    inst = NXinstrument(fermi=NXfermi_chopper(energy=ei))
    inst['name'] = NXfield(value='MAPS', short_name='MAPS')
    m_ = {'units':'metre'}
    # Defines the Fermi chopper
    fermi = inst['fermi']
    fermi_data = {'type':'sloppy', 'distance':NXfield(-1.857, **m_),
                  'rotation_speed':NXfield(600., units='hertz'), 
                  'radius':NXfield(0.049, **m_), 'r_slit':NXfield(1.3, **m_),
                  'slit':NXfield(0.002899, **m_), 'number':9, 'width':NXfield(0.0522, **m_)} 
    for k, v in fermi_data.items():
        fermi[k] = v
    # Defines the Aperture
    d_ap = NXtransformations(AP_AXIS=NXfield(-10.3290, transformation_type='translation',
                                             vector=[0.,0.,1.], depends_on='.', **m_))
    inst['aperture'] = NXslit(x_gap=NXfield(0.0989, **m_), y_gap=NXfield(0.0989, **m_),
                              transforms=d_ap)
    # Defines the Source
    inst['source'] = NXsource(Name='ISIS', type='Spallation Neutron Source',
                              frequency=NXfield(50, units='hertz'), target_material='W')
    # Defines the moderator
    d_mod = NXtransformations(MOD_T_AXIS=NXfield(-12., transformation_type='translation',
                                                 vector=[0.,0.,1.], depends_on='.', **m_),
                              MOD_R_AXIS=NXfield(32., transformation_type='rotation',
                                                 vector=[0.,1.,0.], depends_on='MOD_T_AXIS', units='degree'))
    pulse = NXdata(signal=NXfield(INST_TABLES['maps_mod'][:,1], unit='1/microsecond/meV', name='Intensity'),
                   axes=NXfield(INST_TABLES['maps_mod'][:,0], unit='microsecond', name='Time'))
    inst['moderator'] = NXmoderator(type='H20', temperature=NXfield(300, units='kelvin'),
                                    pulse_shape=pulse, transforms=d_mod)
    return inst

def merlin_instrument(ei):
    INST_TABLES = scipy.io.loadmat('inst_tables.mat')

    inst = NXinstrument(fermi=NXfermi_chopper(energy=ei))
    inst['name'] = NXfield(value='MERLIN', short_name='MER')
    m_ = {'units':'metre'}
    # Defines the Fermi chopper
    inst['name'].replace(NXfield(value='MERLIN', short_name='MER'))
    inst['fermi'].energy = 120.
    fermi = inst['fermi']
    fermi_data = {'type':'g', 'distance':NXfield(-1.8, **m_),
                  'rotation_speed':NXfield(600., units='hertz'), 
                  'radius':NXfield(0.005, **m_), 'r_slit':NXfield(99999., **m_),
                  'slit':NXfield(2.0e-4, **m_), 'number':125, 'width':NXfield(0.05, **m_)} 
    for k, v in fermi_data.items():
        fermi[k] = v
    # Defines the Aperture
    d_ap = NXtransformations(AP_AXIS=NXfield(-10.1570, transformation_type='translation',
                                             vector=[0.,0.,1.], depends_on='.', **m_))
    inst['aperture'] = NXslit(x_gap=NXfield(0.0967, **m_), y_gap=NXfield(0.0967, **m_),
                              transforms=d_ap)
    # Defines the Source
    inst['source'] = NXsource(Name='ISIS', type='Spallation Neutron Source',
                              frequency=NXfield(50, units='hertz'), target_material='W')
    # Defines the moderator
    d_mod = NXtransformations(MOD_T_AXIS=NXfield(-11.837, transformation_type='translation',
                                                 vector=[0.,0.,1.], depends_on='.', **m_),
                              MOD_R_AXIS=NXfield(0., transformation_type='rotation',
                                                 vector=[0.,1.,0.], depends_on='MOD_T_AXIS', units='degree'))
    pulse = NXdata(signal=NXfield(INST_TABLES['merlin_mod'][:,1], unit='1/microsecond/meV', name='Intensity'),
                   axes=NXfield(INST_TABLES['merlin_mod'][:,0], unit='microsecond', name='Time'))
    inst['moderator'] = NXmoderator(type='H20', temperature=NXfield(300, units='kelvin'),
                                    pulse_shape=pulse, transforms=d_mod)
    return inst

if __name__ == '__main__':
    create_inst_nxs('let_inst.nxspe', let_instrument, 3.7, 'mantid_workspace_1')
    create_inst_nxs('maps_inst.nxspe', maps_instrument, 400., 'mantid_workspace_1')
    create_inst_nxs('merlin_inst.nxspe', merlin_instrument, 120., 'mantid_workspace_1', 'detector.dat')
