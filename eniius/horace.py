import numpy as np
import scipy.io
import os

from nexusformat.nexus import *

THISFOLDER = os.path.dirname(os.path.realpath(__file__))
MOD_TABLES = {}
LET_TABLES = None
MU_ = {'units':'metre'}
INST_FACES = {'maps':'TS1_S01_Maps.mcstas', 'merlin':'TS1_S04_Merlin.mcstas', 'let':'TS2.imat'}


def get_let_divergences(ei, version=2):
    global LET_TABLES
    if LET_TABLES is None:
        LET_TABLES = scipy.io.loadmat(f'{THISFOLDER}/instruments/horace_let_tables.mat')
    htab = LET_TABLES[f'ver{version}_horiz_div']
    vtab = LET_TABLES[f'ver{version}_vert_div']
    def get_div(divtab, lam0, typestr):
        angdeg = divtab['angdeg'][0][0].flatten()
        ang = angdeg * np.pi / 180.
        lam = divtab['lam'][0][0].flatten()
        S = divtab['S'][0][0]
        if lam0 < lam[0] or lam0 > lam[-1]:
            raise RuntimeError('The incident neutron wavelength lies outside the range of the divergence lookup table')
        profile = np.array([np.interp(lam0, lam, S[ii,:]) for ii in range(len(ang))])
        profile = (profile / np.sum(profile)) / np.mean(np.diff(ang))
        return NXdata(signal=NXfield(profile, unit='', name='Normalised Beam Profile'),
                      axes=NXfield(angdeg, unit='degree', name=typestr))
    lam = np.sqrt(81.80420126 / ei)
    hdiv = get_div(htab, lam, 'Horizontal Divergence')
    vdiv = get_div(vtab, lam, 'Vertical Divergence')
    return hdiv, vdiv


def load_mcstas_moderator(instrument):
    modfile = os.path.join(THISFOLDER, 'mcstas-comps', 'contrib', 'ISIS_tables', INST_FACES[instrument])
    with open(modfile, 'r') as f:
        dat = f.read().replace('(','').replace(')','').split('\n')
    id0, en, intens = (dat.index(' time '), [], [])
    n = dat.index(' time ', id0 + 1) - id0 - 6
    # time data originally in ns, energy in MeV
    t = np.loadtxt(dat[(id0+1):(id0+n)], usecols=(0)) / 100
    dt = np.diff(t)
    while True:
        e0 = np.loadtxt(dat[id0-2].split()[2:5:2]) * 1.e9
        i0 = np.loadtxt(dat[(id0+1):(id0+n)], usecols=(1))
        en.append(np.mean(e0))
        intens.append((i0[:-1] / dt) / (e0[1] - e0[0]))
        try:
            id0 = dat.index(' time ', id0 + 1)
        except ValueError:
            break
    return {'en':np.log(en), 'intens':np.array(intens), 't':t}


def get_moderator_time_pulse(instrument, ei):
    instrument = instrument.lower()
    global MOD_TABLES
    if instrument not in MOD_TABLES:
        MOD_TABLES[instrument] = load_mcstas_moderator(instrument)
    TMAX = 2000
    t, en, intens = (MOD_TABLES[instrument]['t'], MOD_TABLES[instrument]['en'], MOD_TABLES[instrument]['intens'])
    kp = np.where(t < TMAX)[0]
    ie = np.where(en < np.log(ei))[0][-1]
    frac = (np.log(ei) - en[ie]) / (en[ie+1] - en[ie])
    return intens[ie,kp] + frac * (intens[ie+1,kp] - intens[ie,kp]), t[kp]


def get_fermi_data(instrument, freq, chopper):
    if instrument.lower() == 'maps':
        if chopper.lower().startswith('s'):
            fermi_data = {'type':'sloppy', 'distance':NXfield(-1.857, **MU_),
                          'rotation_speed':NXfield(freq, units='hertz'),
                          'radius':NXfield(0.049, **MU_), 'r_slit':NXfield(1.3, **MU_),
                          'slit':NXfield(0.002899, **MU_), 'number':9, 'width':NXfield(0.0522, **MU_)}
        elif chopper.lower().startswith('a'):
            pass
        else:
            raise RuntimeError(f'Unrecognised chopper "{chopper}" for instrument "{instrument}"')
    elif instrument.lower() == 'merlin':
        if chopper.lower().startswith('s'):
            pass
        elif chopper.lower().startswith('g'):
            fermi_data = {'type':'g', 'distance':NXfield(-1.8, **MU_),
                          'rotation_speed':NXfield(freq, units='hertz'),
                          'radius':NXfield(0.005, **MU_), 'r_slit':NXfield(99999., **MU_),
                          'slit':NXfield(2.0e-4, **MU_), 'number':125, 'width':NXfield(0.05, **MU_)}
        else:
            raise RuntimeError(f'Unrecognised chopper "{chopper}" for instrument "{instrument}"')
    else:
        raise RuntimeError(f'Unrecognised instrument "{instrument}"')
    return fermi_data


def let_instrument(ei, freq=None):
    if freq is None:
        freq = [40., 240.]

    inst = NXinstrument(fermi=NXfermi_chopper(energy=ei))
    inst['name'] = NXfield(value='LET', short_name='LET')
    # Defines the Source
    inst['source'] = NXsource(Name='ISIS', type='Spallation Neutron Source',
                              frequency=NXfield(10, units='hertz'), target_material='W')
    # Defines the moderator
    d_mod = NXtransformations(MOD_T_AXIS=NXfield(-25., transformation_type='translation',
                                                 vector=[0.,0.,1.], depends_on='.', **MU_),
                              MOD_R_AXIS=NXfield(0., transformation_type='rotation',
                                                 vector=[0.,1.,0.], depends_on='MOD_T_AXIS', units='degree'))
    pulse = NXnote(type='ikcarp', data=[42.1304,0,0],
                   description='Empirical Ikeda-Carpenter type moderator pulse model')
    inst['moderator'] = NXmoderator(type='Liquid H2', temperature=NXfield(17., units='kelvin'),
                                    empirical_pulse_shape=pulse, transforms=d_mod)
    # Defines the divergences
    hdiv, vdiv = get_let_divergences(ei)
    inst['horiz_div'] = NXbeam(data=hdiv)
    inst['vert_div'] = NXbeam(data=vdiv)
    # Defines the choppers
    d_ch1 = NXtransformations(CH1_T_AXIS=NXfield(-17.17, transformation_type='translation',
                                                 vector=[0.,0.,1.], depends_on='.', **MU_))
    sl = np.rad2deg(np.arctan2(0.04, 0.28) / 2.)
    inst['shaping_chopper'] = NXdisk_chopper(rotation_speed=NXfield(freq[0], units='hertz'), radius=NXfield(0.28, **MU_),
                                             type='contra_rotating_pair', slit_edges=[-sl, sl], transforms=d_ch1)
    d_ch5 = NXtransformations(CH5_T_AXIS=NXfield(-1.5, transformation_type='translation',
                                                 vector=[0.,0.,1.], depends_on='.', **MU_))
    sl = np.rad2deg(np.arctan2(0.031, 0.28) / 2.)
    inst['mono_chopper'] = NXdisk_chopper(rotation_speed=NXfield(freq[1], units='hertz'), radius=NXfield(0.28, **MU_),
                                          type='contra_rotating_pair', slit_edges=[-sl, sl], transforms=d_ch5)
    return inst


def maps_instrument(ei, freq=None, chopper='S'):
    if freq is None:
        freq = 600.
    inst = NXinstrument(fermi=NXfermi_chopper(energy=ei))
    inst['name'] = NXfield(value='MAPS', short_name='MAPS')
    # Defines the Fermi chopper
    fermi = inst['fermi']
    for k, v in get_fermi_data('maps', freq, chopper).items():
        fermi[k] = v
    # Defines the Aperture
    d_ap = NXtransformations(AP_AXIS=NXfield(-10.3290, transformation_type='translation',
                                             vector=[0.,0.,1.], depends_on='.', **MU_))
    inst['aperture'] = NXslit(x_gap=NXfield(0.0989, **MU_), y_gap=NXfield(0.0989, **MU_),
                              transforms=d_ap)
    # Defines the Source
    inst['source'] = NXsource(Name='ISIS', type='Spallation Neutron Source',
                              frequency=NXfield(50, units='hertz'), target_material='W')
    # Defines the moderator
    d_mod = NXtransformations(MOD_T_AXIS=NXfield(-12., transformation_type='translation',
                                                 vector=[0.,0.,1.], depends_on='.', **MU_),
                              MOD_R_AXIS=NXfield(32., transformation_type='rotation',
                                                 vector=[0.,1.,0.], depends_on='MOD_T_AXIS', units='degree'))
    pulse_signal, pulse_tof = get_moderator_time_pulse('maps', ei)
    pulse = NXdata(signal=NXfield(pulse_signal, unit='1/microsecond/meV', name='Intensity'),
                   axes=NXfield(pulse_tof, unit='microsecond', name='Time'))
    inst['moderator'] = NXmoderator(type='H20', temperature=NXfield(300, units='kelvin'),
                                    pulse_shape=pulse, transforms=d_mod)
    return inst


def merlin_instrument(ei, freq=None, chopper='G'):
    if freq is None:
        freq = 600.
    inst = NXinstrument(fermi=NXfermi_chopper(energy=ei))
    inst['name'] = NXfield(value='MERLIN', short_name='MER')
    # Defines the Fermi chopper
    inst['name'].replace(NXfield(value='MERLIN', short_name='MER'))
    inst['fermi'].energy = ei
    fermi = inst['fermi']
    for k, v in get_fermi_data('merlin', freq, chopper).items():
        fermi[k] = v
    # Defines the Aperture
    d_ap = NXtransformations(AP_AXIS=NXfield(-10.1570, transformation_type='translation',
                                             vector=[0.,0.,1.], depends_on='.', **MU_))
    inst['aperture'] = NXslit(x_gap=NXfield(0.0967, **MU_), y_gap=NXfield(0.0967, **MU_),
                              transforms=d_ap)
    # Defines the Source
    inst['source'] = NXsource(Name='ISIS', type='Spallation Neutron Source',
                              frequency=NXfield(50, units='hertz'), target_material='W')
    # Defines the moderator
    d_mod = NXtransformations(MOD_T_AXIS=NXfield(-11.837, transformation_type='translation',
                                                 vector=[0.,0.,1.], depends_on='.', **MU_),
                              MOD_R_AXIS=NXfield(0., transformation_type='rotation',
                                                 vector=[0.,1.,0.], depends_on='MOD_T_AXIS', units='degree'))
    pulse_signal, pulse_tof = get_moderator_time_pulse('merlin', ei)
    pulse = NXdata(signal=NXfield(pulse_signal, unit='1/microsecond/meV', name='Intensity'),
                   axes=NXfield(pulse_tof, unit='microsecond', name='Time'))
    inst['moderator'] = NXmoderator(type='H20', temperature=NXfield(300, units='kelvin'),
                                    pulse_shape=pulse, transforms=d_mod)
    return inst


