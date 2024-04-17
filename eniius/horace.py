import numpy as np
import scipy.io
import os

from nexusformat.nexus import *
from .pychop.Instruments import Instrument

THISFOLDER = os.path.dirname(os.path.realpath(__file__))
MOD_TABLES = {}
LET_TABLES = None
MU_ = {'units':'metre'}
INST_FACES = {'maps':'TS1_S01_Maps.mcstas', 'merlin':'TS1_S04_Merlin.mcstas', 'let':'TS2.imat',
              'mari':'TS1_S06_Mari.mcstas'}


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
    modfile = os.path.join(THISFOLDER, 'mcstas-comps', 'data', 'ISIS_tables', INST_FACES[instrument])
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
    return intens[ie,kp] + frac * (intens[ie+1,kp] - intens[ie,kp]), (t[kp] + t[kp+1]) / 2.


def maps_flux_gain(ei):
    l = np.sqrt(81.804201263673718 / np.array(ei))
    if l > 4.5:
        gain = 11.9644283524285 + 3.3130483907728 * (l - 4.5)
    else:
        gain = 1 + np.exp(-2.31186757913804 / l) * (23.2965476737752 + \
                - 15.1289612048092 * l + 6.11185114079924 * l * l \
                - 0.863318992430536 * l * l + 0.0439354263500835 * l * l * l)
    return gain


def merlin_flux_gain(ei):
    l = np.sqrt(81.804201263673718 / np.array(ei))
    if l > 3.5:
        gain = 10.92390518545 + 5.0628835969606 * (l - 3.5)
    else:
        gain = 1 + 0.723584927461204 * l - 3.461019858758497 * l**2 \
                 + 6.870176815937414 * l**3 - 3.962250897938358 * l**4 \
                 + 0.960065940459538 * l**5 - 0.084008173502155 * l**6
    return gain


def mari_flux_gain(ei):
    l = np.sqrt(81.804201263673718 / np.array(ei))
    gain = 1 + 2.3607 * l + 0.0652 * l**2 + 0.043 * l**3
    return gain


def get_fermi_data(instrument, freq, chopper):
    if instrument.lower() == 'maps':
        if chopper.lower().startswith('s'):
            typ, rho, dslit = ('sloppy', 1.3, 0.002899)
        elif chopper.lower().startswith('a'):
            typ, rho, dslit = ('a', 1.3, 0.001087)
        elif chopper.lower().startswith('b'):
            typ, rho, dslit = ('b', 0.92, 0.001812)
        else:
            raise RuntimeError(f'Unrecognised chopper "{chopper}" for instrument "{instrument}"')
        fermi_data = {'type':typ, 'distance':NXfield(-1.857, **MU_),
                      'rotation_speed':NXfield(freq, units='hertz'),
                      'radius':NXfield(0.049, **MU_), 'r_slit':NXfield(rho, **MU_),
                      'slit':NXfield(dslit, **MU_), 'number':9, 'width':NXfield(0.0522, **MU_)}
    elif instrument.lower() == 'merlin':
        if chopper.lower().startswith('s'):
            typ, rho, dslit, rad = ('sloppy', 1.3, 0.00228, 0.049)
        elif chopper.lower().startswith('a'):
            typ, rho, dslit, rad = ('a', 1.3, 0.00228, 0.049)
        elif chopper.lower().startswith('b'):
            typ, rho, dslit, rad = ('b', 0.92, 0.00129, 0.049)
        elif chopper.lower().startswith('g'):
            typ, rho, dslit, rad = ('g', 99999., 2.0e-4, 0.005)
        else:
            raise RuntimeError(f'Unrecognised chopper "{chopper}" for instrument "{instrument}"')
        fermi_data = {'type':typ, 'distance':NXfield(-1.8, **MU_),
                      'rotation_speed':NXfield(freq, units='hertz'),
                      'radius':NXfield(rad, **MU_), 'r_slit':NXfield(rho, **MU_),
                      'slit':NXfield(dslit, **MU_), 'number':125, 'width':NXfield(0.05, **MU_)}
    elif instrument.lower() == 'mari':
        if chopper.lower().startswith('s'):
            typ, rho, dslit, rad = ('sloppy', 1.3, 0.00228, 0.049)
        elif chopper.lower().startswith('g'):
            typ, rho, dslit, rad = ('g', 0.8, 2.0e-4, 0.005)
        elif chopper.lower().startswith('a'):
            typ, rho, dslit, rad = ('a', 1.3, 0.00076, 0.049)
        else:
            raise RuntimeError(f'Unrecognised chopper "{chopper}" for instrument "{instrument}"')
        fermi_data = {'type':typ, 'distance':NXfield(-1.7, **MU_),
                      'rotation_speed':NXfield(freq, units='hertz'),
                      'radius':NXfield(rad, **MU_), 'r_slit':NXfield(rho, **MU_),
                      'slit':NXfield(dslit, **MU_), 'number':125, 'width':NXfield(0.05, **MU_)}
    else:
        raise RuntimeError(f'Unrecognised instrument "{instrument}"')
    return fermi_data


def let_instrument(ei, freq=None, variant=None):
    if freq is None:
        freq = [80., 240.]
    if variant is None:
        variant = 'High flux'

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
    inst['shaping_chopper'] = NXdisk_chopper(rotation_speed=NXfield(freq[0] / 2., units='hertz'), radius=NXfield(0.28, **MU_),
                                             type='contra_rotating_pair', slit_edges=[-sl, sl], transforms=d_ch1)
    d_ch5 = NXtransformations(CH5_T_AXIS=NXfield(-1.5, transformation_type='translation',
                                                 vector=[0.,0.,1.], depends_on='.', **MU_))
    slot_widths = {'High flux':0.031, 'Intermediate':0.02, 'High Resolution':0.015}
    sl = np.rad2deg(np.arctan2(slot_widths[variant], 0.28) / 2.)
    inst['mono_chopper'] = NXdisk_chopper(rotation_speed=NXfield(freq[1], units='hertz'), radius=NXfield(0.28, **MU_),
                                          type='contra_rotating_pair', slit_edges=[-sl, sl], transforms=d_ch5)
    # Gets the allowed energies in multirep mode and adds subsidiary reps as NXcollections
    pych = Instrument('LET', variant, freq[::-1])
    other_Eis = []
    for rep in zip(pych.getAllowedEi(ei), pych.getMultiRepFlux(ei)):
        if not np.isnan(rep[1][0]) and np.abs(rep[0] - ei) > 0.01:
            try:
                hdiv, vdiv = get_let_divergences(rep[0])
            except RuntimeError:
                pass
            inst[f'rep_{rep[0]}'] = NXcollection()
            inst[f'rep_{rep[0]}/horiz_div'] = NXbeam(data=hdiv)
            inst[f'rep_{rep[0]}/vert_div'] = NXbeam(data=vdiv)
    return inst


def maps_single(ei):
    inst = {}
    # Defines the Aperture
    d_ap = NXtransformations(AP_AXIS=NXfield(-10.3290, transformation_type='translation',
                                             vector=[0.,0.,1.], depends_on='.', **MU_))
    fac = np.sqrt(maps_flux_gain(ei))
    inst['aperture'] = NXslit(x_gap=NXfield(0.094*fac, **MU_), y_gap=NXfield(0.094*fac, **MU_),
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


def maps_instrument(ei, freq=None, chopper='S', rrm=True):
    if freq is None:
        freq = 600.
    inst = NXinstrument(fermi=NXfermi_chopper(energy=ei))
    inst['name'] = NXfield(value='MAPS', short_name='MAPS')
    # Defines the Fermi chopper
    fermi = inst['fermi']
    for k, v in get_fermi_data('maps', freq, chopper).items():
        fermi[k] = v
    for k, v in maps_single(ei).items():
        inst[k] = v
    # Gets the allowed energies in multirep mode and adds subsidiary reps as NXcollections
    if rrm:
        pych = Instrument('MAPS', chopper, freq)
        other_Eis = []
        for rep in zip(pych.getAllowedEi(ei), pych.getMultiRepFlux(ei)):
            if not np.isnan(rep[1][0]) and np.abs(rep[0] - ei) > 0.01:
                inst[f'rep_{rep[0]}'] = NXcollection()
                for k, v in merlin_single(rep[0]).items():
                    inst[f'rep_{rep[0]}/{k}'] = v
    return inst


def merlin_single(ei):
    inst = {}
    # Defines the Aperture
    d_ap = NXtransformations(AP_AXIS=NXfield(-10.1570, transformation_type='translation',
                                             vector=[0.,0.,1.], depends_on='.', **MU_))
    fac = np.sqrt(merlin_flux_gain(ei))
    inst['aperture'] = NXslit(x_gap=NXfield(0.094*fac, **MU_), y_gap=NXfield(0.094*fac, **MU_),
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


def merlin_instrument(ei, freq=None, chopper='G', rrm=True):
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
    for k, v in merlin_single(ei).items():
        inst[k] = v
    # Gets the allowed energies in multirep mode and adds subsidiary reps as NXcollections
    if rrm:
        pych = Instrument('MERLIN', chopper, freq)
        other_Eis = []
        for rep in zip(pych.getAllowedEi(ei), pych.getMultiRepFlux(ei)):
            if not np.isnan(rep[1][0]) and np.abs(rep[0] - ei) > 0.01:
                inst[f'rep_{rep[0]}'] = NXcollection()
                for k, v in merlin_single(rep[0]).items():
                    inst[f'rep_{rep[0]}/{k}'] = v
    return inst


def mari_single(ei):
    inst = {}
    # Defines the Aperture
    d_ap = NXtransformations(AP_AXIS=NXfield(-10.04, transformation_type='translation',
                                             vector=[0.,0.,1.], depends_on='.', **MU_))
    fac = np.sqrt(mari_flux_gain(ei))
    inst['aperture'] = NXslit(x_gap=NXfield(0.09*fac, **MU_), y_gap=NXfield(0.09*fac, **MU_),
                              transforms=d_ap)
    # Defines the Source
    inst['source'] = NXsource(Name='ISIS', type='Spallation Neutron Source',
                              frequency=NXfield(50, units='hertz'), target_material='W')
    # Defines the moderator
    d_mod = NXtransformations(MOD_T_AXIS=NXfield(-11.7, transformation_type='translation',
                                                 vector=[0.,0.,1.], depends_on='.', **MU_),
                              MOD_R_AXIS=NXfield(0., transformation_type='rotation',
                                                 vector=[0.,1.,0.], depends_on='MOD_T_AXIS', units='degree'))
    pulse_signal, pulse_tof = get_moderator_time_pulse('mari', ei)
    pulse = NXdata(signal=NXfield(pulse_signal, unit='1/microsecond/meV', name='Intensity'),
                   axes=NXfield(pulse_tof, unit='microsecond', name='Time'))
    inst['moderator'] = NXmoderator(type='H20', temperature=NXfield(300, units='kelvin'),
                                    pulse_shape=pulse, transforms=d_mod)
    return inst


def mari_instrument(ei, freq=None, chopper='G', rrm=None):
    if freq is None:
        freq = 600.
    inst = NXinstrument(fermi=NXfermi_chopper(energy=ei))
    inst['name'] = NXfield(value='MARI', short_name='MAR')
    # Defines the Fermi chopper
    inst['name'].replace(NXfield(value='MARI', short_name='MAR'))
    inst['fermi'].energy = ei
    fermi = inst['fermi']
    for k, v in get_fermi_data('mari', freq, chopper).items():
        fermi[k] = v
    for k, v in mari_single(ei).items():
        inst[k] = v
    # Gets the allowed energies in multirep mode and adds subsidiary reps as NXcollections
    if rrm is not None:
        pych = Instrument('MARI', chopper, freq)
        other_Eis = []
        for rep in zip(pych.getAllowedEi(ei), pych.getMultiRepFlux(ei)):
            if not np.isnan(rep[1][0]) and np.abs(rep[0] - ei) > 0.01:
                inst[f'rep_{rep[0]}'] = NXcollection()
                for k, v in mari_single(rep[0]).items():
                    inst[f'rep_{rep[0]}/{k}'] = v
    return inst
