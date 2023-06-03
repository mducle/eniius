import sys, os
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))

import eniius

from nexusformat.nexus import *


def create_inst_nxs(outfile, inst_fun, pars, det_file=None):
    wrapper = eniius.Eniius(inst_fun(*pars), det_file)
    wrapper.to_nxspe(f'{outfile}.nxspe')
    wrapper.to_icp(f'{outfile}.nxs')


if __name__ == '__main__':
    detdir = os.path.join(os.path.dirname(__file__), '..', 'eniius', 'instruments')
    defdet = os.path.join(detdir, 'detector.dat')
    # These dat files can be found at: https://github.com/pace-neutrons/InstrumentFiles/
    for name, pars, det in zip(['let', 'maps', 'merlin'], 
                               [(3.7, [120., 240.]), (400., 600., 's'), (120., 600., 'g')],
                               ['det_LET_cycle222.dat', 'detector_202.dat', 'det_corr_203_process_5_v2.dat']):
        detfile = os.path.join(detdir, det)
        if not os.path.exists(detfile):
            detfile = defdet
        instfun = getattr(eniius.horace, f'{name}_instrument')
        create_inst_nxs(f'horace_{name}_inst', instfun, pars, detfile)
