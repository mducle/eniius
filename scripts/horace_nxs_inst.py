import sys, os
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))

import eniius

from nexusformat.nexus import *


def create_inst_nxs(outfile, inst_fun, ei, det_file=None):
    wrapper = eniius.Eniius(inst_fun(ei), det_file)
    wrapper.to_nxspe(f'{outfile}.nxspe')
    wrapper.to_icp(f'{outfile}.nxs')


if __name__ == '__main__':
    create_inst_nxs('horace_let_inst', eniius.horace.let_instrument, 3.7, 'detector.dat')
    create_inst_nxs('horace_maps_inst', eniius.horace.maps_instrument, 400., 'detector.dat')
    create_inst_nxs('horace_merlin_inst', eniius.horace.merlin_instrument, 120., 'detector.dat')
