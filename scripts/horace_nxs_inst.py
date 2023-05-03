import sys, os
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))

import eniius.writer
import eniius.horace

from nexusformat.nexus import *


def create_inst_nxs(outfile, inst_fun, ei, det_file=None):
    writer = eniius.writer.Writer(inst_fun(ei))
    writer.sample = NXsample()
    writer.to_nxspe(f'{outfile}.nxspe', ei, det_file)
    writer.to_icp(f'{outfile}.nxs', det_file)


if __name__ == '__main__':
    create_inst_nxs('horace_let_inst', eniius.horace.let_instrument, 3.7, 'detector.dat')
    create_inst_nxs('horace_maps_inst', eniius.horace.maps_instrument, 400., 'detector.dat')
    create_inst_nxs('horace_merlin_inst', eniius.horace.merlin_instrument, 120., 'detector.dat')
