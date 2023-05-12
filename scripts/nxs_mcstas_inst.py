import sys, os
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))

import eniius

from nexusformat.nexus import *


def nx2mcstas(nxfile):
    wrapper = eniius.Eniius.from_nxs(nxfile)
    wrapper.name = os.path.basename(nxfile).replace('.', '_')
    wrapper.to_mcstas().show_components()


if __name__ == '__main__':
    nx2mcstas(sys.argv[1])
