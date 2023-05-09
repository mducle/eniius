import sys, os
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))

import eniius


def mcstas2nxs(instrfile):
    detfile = os.path.join(os.path.dirname(eniius.__file__), 'detector.dat')
    wrapper = eniius.Eniius.from_mcstas(instrfile, detfile)
    wrapper.to_json(f'mcstas_{wrapper.name}.json')
    wrapper.to_icp(f'mcstas_{wrapper.name}.nxs')


if __name__ == '__main__':
    mcstas2nxs(sys.argv[1])
