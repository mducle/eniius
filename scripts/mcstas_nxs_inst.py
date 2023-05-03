import sys, os
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))

import eniius


def mcstas2nxs(instrfile):
    wrapper = eniius.Eniius.from_mcstas(instrfile, 'detector.dat')
    wrapper.to_json(f'mcstas_{wrapper.name}.json')
    wrapper.to_icp(f'mcstas_{wrapper.name}.nxs')


if __name__ == '__main__':
    mcstas2nxs(sys.argv[1])
