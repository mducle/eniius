import sys, os
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))

import eniius.writer
import eniius.mcstas
from nexusformat.nexus import NXfield

def mcstas2nxs(instrfile):
    instobj = eniius.mcstas.get_instr(instrfile)
    nxinst = eniius.mcstas.NXMcStas(instobj.component_list).NXinstrument()
    nxinst['name'] = NXfield(value=instobj.name)
    writer = eniius.writer.Writer(nxinst)
    writer.to_json(f'mcstas_{instobj.name}.json')
    writer.to_icp(f'mcstas_{instobj.name}.nxs', 'detector.dat')

if __name__ == '__main__':
    mcstas2nxs(sys.argv[1])
