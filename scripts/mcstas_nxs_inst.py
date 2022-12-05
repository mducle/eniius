from nexusformat.nexus import *
import nexusformat.nexus as nexus
import mcstasscript
import numpy as np
import warnings
import json
import copy
import sys
import os

cur_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
instr_path = os.path.join(cur_path, 'instruments')
comps_path = os.path.join(instr_path, 'mcstas-comps')
sys.path.append(cur_path)
from create_nxs_inst import create_inst_nxs, create_data

def get_instr(instrfile):
    instname = instrfile.replace('.instr', '')
    inst = mcstasscript.McStas_instr(instname, package_path=comps_path)
    reader = mcstasscript.McStas_file(os.path.join(instr_path, instrfile))
    reader.add_to_instr(inst)
    return inst


def to_float(value):
    try:
        return np.array([float(v) for v in value])
    except ValueError:
        raise RuntimeError('Unsupported instr file: Position or rotation not numerical')


def off_geom(**kwargs):
    return None

# Each entry here maps a NeXus component to a McStas component
# The second element is a mapping of NeXus component parameters to McStas parameters
# The third element (if present) is a functor to return the geommetry from the NeXus OFF_GEOMETRY
# The forth element (if present) is a list of NeXus parameters which invalidates the McStas comp
NX2COMP_MAP = dict(
    NXaperture = ['Slit', 
        {'x_gap':'xwidth', 'y_gap':'yheight'},
    ],
    NXcollimator = ['Collimator_linear',
        {'divergence_x':'divergence', 'divergence_y':'divergenceV'},
        off_geom(width='xwidth', height='yheight', length='length', shape='cuboid'),
    ],
    NXdetector = ['Monitor_nD',
        {},
        off_geom(off_file='geometry'),
    ],
    NXdisk_chopper = ['DiskChopper',
        {'slits':'nslit', 'rotation_speed':'nu', 'radius':'radius', 'slit_angle':'theta_0'},
        {},
        [],
        ['slit_edges'],
    ],
    NXfermi_chopper = ['FermiChopper',
        {'rotation_speed':'nu', 'radius':'radius', 'slit':'w', 'r_slit':'curvature',
         'number':'nslit', 'width':'xwidth', 'height':'yheight'},
    ],
    NXguide = ['Guide',
        {'m_value':'m'},
        off_geom(in_width='w1', in_height='h1', out_width='w2', out_height='h2',
                 length='l', shape='trapezoid'),
    ],
    NXslit = ['Slit',
        {'x_gap':'xwidth', 'y_gap':'yheight'},
    ],
    NXsample = ['Incoherent'],
)

class McStasComp2NX():

    COMPGP2NX_MAP = {
        'Guide*':'NXguide',
        'Collimator*':'NXcollimator',
    }

    COMPCAT2NX_MAP = dict(
        sources = 'NXmoderator',
        monitors = 'NXdetector'
    ) 

    COMP2NX_MAP = dict(
        DiskChopper = 'NXdisk_chopper',
        FermiChopper = 'NXfermi_chopper',
        FermiChopper_ILL = 'NXfermi_chopper',
        Fermi_chop2a = 'NXfermi_chopper',
        Filter_gen = 'NXfilter',
        Filter_graphite = 'NXfilter',
        Elliptic_guide_gravity = 'NXguide',
        Mirror = 'NXmirror',
        Monochromator_flat = 'NXmonochromator',
        Monochromator_curved = 'NXmonochromator',
        Monochromator_pol = 'NXpolarizer',
        Pol_SF_ideal = 'NXflipper',
        Pol_bender = 'NXpolarizer',
        Pol_mirror = 'NXpolarizer',
        SNS_source = 'NXmoderator',
        SNS_source_analytic = 'NXmoderator',
        Source_pulsed = 'NXmoderator',
        Selector = 'NXvelocity_selector',
        V_selector = 'NXvelocity_selector',
        ViewModISIS = 'NXmoderator',
    )

    def __init__(self, mcstas_comp, mcstas_order, transforms, **kwargs):
        mcstas_name = mcstas_comp.component_name
        try:
            self.nxobj = getattr(self, mcstas_name)(**kwargs)
        except AttributeError:
            nxtype = self.getNXtype(mcstas_comp)
            ctor = getattr(nexus, nxtype)
            params = {}
            if nxtype in NX2COMP_MAP:
                for nxpar, mcstaspar in NX2COMP_MAP[nxtype][1].items():
                    params[nxpar] = getattr(mcstas_comp, mcstaspar)
            self.nxobj = ctor(**params)
        kwargs['mcstas_component'] = mcstas_comp.component_name
        kwargs['mcstas_order'] = mcstas_order
        self.nxobj['mcstas'] = json.dumps(kwargs)
        self.nxobj['transforms'] = transforms

    @classmethod
    def getNXtype(cls, comp):
        try:
            return cls.COMP2NX_MAP[comp.component_name]
        except KeyError:
            pass
        try:
            return cls.COMPCAT2NX_MAP[comp.category]
        except KeyError:
            pass
        nxtype = [v for k,v in cls.COMPGP2NX_MAP.items() if comp.component_name.startswith(k.split('*')[0])]
        if len(nxtype) == 0:
            return 'NXnote'
        else:
            return nxtype[0]

    @classmethod
    def Slit(cls, **kw):
        params = {}
        params['x_gap'] = kw['xwidth'] if kw['xwidth'] else float(kw['xmax']) - float(kw['xmin'])
        params['y_gap'] = kw['yheight'] if kw['yheight'] else float(kw['ymax']) - float(kw['ymin'])
        return NXslit(**params)

    @classmethod
    def Diaphragm(cls, **kw):
        return cls.Slit(**kw)


class AffineRotate():
    def __init__(self, transformation_matrix, depends_on='.'):
        # Class for a single instance of a transformation (translation + rotation) supported by McStas/NeXus
        self.transform = transformation_matrix
        self.depends_on = depends_on

    @classmethod
    def from_euler_translation(cls, euler_angles, translation_vector, depends_on='.'):
        # Constructor from McStas data using Euler angles and a translation
        transform = np.eye(4)
        transform[:3, :3] = cls.rotmat(euler_angles)
        transform[:3, 3] = translation_vector
        return cls(transform, depends_on)

    @classmethod
    def from_nxfield(cls, nxfield, depends_on='.'):
        # Constructor from NeXus data from a single transformation field
        # The field must have the 'transformation_type', and 'vector' attributes (optionally 'offset')
        assert hasattr(nxfield, 'transformation_type'), "NXfield must have 'transformation_type' attribute"
        assert hasattr(nxfield, 'vector'), "NXfield must have 'vector' attribute"
        transform = np.eye(4)
        offset = np.zeros(3)
        if hasattr(nxfield, 'offset'):
            offset = to_float(nxfield.offset)
        if nxfield.transformation_type == 'translation':
            transform[:3, 3] = to_float(nxfield.vector) * nxfield._value + offset
        elif nxfield.transformation_type == 'rotation':
            transform[:3, :3] = cls.rodrigues(to_float(nx_field.vector), nxfield._value)
            transform[:3, 3] = offset
        else:
            raise RuntimeError('transformation_type must be either "translation" or "rotation"')
        return cls(transform, depends_on)
        
    def axisrot(self):
        # Computes the net rotation axis and rotation angle from the rotation matrix
        # https://en.wikipedia.org/wiki/Rotation_matrix#Conversion_from_rotation_matrix_to_axis%E2%80%93angle
        dd, vv = np.linalg.eig(self.transform[:3, :3])
        idx = np.where(np.abs(dd - 1) < 1.e-6)[0]
        assert len(idx) > 0, "Error: Input transformation is not a chained rotation"
        axis = vv[:, idx[0]]
        rotmat = self.transform[:3, :3]
        # Sign of rotation angle is not defined by the expression below (from matrix trace)
        angle = np.degrees(np.arccos((np.trace(rotmat) - 1.) / 2.))
        # Checks that the computed axis and angles agree with the original matrix
        r2 = self.rodrigues(axis, angle)
        if np.sum(np.abs(rotmat - r2)) > 1.e-5:
            angle = -angle
            r2 = self.rodrigues(axis, angle)
        assert np.sum(np.abs(rotmat - r2)) < 1.e-5, "Error computing the rotation axes and angles"
        return axis, angle

    @staticmethod
    def rotmat(euler):
        # Computes the Euler rotation matrix in the (improper) XYZ convention used by McStas
        # where a positive angle corresponds to a counter-clockwise rotation
        # https://github.com/McStasMcXtrace/McCode/blob/master/common/lib/share/mccode-r.c#L2521
        cc = np.cos(np.radians(euler))
        ss = np.sin(np.radians(euler))
        return np.array([[cc[1]*cc[2], ss[0]*ss[1]*cc[2]+cc[0]*ss[2], ss[0]*ss[2]]-cc[0]*ss[1]*cc[2],
                         [-cc[1]*ss[2], cc[0]*cc[2]-ss[0]*ss[1]*ss[2], ss[0]*cc[2]]+cc[0]*ss[1]*ss[2],
                         [ss[1],       -ss[0]*cc[1],                   cc[0]*cc[1]]])

    @staticmethod
    def get_euler_angles(rotmat):
        # Calculates the euler angles from a rotation matrix under the McStas convention
        assert (np.abs(np.linalg.det(rotmat)) - 1) < 1.e-5, "Error: transformation is not valid"
        return np.array([np.degrees(np.arctan2(-rotmat[2,1], rotmat[2,2])), 
                         np.degrees(np.arctan2(rotmat[2,0], np.sqrt(1 - rotmat[2,0]**2))),
                         np.degrees(np.arctan2(-rotmat[1,0], rotmat[0,0]))])

    @staticmethod
    def rodrigues(axis, angle):
        # Computes the rotation matrix from an axis and an angle using Rodrigues' formula
        kperp = np.array([[0., -axis[2], axis[1]], [axis[2], 0., -axis[0]], [-axis[1], axis[0], 0.]])
        a = np.radians(angle)
        kvec = np.array(axis)[np.newaxis,:]
        # Using the vector version of Rodrigues' rotation formula
        return np.eye(3)*np.cos(a) + (1 - np.cos(a))*(kvec.T * kvec) + np.sin(a)*kperp

    def reverse(self, depends_on=None):
        # Computes the reverse transformation
        transform = np.eye(4)
        transform[:3, :3] = np.transpose(self.transform[:3, :3])
        transform[:3, 3] = -self.transform[:3, 3]
        if depends_on is None:
            depends_on = self.depends_on
        return AffineRotate(transformation_matrix=transform, depends_on=depends_on)

    @property
    def is_translation(self):
        return np.sum(np.abs(self.transform[:3, :3] - np.eye(3))) < 1.e-5

    @property
    def is_rotation(self):
        return not self.is_translation

    def NXfield(self):
        # Returns an NXfield object for this component
        assert np.abs(np.imag(self.transform)).sum() < 1e-5, "Error computing transformation vector"
        self.transform = np.real(self.transform)
        if self.is_translation:
            distance = np.linalg.norm(self.transform[:3, 3])
            vec = copy.deepcopy(self.transform[:3, 3])
            if distance > 1e-5:
                vec /= distance
            return NXfield(distance, vector=vec, depends_on=self.depends_on,
                           transformation_type='translation', units='metre')
        else:
            # If transformation has a rotation, include translation as offset
            axis, angle = self.axisrot()
            assert np.abs(np.imag(axis)).sum() < 1e-5, "Error computing rotation from transform"
            axis = np.real(axis)
            return NXfield(angle, vector=axis, offset=self.transform[:3, 3], depends_on=self.depends_on,
                           transformation_type='rotation', units='degree')


class NXMcStas():
    # Class to convert a McStas instrument definition embodied by a list of components to a NeXus file

    def __init__(self, components_list):
        self.transforms = {}
        self.depends_on = {}
        self.components = []
        self.affinelist = {}
        self.indices = {}
        for ii, comp in enumerate(components_list):
            relate_at = comp.AT_relative.replace('RELATIVE ', '')
            relate_rot = comp.ROTATED_relative.replace('RELATIVE ', '')
            if (relate_at != relate_rot) and (relate_rot != 'ABSOLUTE') and (relate_at != 'ABSOLUTE'):
                raise RuntimeError('Rotation and position relative to different components not supported')
            print(comp.ROTATED_data)
            print(comp.AT_data)
            self.transforms[comp.name] = AffineRotate.from_euler_translation(to_float(comp.ROTATED_data),
                                                                             to_float(comp.AT_data),
                                                                             depends_on=relate_at)
            self.depends_on[comp.name] = relate_at
            self.components.append(comp)
            self.indices[comp.name] = ii
        # Map out the RELATIVE chain of transformations and save as list
        for comp in self.components:
            node = comp.name
            self.affinelist[comp.name] = [self.transforms[node]]
            while self.depends_on[node] != 'ABSOLUTE':
                node = self.depends_on[node]
                self.affinelist[comp.name].append(self.transforms[node])
        # Horace and Mantid sets the origin at the sample position. 
        # For compatibility, we define NeXus files with the origin there if possible
        samp = [comp for comp in components_list if comp.category == 'samples']
        if len(samp) == 0:
            warnings.warn("Instrument does not have a sample. Will use the McStas "
                          "ABSOLUTE positions", warnings.RuntimeWarning)
            self.origin, rev_trans = ('', [])
        else:
            if len(samp) > 1:
                warnings.warn("More than one sample in instrument. Will use the first "
                              "sample position as the origin.", warnings.RuntimeWarning)
            self.origin = samp[0].name
            deps =  [tr.depends_on for tr in self.affinelist[self.origin][::-1]][1:] + ['.']
            rev_trans = [tr.reverse(depends_on=deps[ii]) for ii, tr in enumerate(self.affinelist[self.origin][::-1])]
        for name in [comp.name for comp in self.components]:
            if name == self.origin:
                self.affinelist[name] = [AffineRotate.from_euler_translation([0, 0, 0], [0, 0, 0])]
                continue
            self.affinelist[name] = self._reduce_transforms(self.affinelist[name] + rev_trans)

    def _reduce_transforms(self, affinelist):
        # Checks if successive transformations could be concatenated
        if len(affinelist) < 2:
            return affinelist
        new_list, transform_added = ([], False)
        tr1 = affinelist[0]
        for ii in range(1, len(affinelist)):
            tr2 = affinelist[ii]
            mat = np.matmul(tr1.transform, tr2.transform)
            try:
                tr1 = AffineRotate(transformation_matrix=mat, depends_on=tr2.depends_on)
                transform_added = False
            except AssertionError:
                new_list.append(tr1)
                tr1 = tr2
                transform_added = True
        if not transform_added:
            new_list.append(tr1)
        return new_list

    def NXtransformations(self, name):
        # Returns an NXtransformations group for a component with a name
        transdict = {}
        for idx, trans in enumerate(self.affinelist[name]):
            transdict[f'{name}{idx}'] = trans.NXfield()
        return NXtransformations(**transdict)

    def NXcomponent(self, name, order=0):
        # Returns a NXcomponent corresponding to a McStas component.
        comp = self.components[self.indices[name]]
        mcpars = {p:getattr(comp, p) for p in comp.parameter_names}
        return McStasComp2NX(comp, order, self.NXtransformations(name), **mcpars).nxobj

    def NXinstrument(self):
        nxinst = NXinstrument()
        for comp in self.components:
            nxinst[comp.name] = self.NXcomponent(comp.name)
        return nxinst


def mcstas2nxs(instrfile):
    instobj = get_instr(instrfile)
    nxinst = NXMcStas(instobj.component_list).NXinstrument()
    nxinst['name'] = NXfield(value=instobj.name)
    print(nxinst.tree)
    create_inst_nxs(f'mcstas_{instrfile}.nxs', lambda ei: nxinst, 25)


if __name__ == '__main__':
    mcstas2nxs(sys.argv[1])
