from nexusformat.nexus import *
import nexusformat.nexus as nexus
import mcstasscript
import numpy as np
import warnings
import json
import copy
import sys
import os

cur_path = os.path.abspath(os.path.join(os.path.dirname(__file__)))
instr_path = os.path.join(cur_path, 'instruments')
comps_path = os.path.join(cur_path, 'mcstas-comps')


def get_instr(instrfile):
    instname = os.path.basename(instrfile).replace('.instr', '')
    inst = mcstasscript.McStas_instr(instname, package_path=comps_path)
    if not os.path.exists(instrfile):
        instrfile = os.path.join(instr_path, instrfile)
    reader = mcstasscript.McStas_file(instrfile)
    reader.add_to_instr(inst)
    return inst


def to_float(value):
    try:
        return np.array([float(v) for v in value])
    except ValueError:
        raise RuntimeError('Unsupported instr file: Position or rotation not numerical')


def _sanitize(indict):
    for k in indict:
        if isinstance(indict[k], str):
            indict[k] = indict[k].replace('<nl>', '\n').replace('<tb>', '    ')\
                .replace('<qt>',"'").replace('<bs>','\\')
        elif isinstance(indict[k], dict):
            indict[k] = _sanitize(indict[k])
    return indict


def dict2NXobj(indict):
    outdict = {}
    for k, v in indict.items():
        if isinstance(v, dict) and all([f in v for f in ['type', 'value']]) and v['type'].startswith('NX'):
            if v['type'] == 'NXfield':
                outdict[k] = NXfield(v['value'], **v.pop('attributes', {}))
            else:
                nxobj = getattr(nexus, v.pop('type'))
                outdict[k] = nxobj(**v.pop('value'))
        else:
            outdict[k] = v
    return outdict

# Only in nexusformat >= 1.0.0
if not hasattr(nexus, 'NXoff_geometry'):
    NXoff_geometry = nexus.tree._makeclass('NXoff_geometry')


class NXoff():
    # Class read/create NXoff_geometry fields using an Object File Format (OFF) syntax
    def __init__(self, vertices, faces):
        self.vertices = vertices
        self.faces = faces

    @classmethod
    def from_nexus(cls, nxfield):
        # Create an OFF structure from a NeXus field
        wo = nxfield.winding_order
        fa = list(nxfield.faces) + [len(wo)]
        faces = [wo[fa[ii]:fa[ii+1]] for ii in range(len(fa)-1)]
        return cls(nxfield.vertices, faces)

    @classmethod
    def from_wedge(cls, l, w1, h1, w2=None, h2=None):
        # Create an OFF structure in shape of a wedge (trapezoidal prism)
        if w2 is None:
            w2 = w1
        if h2 is None:
            h2 = h1
        (x1, y1, x2, y2) = tuple([float(v)/2 for v in [w1, h1, w2, h2]])
        # Sets the origin at the centre of the guide entry square face
        vertices = [[-x1, -y1, 0], [-x1, y1, 0], [x1, y1, 0], [x1, -y1, 0],
                    [-x2, -y2, l], [-x2, y2, l], [x2, y2, l], [x2, -y2, l]]
        # Use a clockwise winding, facing out, beam is in +z direction
        faces = [[0, 1, 2, 3], [1, 5, 6, 2], [5, 4, 7, 6],
                 [6, 7, 3, 2], [7, 4, 0, 3], [1, 0, 4, 5]]
        return cls(vertices, faces)

    def to_nexus(self):
        # Returns a NXoff_geometry field
        winding_order = [ww for fa in self.faces for ww in fa]
        faces = [0] + np.cumsum([len(self.faces[ii]) for ii in range(len(self.faces)-1)]).tolist()
        vertices = NXfield(np.array(self.vertices, dtype='float64'), units='metre')
        return NXoff_geometry(vertices=vertices, winding_order=winding_order, faces=faces)

    @staticmethod
    def _get_width_height(pos):
        # Gets the width and height (in x and y) of a set of points
        width = np.max(pos[:,0]) - np.min(pos[:,0])
        height = np.max(pos[:,1]) - np.min(pos[:,1])
        return width, height

    def get_guide_params(self):
        # Gets the guide parameters from the OFF geometry.
        ve = np.array(self.vertices)
        zmean = np.mean(ve[:,2])
        w1, h1 = self._get_width_height(ve[np.where(ve[:,2] < zmean)[0], :])
        w2, h2 = self._get_width_height(ve[np.where(ve[:,2] >= zmean)[0], :])
        l = np.max(ve[:,2]) - np.min(ve[:,2])
        return w1, h1, w2, h2, l


# Each entry here maps a NeXus component to a McStas component
# The second element is a mapping of NeXus component parameters to McStas parameters
# The third element is a mapping of NeXus component parameters to McStas position paramters
NX2COMP_MAP = dict(
    NXaperture = ['Slit',
        {'x_gap':'xwidth', 'y_gap':'yheight'},
    ],
    NXcollimator = ['Collimator_linear',
        {'divergence_x':'divergence', 'divergence_y':'divergenceV'},
    ],
    NXdetector = ['Monitor_nD',
        {},
    ],
    NXdisk_chopper = ['DiskChopper',
        {'slits':'nslit', 'rotation_speed':'nu', 'radius':'radius', 'slit_angle':'theta_0'},
    ],
    NXfermi_chopper = ['FermiChopper',
        {'rotation_speed':'nu', 'radius':'radius', 'slit':'w', 'r_slit':'curvature',
         'number':'nslit', 'width':'xwidth', 'height':'yheight'},
        {'distance': 'set_AT'},
    ],
    NXguide = ['Guide',
        {'m_value':'m'},
    ],
    NXslit = ['Slit',
        {'x_gap':'xwidth', 'y_gap':'yheight'},
    ],
    NXsample = ['Incoherent',
        {},
    ],
    NXmoderator = ['Moderator',
        {},
        {'distance': 'set_AT'},
    ]
)

class McStasComp2NX():
    # Class to map McStas components to NeXus components

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
        id_ext = mcstas_comp.EXTEND.find('eniius_data')
        if id_ext > 0:
            extstr = mcstas_comp.EXTEND[id_ext:]
            id0 = extstr.find('"') + 1
            id1 = extstr[id0:].find('"')
            jsonstr = extstr[id0:(id0+id1)].replace('\n', '').replace('\\"','').replace('\'','"').replace('\\','')
            try:
                extras = json.loads(jsonstr)
            except SyntaxError:
                pass
            else:
                for k, v in dict2NXobj(_sanitize(extras)).items():
                    self.nxobj[k] = v


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

    @classmethod
    def Guide(cls, **kw):
        geometry = NXoff.from_wedge(l=kw['l'], w1=kw['w1'], h1=kw['h1'], w2=kw.pop('w2', None), h2=kw.pop('h2', None))
        params = {'geometry':geometry.to_nexus()}
        if 'm' in kw:
            params['m_value'] = kw['m']
        return NXguide(**params)

    @classmethod
    def Guide_channeled(cls, **kw):
        return cls.Guide(**kw)

    @classmethod
    def Guide_gravity(cls, **kw):
        return cls.Guide(**kw)

    @classmethod
    def Guide_simple(cls, **kw):
        return cls.Guide(**kw)

    @classmethod
    def Guide_wavy(cls, **kw):
        return cls.Guide(**kw)

    @classmethod
    def Collimator_linear(cls, **kw):
        geometry = NXoff.from_wedge(l=kw['length'], w1=kw['xwidth'], h1=kw['yheight'])
        params = {'divergence_x':kw['divergence'], 'divergence_y':kw['divergenceV'], 'geometry':geometry.to_nexus()}
        return NXcollimator(**params)


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
        if hasattr(nxfield, 'depends_on'):
            depends_on = nxfield.depends_on
        transform = np.eye(4)
        offset = np.zeros(3)
        if hasattr(nxfield, 'offset'):
            offset = to_float(nxfield.offset)
        if nxfield.transformation_type == 'translation':
            transform[:3, 3] = to_float(nxfield.vector) * nxfield._value + offset
        elif nxfield.transformation_type == 'rotation':
            transform[:3, :3] = cls.rodrigues(to_float(nxfield.vector), nxfield._value)
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
        for order, comp in enumerate(self.components):
            nxinst[comp.name] = self.NXcomponent(comp.name, order)
        return nxinst


