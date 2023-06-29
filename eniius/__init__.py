import warnings
try:
    from . import mcstas
    from . import horace
    from . import nexus
    from .eniius import Eniius
except ModuleNotFoundError as e:
    import traceback
    warnings.warn('Could not import submodule')
    traceback.print_tb(e.__traceback__)
    print(e)

from . import _version
__version__ = _version.get_versions()['version']

from .pychop import Instruments as PyChop
