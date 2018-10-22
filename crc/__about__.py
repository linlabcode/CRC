"""Central place for package metadata."""

# NOTE: We use __title__ instead of simply __name__ since the latter would
#       interfere with a global variable __name__ denoting object's name.
__title__ = 'crc'
__summary__ = 'Core transcriptional regulatory circuitry analysis'
__url__ = 'https://github.com/linlabcode/CRC'
__download_url__ = 'https://github.com/linlabcode/CRC/tarball/v1.0.0'

# Semantic versioning is used. For more information see:
# https://packaging.python.org/en/latest/distributing/#semantic-versioning-preferred
__version__ = '2.0.0a1'

__maintainer__ = 'Jost Vrabic Koren'

__email__ = 'jost.koren@bcm.edu'

__license__ = 'MIT'

__all__ = (
    '__title__', '__summary__', '__url__', '__version__', '__maintainer__',  '__license__',
    '__email__',
)
