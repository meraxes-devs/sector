import os


__all__ = [
    'RAITER10_logA500_001',
    'RAITER10_logE500_001',
    'RAITER10_Sal500_001',
    'RAITER10_Sal500_050',
    'STARBURST99_Kroupa',
    'STARBURST99_Salpeter',
]


packageDir = os.path.dirname(os.path.abspath(__file__))
RAITER10_logA500_001 = os.path.join(packageDir, 'logA500_001.hdf5')
RAITER10_logE500_001 = os.path.join(packageDir, 'logE500_001.hdf5')
RAITER10_Sal500_001 = os.path.join(packageDir, 'Sal500_001.hdf5')
RAITER10_Sal500_050 = os.path.join(packageDir, 'Sal500_050.hdf5')
STARBURST99_Kroupa = os.path.join(packageDir, 'STARBURST99_Kroupa.hdf5')
STARBURST99_Salpeter = os.path.join(packageDir, 'STARBURST99_Salpeter.hdf5')
