import subprocess
import shutil
import os
import setuptools
from setuptools.extension import Extension
from setuptools.command.build_ext import build_ext

VERSION = 'v1.0.0'

with open("README.md", "r") as fh:
    long_description = fh.read()

class compile_mican(build_ext):
    def build_extension(self, ext):
        # make
        subprocess.run(['make'])
        outpath = os.path.join(self.build_lib, 'gepred/bin/mican_2015.03.09')
        shutil.move('mican_2015.03.09', outpath)

setuptools.setup(
    name="gepred",
    version=VERSION,
    author="Shintaro Minami",
    description="GroE substrate prediction program",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ShintaroMinami/GEpred",
    cmdclass={
        'build_ext': compile_mican,
    },
    ext_modules=[Extension('', [])],
    packages=setuptools.find_packages(),
    include_package_data=True,
    install_requires=[
        'numpy',
        'tqdm'
    ],
    scripts=[
        'scripts/GEpred',
    ],

)
