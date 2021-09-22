from setuptools import setup
import os

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()
install_requires = ['biopython']
long_description = read('README.md')
setup(
    name = "mTAGs",
    version = "1.0.4",
    author = "Hans-Joachim Ruscheweyh",
    author_email = "hansr@ethz.ch",
    description = ("mTAGs - taxonomic profiling using degenerate consensus reference sequences of ribosomal RNA gene."),
    license = "GPLv3",
    include_package_data=True,
    install_requires=install_requires,
    long_description=long_description,
    long_description_content_type='text/markdown',
    keywords = "bioinformatics metagenomics ngs 16S",
    url = "https://github.com/SushiLab/mTAGs",
    packages=['mTAGs'],
    download_url = "https://github.com/SushiLab/mTAGs/archive/refs/tags/1.0.4.tar.gz",
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: OS Independent',
    ],
    entry_points = {
        'console_scripts': ['mtags=mTAGs.mtags:main'],
    }
)
