from setuptools import setup


setup(
    name = "mTAGs",
    version = "1.0.0",
    author = "Hans-Joachim Ruscheweyh",
    author_email = "hansr@ethz.ch",
    description = ("mTAGs - taxonomic profiling using degenerate consensus reference sequences of ribosomal RNA gene."),
    license = "GPLv3",
    keywords = "bioinformatics metagenomics ngs 16S",
    url = "https://github.com/SushiLab/mTAGs",
    packages=['mTAGs'],
    entry_points = {
        'console_scripts': ['mtags=mTAGs.mtags:main'],
    }
)
