from setuptools import setup

setup(
    name = "mTAGs",
    version = "1.0.3",
    author = "Hans-Joachim Ruscheweyh",
    author_email = "hansr@ethz.ch",
    description = ("mTAGs - taxonomic profiling using degenerate consensus reference sequences of ribosomal RNA gene."),
    license = "GPLv3",
    include_package_data=True,
    keywords = "bioinformatics metagenomics ngs 16S",
    url = "https://github.com/SushiLab/mTAGs",
    packages=['mTAGs'],
    download_url = "https://github.com/SushiLab/mTAGs/archive/refs/tags/1.0.3.tar.gz",
    entry_points = {
        'console_scripts': ['mtags=mTAGs.mtags:main'],
    }
)
