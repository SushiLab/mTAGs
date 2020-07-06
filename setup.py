from setuptools import setup


setup(
    name = "mTAGs",
    version = "0.9",
    author = "Hans-Joachim Ruscheweyh",
    author_email = "hansr@ethz.ch",
    description = ("Bioinformatic toolkit for extraction and taxonomic annotation of rRNA reads from metagenomic sequencing data."),
    license = "GPLv3",
    keywords = "bioinformatics metagenomics ngs 16S",
    url = "https://github.com/SushiLab/mTAGs",
    packages=['mTAGs'],
    entry_points = {
        'console_scripts': ['mtags=mTAGs.mtags:main'],
    }
)
