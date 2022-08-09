import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="scRNAvariants",
    version="2.0.0",
    author="Refael Kohen, Yotam Constantini & Elisha Goldstein",
    author_email="elishagoldstein0308@gmail.com",
    description="A script to help you locate cases of RNA variations in single cell RNAseq data.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://bitbucket.org/eligol/scrarevar",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: POSIX",
    ],
    install_requires=['pysam>=0.15', 'samtools', 'bedtools', 'bamtools', 'pandas', 'matplotlib', 'seaborn', 'HTSeq', 'matplotlib_venn'],
    python_requires='>=3.7',
)
# pip packages: HTseq, matplotlib-venn
# seaborn installs pandas
# should we install pysam pybedtools instead of installing bedtools samtools
# pybamtools is not suppute and is not compatibla with python>2.7