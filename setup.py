import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="scRNAvariants",
    version="2.0.0",
    author="Refael Kohen, Yotam Constantini & Elisha Goldstein",
    author_email="elishagoldstein0308gmail.com",
    description="A script to help you locate cases of RNA modifications in single cell RNAseq data.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://bitbucket.org/eligol/scrarevar",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: POSIX",
    ],
    install_requires=["pysam>=0.15"],
    python_requires='>=3.7',
)