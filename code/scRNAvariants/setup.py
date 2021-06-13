
import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="scRNAvariants",
    version="2.0.0",
    author="Refael Kohen & Yotam Constantini",
    author_email="yotamcons@gmail.com",
    description="A script to help you locate cases of RNA modifications in single cell RNAseq data.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    #url="https://github.com/pypa/sampleproject",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    # Project uses pysam, so ensure that the docutils get installed
    install_requires=["pysam>=0.15"],
    python_requires='>=3.7',
)