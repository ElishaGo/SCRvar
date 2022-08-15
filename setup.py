import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="scRNAvariants",
    version="2.0.9",
    author="Refael Kohen, Yotam Constantini & Elisha Goldstein",
    author_email="elishagoldstein0308@gmail.com",
    description="A script to help you locate cases of RNA variations in single cell RNAseq data.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://bitbucket.org/eligol/scrarevar",
    packages=setuptools.find_packages(),
    package_data={"": ["*.sh"]},
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: POSIX",
    ],
    install_requires=['pysam>=0.15', 'pandas', 'seaborn', 'HTSeq', 'matplotlib-venn'],
    python_requires='>=3.7',
    scripts=['sc_rna_variants/helper_scripts/run_steps_1234.py', 'scripts/step0_process_editing_and_snp_DB.py', 'scripts/step5_filtering_and_analysis.py']
)

# not exist on test.pypi: pandas, matplotlib-venn
# not exist on pypi: 'samtools', 'bedtools', 'bamtools',