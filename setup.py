import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="SCRvar",
    version="2.0.6",
    author="Yotam Constantini, Elisha Goldstein, Refael Kohen & Dena Leshkowitz",
    author_email="elishagoldstein0308@gmail.com",
    description="A script to detect RNA variations in single cell RNA-seq data.",
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
    scripts=['scripts/process_aditional_files_step0.py', 'scripts/run_SCRvar_steps_1234.py',
             'scripts/run_filtering_and_analysis_step5.py', 'scripts/run_gene_level_analysis_step6.py']
)

# not exist on test.pypi: pandas, matplotlib-venn
# not exist on pypi: 'samtools', 'bedtools', 'bamtools',
