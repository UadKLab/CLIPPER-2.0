from setuptools import find_packages, setup

setup(
    name="annotator",
    version="0.0.1",
    description="Peptide annotator for Proteome Discoverer \
                and SpectroMine peptide group output",
    packages=find_packages(),
    author="Konstantinos Kalogeropoulos",
    author_email="konskalogero@gmail.com",
    include_package_data=True,
    zip_safe=False,
    install_requires=[
        "statsmodels==0.11.1",
        "numpy==1.20.1",
        "Biopython==1.78",
        "pandas==1.2.4",
        "tqdm==4.60.0",
        "openpyxl==3.0.7",
        "scipy==1.3.1",
        'logomaker==0.8',
        'seaborn==0.11.0',
    ],
)
