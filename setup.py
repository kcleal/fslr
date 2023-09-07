from setuptools import setup, find_packages


setup(
    name="fslr",
    version='0.2.5',
    license='MIT',
    python_requires='>=3.7',
    install_requires=[
            'numpy',
            'pandas',
            'click',
            'matplotlib',
            'networkx',
            'pafpy',
            'pysam',
        ],
    packages=find_packages(where="."),
    package_data={'': ['primers.csv']},
    include_package_data=True,
    entry_points='''
            [console_scripts]
            fslr=fslr.main:pipeline
        ''',
)
