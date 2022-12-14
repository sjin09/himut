# -*- coding: utf-8 -*-

from setuptools import setup
setup(
        name="himut", 
        version='1.0.0',
        project_urls={
            "homepage": "https://github.com/sjin09/himut",
            "repository": "https://github.com/sjin09/himut"
        },
        author='Sangjin Lee',
        author_email='sl17@sanger.ac.uk',
        license='MIT',
        classifiers=[
            'Programming Language :: Python :: 3.6',
            'Programming Language :: Python :: 3.7',
            'Programming Language :: Python :: 3.8',
            'Programming Language :: Python :: 3.9'
        ],
        entry_points={"console_scripts": ["himut = himut.__main__:main"]},
        packages=['himut'],
        package_dir={"": "src"},
        package_data={"himut": ["*.typed"]},
        install_requires=[
            'argparse==1.*,>=1.4.0', 'biopython==1.*,>=1.78.0',
            'click==7.*,>=7.0.0', 'natsort==8.*,>=8.0.0', 'numpy==1.*,>=1.20.2',
            'psutil==5.*,>=5.8.0', 'pysam==0.*,>=0.16.0'
        ],
)
