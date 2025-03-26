from setuptools import setup, find_packages
import os

setup(
  name = "dna-aligner",
  version = "0.1.0",
  author = "Elliott Ou",
  author_email = "elliottou@gmail.com",
  description = "A command-line tool for DNA sequence alignment using global, local, and overlap alignment with affine gap penalties.",
  long_description = open("README.md").read(),
  long_description_content_type = "text/markdown",
  url="https://github.com/ou-elot/DNA_Aligner",
  packages = find_packages(),
  entry_points={
        "console_scripts": [
            "dna-aligner = dna_aligner.cli:main",
        ],
    },
)
  
  
