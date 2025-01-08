import os
from setuptools import setup, find_packages



setup(
      name='functional_eval',
      version='0.0.1',
      description='functional evaluation pipeline for SALT workflow',
      packages=['physiochemical', 'docking', 'immunogenicity', 'docking.utils', 'immunogenicity.utils', 'immunogenicity.utils.bp3', 'immunogenicity.utils.bp3.BP3Models'],
      install_requires=[
            'tqdm',
            'numpy',
            'torch',
            'scipy',
            'wandb',
            'tensorboard',
            'pytorch_lightning'
      ],
)


