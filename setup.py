from setuptools import setup, find_packages

setup(
      name='func_eval',
      version='0.0.1',
      description='functional evaluation pipeline for SALT workflow',
      packages=['physiochemical'],
      install_requires=[
            'tqdm',
            'numpy',
            'torch',
            'scipy',
            'wandb',
            'tensorboard',
            'pytorch_lightning',
      ],
)