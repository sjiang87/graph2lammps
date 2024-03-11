from setuptools import setup, find_packages

setup(
    name='graph2lammps',
    version='0.0.1',
    description='NetworkX graphs to LAMMPS data',
    author='Shengli (Bruce) Jiang',
    author_email='sj0161@princeton.com',
    url='https://github.com/sjiang87/graph2lammps',
    packages=find_packages(),
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],
    python_requires='>=3.6',
    install_requires=[
        'numpy>=1.0',
        'matplotlib>=3.0',
    ],
)
