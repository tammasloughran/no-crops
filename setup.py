import setuptools

with open('requirements.txt', 'r') as requirements_file:
    requirements_list = requirements_file.read().splitlines()

setuptools.setup(
        name='no-crops',
        packages=['no-crops'],
        install_requires=requirements_list)
