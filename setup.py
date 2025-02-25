from setuptools import setup, find_packages

setup(
    name="mobaseq",
    packages=find_packages(include=['mobaseq', 'mobaseq.*']),
    package_dir={'': '.'},
    include_package_data=True
)