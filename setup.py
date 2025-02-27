from setuptools import setup, find_packages

setup(
    name="mobaseq",
    version="0.1.0",
    packages=find_packages(include=['mobaseq', 'mobaseq.*']),
    package_dir={'': '.'},
    include_package_data=True,
    dependency_links=[
        'git+https://github.com/IrenaeusChan/mobaseq_dev.git@main#egg=mobaseq-0.1.0'
    ],
    python_requires='>=3.7'
)