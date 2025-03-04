import pathlib

from setuptools import setup

readme_file = pathlib.Path(__file__).parent.resolve() / 'README.md'
readme_contents = readme_file.read_text()

setup(
    name="bhresist",
    version="0.1",
    packages=['bhr'],
    description="Does great things with boreholes",
    package_data={},
    include_package_data=False,
    long_description=readme_contents,
    long_description_content_type='text/markdown',
    author='Matt Mitchell',
    author_email='mitchute@gmail.com',
    url='https://github.com/NREL/BHResist')
