from distutils.core import setup

setup(
    name='sam2fpkc',
    version='0.0.1',
    author='Michael Imelfort',
    author_email='mike@mikeimelfort.com',
    packages=['sam2fpkc', 'sam2fpkc.test'],
    scripts=['bin/sam2fpkc'],
    url='http://pypi.python.org/pypi/sam2fpkc/',
    license='LICENSE.txt',
    description='sam2fpkc',
    long_description=open('README.txt').read(),
    install_requires=[],
)
