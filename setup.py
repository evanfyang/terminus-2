import setuptools

setuptools.setup(
    name = 'terminus-2',
    version = '0.1.0dev',
    author = 'Evan F. Yang',
    author_email = 'gnay.nave@gmail.com',
    description = 'This is a demo package',
    long_description = open('README.txt').read(),
    url = 'https://github.com/evanfyang/terminus-2',
    project_urls={
        "Bug Tracker": "https://github.com/evanfyang/terminus-2/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    packages = setuptools.find_packages('src'),
    license = 'MIT license',
    python_requires = '>=3.8'
)