##################################################################
#                                                                #
# MultiNEAs: Numerical tools for near-earth asteroid dynamics   #
#            and population                                      #
#                                                                #
##################################################################
# License: GNU Affero General Public License v3 (AGPL-3.0)        #
##################################################################
import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    # ######################################################################
    # BASIC DESCRIPTION
    # ######################################################################
    name='multineas',
    author="Jorge I. Zuluaga, Juanita A. Agudelo",
    author_email="jorge.zuluaga@udea.edu.co, juanita.agudelo@udea.edu.co",
    description="MultiNEAs: Numerical tools for near-earth asteroid dynamics and population",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/seap-udea/MultiNEAs",
    keywords='astronomy asteroids near-earth-asteroids dynamics celestial-mechanics',
    license='AGPL-3.0-only',

    # ######################################################################
    # CLASSIFIER
    # ######################################################################
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Astronomy",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "License :: OSI Approved :: GNU Affero General Public License v3",
        "Operating System :: OS Independent",
    ],
    version='0.3.4',

    # ######################################################################
    # FILES
    # ######################################################################
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    
    # ######################################################################
    # TESTS
    # ######################################################################
    test_suite='pytest',
    tests_require=['pytest'],

    # ######################################################################
    # DEPENDENCIES
    # ######################################################################
    install_requires=[
        'numpy>=1.20.0',
        'scipy>=1.7.0',
        'matplotlib>=3.3.0',
        'spiceypy>=5.0.0',
    ],
    
    python_requires='>=3.8',

    # ######################################################################
    # OPTIONS
    # ######################################################################
    include_package_data=True,
    package_data={"": ["data/*.*", "tests/*.*"]},
)
