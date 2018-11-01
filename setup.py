#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
from setuptools import setup, Extension


class get_pybind_include(object):
    def __init__(self, user=False):
        self.user = user

    def __str__(self):
        import pybind11
        return pybind11.get_include(self.user)


flags = ['-march=native']
if sys.platform == 'darwin':
    flags += ['-mmacosx-version-min=10.9']

ext_modules = [
    Extension(
        'vice.benchmark',
        ['vice/benchmark.cpp'],
        include_dirs=[
            get_pybind_include(),
            get_pybind_include(user=True),
            "vendor/starry",
            "vendor/starry/lib/eigen_3.3.3",
            "vendor/starry/lib/boost_1_66_0",
        ],
        extra_compile_args=[
            '-O2',
            '-DNDEBUG',
            '-march=native',
            '-std=c++14',
            '-mmacosx-version-min=10.9',
        ],
        extra_link_args=[
            '-march=native',
            '-mmacosx-version-min=10.9',
        ],
    )
]

setup(
    name='vice',
    version='0',
    packages=['vice'],
    ext_modules=ext_modules,
    install_requires=['pybind11>=2.2'],
    include_package_data=True,
    zip_safe=False,
)
