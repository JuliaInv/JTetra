[![Build Status](https://travis-ci.org/cschwarzbach/JTetra.svg?branch=master)](https://travis-ci.org/cschwarzbach/JTetra)
[![Coverage Status](https://coveralls.io/repos/github/cschwarzbach/JTetra/badge.svg?branch=master)](https://coveralls.io/github/cschwarzbach/JTetra?branch=master)
[![Build status](https://ci.appveyor.com/api/projects/status/u4gaso4i7kx3iivo/branch/master?svg=true)](https://ci.appveyor.com/project/cschwarzbach/jtetra/branch/master)

# JTetra
Tetrahedral mesh support for jInv

With the module `JTetra` tetrahedral meshes can be used in the [`jInv`](https://github.com/JuliaInv/jInv.jl) package.

# Requirements

`JTetra` is intended for use with julia 0.5.x. [`jInv`](https://github.com/JuliaInv/jInv.jl) is the only required dependency.

# Examples

A forward modeling example for electromagnetics on a tetrahedral mesh can be found in the `examples` folder.

# Installation

In julia type:
```
Pkg.clone("https://github.com/jInv/JTetra.git")
```

# Tests

From the julia prompt run
```
Pkg.test("JTetra")
```
