[![Build Status](https://travis-ci.org/JuliaInv/JTetra.svg?branch=master)](https://travis-ci.org/JuliaInv/JTetra)
[![Coverage Status](https://coveralls.io/repos/github/JuliaInv/JTetra/badge.svg?branch=master)](https://coveralls.io/github/JuliaInv/JTetra?branch=master)
[![Build status](https://ci.appveyor.com/api/projects/status/dbdcs1u4wv6255r9/branch/master?svg=true)](https://ci.appveyor.com/project/cschwarzbach/jtetra-yqlgs/branch/master)

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
Pkg.clone("https://github.com/JuliaInv/JTetra.git")
```

# Tests

From the julia prompt run
```
Pkg.test("JTetra")
```
