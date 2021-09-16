# RationalFunctionFields

The package provides convenenient interface for fields of rational polynomial functions via the main struct RationalFunctionField. Currently, the main contribution is the assession of field membership.

Assets are freely available here: 
- source: https://github.com/sumiya11/RationalFunctionFields.jl
- documentation: ???

## Example

Installing the package via Pkg from Github

```
import Pkg;
Pkg.add(url="https://github.com/sumiya11/RationalFunctionFields.jl")

using RationalFunctionFields
```

Creating generators of the desired field 

```
using Nemo

R, (x, y) = PolynomialRing(QQ, ["x", "y"])

generators = [
            (x^2 + y^2) // R(1),
            1 // (x^2)
]
```

Initializing the field of rational functions

```
rff = RationalFunctionField(generators)
```

Should print `true` with the probability of at least *99.9%*

```
contains(rff, y^2 // R(1), p=0.999)
```

The usage of deterministic algorithm
```
contains(rff, y^2 // R(1), p=1)
contains(rff, x*y // R(1), p=1)
```

