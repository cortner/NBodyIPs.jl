# NBodyIPs.jl

<!--
[![Build Status](https://travis-ci.org/cortner/NBodyIPs.jl.svg?branch=master)](https://travis-ci.org/cortner/NBodyIPs.jl)

[![Coverage Status](https://coveralls.io/repos/cortner/NBodyIPs.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/cortner/NBodyIPs.jl?branch=master)

[![codecov.io](http://codecov.io/github/cortner/NBodyIPs.jl/coverage.svg?branch=master)](http://codecov.io/github/cortner/NBodyIPs.jl?branch=master)
-->


## Implementation Notes

As I am rewriting / redesigning some aspects of this package, I will slowly
update some thoughts here about how NBodyIPs works internally.


### `hash`

In `NBodyIPs` there is a dummy-type `BASIS` defined. This is primarily used
to dispatch the `hash` function as follows: Each basis function implements
```
hash(::BASIS, V::T) = ...
```
which returns a hash key that describes only the basis function `V`
but **disregards** its parameters. For example, `OneBody(1.234)` has a
hash that includes the type information as well as the value `1.234` while
```
hash(::BASIS, V::OneBody) = hash(OneBody)
```
i.e. this ignores the parameter `E0 = 1.234`.

These basis-hashs are used to group and combine the basis functions to
allow for faster assembly.


### Reference potential

The `lsqfit` function allows the kw argument `Vref` which can be an arbitrary
`AbstractCalculator` that can be used as a reference potential from which to
start the fit. I.e. the final fit will be `Vref + Vfit`. If the
kwarg `E0` is provided, then `Vref = OneBody(E0)` is used as default.
Providing both `E0` and `Vref` means that `E0` will be ignored.

## Descriptors

The two main descriptors in `NBodyIPs` are `BondAngleDesc` and `BondLengthDesc`.
Eacha are defined via a space transform and a cutoff function. Please
see inline documentation for these.

### Cutoff functions

* CosCut
* CosCut2s
* PolyCut
* PolyCutSym
* PolyCut2sA

### Space transforms

* ExpTransform : x = exp(- A*(r/rnb - 1)
* PolyTransform : x = (r0/r)^p
* IdTransform : x = r
* CoulombTransform : x = 1/r
* MorseTransform : c = exp(-r)
* AnalyticTransform: this used to be called SpaceTransform, use only for experimenting. If a good new transform can be found then it can be implemented into `NBodyIPs` as a separate type.

## Regularisation

Regularisers


## More on Descriptors and Basis sets

### Bond-length descriptor

### Bond-angle descriptor

### Environment-dependent N-Body potentials




<!--

TODO
 * regularisation
 * partial data
    - unclear this is a good idea, leave for now
 * LASSO              >>> Genevieve
 * sparse grid basis  >>> with Genevieve next week??

Mess around at night:
 * W fit
 * Al fit
 * Si fit => or not

Big ones
 * environment dependence for fullsymm
 * bond-angle potentials + environment dependence

-->
