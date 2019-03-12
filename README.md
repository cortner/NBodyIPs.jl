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
