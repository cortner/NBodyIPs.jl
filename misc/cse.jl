
using CommonSubexpressions

ex1a = :( x[1]*x[2]*x[3] + x[1]*x[2]*x[4] + x[1]*x[2]*x[5] )
cse(ex1a)

ex1b = :( x1*x2*x3 + x1*x2*x4 + x1*x2*x5 )
cse(ex1b)

display(@macroexpand(@cse x[1]*x[2]*x[3] + x[1]*x[2]*x[4] + x[1]*x[2]*x[5]))
display(@macroexpand(@cse x1*x2*x3 + x1*x2*x4 + x1*x2*x5))

ex = :(
   x1*x2*x3*x4*x5 +
   x1*x2*x3*x4*x6 +
   x1*x2*x3*x4*x8 +
   x1*x2*x5*x6*x7 +
   x1*x3*x5*x6*x7 +
   x1*x5*x6*x7*x8 +
   x1*x2*x3*x4*x7 +
   x1*x4*x5*x6*x7 +
   x1*x2*x5*x8*x9 +
   x2*x3*x5*x8*x9 +
   x2*x5*x6*x8*x9 +
   x1*x2*x3*x4*x9 +
   x2*x4*x5*x8*x9 +
   x1*x5*x6*x7*x9 +
   x2*x5*x7*x8*x9 +
   x1*x3*x6*x8*x10 +
   x2*x3*x6*x8*x10 +
   x3*x5*x6*x8*x10 +
   x1*x2*x3*x4*x10 +
   x3*x4*x6*x8*x10 +
   x1*x5*x6*x7*x10 +
   x3*x6*x7*x8*x10 +
   x2*x5*x8*x9*x10 +
   x3*x6*x8*x9*x10 +
   x1*x4*x7*x9*x10 +
   x2*x4*x7*x9*x10 +
   x4*x5*x7*x9*x10 +
   x3*x4*x7*x9*x10 +
   x4*x6*x7*x9*x10 +
   x4*x7*x8*x9*x10
)

const ex_cse = cse(ex)
