# ngeom

ngeom is a library for doing geometry in N dimensions.

ngeom ships with modules for rigid Euclidean 2D and 3D geometry.
The built-in modules provide primitives for points, lines, and planes,
along with rigid transformations (rotations, translations, and reflections.)

If you prefer to own your geometric primitives,
or if you desire a space with different dimensionality, metric, handedness, etc.
then you can use the provided macro
to generate all the necessary geometric operations on your own `struct`s.

ngeom uses [homogeneous coordinates](https://en.wikipedia.org/wiki/Homogeneous_coordinates) to express ideal/infinite points,
ideal/infinite lines, etc. and to provide for exception-free meet & join operations.

ngeom is generic over the [scalar] datatype, and can be used with `f32`, `f64`, or custom datatypes.
It can even use integers, with some restrictions--
most functionality requires only scalar addition and multiplication.

With some effort, usage code can be written to be generic
over the dimensionality or even the metric of the space.

ngeom does not use SIMD intrinsics,
but rather emits compact expressions and structs
that are easy for the compiler to optimize into SIMD instructions
in an architecture-independent way.

ngeom is dependency-free.
ngeom is `no_std`-compatible, with most functionality available,
and the option to implement the rest.

Under the hood, ngeom is able to make these generalizations by using [geometric algebra](https://en.wikipedia.org/wiki/Geometric_algebra).
Understanding geometric algebra is helpful but optional,
as the functions are named for their geometric interpretation.
Lower-level functions named after their geometric algebra expressions are available if needed.
