# Operations

TAME takes advantage of established libraries for heavy numerical tasks.
Currently, the computations are performed with the Numpy library. Other
computation backends might be implemented in hte future.

!!! danger "documentation in construction"

    The documentation now contains an **incomplete** list of available
    operations in TAME, with very brief introductions of what they do,
    a more detailed API documentation should be available later.

## Time-dependent Operations

Below lists operations that are time-dependent, which means the value of each
frame depend on previous frames.

| function       | operation                        |
| -------------- | -------------------------------- |
| `tm.op.unwrap` | unwraps the coordinate [^unwrap] |

[^unwrap]:
    TAME unwraps the trajectory using the minimal displacement vector at each
    time step, to avoid possible artifact when unwrapping constant-pressure
    simulations, which is discussed by von Bülow _et al._ in
    detail.[@2020_BuelowBullerjahnEtAl]

## Basic Functions

The below functions are basic mathematical operations, those functions are
also defined as python operators between `tm.FrameArray`s.

| function     | syntax | operation      |
| ------------ | ------ | -------------- |
| `tm.plus`    | `a+b`  | addition       |
| `tm.minus`   | `a-b`  | subtraction    |
| `tm.product` | `a*b`  | multiplication |

\bibliography
