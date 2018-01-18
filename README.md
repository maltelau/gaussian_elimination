

# Gaussian elimination
A python implementation of the algorithm for solving systems of linear equations taught in MAT121: Linear algebra at University of Bergen.

For example, we start out with a matrix representing the equation system:

```
|     2,     -1,      1 |     0 |
|     2,     -1,      1 |     0 |
|     4,      1,      1 |     2 |
```

Specified as 

```python
eq1 = equation_system([[2, -1, 1],  # lhs
                       [2, -1, 1],
                       [4, 1, 1]],
                       [0, 0, 2])   # rhs
```

We can run through the algorithm and see what steps are necessary to solve the system:

```python
eq1.to_reduced_triangular()
```
```
============================ LOG =================================

|     2,     -1,      1 |     0 |
|     2,     -1,      1 |     0 |
|     4,      1,      1 |     2 |

Row 2 -> Row 2 - 1.000 * Row 1
Row 3 -> Row 3 - 2.000 * Row 1
Row 2 <-> Row 3
Now in triangular form.

|     2,     -1,      1 |     0 |
|   0.0,    3.0,   -1.0 |   2.0 |
|   0.0,    0.0,    0.0 |   0.0 |

Row 2 -> Row 2 * 0.333
Row 1 -> Row 1 * 0.500
Row 1 -> Row 1 + -0.500 * Row 2
Now in reduced triangular.

|   1.0,    0.0,  0.333 | 0.333 |
|   0.0,    1.0, -0.333 | 0.667 |
|   0.0,    0.0,    0.0 |   0.0 |

```

It can circle the pivot points for us:
```

| (1.0),    0.0,  0.333 | 0.333 |
|   0.0,  (1.0), -0.333 | 0.667 |
|   0.0,    0.0,    0.0 |   0.0 |

```

And it can print a little more information if we want it.

```python
eq2 = equation_system([[2, -1, 1],  # lhs
                       [2, -1, 1],
                       [4, 1, 1]],
                       [0, 0, 2],   # rhs
                      debug = True, full_log = True)
```

```
Not all coefficients below the leading coefficient [1, 1] are zero.

|     2,     -1,      1 |     0 |
|     2,     -1,      1 |     0 |
|     4,      1,      1 |     2 |
Has to be in triangular before it can be reduced.

============== Beginning algorithm to get triangular  ===========

|     2,     -1,      1 |     0 |
|     2,     -1,      1 |     0 |
|     4,      1,      1 |     2 |
Improving row 2 by subtracting 1.0 times row 1

|     2,     -1,      1 |     0 |
|   0.0,    0.0,    0.0 |   0.0 |
|     4,      1,      1 |     2 |
Improving row 3 by subtracting 2.0 times row 1

|     2,     -1,      1 |     0 |
|   0.0,    0.0,    0.0 |   0.0 |
|   0.0,    3.0,   -1.0 |   2.0 |
There's a nonzero row below a zero row.
The pivotposition is zero, swapping rows.

|     2,     -1,      1 |     0 |
|   0.0,    3.0,   -1.0 |   2.0 |
|   0.0,    0.0,    0.0 |   0.0 |

|     2,     -1,      1 |     0 |
|   0.0,    3.0,   -1.0 |   2.0 |
|   0.0,    0.0,    0.0 |   0.0 |

======== Beginning algorithm to get reduced triangular ========

|     2,     -1,      1 |     0 |
|   0.0,    3.0,   -1.0 |   2.0 |
|   0.0,    0.0,    0.0 |   0.0 |
The leading coefficient of Row 2 is not 1

|     2,     -1,      1 |     0 |
|   0.0,    1.0, -0.333 | 0.667 |
|   0.0,    0.0,    0.0 |   0.0 |
The leading coefficient of Row 1 is not 1

|   1.0,   -0.5,    0.5 |   0.0 |
|   0.0,    1.0, -0.333 | 0.667 |
|   0.0,    0.0,    0.0 |   0.0 |
The leading coefficient at [2, 2] is not the only one in its column! (Row 1).
Improving row 1 by subtracting -0.5 times row 2

|   1.0,    0.0,  0.333 | 0.333 |
|   0.0,    1.0, -0.333 | 0.667 |
|   0.0,    0.0,    0.0 |   0.0 |

============================ LOG =================================

|     2,     -1,      1 |     0 |
|     2,     -1,      1 |     0 |
|     4,      1,      1 |     2 |

Row 2 -> Row 2 - 1.000 * Row 1
|     2,     -1,      1 |     0 |
|   0.0,    0.0,    0.0 |   0.0 |
|     4,      1,      1 |     2 |

Row 3 -> Row 3 - 2.000 * Row 1
|     2,     -1,      1 |     0 |
|   0.0,    0.0,    0.0 |   0.0 |
|   0.0,    3.0,   -1.0 |   2.0 |

Row 2 <-> Row 3
|     2,     -1,      1 |     0 |
|   0.0,    3.0,   -1.0 |   2.0 |
|   0.0,    0.0,    0.0 |   0.0 |
Now in triangular form.

|     2,     -1,      1 |     0 |
|   0.0,    3.0,   -1.0 |   2.0 |
|   0.0,    0.0,    0.0 |   0.0 |

Row 2 -> Row 2 * 0.333
|     2,     -1,      1 |     0 |
|   0.0,    1.0, -0.333 | 0.667 |
|   0.0,    0.0,    0.0 |   0.0 |

Row 1 -> Row 1 * 0.500
|   1.0,   -0.5,    0.5 |   0.0 |
|   0.0,    1.0, -0.333 | 0.667 |
|   0.0,    0.0,    0.0 |   0.0 |

Row 1 -> Row 1 + -0.500 * Row 2
|   1.0,    0.0,  0.333 | 0.333 |
|   0.0,    1.0, -0.333 | 0.667 |
|   0.0,    0.0,    0.0 |   0.0 |

Now in reduced triangular.

|   1.0,    0.0,  0.333 | 0.333 |
|   0.0,    1.0, -0.333 | 0.667 |
|   0.0,    0.0,    0.0 |   0.0 |
```



Use responsibly.

