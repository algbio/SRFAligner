# efg-simplify
Program that checks whether a given (Elastic) Founder Graph (in xGFA format) respects the repeat-free property or the semi-repeat-free property.

```
Usage: efg-simplify inputgraph.xgfa simplifiedgraph.xgfa
Program to transform and simplify an Elastic Founder Graph given in xGFA
format.

The program takes an Elastic Founder Graph in xGFA format and merges adjacent
blocks that only contain parallel paths.

  -h, --help     Print help and exit
  -V, --version  Print version and exit
```

## GFA format (xGFA)
See [here](https://github.com/algbio/founderblockgraphs/blob/master/xGFAspec.md).

## known issues

- as we merge blocks, paths warp a little bit: the program extends the original paths up to the first non-simplified node at the beginning and end

## todo

- tests
- solve the paths extension
