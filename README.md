# Fractiles

Computes hazard fractiles and the mean hazard by Monte Carlo sampling of the
epistemic branches of a HAZ run (versions 45.1 to 45.4 fault files; out1 from
HAZ 45.3 and 45.4). It reads the HAZ run file (GMM weights), the fault file
(seismic source branches and weights) and out1 (each source's branch curves).

## Build

```
mkdir build
cd build
cmake ..
make
```

## Run

```
Fractiles
```

and enter the fractile run file name. See `Reference_Test/Fractiles/Input/Run_Fractiles.txt`.

In the UHS table, a negative value flags a fractile curve that does not cross
the hazard level within the test ground motions: minus the highest test
ground motion if the curve is still above the hazard level there, minus the
last ground motion with non-zero hazard if the curve drops to zero first, or
minus the lowest test ground motion if the curve is already below the hazard
level there.

## Clustered sequences (source type 8)

A cluster member is read as source type 1. In out1, HAZ writes zeros in the
members' blocks and, after the last member, the cluster's branch curves for
each segmentation model of its fault system (for example each combination of
alternative geometries), numbered after the last source. The program samples
the cluster as one source: it draws a segmentation model and one branch (the
members' branches are matched by index), and adds that curve.

## Tests

`Reference_Test`, `PEER_Verification_Test`, `Cluster_Test`,
`Cluster_Geometries_Test` and `Source_Types_Test` each hold a Hazard case (HAZ inputs and outputs) and a
Fractiles case. From `scripts/`:

```
FRACT_BIN=../build/Fractiles python -m pytest test_reference.py
```

Setting `HAZ_BIN` to the HAZ executable also reruns the Hazard cases (a few
minutes) and compares their out1 and out3 files.
