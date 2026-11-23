Clustered earthquake sequence with alternative geometries and epistemic
branches: test C of the HAZ cluster tests (HAZ tests/cluster/TC_Geometries,
Cluster run, site 1, Memphis).

Six segments (two geometries of each of NMS, RFT and NMN) and eight
segmentation models, each picking one geometry of each fault. Every member
has the same branches, matched by index: three thicknesses, two rates times
two rupture sets, two mechanisms; and two GMMs.

In out1, HAZ writes zeros in the members' blocks and, after the last member,
the cluster's branch curves for each segmentation model as sources 7 to 14.
The fractile program draws a segmentation model and one branch, and adds that
curve. The hazard is therefore one of 384 branch curves with known weights;
scripts/test_reference.py checks the fractiles against the exact quantiles of
that distribution and the mean against its exact mean.
