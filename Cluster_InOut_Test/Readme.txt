In- and out-of-cluster branches in one fault system: test D of the HAZ
cluster tests (HAZ tests/cluster/TD_InOut, Cluster run, site 1, Memphis).

The fault system holds the six members of a clustered sequence (two
geometries of each of NMS, RFT and NMN) and RFT-S and RFT-L as ordinary
segments with out-of-cluster rates. Its eleven segmentation models follow
HID Figure H-5.5-1: the eight geometry combinations of the cluster (0.9 x the
geometry weight), RFT-S or RFT-L alone (0.05 x 0.7 / 0.3), and no segment at
all (0.05).

In out1, HAZ writes zeros in the members' blocks, then the cluster's curves
for each of the eleven models as sources 9 to 19 (zeros for the three models
with no member), then the ordinary segments' blocks. The fractile program
draws a model, adds the cluster's curve for it on one drawn branch, and adds
the ordinary segments the model includes as usual. Every fractile below 0.05
is zero (no activity). scripts/test_reference.py checks the fractiles
against the exact quantiles of the mixture and the mean against its exact
mean.
