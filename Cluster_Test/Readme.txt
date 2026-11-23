Clustered earthquake sequence (source type 8): test B of the HAZ cluster tests
(HAZ tests/cluster/TB_NMFS, Cluster run, site 1, Memphis).

The NMS, RFT and NMN faults of the New Madrid Fault System are the members of
one cluster. In out1, HAZ writes the cluster's branch curves in the first
member's block and zeros in the other members' blocks, so the fractile
program samples the cluster as one source.

The cluster has two GMMs weighted 0.5 / 0.5 and no other epistemic branch, so
each Monte Carlo sample is one of two curves. scripts/test_reference.py checks
that the fractiles below 0.5 are the lower curve and those above 0.5 the
higher one, that the mean is within Monte Carlo error of the average of the
two, and that the 5th and 95th UHS values are the two curves at each hazard
level. The two periods (PGA and 1.0 s) have different test ground motions.
