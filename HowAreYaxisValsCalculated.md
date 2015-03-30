# This is how the Y-axis values are calculated (Step by step) #

  1. Each genomic region (e.g. TSS+/-2Kb) is first extended by the fragment length on both sides (to create a "buffer" in the flanking regions);
  1. The short reads that are overlapping the extended regions are retrieved from BAM files;
  1. Coverage (or depth) at single base resolution are calculated by extending each short read to the length of the fragment;
  1. If the algorithm is spline: the coverage vector is fit to a spline and then 101 points are sampled at equal interval;
  1. If the algorithm is bin: the coverage vector is separated into 101 equal-sized bins and the averaged value for each bin is calculated;
  1. The negative values in the coverage vectors are forced to zeros (This is caused by spline fit);
  1. Each value is normalized by library size using: val = val / libsize `*` 1e6;
  1. If the setting is bam-pair, repeat the above procedure for the background bam file; The final value = log2(foreground value / background value) with pseudo count used to avoid division by zeros;
  1. The coverage (or log2 ratio) vectors for each gene list are averaged to produce the average profiles.
**DONE**