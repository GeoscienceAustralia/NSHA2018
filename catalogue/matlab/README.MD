Steps to generate matlab mat file, to be used in magnitude revision and conversion routines

1) Parse Gary Gibson's earthquake catalogue using "ggcat2mat.m"
2) Merge GA magnitudes for larger earthquakes using "mergeGA.m"
3) Append alternate ML values using "merge_alt_mags.m"
4) Append update MW values using "append_mw/append_mw.m"