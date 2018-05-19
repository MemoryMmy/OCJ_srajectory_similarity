## OCJ_srajectory_similarity
# An Efficient Query Algorithm for Trajectory Similarity Based on MPI using the Fréchet Distance Threshold 

This is a similar trajectories finding program we submitted to the ACM SIGSPATIAL Cup 2017 (http://sigspatial2017.sigspatial.org/giscup2017/). In this competition, the goal is to find the similar trajectories of given trajectories by using Fréchet distance as similarity measurement. In our method, we create spatial indexes for the first and last points of the trajectories in the dataset.txt separately, which are used to filter out most of the dissimilar trajectories and generate a much smaller candidate set. Then an Ordered-Coverage-Judge Fréchet distance algorithm is presented to search the similar trajectories form the candidate set.
