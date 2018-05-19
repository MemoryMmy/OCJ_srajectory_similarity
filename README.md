# An Efficient Query Algorithm for Trajectory Similarity Based on MPI using the Fréchet Distance Threshold 

This is a similar trajectories finding program we submitted to the ACM SIGSPATIAL Cup 2017 (http://sigspatial2017.sigspatial.org/giscup2017/). In this competition, the goal is to find the similar trajectories of given trajectories by using Fréchet distance as similarity measurement. In our method, we create spatial indexes for the first and last points of the trajectories in the dataset.txt separately, which are used to filter out most of the dissimilar trajectories and generate a much smaller candidate set. Then an Ordered-Coverage-Judge Fréchet distance algorithm is presented to search the similar trajectories form the candidate set.

# Software dependencies:
	
	MPICH2
	Boost 1.62


# Compile:
	
	Change the path of Boost 1.62 in Makefile, then run the following command to generate the executable program:
	> make

# Run:
	Note: When the number of trajectories in the dataset.txt is more than 50000, the minimun RAM required is 32GB.
	Use the following command to run the program.
	> mpirun -np ${cpu_cores} Trajquery datasetfile queryfile outputdir
	${cpu_cores} represents the number of cores of the computer, it can be gotten through the following command:
	> grep 'core id' /proc/cpuinfo | sort -u | wc -l
	Example:	
	> mpirun -np 4 ./Trajquery ./data/dataset.txt ./data/queries.txt ./result/
	
# Central idea:

In our method, we create R-tree indexes for the first and last points of the trajectories in the dataset.txt separately, which could filter out most of the dissimilar trajectories and generate a much smaller candidate trajectories set. Then an Ordered-Coverage-Judge Fréchet distance algorithm is presented to find the similar trajectories from the candidate set. The Ordered-Coverage-Judge Fréchet distance algorithm is a depth-first heuristic search algorithm. This algorithm can search out whether there exists a match distance of two discrete trajectories within the given Fréchet distance using fewer searches.

# Paper:

Guo N, Ma M, Xiong W, et al. An Efficient Query Algorithm for Trajectory Similarity Based on Fréchet Distance Threshold[J]. International Journal of Geo-Information, 2017, 6(11):326.

# Contact:

	Mengyu Ma @ National University of Defense Technology
	Email: mamengyu10@nudt.edu.cn

