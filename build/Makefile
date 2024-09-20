CC = g++
FLANN_ROOT = /home/juliane/bienen_bahnen/src/track_segmenter/flann/src/cpp
EIGEN_ROOT = /usr/include/eigen3
LOWESS_ROOT = /home/juliane/bienen_bahnen/src/track_segmenter/loess/CppWeightedLowess/include/
FAST_MST_ROOT = /home/juliane/bienen_bahnen/src/track_segmenter/EuclideanMST/standalone_header/
CFLAGS = -Wall -llz4 -O3 -std=c++17 -Wfloat-equal -ggdb3 -ltbb -fopenmp 

all: out

#out: main.o option.o fileio.o pointcloud.o mathstuff.o octree.o space.o dbscan.o kdtree.o mst.o assign.o loess.o postprocess.o
#	${CC} -o out main.o option.o fileio.o pointcloud.o octree.o mathstuff.o space.o dbscan.o kdtree.o mst.o assign.o loess.o postprocess.o ${CFLAGS} 

out: main.o option.o fileio.o pointcloud.o octree.o dbscan.o mathstuff.o assign.o mst.o loess.o postprocess.o
	${CC} -o out main.o option.o fileio.o pointcloud.o octree.o dbscan.o mathstuff.o assign.o mst.o loess.o postprocess.o ${CFLAGS} 

main.o: ../main.cpp
	${CC} -c ../main.cpp ${CFLAGS} -I ${FLANN_ROOT}

loess.o: ../loess/loess.cpp ../loess/loess.hpp 
	${CC} -c ../loess/loess.cpp -I ${FLANN_ROOT} -I ${LOWESS_ROOT} ${CFLAGS}

option.o: ../utils/option.cpp ../utils/option.hpp
	${CC} -c ../utils/option.cpp ${CFLAGS}

fileio.o: ../utils/fileio.cpp ../utils/fileio.hpp
	${CC} -c ../utils/fileio.cpp ${CFLAGS}

postprocess.o: ../postprocess/postprocess.cpp ../postprocess/postprocess.hpp
	${CC} -c ../postprocess/postprocess.cpp ${CFLAGS}

mathstuff.o: ../utils/mathstuff.cpp ../utils/mathstuff.hpp
	${CC} -c ../utils/mathstuff.cpp ${CFLAGS} -I ${EIGEN_ROOT}	

dbscan.o: ../dbscan/dbscan.cpp ../dbscan/dbscan.hpp
	${CC} -c ../dbscan/dbscan.cpp ${CFLAGS} -I ${FLANN_ROOT} -I ${EIGEN_ROOT}

kdtree.o: ../kdtree/kdtree.cpp ../kdtree/kdtree.hpp
	${CC} -c ../kdtree/kdtree.cpp ${CFLAGS}

mst.o: ../mst/mst.cpp ../mst/mst.hpp
	${CC} -c  ../mst/mst.hpp  ../mst/mst.cpp ${CFLAGS} -I ${FAST_MST_ROOT}

pointcloud.o: ../utils/pointcloud.cpp ../utils/pointcloud.hpp
	${CC} -c ../utils/pointcloud.cpp ${CFLAGS} -I ${EIGEN_ROOT}

octree.o: ../preprocessing/octree.cpp ../preprocessing/octree.hpp
	${CC} -c ../preprocessing/octree.cpp ${CFLAGS}

space.o: ../preprocessing/space.cpp ../preprocessing/space.hpp
	${CC} -c ../preprocessing/space.cpp ${CFLAGS} -I ${FLANN_ROOT}

assign.o: ../assign_to_cluster/assign.cpp ../assign_to_cluster/assign.hpp
	${CC} -c ../assign_to_cluster/assign.cpp ${CFLAGS} -I ${FLANN_ROOT}


clean: 
	rm -f *.o
