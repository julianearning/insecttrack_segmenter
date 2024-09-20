# Insecttrack Segmenter   

This C++-project is the result of my master's thesis. With it, a DVS-Pointcloud can be read and the instances of insect-tracks can be segmented. Additionally, it calculates 2D-curves for each instance. 

# Compilation  

There are several requirements. 
Please install the following packages: Boost, tbb, lz4, Eigen3
To install the software please change directory to the build directory and make
```
cd build
make
```


This generates the executable 'out'.

# Usage  

Input a CSV, that is comma-seperated and has no header and has at leas three columns: (x,y,t). You can only input 4D data like this: (x,y,z,t) if you add the parameter -4d. 
Always add the path to the input file as the second argument ( ./out [path] { arguments ... }  ). This, by default, outputs debug information, to shut it off add -silent as a parameter. Further help of the arguments can be accessed by 
adding a parameter that isn't known (like -help). 

The output has this format: (x,y,t,class,r,g,b,scalar) and the default output filename is test.csv to change this add -out_filename <filename> to the parameters. The class information is in the fourth column, the r, g, b columns are useful if you want to visualize the data in CloudCompare and the last one is for debugging. 

There are four example datasets, two that are small and two that are large. To visualize the small datasets, you can use the python script csv2gnuplot.py, but for the big ones you should use CloudCompare:

```
./out ../data/example_crossing.csv
python3 ../csv2gnuplot.py test.csv | gnuplot
./out ../data/example_big_track.csv
python3 ../csv2gnuplot.py test.csv | gnuplot
```

This should plot the clusterings. To plot the LOWESS lines, there is a R-Scipt included. For this the 'scales' package has to be installed, because the points in the background have transparency. 

```
cd ..
R
> source("plot_lowess.r")
```

There is one example dataset that doesn't have noise and one with noise:

```
cd build
./out ../data/only_insecttracks.csv -silent
./out ../data/noisy.csv -silent -octree_denoise 3 -mst_max_mean_edge_len 10
```
these can be visualized in CloudCompare. 


# Authors and Copyright  

Juliane Arning (2024), the license can be found in the file LICENSE


Also, this project includes several external libraries: EuclidianMST (https://github.com/AndrewB330/EuclideanMST by Andrii Borziak ), FLANN (https://github.com/flann-lib/flann by Marius Muja et al), CppWeightedLowess (https://github.com/LTLA/CppWeightedLowess/ by Aaron Lun )
