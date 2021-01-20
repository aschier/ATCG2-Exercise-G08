# ATCG1 Exercise
Welcome to ATCG1 Exercise

This file explains everything you need to know on how to build the project.

We use:
* [CMake](https://cmake.org/) - to build the project
* [libigl](https://github.com/libigl/libigl) - to visualize point clouds and meshes

## Folder hierarchy

We will explain the folder hierarchy to make sure you know where everything is located.

- **.** is the root of our project.
    ```
    It contains the setup.sh, and the download_*.sh and all the other files and folders
    ```
- **src** contains all the files you need to edit to complete the assignments
    ```
    It contains the main.cpp, lib_*.cpp
    ```
- **data** contains all the data for the assignment
    ```
    It contains the *.xyz, *.txt, *.off, *.obj etc. data files
    ```
- **include** contains all the headers for the files in src
    ```
    It contains headers to the files in src
    ```
- **lib** contains all the files we provide for the assignment
    ```
    It contains mostly headers you can use
    ```
- **externals** contains all frameworks from external sources
    ```
    It contains libigl
    ```
- **cmake** contains all files for CMake to find the files in externals
    ```
    It contains the cmake file for finding libigl
    ```
The following two folders will be re-/generated by setup.sh
- **bin** contains the binary assignmentsheet1
    ```
    It contains the executable
    ```
- **build** contains the CMake cache and all compiled libraries
    ```
    It contains the compiled stuff of libigl
    ```
## You need to make the externals folder manually and then run download_*.sh

## Why is there a setup.sh?
We want a clean build of our program and put all generated files into specific folders.

CMake doesn't know how to clean the project properly, especially to remove the CMakeCache.txt from **build** when you want a clean build of your program.

Even though everything has to be compiled again, you will avoid a lot of problems this way.
## We only guarantee support for Linux
## Getting Started on Linux

1. Make sure you are using at least CMake 2.8.12
2. Download git
3. Download libigl by running in your shell
    ```
    sh download_libigl.sh
    ```

4. Setup your favorite IDE for C/C++ to run the binary assignmentsheet1 with 2 parameters:
    ```
    data/figurine.xyz data/measurements.txt data/maxear.off data/maxsimple.obj
    ```

5. Build and run the project with:
    * the **IDE** 
    * or from the command line with:
        1. Compile project by running **once**
        
            ```sh setup.sh```
            
            then run
        2. ```sh build.sh```
        3. ```./bin/assignmentsheet1 data/figurine.xyz data/measurements.txt data/maxear.off data/maxsimple.obj```
## Getting Started on Windows

1. Make sure you are using at least CMake 2.8.12
2. Download git for windows, you need the git bash console
3. Download libigl by running in your git bash
    ```
    sh download_libigl.sh
    ```
4. open build.sh with your favorite texteditor and change this line:

	```
	cmake --build build -- -j 13
	```

	to

	```
	cmake --build build #-- -j 13
	```
5. open setup.sh with your favorite texteditor and add this line:

 	```
	-DCMAKE_GENERATOR_PLATFORM=x64
	```
	to the line starting with

 	```
	cmake -E chdir build cmake
	```	

6. Setup your favorite IDE for C/C++ to run the binary assignmentsheet1 with 2 parameters:
    ```
    data/figurine.xyz data/measurements.txt data/maxear.off data/maxsimple.obj
    ```
7. Build and run the project with:
    * the **IDE** 
    * or from the git bash with:
        1. Compile project by running **once**
		
            ```sh setup.sh```
			
            then run
			
        2. ```sh build.sh```
        3. ```./bin/assignmentsheet1 data/figurine.xyz data/measurements.txt data/maxear.off data/maxsimple.obj```
## Getting Started on Mac

1. Make sure you are using at least CMake 2.8.12
2. Download git
3. Download libigl by running in your git bash
    ```
    sh download_libigl.sh
    ```
4. Setup your favorite IDE for C/C++ to run the binary assignmentsheet1 with 2 parameters:
    ```
    data/figurine.xyz data/measurements.txt data/maxear.off data/maxsimple.obj
    ```
5. Build and run the project with:
    * the **IDE** 
    * or from the git bash with:
        1. Compile project by running **once**
		
            ```sh setup.sh```
			
            then run
			
        2. ```sh build.sh```
        3. ```./bin/assignmentsheet1 data/figurine.xyz data/measurements.txt data/maxear.off data/maxsimple.obj```

## Hints
- **you can use the setup.sh in your ide to build the whole framework**
- **you can use the build.sh in your ide to build the only your changes**
                           
## Troubleshooting

- **you did not run setup.sh first** -> run setup.sh
- **your ide is not working in the correct working path** -> set the working path to the directory containing this README.md
- **your ide does generate a lot of linking errors** -> delete the CMakeCache.txt in build once before building with your ide
- **your ide tells you that you need some arguments** -> add the arguments to your ide binary parameters
- **you still need help** -> write to the mailing list, other might have the same problem and know a solution
## Resources which might be helpful

* [libigl-tutorial](http://libigl.github.io/libigl/tutorial/tutorial.html#meshrepresentation)
* [CMake-tutorial](https://cmake.org/cmake-tutorial/)
* [CLion](https://www.jetbrains.com/clion/) - a cross platform IDE with students-license


