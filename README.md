# novaSTA

C++ version of subtomogram averaging (SA) scripts from TOM/AV3 package [https://doi.org/10.1073/pnas.0409178102](https://doi.org/10.1073/pnas.0409178102). Both CPU and GPU parallelization is supported although the latter performs significantly worse in terms of processing time (the code is not well optimized) and is thus not recommended for larger datasets.

The program requires both CUDA and OpenMP libraries and was tested with foss 2017b and CUDA 9.1.85. 

Compilation:

`make includepath="path_include_files" libpath="path_to_libraries" cudapath="path_to_cuda_include_files"`

Example:

`make includepath="/path/to/openmpi/OpenMPI/2.1.1/include /path/to/fftw/fftw-3.3.4/include" libpath="/usr/lib /path/to/fftw/fftw-3.3.4/lib /path/to/fftw/fftw-3.1.2/lib64" cudapath="/path/to/cuda/CUDA/9.1.85"`

Run with CPU parallelization:

`mpirun -n #numbeOfCores novaSTA -param parameter_file.txt`

Run on GPU(s):

`./novaSTA -param parameter_file.txt -useGPU 1`

The source code for novaSTA is distributed under an GPL v.3 license. The code can be redistributed and/or modified under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or any later version.

The code is distributed WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
