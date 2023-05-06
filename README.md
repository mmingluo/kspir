Private Information Retrieval (PIR)
 Protocol kspir
=====

***

 
The [kspir](https://github.com/parsear/kspir) is the implement of the paper [Faster FHE-based Single-Server PIR].

## Requirements
The runnig of this code requires a basic c++ toolchain, including

* g++ compiler (clang++12 or newer is optimal)
* make
* cmake


## Build and Run
We adapt [Intel HEXL](https://github.com/intel/hexl) library to implement the NTTs.
Therefore, please install it firstly.
We recommend that install it to your default user's program path `/usr/local`, and the suggeted version of them are showed in following table. What's more, using `clang++12` or newer compiler can achieve better performance.

| Libraries | Version | Website |
| ---- | ---- | ---- |
| HEXL  | v1.2.5 | https://github.com/intel/hexl |


We recommend using the following command to install HEXL. You can also refer to the website for more detailed information.

```
git clone https://github.com/intel/hexl.git -b v1.2.5 --single-branch 
cd hexl
mkdir build
cd build
cmake ..
make
sudo make install
```

After install [Intel HEXL](https://github.com/intel/hexl), you can build and run our tests by

```
cd kspir
mkdir build
cd build
cmake ..
make
./test-first-dimension
```
