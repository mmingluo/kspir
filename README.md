Private Information Retrieval (PIR)
 Protocol kspir
=====

***

 
The [kspir](https://github.com/parsear/kspir) is the implement of the paper [Faster FHE-based Single-Server Private Information Retrieval]. This library will be continuously updated to support more features.

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

After install [Intel HEXL](https://github.com/intel/hexl), you can build and run our PIR protocol by

```
cd kspir
mkdir build
cd build
cmake ..
make
./tests/test-pir
```


An example output:
```
Packing number: 16
Database configuration: 32768 * 8 KB, total database size 256 MB
BSGS parameters: (N1: 128, N2: 16)

target_col: 177, target_packing: 4

 server preprocessing costs 4144 ms.
 query costs 7873 us.
 online server response costs 236635 us.
 decrypt costs 284 us.

the recovered result = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, ..., 4096]
```
