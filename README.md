Installing
==========

This will currently not compiled without supporting fortran library code, which has, as yet, not been made available. With that in mind ...

Set the location for your module directory, which contains all the supporting fortran library code:

```
export FORTRANMODULEDIR='/path/to/directory/'
```

Clone the repo and then:

```
cd oceanbot
mkdir build
cd build
cmake ..
make
```

