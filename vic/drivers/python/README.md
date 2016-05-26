# Python Driver Readme

### Requirements
- [CFFI](http://cffi.readthedocs.org/en/latest/index.html) version 1.2 or greater

### Installing
run `python setup.py install` from the `vic/drivers/python` directory. `setup.py` will automatically generate the headers (`vic_headers.py`) file that `CFFI` requires for the C-Python bindings.

### Usage
```python
from vic import lib as vic_lib

vic_lib.print_license()
```
