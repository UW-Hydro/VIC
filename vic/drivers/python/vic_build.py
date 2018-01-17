from cffi import FFI
from vic_headers import headers

ffi = FFI()
ffi.cdef(headers)
ffi.set_source('vic._vic', None)

if __name__ == '__main__':
    ffi.compile()
