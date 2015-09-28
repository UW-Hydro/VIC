#!/usr/bin/env bash

# Run uncrustify on all VIC source files

files=`find ../../vic -name '*.c' -o -name '*.h'`
uncrustify -c uncrustify_VIC_c.cfg -l C --replace --no-backup $files
