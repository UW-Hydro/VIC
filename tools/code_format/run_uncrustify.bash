#!/usr/bin/env bash

# Run uncrustify on all VIC source files

files=`find ../../drivers ../../vic_run -name '*.c' -o -name '*.h'`
uncrustify -c uncrustify_VIC_c.cfg -l C --replace --no-backup $files
