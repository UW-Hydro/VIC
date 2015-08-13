# Other VIC Model Output Files

In addition to the default output files created by VIC, there are a few other formats that have historically been used in various projects. They are documented here.

## LDAS Project Output Files

The LDAS project used a specific format designed to minimize file size by writing a limited number of variables and storing them in the most compressed file format available. Most of the output variables are stored as short int values in a binary file. Short int values require only 2 bytes for storage but do not allow the range of values normally stored so most variables are modified by a constant before storage.

For LDAS Project output files, [click here](LDASOutputFile.md).
