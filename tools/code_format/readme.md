Readme file for VIC source code formating using uncrustify (http://uncrustify.sourceforge.net/)

Code Format Guidelines and Best Practices:
1.  Run uncrustify before commiting to git.
2.  Add space after typecast e.g. `(double) myvar / other_var`
3.  Be explicit about integer types, e.g. declarations: `unsigned int myvar`.and typecasting: `(unsigned short int) myvar`.
4.  Use integer and float printf format specifiers for all integer and floating point types, regarless of signedness or precision.

Purpose:
To create easy to read, constistantly formatted source code.

To format a VIC source code file:

    uncrustify -c uncrustify_VIC_c.cfg [files ...]

or with replacement and without backup files:

    uncrustify -c uncrustify_VIC_c.cfg --replace --no-backup [files ...]

There are other command line options available:
man uncrustify

