Readme file for VIC source code formating using uncrustify (http://uncrustify.sourceforge.net/)

Purpose:
To create easy to read, constistantly formatted source code.

To format a VIC source code file:

    uncrustify -c uncrustify_VIC_c.cfg [files ...]

or with replacement and without backup files:

    uncrustify -c uncrustify_VIC_c.cfg --replace --no-backup [files ...]

There are other command line options available:
man uncrustify

