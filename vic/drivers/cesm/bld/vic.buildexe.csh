#! /bin/csh -f

cd $OBJROOT/lnd/obj

set comp = 'unknown'
if ($COMP_INTERFACE == 'MCT' ) set comp = mct
if ($COMP_INTERFACE == 'ESMF') set comp = esmf

cat >! Filepath << EOF
$CASEROOT/SourceMods/src.vic
$CODEROOT/lnd/vic/vic
$CODEROOT/lnd/vic/vic/vic_run/src
$CODEROOT/lnd/vic/vic/vic_run/include
$CODEROOT/lnd/vic/vic/drivers/shared/src
$CODEROOT/lnd/vic/vic/drivers/shared/include
$CODEROOT/lnd/vic/vic/drivers/rasm/src
$CODEROOT/lnd/vic/vic/drivers/rasm/include
EOF

set vicdefs = "`cat $CASEBUILD/vicconf/CESM_cppdefs`"

gmake complib -j $GMAKE_J MODEL=vic COMPLIB=$LIBROOT/liblnd.a USER_CPPDEFS="$vicdefs" -f $CASETOOLS/Makefile MACFILE=$CASEROOT/Macros.$MACH || exit 2
