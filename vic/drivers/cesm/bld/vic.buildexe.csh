#! /bin/csh -f

cd $OBJROOT/lnd/obj

set comp = 'unknown'
set ROUT = 'rout_stub'

if ($COMP_INTERFACE == 'MCT' ) set comp = mct
if ($COMP_INTERFACE == 'ESMF') set comp = esmf

cat >! Filepath << EOF
$CASEROOT/SourceMods/src.vic
$CODEROOT/lnd/vic/vic
$CODEROOT/lnd/vic/vic/vic_run/src
$CODEROOT/lnd/vic/vic/vic_run/include
$CODEROOT/lnd/vic/vic/drivers/shared_all/src
$CODEROOT/lnd/vic/vic/drivers/shared_all/include
$CODEROOT/lnd/vic/vic/drivers/shared_image/src
$CODEROOT/lnd/vic/vic/drivers/shared_image/include
$CODEROOT/lnd/vic/vic/drivers/cesm/src
$CODEROOT/lnd/vic/vic/drivers/cesm/include
$CODEROOT/lnd/vic/vic/drivers/cesm/cpl_$comp/
$CODEROOT/lnd/vic/vic/extensions/rout_stub/src
$CODEROOT/lnd/vic/vic/extensions/rout_stub/include
EOF

set vicdefs = ""

gmake complib -j $GMAKE_J MODEL=vic COMPLIB=$LIBROOT/liblnd.a USER_CPPDEFS="$vicdefs" \
    -f $CASETOOLS/Makefile MACFILE=$CASEROOT/Macros.$MACH || exit 2
