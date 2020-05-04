#!/bin/bash

set -e

# Extract & adapt parts from MS-MPI SDK & runtime
# MS-MPI 10.1.2

source=(
  'https://download.microsoft.com/download/a/5/2/a5207ca5-1203-491a-8fb8-906fd68ae623/msmpisdk.msi'
  'https://download.microsoft.com/download/a/5/2/a5207ca5-1203-491a-8fb8-906fd68ae623/msmpisetup.exe'
)
sha256sums=(

)

if [[ $MSYSTEM_CARCH != "x86_64" ]]; then 
	echo "ERROR: this scripts requires 64-bit MSYS2 environment"
	exit 1
fi

# exe
mpiexec=`cygpath -m "$MSMPI_BIN/mpiexec.exe"`
if [[ ! -f "$mpiexec" ]]; then
	echo "ERROR: mpiexec.exe is missing. Is MS-MPI installed?"
	exit 1
fi
cp "$mpiexec" .

# x64
msmpi=`cygpath -m "$WINDIR/system32/msmpi.dll"`
if [[ ! -f "$msmpi" ]]; then
	echo "ERROR: 64-bit msmpi.dll is missing. Is MS-MPI runtime installed?"
	exit 1
fi
/mingw64/bin/gendef "$msmpi" 2> /dev/null
cat msmpi.def | tr -d '\015' > msmpi.def.x86_64
rm msmpi.def

## x32
#msmpi=`cygpath -m "$WINDIR/syswow64/msmpi.dll"`
#if [[ ! -f "$msmpi" ]]; then
#	echo "ERROR: 32-bit msmpi.dll is missing. Is MS-MPI runtime installed?"
#	exit 1
#fi
#/mingw32/bin/gendef "$msmpi" 2> /dev/null
#cp msmpi.def msmpi.def~
#ruby -n -e '$_.strip!; puts "#{$1}@0=#{$1}\n" if ($_ != $_.upcase) && ($_ != $_.downcase) && /^([a-zA-Z0-9_]*)$/.match($_)' msmpi.def >> msmpi.def~
#cat msmpi.def~ | tr -d '\015' > msmpi.def.i686
#rm msmpi.def msmpi.def~

rootdir=`pwd`

(
	cd "`cygpath -m "${MSMPI_INC}"`"
	if [[ ! -f mpi.h ]]; then
		echo "ERROR: \$MSMPI_INC directory does not exist or contains no mpi.h file. Is MS-MPI SDK installed?"
		exit 1
	fi
	for f in mpi.h mpif.h mpi.f90; do
		cat $f | tr -d '\015' > "$rootdir/$f"
	done
	cat x64/mpifptr.h | tr -d '\015' > "$rootdir/mpifptr.h.x86_64"
	#cat x86/mpifptr.h | tr -d '\015' > "$rootdir/mpifptr.h.i686"
)

# Export file signatures to be embed into PKGBUILD in order to ensure build integrity
ruby export.rb mpi.c mpi.h mpif.h mpi.f90 mpifptr.h.x86_64 msmpi.def.x86_64 mpiexec.exe
#ruby export.rb mpi.c mpi.h mpif.h mpi.f90 mpifptr.h.{x86_64,i686} msmpi.def.{x86_64,i686} mpiexec.exe
