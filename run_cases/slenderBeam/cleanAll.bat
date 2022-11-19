@echo off

set CURRENT=%CD%

cd ../../source/preprocessing/

make -f Makefile.mak clean

cd %CURRENT%

cd ../../source/processing/

make -f Makefile.mak clean

cd %CURRENT%

cd ../../source/postprocessing

make -f Makefile.mak clean