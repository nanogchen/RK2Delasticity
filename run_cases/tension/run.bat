@echo off

del *.exe

set CURRENT=%CD%

if exist "data_post" (
echo folder existed
) else (
mkdir data_post
)

if exist "data_pre" (
echo folder existed
) else (
mkdir data_pre
)

if exist "data_pro" (
echo folder existed
) else (
mkdir data_pro
)

cd ../../source/preprocessing/

make -f Makefile.mak clean

make -f Makefile.mak

cp RKpre.exe %CURRENT%

cd %CURRENT%

RKpre.exe

cd ../../source/processing/

make -f Makefile.mak clean

make -f Makefile.mak

cp RK.exe %CURRENT%

cd %CURRENT%

RK.exe

cd ../../source/postprocessing

make -f Makefile.mak clean

make -f Makefile.mak

cp RKpost.exe %CURRENT%

cd %CURRENT%

RKpost.exe

cd data_post/

if exist "postProcessing.m" (
 echo MATLAB files exist
) else (
cp -r ../../../MATLAB-postprocessing/. .
)

cp ../data_pre/xynodes.txt .

rem matlab postProcessing.m
