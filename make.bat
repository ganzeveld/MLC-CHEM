echo off
echo CLEANING....
del messy_emdep_emis.o messy_emdep_emis_mem.o messy_emdep.o messy_emdep_mem.o messy_emdep_xtsurf_box.o messy_emdep_xtsurf.o messy_main_constants_mem.o messy_main_tools.o *.mod *.old emdep_xtsurf.exe 

echo COMPILING...
g95 -cpp -O0 -c -ftrace=full -fbounds-check messy_main_constants_mem.f90
g95 -cpp -O0 -c -ftrace=full -fbounds-check messy_main_tools.f90
g95 -cpp -O0 -c -ftrace=full -fbounds-check messy_emdep.f90
g95 -cpp -O0 -c -ftrace=full -fbounds-check messy_emdep_emis_mem.f90
g95 -cpp -O0 -c -ftrace=full -fbounds-check messy_emdep_emis.f90
g95 -cpp -O0 -c -ftrace=full -fbounds-check messy_emdep_mem.f90
g95 -cpp -O0 -c -ftrace=full -fbounds-check messy_emdep_xtsurf.f90
g95 -cpp -O0 -c -ftrace=full -fbounds-check messy_emdep_xtsurf_box.f90

echo LINKING...
g95 -cpp -O0 messy_emdep_emis.o messy_emdep_emis_mem.o messy_emdep.o messy_emdep_mem.o messy_emdep_xtsurf_box.o messy_emdep_xtsurf.o messy_main_constants_mem.o messy_main_tools.o -o emdep_xtsurf.exe
