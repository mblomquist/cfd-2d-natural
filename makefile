# makefile for 2D CFD Repository
#
# Written by Matt Blomquist
# Last Update: 2018-02-25 (YYYY-MM-DD)
#
# This file compiles and links the cfd-2d repository source code into the
# executable file main2d.out
#
main2d.out : convergence2d.o solver2d_gmres.o solver2d_bicg.o solver2d_bicgstab.o solver2d_bicgstab2.o solver2d_tdma.o geometry2d.o initialize2d.o output_results2d.o pressure_solve2d.o pressure_correct2d.o pressure2d.o simpler2d.o temperature_solve2d.o temperature_source2d.o temperature2d.o velocity_correct2d.o pseudo_solve2d.o velocity_solve2d.o velocity_source2d.o velocity2d.o main2d.o
	ifort -o build/main2d.out -mkl build/solver2d_tdma.o build/solver2d_gmres.o build/solver2d_bicg.o build/solver2d_bicgstab.o build/solver2d_bicgstab2.o build/convergence2d.o build/geometry2d.o build/initialize2d.o build/output_results2d.o build/pressure_solve2d.o build/pressure_correct2d.o build/pressure2d.o build/pseudo_solve2d.o build/simpler2d.o build/temperature_solve2d.o build/temperature_source2d.o build/temperature2d.o build/velocity_correct2d.o build/velocity_solve2d.o build/velocity_source2d.o build/velocity2d.o build/main2d.o

main2d.o : src/main2d.f90
	ifort -o build/main2d.o -c src/main2d.f90

convergence2d.o : src/convergence2d.f90
	ifort -o build/convergence2d.o -c src/convergence2d.f90

geometry2d.o : src/geometry2d.f90
	ifort -o build/geometry2d.o -c src/geometry2d.f90

initialize2d.o : src/initialize2d.f90
	ifort -o build/initialize2d.o -c src/initialize2d.f90

output_results2d.o : src/output_results2d.f90
	ifort -o build/output_results2d.o -c src/output_results2d.f90

pressure_solve2d.o : src/pressure_solve2d.f90
	ifort -o build/pressure_solve2d.o -c src/pressure_solve2d.f90

pressure_correct2d.o : src/pressure_correct2d.f90
	ifort -o build/pressure_correct2d.o -c src/pressure_correct2d.f90

pressure2d.o : src/pressure2d.f90
	ifort -o build/pressure2d.o -c src/pressure2d.f90

pseudo_solve2d.o : src/pseudo_solve2d.f90
	ifort -o build/pseudo_solve2d.o -c src/pseudo_solve2d.f90

simpler2d.o : src/simpler2d.f90
	ifort -o build/simpler2d.o -c src/simpler2d.f90

solver2d_bicg.o : src/solver2d_bicg.f90
	ifort -o build/solver2d_bicg.o -c src/solver2d_bicg.f90

solver2d_bicgstab.o : src/solver2d_bicgstab.f90
	ifort -o build/solver2d_bicgstab.o -c src/solver2d_bicgstab.f90

solver2d_bicgstab2.o : src/solver2d_bicgstab2.f90
	ifort -o build/solver2d_bicgstab2.o -c src/solver2d_bicgstab2.f90

solver2d_gmres.o : src/solver2d_gmres.f90
	ifort -o build/solver2d_gmres.o -c src/solver2d_gmres.f90

solver2d_tdma.o : src/solver2d_tdma.f90
	ifort -o build/solver2d_tdma.o -c src/solver2d_tdma.f90

temperature_solve2d.o : src/temperature_solve2d.f90
	ifort -o build/temperature_solve2d.o -c src/temperature_solve2d.f90

temperature_source2d.o : src/temperature_source2d.f90
	ifort -o build/temperature_source2d.o -c src/temperature_source2d.f90

temperature2d.o : src/temperature2d.f90
	ifort -o build/temperature2d.o -c src/temperature2d.f90

velocity_correct2d.o : src/velocity_correct2d.f90
	ifort -o build/velocity_correct2d.o -c src/velocity_correct2d.f90

velocity_solve2d.o : src/velocity_solve2d.f90
	ifort -o build/velocity_solve2d.o -c src/velocity_solve2d.f90

velocity_source2d.o : src/velocity_source2d.f90
	ifort -o build/velocity_source2d.o -c src/velocity_source2d.f90

velocity2d.o : src/velocity2d.f90
	ifort -o build/velocity2d.o -c src/velocity2d.f90

clean :
	rm build/convergence2d.o build/geometry2d.o build/solver2d_bicgstab.o build/solver2d_bicgstab2.o build/solver2d_gmres.o build/solver2d_bicg.o build/solver2d_tdma.o build/initialize2d.o build/output_results2d.o build/pressure_solve2d.o build/pressure_correct2d.o build/pressure2d.o build/pseudo_solve2d.o build/simpler2d.o build/temperature_solve2d.o build/temperature_source2d.o build/temperature2d.o build/velocity_correct2d.o build/velocity_solve2d.o build/velocity_source2d.o build/velocity2d.o build/main2d.o
