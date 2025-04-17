cd demoni
make
cp demoni.out ../run_mond/demoni.out
cd ..
cd newhexi
make
cp newhexi.out ../run_newton/newhexi.out
cd ..
cd genetic
ifort geneticalgorithm.f90 -o genetic.out
cp genetic.out ../run_mond/genetic.out
cp genetic.out ../run_newton/genetic.out
cd ..
cd mkinput
ifort mkinput.f90 -o mkinput.out
cp mkinput.out ../run_mond/mkinput.out
cp mkinput.out ../run_newton/mkinput.out
cd ..
cd model
ifort mkmodel.f90 -o mkmodel.out
cp mkmodel.out ../run_mond/mkmodel.out
cp mkmodel.out ../run_newton/mkmodel.out
cd ..
cd outputbest
ifort outputbest.f90 -o outputbest.out
cp outputbest.out ../run_mond/outputbest.out
cp outputbest.out ../run_newton/outputbest.out
cd ..
cd cutter
ifort cutter.f90 -o cutter.out
cp cutter.out ../run_mond/cutter.out
cp cutter.out ../run_newton/cutter.out
cd ..
cd real_distribution
ifort realdist.f90 -o real.out
cp real.out ../run_mond/real.out
cp real.out ../run_newton/real.out
cd ..
 
