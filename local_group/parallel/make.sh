cd demoni
make
cp demoni.out ../run_mond_main/demoni.out
cp demoni.out ../run_mond_new/demoni.out
cp demoni.out ../run_mond_heavy/demoni.out
cp demoni.out ../run_mond_4m/demoni.out
cp demoni.out ../run_mond_2m/demoni.out
cd ..
cd newhexi
make
cp newhexi.out ../run_newton_main/newhexi.out
cp newhexi.out ../run_newton_new/newhexi.out
cp newhexi.out ../run_newton_4m/newhexi.out
cp newhexi.out ../run_newton_2m/newhexi.out
cp newhexi.out ../run_newton_bighalo/newhexi.out
cp newhexi.out ../run_newton_smallhalo/newhexi.out
cd ..
cd genetic
mpif90 geneticalgorithm.f90 -o genetic.out
cp genetic.out ../run_mond_main/genetic.out
cp genetic.out ../run_mond_new/genetic.out
cp genetic.out ../run_mond_heavy/genetic.out
cp genetic.out ../run_mond_4m/genetic.out
cp genetic.out ../run_mond_2m/genetic.out
cp genetic.out ../run_newton_main/genetic.out
cp genetic.out ../run_newton_new/genetic.out
cp genetic.out ../run_newton_4m/genetic.out
cp genetic.out ../run_newton_2m/genetic.out
cp genetic.out ../run_newton_bighalo/genetic.out
cp genetic.out ../run_newton_smallhalo/genetic.out
cd ..
cd mkinput
mpif90 mkinput.f90 -o mkinput.out
cp mkinput.out ../run_mond_main/mkinput.out
cp mkinput.out ../run_mond_new/mkinput.out
cp mkinput.out ../run_mond_heavy/mkinput.out
cp mkinput.out ../run_mond_4m/mkinput.out
cp mkinput.out ../run_mond_2m/mkinput.out
cp mkinput.out ../run_newton_main/mkinput.out
cp mkinput.out ../run_newton_new/mkinput.out
cp mkinput.out ../run_newton_4m/mkinput.out
cp mkinput.out ../run_newton_2m/mkinput.out
cp mkinput.out ../run_newton_bighalo/mkinput.out
cp mkinput.out ../run_newton_smallhalo/mkinput.out
cd ..
cd model
mpif90 mkmodel.f90 -o mkmodel.out
cp mkmodel.out ../run_mond_main/mkmodel.out
cp mkmodel.out ../run_mond_new/mkmodel.out
cp mkmodel.out ../run_mond_heavy/mkmodel.out
cp mkmodel.out ../run_mond_4m/mkmodel.out
cp mkmodel.out ../run_mond_2m/mkmodel.out
cp mkmodel.out ../run_newton_main/mkmodel.out
cp mkmodel.out ../run_newton_new/mkmodel.out
cp mkmodel.out ../run_newton_4m/mkmodel.out
cp mkmodel.out ../run_newton_2m/mkmodel.out
cp mkmodel.out ../run_newton_bighalo/mkmodel.out
cp mkmodel.out ../run_newton_smallhalo/mkmodel.out
cd ..
cd outputbest
mpif90 outputbest.f90 -o outputbest.out
cp outputbest.out ../run_mond_main/outputbest.out
cp outputbest.out ../run_mond_new/outputbest.out
cp outputbest.out ../run_mond_heavy/outputbest.out
cp outputbest.out ../run_mond_4m/outputbest.out
cp outputbest.out ../run_mond_2m/outputbest.out
cp outputbest.out ../run_newton_main/outputbest.out
cp outputbest.out ../run_newton_new/outputbest.out
cp outputbest.out ../run_newton_4m/outputbest.out
cp outputbest.out ../run_newton_2m/outputbest.out
cp outputbest.out ../run_newton_bighalo/outputbest.out
cp outputbest.out ../run_newton_smallhalo/outputbest.out
cd ..
cd cutter
mpif90 cutter.f90 -o cutter.out
cp cutter.out ../run_mond_main/cutter.out
cp cutter.out ../run_mond_new/cutter.out
cp cutter.out ../run_mond_heavy/cutter.out
cp cutter.out ../run_mond_4m/cutter.out
cp cutter.out ../run_mond_2m/cutter.out
cp cutter.out ../run_newton_main/cutter.out
cp cutter.out ../run_newton_new/cutter.out
cp cutter.out ../run_newton_4m/cutter.out
cp cutter.out ../run_newton_2m/cutter.out
cp cutter.out ../run_newton_bighalo/cutter.out
cp cutter.out ../run_newton_smallhalo/cutter.out
cd ..
cd real_distribution
mpif90 realdist.f90 -o real.out
cp real.out ../run_mond_main/real.out
cp real.out ../run_mond_new/real.out
cp real.out ../run_mond_heavy/real.out
cp real.out ../run_mond_4m/real.out
cp real.out ../run_mond_2m/real.out
cp real.out ../run_newton_main/real.out
cp real.out ../run_newton_new/real.out
cp real.out ../run_newton_4m/real.out
cp real.out ../run_newton_2m/real.out
cp real.out ../run_newton_bighalo/real.out
cp real.out ../run_newton_smallhalo/real.out
cd ..
cd run_newton_main
sh prepare.sh
cd ..
cd run_newton_new
sh prepare.sh
cd ..
cd run_mond_heavy
sh prepare.sh
cd ..
cd run_newton_4m
sh prepare.sh
cd ..
cd run_newton_2m
sh prepare.sh
cd ..
cd run_newton_bighalo
sh prepare.sh
cd ..
cd run_newton_smallhalo
sh prepare.sh
cd ..
cd run_mond_main
sh prepare.sh
cd ..
cd run_mond_new
sh prepare.sh
cd ..
cd run_mond_4m
sh prepare.sh
cd ..
cd run_mond_2m
sh prepare.sh
cd ..
