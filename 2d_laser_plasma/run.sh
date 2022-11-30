cd /Users/kevinambrogioni/Thesis/Examples/pic_tests/2d_laser_plasma/smilei
export OMP_NUM_THREADS=1 
mpirun -np 2 ./smilei input.py 

cd /Users/kevinambrogioni/Thesis/Examples/pic_tests/2d_laser_plasma/warpx
mpirun -np 2 ./warpx.2d.MPI.OMP.DP.PDP.OPMD.QED input.txt

cd /Users/kevinambrogioni/Thesis/Examples/pic_tests/2d_laser_plasma/epoch
echo Data | mpirun -np 2 ./epoch2d 