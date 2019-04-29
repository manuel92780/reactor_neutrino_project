for dist in 50 500 #10 100 1000
do
    python neutrino_injector.py -num 100000 -dist ${dist} -nu_type solar
    python neutrino_injector.py -num 100000 -dist ${dist} -nu_type reactor_Ginna
done