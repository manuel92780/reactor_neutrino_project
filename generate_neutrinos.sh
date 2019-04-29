for dist in 10 100 1000 10000
do
    python neutrino_injector.py -num 100000 -dist ${dist} -nu_type solar
    python neutrino_injector.py -num 100000 -dist ${dist} -nu_type reactor_Ginna
done