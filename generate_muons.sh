for depth in 10 100 1000 #10000
do
    python muon_injector.py -num 100000 -depth ${depth}
done