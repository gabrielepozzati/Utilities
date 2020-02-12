while read p; do
    wget http://www.rcsb.org/pdb/files/$p.pdb.gz
    #wget https://files.rcsb.org/download/$p.pdb1
done<$1

