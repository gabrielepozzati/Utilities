#!/usr/bin/env python3
import argparse


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = ' ')
    parser.add_argument('-l', required= True, help='list')
    parser.add_argument('-f', required= True, help='msa folder')
    parser.add_argument('-o', required= True, help='out folder')
    ns = parser.parse_args()

    idlist = open(ns.l, 'r')

    for msa_id in idlist: 
        msa = open(ns.f+msa_id.rstrip(), 'r')
        count = 0
        cutted_msa = open(ns.o+msa_id.rstrip().rstrip('?'), 'w')
        for line in msa:

            if line.startswith('>'): count += 1
            if count <= 1000: cutted_msa.write(line)

            if line.startswith('>'):
                lastseq = []
                lastseq.append(line)
            else: lastseq.append(line)

        for line in lastseq: cutted_msa.write(line)
        cutted_msa.close()
        msa.close()
    idlist.close()



