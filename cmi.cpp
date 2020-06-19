#include <unordered_map>
#include <stdexcept> 
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cctype>
#include <math.h>
#include <cmath>
#include <array>
#include <map>

#include"cnpy.h"

std::string split(std::string str, const char delim, int pos){
    int count = 0, pointer = 0;
    std::string token, prev;
    std::stringstream ss(str);

    while (std::getline (ss, token, delim)) {count+=1;}
    ss.str(""); ss.clear(); ss.str(str);

    std::vector< std::string > splitted;

    splitted.resize(count);

    while (std::getline (ss, token, delim)) {splitted[pointer]=token; pointer+=1;}

    if (pos >= 0){
	if (pos >= count) throw std::runtime_error("Out of splitting range!");
        token=splitted[pos]; 
    }
    else {
	if (pos < -count) throw std::runtime_error("Out of splitting range!");
        token=splitted[count+pos];
    }

    return token;
}

double MI(std::vector< std::vector < int > > parsedmsa,
          std::vector< std::unordered_map < int, double > > msafreqs,
	  int i, int j) {

    int x, y;
    double mi, Pij;
    unsigned int seq, count;
    //std::unordered_map < int, double > msaposcount;
    std::array< int, 2 > couple;
    std::map < std::array< int, 2 >, double > jointfreqs;

    std::unordered_map < int, double > Pi = msafreqs[i], Pj = msafreqs[j];

    for (seq=0; seq<parsedmsa.size(); seq++) {
	if (parsedmsa[seq][i]==21 || parsedmsa[seq][j]==21) {continue;}
    	couple[0]=parsedmsa[seq][i];
	couple[1]=parsedmsa[seq][j];
	jointfreqs[couple] ++;
    }

//    for (auto coupleidx=jointfreqs.begin(); coupleidx!=jointfreqs.end(); coupleidx++) {
//        count += coupleidx->second;
//    }
//    for (auto coupleidx=jointfreqs.begin(); coupleidx!=jointfreqs.end(); coupleidx++) {
//        jointfreqs[coupleidx->first] = coupleidx->second/count;
//    }
    for (auto coupleidx=jointfreqs.begin(); coupleidx!=jointfreqs.end(); coupleidx++) {
	Pij = coupleidx->second;
        x = coupleidx->first[0];
	y = coupleidx->first[1];
        mi += Pij*log10(Pij/(Pi[x]*Pj[y]));
    }
    
    return mi;
}

int main (int argc, char** argv){

    if (argc != 4) {
    std::cout << "3 command line arguments required (in this order): " << std::endl;
    std::cout << "	path of input MSA file in a3m format;" << std::endl;
    std::cout << "      path of file with sequence lenghts for APC" << std::endl;
    std::cout << "      path to save output" << std::endl;
    return 0;
    }

    std::string msapath = argv[1], apcpath = argv[2], outpath = argv[3];

    double mi;
    unsigned int i = 0, j = 0, minlength = 0, msalength = 0, gapcount = 0, hitcount = 0, count = 0;
    std::string msaline, apcline, msahit = "", residue, names, id1, id2;
    std::ifstream msafile, apcfile;
    std::vector< int > parsedmsaline;
    std::vector< std::vector < int > > parsedmsa;

    std::vector< double > mi_line;
    std::vector< std::vector < double > > mi_matrix;

    std::vector< std::unordered_map < int, double > > msafreqs;

    std::unordered_map< std::string, int > aamap;

    aamap["A"] = 1, aamap["C"] = 2, aamap["D"] = 3, aamap["E"] = 4, aamap["F"] = 5;
    aamap["G"] = 6, aamap["H"] = 7, aamap["I"] = 8, aamap["K"] = 9, aamap["L"] = 10;
    aamap["M"] = 11, aamap["N"] = 12, aamap["P"] = 13, aamap["Q"] = 14, aamap["R"] = 15;
    aamap["S"] = 16, aamap["T"] = 17, aamap["V"] = 18, aamap["W"] = 19, aamap["Y"] = 20;
    aamap["U"] = 21, aamap["Z"] = 21, aamap["X"] = 21, aamap["J"] = 21, aamap["B"] = 21; 
    aamap["O"] = 21, aamap["-"] = 21;

    names = split(msapath, '/', -1);
    names = split(msapath, '.', 0);
    id1 = split(names, '_', 0);
    id2 = split(names, '_', 1);

    // checking MSA dimensions //
    msafile.open(msapath); apcfile.open(apcpath);
    if(!msafile) throw std::runtime_error("Unable to open msa file!");
    if(!apcfile) throw std::runtime_error("Unable to open APC file!");

    while (std::getline(msafile, msaline)) {
        if (msaline.at(0) == '>') {
	    if (minlength == 0 || minlength < msalength) {minlength = msalength;}
	    msalength = 0; hitcount += 1;
	}
        else {msalength += msaline.length();}
    }

    std::cout << minlength << " " << hitcount << std::endl;

    msafile.clear();
    msafile.seekg(0, std::ios::beg);

    parsedmsaline.resize(minlength, 21);
    parsedmsa.resize(hitcount, parsedmsaline);
    msafreqs.resize(minlength);

    // MSA parsing, removing lowcase residues and seqs with more than 80% gaps  //
    j = 0;
    while (std::getline(msafile, msaline)) {
        std::cout << msahit.length() << " " << minlength << " " << gapcount << std::endl;

        if (msahit.length() < minlength) {
            for (i=0; i<msaline.length(); i++) {if (msaline.at(i) == '-') {gapcount += 1;}}
            msahit += msaline;
        }

        if (msahit.length() >= minlength && gapcount/msahit.length() <= 0.9) {
            std::cout << "Wrinting!" << std::endl;
            for (i=0; i<msahit.length(); i++) {
		residue.push_back(msahit.at(i));
                if (aamap.find(residue) != aamap.end()) {parsedmsa[j][i] = aamap.at(residue); residue="";}
	    }

            j += 1;
        }

        if (msaline.at(0) == '>') {msahit = ""; gapcount = 0; continue;}

    }
    parsedmsa.resize(j);
       
    std::cout << parsedmsa.size() << " " << j << std::endl;
    std::cout << parsedmsa[0].size() << std::endl;
    for (count=0; count<minlength; count++) { std::cout << parsedmsa[5][count] << " ";}
    std::cout << std::endl;
    

    // AA frequency for each position of the MSA //
    for (j=0; j<parsedmsa.size(); j++) {
        for (i=0; i<minlength; i++) {msafreqs[i][parsedmsa[j][i]]++;}
    }

//    for (i=0; i<msafreqs.size(); i++) {
//	for (j=1; j==20; j++) {if (msafreqs[i].find(j) != msafreqs[i].end()) {count += msafreqs[i][j];}}
//	for (j=1; j==20; j++) {if (msafreqs[i].find(j) != msafreqs[i].end()) {msafreqs[i][j] = msafreqs[i][j]/count;}}
//    }

    // MI calculation //
    mi_line.resize(minlength, 0.0);
    mi_matrix.resize(minlength, mi_line);
    for (i=0; i<=msalength; i++) {
        for (j=i+1; j<=msalength; j++) {
            mi = MI(parsedmsa, msafreqs, i, j);
            mi_matrix[i][j] = mi;
	    mi_matrix[j][i] = mi;
        }
    }

    cnpy::npy_save("array_cpp.npy", &mi_matrix[0][0], {minlength, minlength}, "w"); 
}
