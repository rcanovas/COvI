/* covil - Compressed Overlap Index Lib
    Copyright (C)2016-2017 Rodrigo Canovas
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/ .
*/

//Example test of how to use a COvI file

#include <random>
#include <iostream>

#include "../include/covi/covi.h"
#include "../include/covi/fullAC.h"


using namespace std;


template<class idx_type>
void
test_overlap(string file) {
    using timer = std::chrono::high_resolution_clock;
    auto start = timer::now();
    idx_type idx;
    std::ifstream f_in(file, std::ios::in | std::ios::binary);
    if(!f_in) {
        std::cerr << "Failed to open file " << file;
        exit(1);
    }
    idx.load(f_in);
    auto stop = timer::now();
    auto elapsed = stop - start;
    cout << "Loading Time: " << (double)((chrono::duration_cast<chrono::milliseconds>(elapsed).count() * 1.0)) << " ms" << endl;

    std::random_device rd;
    std::mt19937 mt(rd());
    //if you want to always generate the same random numbers whenever the test is run use this:
    //std::default_random_engine generator;

    std::uniform_real_distribution<double> dist(1.0, 1.0 * idx.number_of_words() -1);

    start = timer::now();
    srand (time(NULL)); // initialize
    size_t x, y, r;
    for (int j = 0; j < 10000; ++j) {
        x = (uint64_t)dist(mt) + 1;
        y = (uint64_t)dist(mt) + 1;
        //x = (uint64_t)dist(generator) + 1;
        //y = (uint64_t)dist(generator) + 1;
        r = idx.max_ov(x,y);
    }
    stop = timer::now();
    elapsed = stop - start;
    cout << "Max-ov(x,y) Time: " << (double)((chrono::duration_cast<chrono::microseconds>(elapsed).count() * 1.0))/10000.0 << " microsec" << endl;

    start = timer::now();
    for (int j = 0; j < 10000; ++j) {
        x = (uint64_t)dist(mt) + 1;
        y = (uint64_t)dist(mt) + 1;
        //x = (uint64_t)dist(generator) + 1;
        //y = (uint64_t)dist(generator) + 1;
        auto vec = idx.correlation(x, y);
    }
    stop = timer::now();
    elapsed = stop - start;
    cout << "Correlation(x,y) Time: " << (double)((chrono::duration_cast<chrono::microseconds>(elapsed).count() * 1.0))/10000.0 << " microsec" << endl;

    start = timer::now();
    for (int j = 0; j < 10; ++j) {  //change to 10 or more
        x = (uint64_t)dist(mt) + 1;
        //x = (uint64_t)dist(generator) + 1;
        auto vec = idx.all_right_overlaps(x);
    }
    stop = timer::now();
    elapsed = stop - start;
    cout << "all_right_overlaps(x) Time: " << (double)((chrono::duration_cast<chrono::milliseconds>(elapsed).count() * 1.0))/10.0 << " millis" << endl;

    start = timer::now();
    for (int j = 0; j < 10; ++j) {
        x = (uint64_t)dist(mt) + 1;
        //x = (uint64_t)dist(generator) + 1;
        auto vec = idx.all_left_overlaps(x);
    }
    stop = timer::now();
    elapsed = stop - start;
    cout << "all_left_overlaps(x) Time: " << (double)((chrono::duration_cast<chrono::milliseconds>(elapsed).count() * 1.0))/10.0 << " millis" << endl;

    start = timer::now();
    size_t d;
    for (int j = 0; j < 2; ++j) {
        std::cout << j << std::endl;
        auto vec = idx.global_max_overlap(d);
        cout << "max_ov size: " << vec.size() << "  d: " << d << endl;
    }
    stop = timer::now();
    elapsed = stop - start;
    cout << "global_max_overlap() Time: " << ((double)((chrono::duration_cast<chrono::milliseconds>(elapsed).count() * 1.0))/2.0) << " millis" << endl;

}


int main(int argc, char* argv[]) {
    if (argc < 2) {
        cout << "Usage: " << argv[0] << " covi_file_name <opt>" << endl;
        cout << "opt: " << endl;
        cout << "-w Index_type. Default = 0" << endl;
        cout << " # | Index" << endl;
        cout << "---+--------------------" << endl;
        cout << " 0 | COvI" << endl;
        cout << " 1 | fullAC" << endl;
        return 0;
    }

    string file = argv[1];
    size_t op = 0;
    uint8_t w = 0;

    int c;
    while ((c = getopt(argc, argv, "w:")) != -1) {
        switch (c) {
            case 'w': w = atoi(optarg); break;
            case '?':
                if (optopt == 'w' )
                    fprintf(stderr, "Option -%c requires an argument.\n", optopt);
                else
                    fprintf(stderr, "Unknown option character `\\x%x'.\n", optopt);
                return 1;
            default:
                abort();
        }
    }


    switch (w) {
        case 0:
            test_overlap<covil::covi<>>(file);
            break;
        case 1:
            test_overlap<covil::fullAC<>>(file);
            break;
        default:
            cout << "index_type must be a value in [0,1]" << endl;
            break;
    }
    return 0;
}
