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


#include <iostream>

#include "../include/covi/covi.h"

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

    start = timer::now();
    srand (time(NULL)); // initialize
    size_t x, y, r;
    for (int j = 0; j < 10000; ++j) {
        x = rand() % (idx.number_of_words()  - 1) + 1;
        y = rand() % (idx.number_of_words()  - 1) + 1;
        r = idx.overlap(x,y);
    }
    stop = timer::now();
    elapsed = stop - start;
    cout << "Overlap(x,y) Time: " << (double)((chrono::duration_cast<chrono::microseconds>(elapsed).count() * 1.0))/10000.0 << " microsec" << endl;

    start = timer::now();
    for (int j = 0; j < 10000; ++j) {
        x = rand() % (idx.number_of_words()  - 1) + 1;
        y = rand() % (idx.number_of_words()  - 1) + 1;
        auto vec = idx.all_overlaps(x, y);
    }
    stop = timer::now();
    elapsed = stop - start;
    cout << "Overlap_all(x,y) Time: " << (double)((chrono::duration_cast<chrono::microseconds>(elapsed).count() * 1.0))/10000.0 << " microsec" << endl;

    start = timer::now();
    for (int j = 0; j < 10000; ++j) {
        x = rand() % (idx.number_of_words()  - 1) + 1;
        auto vec = idx.all_right_overlaps(x);
    }
    stop = timer::now();
    elapsed = stop - start;
    cout << "all_right_overlaps(x) Time: " << (double)((chrono::duration_cast<chrono::milliseconds>(elapsed).count() * 1.0))/10.0 << " millis" << endl;

    start = timer::now();
    for (int j = 0; j < 10000; ++j) {
        x = rand() % (idx.number_of_words()  - 1) + 1;
        auto vec = idx.all_left_overlaps(x);
    }
    stop = timer::now();
    elapsed = stop - start;
    cout << "all_left_overlaps(x) Time: " << (double)((chrono::duration_cast<chrono::milliseconds>(elapsed).count() * 1.0))/10.0 << " millis" << endl;

    start = timer::now();
    size_t d;
    for (int j = 0; j < 100; ++j) {
        auto vec = idx.max_overlap(d);
    }
    stop = timer::now();
    elapsed = stop - start;
    cout << "max_ovelap() Time: " << ((double)((chrono::duration_cast<chrono::milliseconds>(elapsed).count() * 1.0))/5.0) << " millis" << endl;
}


int main(int argc, char* argv[]) {
    if (argc != 2) {
        cout << "Usage: " << argv[0] << " covi_file_name" << endl;
        return 0;
    }
    string file = argv[1];
    test_overlap<covil::covi<>>(file);
    return 0;
}