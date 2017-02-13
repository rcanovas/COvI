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

//Example of how to create the COvI of a file

#include <iostream>

#include "../include/covi/covi.h"


using namespace std;


template<class idx_type>
void
create_index(string file, string out_file, size_t separator) {
    using timer = std::chrono::high_resolution_clock;
    auto start = timer::now();
    auto idx = idx_type(file, separator);
    ofstream out(out_file);
    idx.serialize(out);
    auto stop = timer::now();
    auto elapsed = stop - start;
    cout << "Construction Time: " << (double)((chrono::duration_cast<chrono::milliseconds>(elapsed).count() * 1.0)) << " ms" << endl;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cout << "Usage: " << argv[0] << " file_name <opt>" << endl;
        cout << "opt: " << endl;
        cout << "-o output_name:  String containing the name of the output file. Default file_name.covi" << endl;
        cout << "-s separator:  (ASCII) Value used as separator between words. Default separator = 0" << endl;
        return 0;
    }

    string file = argv[1];
    string out_file = file;
    size_t separator = '\0';

    int c;
	while((c = getopt (argc, argv, "o:w:s:")) != -1){
		switch (c) {
            case 'o': out_file = optarg;  break;
            case 's': separator = atoi(optarg); break;
			case '?': if(optopt == 'o' || optopt == 's')
                        fprintf (stderr, "Option -%c requires an argument.\n", optopt);
                      else
                        fprintf(stderr,"Unknown option character `\\x%x'.\n",	optopt);
                      return 1;
            default:  abort ();
		}
	}

    //create index
    out_file += ".covi";
    create_index<covil::covi<>>(file, out_file, separator);
    return 0;
}