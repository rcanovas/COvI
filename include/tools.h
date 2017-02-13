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

/*! \file tools.h
    \brief tools.h contains implementations of different stats
            used in our programs.
    \author Rodrigo Canovas
*/

#ifndef C_TOOLS_H
#define C_TOOLS_H

#include <time.h>
#include <sdsl/int_vector.hpp>

namespace covil {

//STATS Functions

    //Get memory usage by the program so far
    void
    print_proc_statm(int pid) {
        int ret;
        std::stringstream ss;
        ss << pid;
        std::string operation = "cat /proc/" + ss.str() + "/statm";
        ret = system(operation.c_str());
        // 1 page is equal to 4 kB
    }


    //Prints length, entropy, and number of words of the text
    void
    get_stats(sdsl::int_vector<> text, size_t separator) {
        double ent = 0;
        typedef typename sdsl::int_vector<>::size_type size_type;
        std::map<size_type, size_type> t_alphabet;
        size_type n = text.size(), n_alph = 0;
        size_type num_words = 0, length = 0;
        size_type min = 300, max = 0, w_s = 0;
        size_type occ = 0;
        for (size_type i = 0; i < n; ++i) {
            occ = t_alphabet[text[i]] + 1;
            t_alphabet[text[i]] = occ;
            if (text[i] == separator){
                if (w_s > max)
                    max = w_s;
                if (w_s < min)
                    min = w_s;
                ++num_words;
                w_s = 0;
            }
            else {
                ++w_s;
                ++length;
            }
        }
        n_alph = t_alphabet.size();
        for (auto it = t_alphabet.begin(); it != t_alphabet.end(); ++it)
            ent += ((it->second * 1.0) / (n * 1.0)) * log2((n * 1.0) / (it->second * 1.0));
        std::cout << "Text Length: " << text.size() << std::endl;
        std::cout << "H0: " << ent << std::endl;
        std::cout << "Number of words: " << num_words << "  Average size: " << (length * 1.0) / num_words
                  << std::endl;
        std::cout << "Min: " << min << "  Max: " << max << std::endl;
    }
}

#endif //C_TOOLS_H
