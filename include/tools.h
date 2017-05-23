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

#ifndef COV_TOOLS_H
#define COV_TOOLS_H

#include <time.h>
#include <sdsl/int_vector.hpp>

#ifndef FACT_RANK
#define FACT_RANK 20
#endif

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

    //computes the accumulate frequency of the values in the set of size n
    //we assumed that the values can not be bigger than n
    sdsl::int_vector<>
    accumulateFre(sdsl::int_vector<> &set) {
        size_t max_value = 0, n = set.size();
        sdsl::int_vector<> fre_ac(n, 0, (uint8_t)(sdsl::bits::hi(n) + 1));
        for (size_t i = 0; i < n; ++i) {
            fre_ac[set[i]]++;
            if(max_value < set[i])
                max_value = set[i];
        }
        for(size_t i = 1; i <= max_value; i++)
            fre_ac[i] = fre_ac[i] + fre_ac[i - 1];
        fre_ac.resize(max_value + 1);
        return fre_ac;
    }

    //computes the optimal number of bits used to represent
    // each level of a DAC (direct access code) given a set
    std::vector<uint8_t>
    bits_per_level(sdsl::int_vector<> &set) {
        std::vector<uint8_t> kvalues;
        sdsl::int_vector<> fre_ac = accumulateFre(set);
        size_t max_value = fre_ac.size() - 1;
        size_t list_length = fre_ac[max_value];
        size_t nBits = sdsl::bits::hi(max_value) + 1;
        //This table will contain the size of the best option for store the first x bits
        std::vector<size_t> tableSize(nBits+1);
        std::vector<size_t> tableNLevels(nBits+1);
        std::vector<std::vector<size_t>> tableKvalues(nBits+1);

        size_t maxSize = 0, maxPos = 0;
        size_t posVocInf, posVocSup, currentSize;
        tableSize[0] = 0;
        tableNLevels[0] = 0;
        for (size_t i = 1; i <= nBits; i++) {
            maxSize = (size_t)-1;
            maxPos = 0;
            for (size_t j = 0; j < i; j++) {
                if (i == nBits)
                    posVocInf = 0;
                else
                    posVocInf = 1 << (nBits - i);
                posVocSup = (1 << (nBits - j));
                if (posVocSup >= max_value)
                    posVocSup = max_value;
                if (j == 0)
                    currentSize = tableSize[j] +
                            ((size_t) (fre_ac[max_value] - fre_ac[posVocInf])) * ((i - j));
                else
                    currentSize = tableSize[j] +
                            ((size_t) (fre_ac[max_value] - fre_ac[posVocInf])) * ((i - j) + 1) +
                            (fre_ac[max_value] - fre_ac[posVocInf]) / FACT_RANK;
                if (maxSize > currentSize) {
                    maxSize = currentSize;
                    maxPos = j;
                }
            }
            tableSize[i] = maxSize;
            tableNLevels[i] = tableNLevels[maxPos] + 1;
            tableKvalues[i] = std::vector<size_t>(tableNLevels[i]);
            for (size_t j = 0; j < tableNLevels[i] - 1; j++)
                tableKvalues[i][j] = tableKvalues[maxPos][j];
            tableKvalues[i][tableNLevels[i] - 1] = i - maxPos;
        }

        return kvalues;
    }

}

#endif //COV_TOOLS_H
