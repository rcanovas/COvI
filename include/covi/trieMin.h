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

#ifndef COVIL_TRIEMIN_H
#define COVIL_TRIEMIN_H


#include "../tools.h"
#include <vector>
#include <stack>
#include <map>
#include <queue>


namespace covil {


    template<uint8_t t_width = 8>
    class trieMin {
    private:
        static_assert(t_width <= 64, "width of elements of the trie must be at most 64bits");

    public:
        typedef typename sdsl::int_vector<>::size_type                  size_type;
        typedef typename sdsl::int_vector_trait<t_width>::value_type    value_type;


    private:
        sdsl::int_vector<> m_leafs;
        sdsl::int_vector<> m_neighbour;
        sdsl::int_vector<> m_alphabet;
        size_type sigma; //size of the alphabet
        size_type num_words;
        size_type max_depth;
        value_type m_separator;
        std::vector<size_type> m_permutation;


    public:
        const sdsl::int_vector<> &leafs = m_leafs;
        const sdsl::int_vector<> &neighbour = m_neighbour;
        const size_type &Mdepth = max_depth;
        const value_type &Mseparator = m_separator;
        const std::vector<size_type>  &permutation = m_permutation;


        //! Empty Constructor
        trieMin() {}

        //! Constructor
        trieMin(std::string file_name, size_type separator, bool perm = false) {
            std::ifstream f_in(file_name, std::ios::in | std::ios::binary);
            if(!f_in) {
                std::cerr << "Failed to open file " << file_name;
                exit(1);
            }
            m_separator = separator;
            create_trie(f_in, perm);
        }

        bool
        is_leaf(size_type node) {
            return (bool)m_leafs[node];
        }

        size_type
        first_child(size_type node) {
            if (is_leaf(node))
                return 0;
            else
                return node + 1;
        }

        size_type
        get_neighbour(size_type node) {
            return m_neighbour[node];
        }

         //! Returns the position of the child of "node" whose edge starts with the value symb
        size_type
        get_child(value_type symb, size_type node) {
            size_type child;
            if (m_leafs[node])
                return 0; //is a leaf
            else {
                child = node + 1; //first child
                while (child != 0) {
                    if (m_alphabet[child - 1] == symb)
                        return child;
                    child = m_neighbour[child];
                }
            }
            return 0;
        }

        void
        get_failure_links_and_ehog_nodes(sdsl::int_vector<> &failure_link,
                                         sdsl::int_vector<> &marked_nodes, size_type &num_marked) {
            sdsl::int_vector<> f_links(m_neighbour.size(), 0, (uint8_t)(sdsl::bits::hi(m_neighbour.size()) + 1));
            sdsl::int_vector<> mark(m_neighbour.size(), 0, 1);
            std::queue<size_type> to_check;
            std::queue<size_type> parent;
            size_type current_node, current_parent, pos_label = 0, child, fl;
            value_type symb = 0;
            mark[0] = 1;
            num_marked = 1;
            to_check.push(0);
            parent.push(0);
            while (!to_check.empty()) {
                current_node = to_check.front();
                current_parent = parent.front();
                to_check.pop();
                parent.pop();
                //compute the failure link of the current node
                if (current_parent == 0)
                        f_links[current_node] = 0;
                else {
                    symb = (value_type)m_alphabet[current_node - 1];
                    fl = f_links[current_parent];
                    while (fl != 0) {
                        child = get_child(symb, fl);
                        if (child != 0) {
                            fl = child;
                            break;
                        }
                        fl = f_links[fl];
                    }
                    if (fl == 0) //check the root
                        fl = get_child(symb, 0);
                    f_links[current_node] = fl;
                }
                if (!is_leaf(current_node)) { //is not a leaf -> push children info to the queue
                    child = current_node + 1; // first child
                    while (child != 0) {
                        to_check.push(child);
                        parent.push(current_node);
                        child = m_neighbour[child];
                    }
                }
                else {  // current node is a leaf--> mark nodes now
                    fl = current_parent; //given how we are storing the words in the trie
                    while (mark[fl] != 1) {
                        mark[fl] = 1;
                        ++num_marked;
                        fl = f_links[fl];
                    }

                }
            }
            f_links.swap(failure_link);
            mark.swap(marked_nodes);

        }

        void
        get_louds(sdsl::int_vector<> &louds_bitvector, sdsl::int_vector<> &alphabet) {
            sdsl::int_vector<> trie_louds(2 * m_neighbour.size() + 1, 0, 1);
            sdsl::int_vector<> trie_alphabet(m_neighbour.size() - 1, 0, (int8_t)(sdsl::bits::hi(sigma) + 1));
            trie_louds[1] = 1;
            std::queue<size_type> to_check;
            std::queue<value_type> label;
            size_type current_node, pos_louds = 2, pos_label = 0;
            to_check.push(0);
            while (!to_check.empty()) {
                current_node = to_check.front();
                if (m_leafs[current_node]) { //is a leaf
                    trie_alphabet[pos_label] = label.front();
                    pos_label++;
                    label.pop();
                }
                else {
                    current_node += 1; // first child
                    pos_louds ++;
                    to_check.push(current_node);
                    label.push((value_type)m_alphabet[current_node - 1]);
                    while (m_neighbour[current_node] != 0) {
                        current_node = m_neighbour[current_node];
                        to_check.push(current_node);
                        label.push((value_type)m_alphabet[current_node - 1]);
                        pos_louds ++;
                    }

                }
                trie_louds[pos_louds] = 1;
                pos_louds ++;
                to_check.pop();
            }
            alphabet.swap(trie_alphabet);
            louds_bitvector.swap(trie_louds);
        }


    private:

        void
        create_trie(std::istream& in, bool perm) {
            size_type last_pos = 0, current_pos = 0, final_size = 0;
            std::vector<size_type> tmp_neigh;
            std::vector<bool> tmp_leafs;
            std::vector<value_type> alphabet;
            std::map<value_type, bool> letters;
            value_type tmp_value;
            size_type child;
            tmp_neigh.push_back(0);
            tmp_leafs.push_back(0);
            num_words = max_depth = 0;
            size_type d = 0;
            char c;
            while (in.get(c)) {
                letters[c] = true;
                ++d;
                if (current_pos == tmp_neigh.size() - 1) { //no children to check
                    tmp_neigh.push_back(0);
                    tmp_leafs.push_back(0);
                    alphabet.push_back((value_type)c);
                    ++current_pos;
                }
                else { //check children
                    child = current_pos + 1;
                    tmp_value = alphabet[child - 1];
                    while (tmp_value != (value_type)c) {
                        if (tmp_neigh[child] == 0) { //not more child
                            tmp_neigh.push_back(0);
                            tmp_leafs.push_back(0);
                            tmp_neigh[child] = tmp_neigh.size() - 1;
                            child = tmp_neigh[child];
                            alphabet.push_back((value_type)c);
                            break;
                        }
                        else {
                            child = tmp_neigh[child];
                            tmp_value = alphabet[child - 1];
                        }
                    }
                    current_pos = child;
                }
                if ((value_type)c == m_separator) { // check  word condition
                    tmp_leafs[current_pos] = 1;
                    if (perm)
                        m_permutation.push_back(last_pos); //the parent of this leaf is a word
                    ++num_words;
                    if (d > max_depth)
                        max_depth = d;
                    d = 0;
                    current_pos = 0;
                }
                last_pos = current_pos;
            }
            sigma = letters.size();
            letters.clear();
            std::cout << "Number of Words: " << num_words << std::endl;
            final_size = tmp_neigh.size();
            std::cout << "Number of nodes in the Trie: " << final_size << std::endl;
            m_leafs = sdsl::int_vector<>(final_size, 0, 1);
            for (size_type i = 0; i < final_size; ++i)
                m_leafs[i] = tmp_leafs[i];
            tmp_leafs.clear();
            m_neighbour = sdsl::int_vector<>(final_size, 0, (uint8_t)(sdsl::bits::hi(final_size) + 1));
            for (size_type i = 0; i < final_size; ++i)
                m_neighbour[i] = tmp_neigh[i];
            tmp_neigh.clear();
            m_alphabet = sdsl::int_vector<>(final_size - 1, 0, t_width);
            for (size_type i = 0; i < final_size - 1; ++i)
                m_alphabet[i] = alphabet[i];
            alphabet.clear();
        }

    };

}

#endif //COVIL_TRIEMIN_H
