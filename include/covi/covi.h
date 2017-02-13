/* covil - Compressed Overlap Index Lib
    Copyright (C)2017 Rodrigo Canovas
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

//For the moment only works with input text of size 8 bits per value

#ifndef COVIL_COVI_H
#define COVIL_COVI_H

#include "trieMin.h"
#include "./../treetop/tree_topology.hpp"

namespace covil {


    template<uint8_t t_width = 8>
    class covi {

    private:
        static_assert(t_width <= 64, "width of elements of the covi must be at most 64bits.");

    public:
        typedef typename sdsl::int_vector<>::size_type                  size_type;
        typedef typename sdsl::int_vector_trait<t_width>::value_type    value_type;

    private:
        tree_louds<>       m_tree; //louds topology of the final tree
        sdsl::int_vector<> permutation;
        sdsl::int_vector<> failure_links;
        sdsl::int_vector<> depths;

        size_type num_nodes;
        size_type max_depth;

    public:

        covi() {
            num_nodes =  max_depth = 0;
        }

        covi(std::string file_name, size_type separator) {
            using timer = std::chrono::high_resolution_clock;
            auto start = timer::now();
            trieMin<t_width> Trie(file_name, separator, true);
            max_depth = Trie.Mdepth - 1; //erase the separator symbol
            auto stop = timer::now();
            auto elapsed = stop - start;
            std::cout << "Time to create Trie: " << (double)((std::chrono::duration_cast<std::chrono::seconds>(elapsed).count() * 1.0)) << " sec" << std::endl;
            std::cout << "trie was created" << std::endl;
            start = timer::now();
            std::map<size_type, size_type> end_word_to_louds;
            {
                sdsl::int_vector<> flink;
                sdsl::int_vector<> mark_nodes;
                Trie.get_failure_links_and_ehog_nodes(flink, mark_nodes, num_nodes);
                stop = timer::now();
                elapsed = stop - start;
                std::cout << "Time to compute Failure Links: " << (double)((std::chrono::duration_cast<std::chrono::seconds>(elapsed).count() * 1.0)) << " sec" << std::endl;
                std::cout << "failure link and mark nodes computed" << std::endl;
                start = timer::now();
                create_components(Trie, mark_nodes, flink, end_word_to_louds);
                std::cout << "components finished" << std::endl;
            }
            size_type number_of_words = Trie.permutation.size(), position_trie;
            permutation = sdsl::int_vector<>(number_of_words, 0, (uint8_t)(sdsl::bits::hi(2 * num_nodes + 1) + 1));
            for (size_type i = 0; i < number_of_words; ++i) {
                position_trie = Trie.permutation[i];
                permutation[i] = end_word_to_louds[position_trie];
            }
            stop = timer::now();
            elapsed = stop - start;
            std::cout << "Time to create COvI: " << (double)((std::chrono::duration_cast<std::chrono::seconds>(elapsed).count() * 1.0)) << " sec" << std::endl;
        }

        size_type
        serialize(std::ostream& out, sdsl::structure_tree_node* v=nullptr, std::string name="") const {
            sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += write_member(num_nodes, out, child, "number of nodes");
            written_bytes += permutation.serialize(out, child, "permutation");
            written_bytes += failure_links.serialize(out, child, "failure links");
            written_bytes += depths.serialize(out, child, "depths");
            written_bytes += m_tree.serialize(out, child, "louds tree");
            sdsl::structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Load from a stream.
        void
        load(std::istream& in) {
            sdsl::read_member(num_nodes, in);
            permutation.load(in);
            failure_links.load(in);
            depths.load(in);
            m_tree.load(in);
        }

        size_type
        number_of_words() {
            return permutation.size();
        }

        //! Returns the maximum right overlap between x and y (x->y)
        // x and y are the positions of the words in the original file
        size_type
        overlap(size_type x, size_type y) {
            typename tree_louds<>::node_type tree_node;
            tree_node = m_tree.node_select(permutation[y]);
            tree_node = m_tree.parent(tree_node);
            size_type t_y = m_tree.id(tree_node);
            size_type t_x = permutation[x];
            t_x = failure_links[t_x];
            size_type t_x_depth = depths[t_x];
            size_type t_y_depth = depths[t_y];
            while (t_x_depth != 0) {
                if (t_x == t_y)
                    return t_x_depth;
                else if (t_x_depth >= t_y_depth) {
                    t_x = failure_links[t_x];
                    t_x_depth = depths[t_x];
                }
                else {
                    tree_node = m_tree.parent(tree_node);
                    t_y = m_tree.id(tree_node);
                    t_y_depth = depths[t_y];
                    if (t_y == 0) // pointing to the root
                        break;
                }
            }
            return 0;
        }

         //! Returns all the right overlaps between x and y (x->y)
        // x and y are the positions of the words in the original file
        std::vector<size_type>
        all_overlaps(size_type x, size_type y) {
            std::vector<size_type> ov;
            typename tree_louds<>::node_type tree_node;
            tree_node = m_tree.node_select(permutation[y]);
            tree_node = m_tree.parent(tree_node);
            size_type t_y = m_tree.id(tree_node);
            size_type t_x = permutation[x];
            t_x = failure_links[t_x];
            size_type t_x_depth = depths[t_x];
            size_type t_y_depth = depths[t_y];
            while (t_x_depth != 0) {
                if (t_x == t_y)
                    ov.push_back(t_x_depth);
                if (t_x_depth >= t_y_depth) {
                    t_x = failure_links[t_x];
                    t_x_depth = depths[t_x];
                }
                else {
                    tree_node = m_tree.parent(tree_node);
                    t_y = m_tree.id(tree_node);
                    t_y_depth = depths[t_y];
                    if (t_y == 0) // pointing to the root, no more overlaps
                        break;
                }
            }
            return ov;
        }


        std::vector<std::pair<size_type, size_type>>
        right_overlaps(size_type x, size_type k) {
            typename tree_louds<>::node_type tree_node;
            size_type n = permutation.size();
            std::vector<std::pair<size_type, size_type>> ov;
            sdsl::int_vector<> num(num_nodes, (size_type)-1, (uint8_t)(sdsl::bits::hi(max_depth) + 1));
            num[0] = 0;
            size_type t_y, length;
            size_type t_x = permutation[x];
            t_x = failure_links[t_x];
            while (t_x != 0) {
                num[t_x] = depths[t_x];
                t_x = failure_links[t_x];
            }
            for (size_type i = 0; i < n; i++) {
                tree_node = m_tree.node_select(permutation[i]);
                tree_node = m_tree.parent(tree_node);
                t_y = m_tree.id(tree_node);
                length = ove_right_num(t_y, num, k);
                if (length >= k)
                    ov.push_back(std::make_pair(i, length));
            }
            return ov;
        }

        std::vector<std::pair<size_type, size_type>>
        left_overlaps(size_type x, size_type k) {
            typename tree_louds<>::node_type tree_node;
            size_type n = permutation.size();
            std::vector<std::pair<size_type, size_type>> ov;
            sdsl::int_vector<> num(num_nodes, (size_type)-1, (uint8_t)(sdsl::bits::hi(max_depth) + 1));
            num[0] = 0;
            size_type t_y, t_x, length;
            t_x = m_tree.id(m_tree.parent(m_tree.node_select(permutation[x]))); //parent of x
            while (t_x != 0) {
                num[t_x] = depths[t_x];
                t_x = m_tree.id(m_tree.parent(m_tree.node_select(t_x))); //parent of x
            }
            for (size_type i = 0; i < n; i++) {
                t_y = permutation[i];
                t_y = failure_links[t_y];
                length = ove_left_num(t_y, num, k);
                if (length >= k)
                    ov.push_back(std::make_pair(i, length));
            }
            return ov;
        }

        //! Returns an array containing the length of the maximum right
        //overlap of each word with the x-th word (of the form x->y)
        std::vector<size_type>
        all_right_overlaps(size_type x) {
            typename tree_louds<>::node_type tree_node;
            size_type n = permutation.size();
            std::vector<size_type> ov(n);
            sdsl::int_vector<> num(num_nodes, (size_type) -1, (uint8_t) (sdsl::bits::hi(max_depth) + 1));
            num[0] = 0;
            size_type t_y;
            size_type t_x = permutation[x];
            t_x = failure_links[t_x];
            while (t_x != 0) {
                num[t_x] = depths[t_x];
                t_x = failure_links[t_x];
            }
            for (size_type i = 0; i < n; i++) {
                tree_node = m_tree.node_select(permutation[i]);
                tree_node = m_tree.parent(tree_node);
                t_y = m_tree.id(tree_node);
                ov[i] = ove_right_num(t_y, num);
            }
            return ov;
        }

        //! Returns an aray containing the length of the maximum left
        //overlap of each word with the x-th word
        std::vector<size_type>
        all_left_overlaps(size_type x) {
            typename tree_louds<>::node_type tree_node;
            size_type n = permutation.size();
            std::vector<size_type> ov (n);
            sdsl::int_vector<> num(num_nodes, (size_type)-1, (uint8_t)(sdsl::bits::hi(max_depth) + 1));
            num[0] = 0;
            size_type t_y, t_x;
            t_x = m_tree.id(m_tree.parent(m_tree.node_select(permutation[x]))); //parent of x
            while (t_x != 0) {
                num[t_x] =  depths[t_x];
                t_x = m_tree.id(m_tree.parent(m_tree.node_select(t_x)));
            }
            for (size_type i = 0; i < n; i++) {
                t_y = permutation[i];
                t_y = failure_links[t_y];
                ov[i] = ove_left_num(t_y, num);
            }
            return ov;
        }

        std::vector<size_type>
        max_overlap(size_type &max_d) {
            typename tree_louds<>::node_type tree_node;
            sdsl::int_vector<> candidate(num_nodes, 0, 1);
            size_type n = permutation.size(), d;
            max_d = 0;
            std::vector<size_type> max_ov_nodes;
            size_type t_x;
            for (size_type i = 0; i < n; i++) {
                t_x = m_tree.id(m_tree.parent(m_tree.node_select(permutation[i])));
                candidate[t_x] = 1;
            }
            for (size_type i = 0; i < n; i++) {
                t_x = permutation[i];
                t_x = failure_links[t_x];
                if (candidate[t_x]) { //is the maximum overlap with at least one node
                    d = depths[t_x];
                    if (d >= max_d) {
                        if (d > max_d) {
                            max_d = d;
                            max_ov_nodes.clear();
                        }
                        max_ov_nodes.push_back(i);
                    }
                }
            }
            return max_ov_nodes;
        }

    private:

        void
        create_components(trieMin<t_width> Trie,
                          sdsl::int_vector<> &mark_nodes, sdsl::int_vector<> &flink,
                          std::map<size_type, size_type> &end_word_to_louds){
            std::cout << "Number of nodes of COvI: " << num_nodes << std::endl;
            sdsl::bit_vector bit_sequence(2 * num_nodes + 1, 0);
            depths = sdsl::int_vector<>(num_nodes, 0, (uint8_t)(sdsl::bits::hi(Trie.Mdepth) + 1));
            failure_links = sdsl::int_vector<>(num_nodes, 0, (uint8_t)(sdsl::bits::hi(flink.size()) + 1)); //could be improved
            bit_sequence[1] = 1;
            std::queue<size_type> to_check;
            std::queue<size_type> tmp_depth;
            size_type current_node, pos_sequence = 2, pos_depth = 0, tmp_fl;
            to_check.push(0);
            tmp_depth.push(0);
            while (!to_check.empty()) {
                current_node = to_check.front();
                depths[pos_depth] = tmp_depth.front();
                tmp_depth.pop();
                to_check.pop();
                failure_links[pos_depth] = flink[current_node];
                flink[current_node] = pos_depth;
                if (Trie.get_child(Trie.Mseparator, current_node) != 0)  //this node represent a work
                    end_word_to_louds[current_node] = pos_depth; //id are given by louds order
                if (!Trie.is_leaf(current_node)) {   //is not a leaf
                    ++current_node;
                    while (current_node != 0) {
                        pileup_children(current_node, depths[pos_depth] + 1, Trie, pos_sequence, mark_nodes,
                                        to_check, tmp_depth);
                        current_node = Trie.get_neighbour(current_node);
                    }
                }
                bit_sequence[pos_sequence] = 1;
                ++pos_sequence;
                ++pos_depth;
            }
            for (size_type i = 0; i < pos_depth; ++i)   //change failure link pointers
                failure_links[i] = flink[failure_links[i]];
            tree_louds<> louds(bit_sequence, sdsl::int_vector<>());
            m_tree.swap(louds);
        }

        void
        pileup_children(size_type current_node, size_type d,
                        trieMin<t_width> &Trie,
                        size_type &pos_sequence, sdsl::int_vector<> &mark_nodes,
                        std::queue<size_type> &to_check, std::queue<size_type> &tmp_depth) {
            if (mark_nodes[current_node]) {
                to_check.push(current_node);
                tmp_depth.push(d);
                pos_sequence = pos_sequence + 1;
                return;
            }
            else {
                if (Trie.is_leaf(current_node))
                    return;
                current_node = current_node + 1;
                while (current_node != 0) {
                    pileup_children(current_node, d+1, Trie, pos_sequence,
                                    mark_nodes, to_check, tmp_depth);
                    current_node = Trie.get_neighbour(current_node);
                }
            }
        }

        //!Returns the maximum right overlap (x->y) and update the num value of
        //each node visited
        size_type
        ove_right_num(size_type y, sdsl::int_vector<> &num, size_type k = 0) {
            if (num[y] < max_depth) //had been computed before
                return num[y];
            else {
                y = m_tree.id(m_tree.parent(m_tree.node_select(y))); //get the parent of y
                if (y == 0 or depths[y] < k)
                    num[y] = 0;
                else
                    num[y] = ove_right_num(y, num);
            }
            return num[y];
        }

        //!Returns the maximum right overlap (x->y) and update the num value of
        //each node visited
        size_type
        ove_left_num(size_type y, sdsl::int_vector<> &num, size_type k = 0) {
            if (num[y] < max_depth) //had been computed before
                return num[y];
            else {
                if (failure_links[y] == 0 or depths[failure_links[y]] < k) //points to the root
                    num[y] = 0;
                else
                    num[y] = ove_left_num(failure_links[y], num);
            }
            return num[y];
        }


    };
}

#endif
