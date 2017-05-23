/* covi - Compressed Overlap Index Lib
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

#ifndef COVIL_AHOCO_H
#define COVIL_AHOCO_H

#include "trieMin.h"
#include "./../treetop/tree_topology.hpp"

namespace covil {


    template<uint8_t t_width = 8>
    class fullAC {

    private:
        static_assert(t_width <= 64, "width of elements of the covi must be at most 64bits.");

    public:
        typedef typename sdsl::int_vector<>::size_type                  size_type;
        typedef typename sdsl::int_vector_trait<t_width>::value_type    value_type;

    private:
        tree_bp<>       m_tree; //should use louds topology of the final tree
        sdsl::int_vector<> permutation;
        sdsl::int_vector<> failure_links;

        size_type num_nodes;
        size_type max_depth;

    public:

        fullAC() {
            num_nodes =  max_depth = 0;
        }

        fullAC(std::string file_name, size_type separator) {
            using timer = std::chrono::high_resolution_clock;
            auto start = timer::now();
            trieMin<t_width> Trie(file_name, separator, true);
            max_depth = Trie.Mdepth - 1; //erase the separator symbol
            auto stop = timer::now();
            auto elapsed = stop - start;
            std::cout << "Time to create Trie: " << (double)((std::chrono::duration_cast<std::chrono::seconds>(elapsed).count() * 1.0)) << " sec" << std::endl;
            std::cout << "trie was created" << std::endl;
            start = timer::now();
            std::map<size_type, size_type> end_word_to_bp;
            {
                sdsl::int_vector<> flink;
                num_nodes = Trie.num_nodes();
                Trie.get_failure_links(flink);
                stop = timer::now();
                elapsed = stop - start;
                std::cout << "Time to compute Failure Links: " << (double)((std::chrono::duration_cast<std::chrono::seconds>(elapsed).count() * 1.0)) << " sec" << std::endl;
                std::cout << "failure links computed" << std::endl;
                start = timer::now();
                create_components(Trie, flink, end_word_to_bp);
                std::cout << "components finished" << std::endl;
            }
            size_type number_of_words = Trie.permutation.size(), position_trie;
            permutation = sdsl::int_vector<>(number_of_words, 0, (uint8_t)(sdsl::bits::hi(2 * num_nodes) + 1));
            for (size_type i = 0; i < number_of_words; ++i) {
                position_trie = Trie.permutation[i];
                permutation[i] = end_word_to_bp[position_trie];
            }
            stop = timer::now();
            elapsed = stop - start;
            std::cout << "Time to create Aho-Corasick: " << (double)((std::chrono::duration_cast<std::chrono::seconds>(elapsed).count() * 1.0)) << " sec" << std::endl;
        }

        size_type
        serialize(std::ostream& out, sdsl::structure_tree_node* v=nullptr, std::string name="") const {
            sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
            size_type written_bytes = 0, aux_size;
            written_bytes += write_member(num_nodes, out, child, "number of nodes");
            written_bytes += write_member(max_depth, out, child, "max depth");
            aux_size = permutation.serialize(out, child, "permutation");
            std::cout << "Permutation: " << aux_size << " bytes" << std::endl;
            written_bytes += aux_size;
            aux_size = failure_links.serialize(out, child, "failure links");
            std::cout << "FLinks: " << aux_size << " bytes  <--> " << aux_size * 1.0 / num_nodes <<  "n bytes" << std::endl;
            written_bytes += aux_size;
            aux_size = m_tree.serialize(out, child, "tree bp");
            std::cout << "Tree: " << aux_size << " bytes  <--> " << aux_size * 1.0 / num_nodes <<  "n bytes" << std::endl;
            written_bytes += aux_size;
            sdsl::structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Load from a stream.
        void
        load(std::istream& in) {
            sdsl::read_member(num_nodes, in);
            sdsl::read_member(max_depth, in);
            permutation.load(in);
            failure_links.load(in);
            m_tree.load(in);
        }

        size_type
        number_of_words() {
            return permutation.size();
        }

        //! Returns the maximum right overlap between x and y (x->y)
        // x and y are the positions of the words in the original file
        size_type
        max_ov(size_type x, size_type y) {
            typename tree_bp<>::node_type node_x, node_y;
            size_type t_y, t_y_depth, t_x, t_x_depth;
            node_y = m_tree.node_select(permutation[y]);
            node_y = m_tree.parent(node_y);
            t_y = m_tree.id(node_y);
            t_y_depth = m_tree.depth(node_y);
            t_x = permutation[x];
            t_x = failure_links[t_x];
            node_x = m_tree.node_select(t_x);
            t_x_depth = m_tree.depth(node_x);
            while (t_x_depth != 0) {
                if (t_x == t_y)
                    return t_x_depth;
                else if (t_x_depth >= t_y_depth) {
                    t_x = failure_links[t_x];
                    node_x = m_tree.node_select(t_x);
                    t_x_depth = m_tree.depth(node_x);
                }
                else {
                    node_y = m_tree.parent(node_y);
                    t_y = m_tree.id(node_y);
                    t_y_depth = t_y_depth - 1;
                    if (t_y == 0) // pointing to the root
                        break;
                }
            }
            return 0;
        }

        //! Returns all the right overlaps between x and y (x->y)
        // x and y are the positions of the words in the original file
        std::vector<size_type>
        correlation(size_type x, size_type y) {
            std::vector<size_type> ov;
            typename tree_bp<>::node_type node_x, node_y;
            size_type t_y, t_y_depth, t_x, t_x_depth;
            node_y = m_tree.node_select(permutation[y]);
            node_y = m_tree.parent(node_y);
            t_y = m_tree.id(node_y);
            t_y_depth = m_tree.depth(node_y);
            t_x = permutation[x];
            t_x = failure_links[t_x];
            node_x = m_tree.node_select(t_x);
            t_x_depth = m_tree.depth(node_x);
            while (t_x_depth != 0) {
                if (t_x == t_y)
                    ov.push_back(t_x_depth);
                if (t_x_depth >= t_y_depth) {
                    t_x = failure_links[t_x];
                    node_x = m_tree.node_select(t_x);
                    t_x_depth = m_tree.depth(node_x);
                }
                else {
                    node_y = m_tree.parent(node_y);
                    t_y = m_tree.id(node_y);
                    t_y_depth = t_y_depth - 1;
                    if (t_y == 0) // pointing to the root, no more overlaps
                        break;
                }
            }
            return ov;
        }

        //! Returns an array containing all the pairs (y, length) such that
        // the maximum right overlap of the form x->y is greater or equal
        // to k
        std::vector<std::pair<size_type, size_type>>
        right_overlaps(size_type x, size_type k) {
            typename tree_louds<>::node_type tree_node;
            size_type n = permutation.size();
            std::vector<std::pair<size_type, size_type>> ov;
            sdsl::int_vector<> num(num_nodes, (size_type) -1, (uint8_t) (sdsl::bits::hi(max_depth) + 1));
            num[0] = 0;
            size_type t_y, length;
            size_type t_x = permutation[x], d = 0;
            t_x = failure_links[t_x];
            while (t_x != 0) {
                tree_node = m_tree.node_select(t_x);
                d = m_tree.depth(tree_node);
                num[t_x] = d;
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

        //! Returns an array containing all the pairs (y, length) such that
        // the maximum left overlap of the form y->x is greater or equal
        // to k
        std::vector<std::pair<size_type, size_type>>
        left_overlaps(size_type x, size_type k) {
            typename tree_bp<>::node_type tree_node;
            size_type n = permutation.size();
            std::vector<std::pair<size_type, size_type>> ov;
            sdsl::int_vector<> num(num_nodes, (size_type) -1, (uint8_t) (sdsl::bits::hi(max_depth) + 1));
            num[0] = 0;
            size_type t_y, t_x, length;
            tree_node = m_tree.parent(m_tree.node_select(permutation[x]));
            t_x = m_tree.id(tree_node); //parent of x
            while (t_x != 0) {
                num[t_x] = m_tree.depth(tree_node);
                tree_node = m_tree.parent(tree_node);
                t_x = m_tree.id(tree_node); //parent of x
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
            typename tree_bp<>::node_type tree_node;
            size_type n = permutation.size(), d = 0;
            std::vector<size_type> ov(n);
            sdsl::int_vector<> num(num_nodes, (size_type) -1, (uint8_t) (sdsl::bits::hi(max_depth) + 1));
            num[0] = 0;
            size_type t_y;
            size_type t_x = permutation[x];
            t_x = failure_links[t_x];
            while (t_x != 0) {
                tree_node = m_tree.node_select(t_x);
                d = m_tree.depth(tree_node);
                num[t_x] = d;
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


        //! Returns an array containing the length of the maximum left
        //overlap of each word with the x-th word
        std::vector<size_type>
        all_left_overlaps(size_type x) {
            typename tree_bp<>::node_type tree_node;
            size_type n = permutation.size();
            std::vector<size_type> ov(n);
            sdsl::int_vector<> num(num_nodes, (size_type) -1, (uint8_t) (sdsl::bits::hi(max_depth) + 1));
            num[0] = 0;
            size_type t_y, t_x;
            tree_node = m_tree.parent(m_tree.node_select(permutation[x]));
            t_x = m_tree.id(tree_node); //parent of x
            while (t_x != 0) {
                num[t_x] = m_tree.depth(tree_node);
                tree_node = m_tree.parent(tree_node);
                t_x = m_tree.id(tree_node);
            }
            for (size_type i = 0; i < n; i++) {
                t_y = permutation[i];
                t_y = failure_links[t_y];
                ov[i] = ove_left_num(t_y, num);
            }
            return ov;
        }

        //! Returns the set of word index such that the overlap between them is equal
        // to the maximum possible overlap between all the words of the set.
        std::vector<size_type>
        global_max_overlap(size_type &max_d) {
            std::vector<size_type> max_ov_nodes;
            typename tree_bp<>::node_type node_x;
            sdsl::int_vector<> candidate(num_nodes, 0, 1);
            size_type n = permutation.size();
            size_type t_x, d, min_d = (size_type)-1;
            max_d = 0;
            for (size_type i = 0; i < n; i++) {
                t_x = failure_links[permutation[i]];
                candidate[t_x] = 1;
                node_x = m_tree.node_select(t_x);
                d = m_tree.depth(node_x);
                if (d <  min_d)
                    min_d = d;
            }
            for (size_type i = 0; i < n; i++) {
                node_x = m_tree.node_select(permutation[i]);
                node_x = m_tree.parent(node_x);
                d = m_tree.depth(node_x);
                while (d >= min_d) {
                    t_x = m_tree.id(node_x);
                    if (candidate[t_x]) { //is the maximum overlap with at least one node
                        if (d >= max_d) {
                            if (d > max_d) {
                                max_d = d;
                                min_d = d;
                                max_ov_nodes.clear();
                            }
                            max_ov_nodes.push_back(i);
                        }
                    }
                    node_x = m_tree.parent(node_x);
                    --d;
                }
            }
            return max_ov_nodes;
        }

    private:

        void
        create_components(trieMin<t_width> &Trie, sdsl::int_vector<> &flink,
                          std::map<size_type, size_type> &end_word_to_bp) {
            sdsl::bit_vector bit_sequence(2 * num_nodes, 0);
            failure_links = sdsl::int_vector<>(num_nodes, 0, (uint8_t) (sdsl::bits::hi(flink.size()) + 1));
            size_type pos_sequence = 0, current_node = 0, id = 0;
            in_order_check(Trie, current_node, pos_sequence, id,
                           flink, bit_sequence, end_word_to_bp);
            std::cout << "id: " << id << "  pos_seq: " << pos_sequence << std::endl;
            for (size_type i = 0; i <= id; ++i)   //change failure link pointers
                failure_links[i] = flink[failure_links[i]];
            tree_bp<> bp(bit_sequence, sdsl::int_vector<>()); //empty labels
            m_tree.swap(bp);
        }

        void
        in_order_check(trieMin<t_width> &Trie, size_type current_node, size_type &pos_sequence,
                       size_type &id, sdsl::int_vector<> &flink,
                       sdsl::bit_vector &bit_sequence,
                       std::map<size_type, size_type> &end_word_to_bp) {
            bit_sequence[pos_sequence] = 1; //open parenthesis
            ++ pos_sequence;
            failure_links[id] = flink[current_node];
            flink[current_node] = id;
            if (Trie.get_child(Trie.Mseparator, current_node) != 0)  //this node represent a word
                end_word_to_bp[current_node] = id; //id given by bp
            if (!Trie.is_leaf(current_node)) { //check the children
                ++ current_node;
                while (current_node != 0) {
                    ++ id;
                    in_order_check(Trie, current_node, pos_sequence, id,
                                   flink, bit_sequence, end_word_to_bp);
                    current_node = Trie.get_neighbour(current_node);
                }
            }
            ++ pos_sequence; //close parenthesis
        }

        size_type
        ove_right_num(size_type y, sdsl::int_vector<> &num, size_type k = 0) {
            typename tree_bp<>::node_type tree_node;
            if (num[y] < max_depth) //had been computed before
                return num[y];
            else {
                tree_node = m_tree.parent(m_tree.node_select(y));
                y = m_tree.id(tree_node); //get the parent of y
                if (y == 0 or m_tree.depth(tree_node) < k)
                    num[y] = 0;
                else
                    num[y] = ove_right_num(y, num);
            }
            return num[y];
        }


        size_type
        ove_left_num(size_type y, sdsl::int_vector<> &num, size_type k = 0) {
            typename tree_bp<>::node_type tree_node;
            if (num[y] < max_depth) //had been computed before
                return num[y];
            else {
                tree_node = m_tree.node_select(failure_links[y]);
                if (failure_links[y] == 0 or m_tree.depth(tree_node) < k) //points to the root
                    num[y] = 0;
                else
                    num[y] = ove_left_num(failure_links[y], num);
            }
            return num[y];
        }

    };
}

#endif
