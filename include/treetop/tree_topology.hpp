/*
    Copyright (C)2016 Rodrigo Canovas

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


#ifndef DISD_TREE_TOPOLOGY_HPP
#define DISD_TREE_TOPOLOGY_HPP

namespace covil {
//! A class for the node representation of tree
    struct tree_node {
    public:
        typedef sdsl::bit_vector::size_type size_type;
        size_type nr;             // node number
        size_type pos;            // position in the bitmap

        tree_node(size_type f_nr = 0, size_type f_pos = 0) : nr(f_nr), pos(f_pos) { }

        //! Copy constructor
        tree_node(const tree_node &v) = default;

        //! Move copy constructor
        tree_node(tree_node &&v) = default;

        //! Equality operator.
        bool operator==(const tree_node &v) const {
            return nr == v.nr and pos == v.pos;
        }

        //! Inequality operator.
        bool operator!=(const tree_node &v) const {
            return !(v == *this);
        }

        //! Assignment operator.
        tree_node &operator=(const tree_node &v) = default;

        //! Move assignment
        tree_node &operator=(tree_node &&v) = default;
    };

}

#include "tree_louds.hpp"


#endif //DISD_TREE_TOPOLOGY_HPP
