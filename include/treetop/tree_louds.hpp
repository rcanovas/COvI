/* dictionary succinct data library
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


#ifndef DISD_TREE_LOUDS_HPP
#define DISD_TREE_LOUDS_HPP


#include <sdsl/bit_vectors.hpp>
#include <sdsl/util.hpp>

namespace covil {

    //! A tree class based on the level order unary degree sequence (LOUDS) representation.
    /*!
     * \tparam bit_vec_t  The bit vector representation used for LOUDS.
     * \tparam select_1_t A select_support on 1-bits required for the child(v,i) operation.
     * \tparam select_0_t A select_support on 0-bits required for the parent operation.
     * \tparam t_width    Number of bits used per symbol in the labels
     */
     template<class bit_vec_t = sdsl::bit_vector,
             class select_1_t = typename bit_vec_t::select_1_type,
             class select_0_t = typename bit_vec_t::select_0_type,
             uint8_t t_width = 8>
    class tree_louds {

    public:
        typedef sdsl::bit_vector::size_type                             size_type;
        typedef tree_node                                               node_type;
        typedef bit_vec_t                                               bit_vector_type;
        typedef select_1_t                                              select_1_type;
        typedef select_0_t                                              select_0_type;
        typedef typename sdsl::int_vector_trait<t_width>::value_type    value_type;

    private:
        sdsl::int_vector<>          m_labels;
        bit_vector_type             m_bv;         // bit vector for the LOUDS sequence
        select_1_type               m_bv_select1; // select support for 1-bits on m_bv
        select_0_type               m_bv_select0; // select support for 0-bits on m_bv

    public:
        const bit_vector_type& bv = m_bv;    // const reference to the LOUDS sequence

        //! Default constructor
        tree_louds() {}

        //! Constructor given a bit sequence as parameter.
        // The sequence must be a valid loud representation. That is,
        //  for example, a node with 3 children is represented by 0001. Also
        // the sequence must prepend a 01  which will help to avoid several
        // special cases
        tree_louds(const sdsl::bit_vector sequence, const sdsl::int_vector<> labels);

        tree_louds(const tree_louds& lt) {
            *this = lt;
        }

        tree_louds(tree_louds&& lt) {
            *this = std::move(lt);
        }

        tree_louds& operator=(const tree_louds& lt) {
            if (this != &lt) {
                m_labels = lt.m_labels;
                m_bv = lt.m_bv;
                m_bv_select1 = lt.m_bv_select1;
                m_bv_select1.set_vector(&m_bv);
                m_bv_select0 = lt.m_bv_select0;
                m_bv_select0.set_vector(&m_bv);
            }
            return *this;
        }

        tree_louds& operator=(tree_louds&& lt) {
            if (this != &lt) {
                m_labels = std::move(lt.m_labels);
                m_bv = std::move(lt.m_bv);
                m_bv_select1 = std::move(lt.m_bv_select1);
                m_bv_select1.set_vector(&m_bv);
                m_bv_select0 = std::move(lt.m_bv_select0);
                m_bv_select0.set_vector(&m_bv);
            }
            return *this;
        }

        void
        swap(tree_louds& tl) {
            if (this != &tl) {
                m_labels.swap(tl.m_labels);
                m_bv.swap(tl.m_bv);
                sdsl::util::swap_support(m_bv_select0, tl.m_bv_select0, &m_bv, &(tl.m_bv));
                sdsl::util::swap_support(m_bv_select1, tl.m_bv_select1, &m_bv, &(tl.m_bv));
            }
        }

        size_type
        serialize(std::ostream& out, sdsl::structure_tree_node* v=nullptr, std::string name="") const {
            sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += m_labels.serialize(out, child, "labels");
            written_bytes += m_bv.serialize(out, child, "bitvector");
            written_bytes += m_bv_select1.serialize(out, child, "select1");
            written_bytes += m_bv_select0.serialize(out, child, "select0");
            sdsl::structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Load from a stream.
        void
        load(std::istream& in) {
            m_labels.load(in);
            m_bv.load(in);
            m_bv_select1.load(in, &m_bv);
            m_bv_select0.load(in, &m_bv);
        }


        //! Returns size of the LOUDS bitmap
        size_type
        size() const {
            return m_bv.size();
        }

        //! Returns the number of nodes in the tree.
        size_type
        nodes() const {
            return  (m_bv.size() - 1)/2;
        }

        //! Returns an unique id for each node in [0..size()-1]
        size_type
        id(const node_type& v)const {
            return v.nr;
        }

        //! Returns the root node
        node_type
        root() const {
            return node_type(0, 2);
        }

        //! Returns the parent of a node v or root() if v==root().
        node_type
        parent(const node_type& v) const {
            if (v == root()) {
                return root();
            }
            //v.nr + 1 = rank_1(v.pos - 1) = number of nodes before v
            size_type zero_pos = m_bv_select0(v.nr + 1); //select where the node is represented as a child
            size_type parent_nr = (zero_pos + 1) - v.nr - 2;
            size_type parent_pos = sdsl::util::prev_bit(m_bv ,zero_pos) + 1; //finally compute the prev_1 + 1
            return node_type(parent_nr, parent_pos);
        }

        //! Returns the parent of a node v and the label that get to it, or root() if v==root().
        node_type
        parent_label(const node_type& v, value_type &c) const {
            if (v == root()) {
                return root();
            }
            //v.nr + 1 = rank_1(v.pos - 1) = number of nodes before v
            size_type zero_pos = m_bv_select0(v.nr + 1); //select where the node is represented as a child
            size_type parent_nr = (zero_pos + 1) - v.nr - 2;
            size_type parent_pos = sdsl::util::prev_bit(m_bv, zero_pos) + 1; //finally compute the prev_1 + 1
            c = m_labels[zero_pos - parent_nr - 2];
            return node_type(parent_nr, parent_pos);
        }

        //! Returns the i-child of a node. Otherwise returns the root
        /*!
         * \param v    The parent node.
         * \param i    Index of the child. Indexing starts at 1.
        */
        node_type
        child(const node_type& v, size_type i)const {
            if (degree(v) < i)
                return root();
            size_type pos =  v.pos - v.nr + i - 1; // rank_0(v_pos) + i - 1
            return node_type(pos - 1, m_bv_select1(pos) + 1);
        }

        //!Return the ancestor w of v such that depth(w) = depth(v) − d, or the root
        node_type
        level_ancestor(const node_type& v, size_type delta){
            size_type d = delta;
            node_type an = v;
            while (an != root() and d > 0) {
                an = parent(an);
                --d;
            }
            return an;
        }

        //! Returns the tree depth of the node v
        size_type
        depth(const node_type& v) {
            node_type w = v;
            size_type d = 0;
            while (v != root()) {
                d++;
                w = parent(w);
            }
            return d;
        }

        //! Returns the number of children of v
        size_type
        degree(const node_type& v) const {
            return sdsl::util::next_bit(m_bv, v.pos) - v.pos;
        }

        size_type
        isLeaf(const node_type& v) const{
            return m_bv[v.pos];
        }

        //! Returns the labeled child of a node. Otherwise returns the root
        /*!
         * \param v    The parent node.
         * \param s     Label that lead to the child node.
        */
        node_type
        childLabel(const node_type& v, value_type s) const {
            if (isLeaf(v))
                return root();
            size_type d = degree(v);
            if(d == 0)
                return root();
            //find the first label child
            size_type ini = v.pos - v.nr - 2; //position in the label array (number of zeros - 2)
            size_type j = 0, mid;
            value_type c = 0;
            while (d - j > 5) { //do binary search
                mid = (d + j) / 2;
                c = m_labels[mid + ini];
                if (c == s)
                    return child(v, mid + 1);
                else if (c > s)
                    d = mid;
                else
                    j = mid;
            }
            //for the last 5 is better just do sequential search
            while (j < d) {
                c = m_labels[j + ini];
                if (c == s)
                    return child(v, j + 1);
                if (c > s)
                    break;
                j++;
            }
            return root();
        }

        //! Returns the number of siblings to the left plus 1. That is, which
        // child number is v of its parent. Otherwise, returns 0
        size_type
        child_rank(const node_type& v){
            if (v == root())
                return 0;
            size_type zero_pos = m_bv_select0(v.nr + 1); //select where the node is represented as a child
            zero_pos = zero_pos -  sdsl::util::prev_bit(m_bv ,zero_pos); // compute its local position
            return zero_pos;
        }

        //!Return the node with id i. Otherwise return the root
        node_type
        node_select(size_type i) const {
            if (i >= nodes() or i == 0)
                return root();
            size_type one_pos = m_bv_select1(i + 1); //select the previous node ends
            return node_type(i, one_pos + 1);
        }

        //!Returns the depth of the node in the tree
        size_type
        depth(const node_type& v) const {
            node_type x = v;
            size_type d = 0;
            while (x != root()) {
                x = parent(x);
                ++d;
            }
            return d;
        }

        //! Returns true if v is an ancestor of w
        bool
        ancestor(const node_type& v, const node_type& w) {
            if (v == root())
                return true;
            node_type x = w;
            while (x != root()) {
                if (x == v)
                    return true;
                else if(x.pos < v.pos)
                    return false;
                x = parent(x);
            }
            return false;
        }

        //! Returns the lower common ancestor between v and w;
        node_type
        lca(const node_type& v, const node_type& w) {
            node_type x = v;
            node_type y = w;
            if (x == root() or y == root())
                return root();
            while (x != y) {
                if (x.pos < y.pos)
                    y = parent(y);
                else
                    x = parent(x);
            }
            return x;
        }

        //! Returns the next sibling of v. Otherwise, returns the root
        node_type
        next_sibling(const node_type& v) {
            if (v == root())
                return root();
            size_type zero_pos = m_bv_select0(v.nr + 1); //select where the node is represented as a child
            if (m_bv[zero_pos + 1] == 1) //no more siblings
                return root();
            return node_type(v.nr + 1, m_bv_select1(v.nr + 2) + 1);
        }

        //! Returns the previous sibling of v. Otherwise, returns the root
        node_type
        prev_sibling(const node_type& v) {
            if (v == root())
                return root();
            size_type zero_pos = m_bv_select0(v.nr + 1); //select where the node is represented as a child
            if (m_bv[zero_pos - 1] == 1) //no more siblings
                return root();
            return node_type(v.nr - 1, m_bv_select1(v.nr) + 1);
        }

        size_type
        subtree_size(const node_type& v) {
            return 0;  //no implemented. Could be done using brute force.
        }

    private:


    }; //end class tree_louds

    template<class bit_vec_t, class select_1_t, class select_0_t, uint8_t t_width>
    tree_louds<bit_vec_t, select_1_t, select_0_t, t_width>::tree_louds(const sdsl::bit_vector sequence, const sdsl::int_vector<> labels) {
        m_labels = labels;
        m_bv = sequence;
        sdsl::util::init_support(m_bv_select1, &m_bv);
        sdsl::util::init_support(m_bv_select0, &m_bv);
    }


}// end namespace


#endif //DISD_TREE_LOUDS_HPP
