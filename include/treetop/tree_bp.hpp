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


#ifndef DISD_TREE_BP_HPP
#define DISD_TREE_BP_HPP


#include <sdsl/bit_vectors.hpp>
#include <sdsl/bp_support.hpp>
#include <sdsl/util.hpp>

namespace covil {

    //! A tree class based on the balanced parentheses (BP) representation.
    //An opening parenthesis in the balanced parentheses sequence is represented
    // by a 1 in the bit_vector and a closing parenthesis by a 0.
    /*!
     * \tparam bp_support  The balanced parentheses support used
     * \tparam t_width    Number of bits used per symbol in the labels
     */
     template<class bp_support = sdsl::bp_support_sada<>, uint8_t t_width = 8>
    class tree_bp {

    public:
        typedef sdsl::bit_vector::size_type                             size_type;
        typedef tree_node                                               node_type;
        typedef bp_support                                              bp_type;
        typedef typename sdsl::int_vector_trait<t_width>::value_type    value_type;

    private:
        sdsl::int_vector<>          m_labels;
        sdsl::bit_vector            m_sequence;         //  The BP sequence
        bp_type                     m_BP;


    public:
        const sdsl::bit_vector& sequence = m_sequence;    // const reference to the sequence

        //! Default constructor
        tree_bp() {}

        //! Constructor given a bit sequence as parameter.
        tree_bp(const sdsl::bit_vector s, const sdsl::int_vector<> labels);

        tree_bp(const tree_bp& lt) {
            *this = lt;
        }

        tree_bp(tree_bp&& lt) {
            *this = std::move(lt);
        }

        tree_bp& operator=(const tree_bp& lt) {
            if (this != &lt) {
                m_labels = lt.m_labels;
                m_sequence = lt.m_sequence;
                m_BP = lt.m_BP;
                m_BP.set_vector(&m_sequence);
            }
            return *this;
        }

        tree_bp& operator=(tree_bp&& lt) {
            if (this != &lt) {
                m_labels = std::move(lt.m_labels);
                m_sequence = std::move(lt.m_sequence);
                m_BP = std::move(lt.m_BP);
                m_BP.set_vector(&m_sequence);
            }
            return *this;
        }

        void
        swap(tree_bp& tl) {
            if (this != &tl) {
                m_labels.swap(tl.m_labels);
                m_sequence.swap(tl.m_sequence);
                sdsl::util::swap_support(m_BP, tl.m_BP, &m_sequence, &(tl.m_sequence));
            }
        }

        size_type
        serialize(std::ostream& out, sdsl::structure_tree_node* v=nullptr, std::string name="") const {
            sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += m_labels.serialize(out, child, "labels");
            written_bytes += m_sequence.serialize(out, child, "bitvector");
            written_bytes += m_BP.serialize(out, child, "balanced parentheses");
            sdsl::structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Load from a stream.
        void
        load(std::istream& in) {
            m_labels.load(in);
            m_sequence.load(in);
            m_BP.load(in, &m_sequence);
        }


        //! Returns size of the BP bitmap
        size_type
        size() const {
            return m_sequence.size();
        }

        //////////////////////////////////////////////HERE

        //! Returns the number of nodes in the tree.
        size_type
        nodes() const {
            return  m_sequence.size()/2;
        }

        //! Returns an unique id for each node in [0..size()-1]
        size_type
        id(const node_type& v)const {
            return v.nr;
        }

        //! Returns the root node
        node_type
        root() const {
            return node_type(0, 0);
        }

        bool
        isRoot(const node_type& v) const {
            return v == root();
        }

        //! Returns the parent of a node v or root() if v==root().
        node_type
        parent(const node_type& v) const {
            if (v == root()) {
                return root();
            }
            size_type parent_pos = m_BP.enclose(v.pos);
            size_type parent_nr = m_BP.rank(parent_pos) - 1; //rank of 1's until (including) parent_pos
            return node_type(parent_nr, parent_pos);
        }


        //! Returns the parent of a node v and the label that get to it, or root() if v==root().
        node_type
        parent_label(const node_type& v, value_type &c) const {
            if (v == root()) {
                return root();
            }
            size_type parent_pos = m_BP.enclose(v.pos);
            size_type parent_nr = m_BP.rank(parent_pos) - 1; //rank of 0's until parent_pos
            c = m_labels[v.nr - 1];
            return node_type(parent_nr, parent_pos);
        }


        //! Returns the i-child of a node. Otherwise returns the root
        node_type
        child(const node_type& v, size_type i)const {
            if (isLeaf(v))
                return root();
            //get first child
            size_type pos = v.pos + 1;
            node_type child_node(v.nr + 1, pos);
            while (i > 1 and child_node != root()) { //while i>1 compute next sibling
                child_node = next_sibling(child_node);
                -- i;
            }
            return child_node;
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
            return m_BP.excess(v.pos) - 1;
        }

        //! Returns the number of children of v
        size_type
        degree(const node_type& v) const {
            if (isLeaf(v))
                return 0;
            //get first child
            size_type d = 0;
            size_type pos = v.pos + 1;
            node_type child_node(v.nr + 1, pos);
            while (child_node != root()) { //while i>1 compute next sibling
                ++ d;
                child_node = next_sibling(child_node);
            }
            return d;
        }

        size_type
        isLeaf(const node_type& v) const{
            return !(m_sequence[v.pos + 1]);
        }

        //! Returns the labeled child of a node. Otherwise returns the root
        node_type
        childLabel(const node_type& v, value_type s) const {
            if (isLeaf(v))
                return root();
            //find the first label child
            size_type pos = v.pos + 1;
            value_type c = 0;
            node_type child_node(v.nr + 1, pos);
            while (child_node != root()) {
                c = m_labels[child_node.nr - 1];
                if (c == s)
                    return child_node;
                child_node = next_sibling(child_node);
            }
            return child_node; //root()
        }

        //! Returns the number of siblings to the left plus 1. That is, which
        // child number is v of its parent. Otherwise, returns 0
        size_type
        child_rank(const node_type& v){
            if (v == root())
                return 0;
            size_type r = 0;
            node_type x = v;
            while (x != root()) {
                ++ r;
                x = prev_sibling(x);
            }
            return r;
        }

        //!Return the node with id i. Otherwise return the root
        node_type
        node_select(size_type i) const {
            if (i >= nodes() or i == 0)
                return root();
            size_type pos = m_BP.select(i + 1);
            return node_type(i, pos);
        }

        //! Returns true if v is an ancestor of w
        bool
        ancestor(const node_type& v, const node_type& w) {
            if (v == root())
                return true;
            size_type v_close = m_BP.find_close(v.pos);
            if (w.pos >= v.pos and w.pos < v_close)
                return true;
            return false;
        }

        //! Returns the lower common ancestor between v and w;
        node_type
        lca(const node_type& v, const node_type& w) {
            size_type v_pos = v.pos, w_pos = w.pos, r;
            if (v_pos > w_pos)
                std::swap(v_pos, w_pos);
            else if (v_pos == w_pos)
                return v;
            if (v == root())
                return root();
            r = m_BP.double_enclose(v_pos, w_pos);
            return node_type(m_BP.rank(r) - 1, r);
        }

        //! Returns the next sibling of v. Otherwise, returns the root
        node_type
        next_sibling(const node_type& v) {
            if (v == root())
                return root();
            size_type sibling = m_BP.find_close(v.pos) + 1;
            if (!(m_sequence[sibling])) //no more siblings
                return root();
            else {
                size_type id = m_BP.rank(sibling) - 1;
                return node_type(id, sibling);
            }
        }

        //! Returns the previous sibling of v. Otherwise, returns the root
        node_type
        prev_sibling(const node_type& v) {
            if (v == root())
                return root();
            size_type sibling = v.pos - 1;
            if (!(m_sequence[sibling])) {
                sibling = m_BP.find_open(sibling);
                size_type id = m_BP.rank(sibling) - 1;
                return node_type(id, sibling);
            }
            else //no more siblings
                return root();
        }

        size_type
        subtree_size(const node_type& v) {
            return (m_BP.find_close(v.pos) - v.pos + 1) / 2;
        }

    private:


    }; //end class tree_bp

    template<class bp_support, uint8_t t_width>
    tree_bp<bp_support, t_width>::tree_bp(const sdsl::bit_vector s, const sdsl::int_vector<> labels) {
        m_labels = labels;
        m_sequence = s;
        sdsl::util::init_support(m_BP, &m_sequence);
    }


}// end namespace


#endif //DISD_TREE_LOUDS_HPP
