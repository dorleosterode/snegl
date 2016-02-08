/*
    Copyright (c) 2015, 2016 Dorle Osterode

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

#ifndef DNA_ALGORITHMS_H
#define DNA_ALGORITHMS_H

#include "sdsl_types.hpp"
#include "utility.hpp"

// header file for the dna-algorithms

class mem_iterator
{
  private:
    // arguments from the caller
    cst_t& p_cst;
    std::string& p_pattern;
    unsigned long p_l;
    bool p_mums;
    bool p_smems;
    unsigned long p_max_val;
    bool p_rev;
    // variables for the internal use
    unsigned long p_plen;
    unsigned long p_curr_pos;
    unsigned long p_left, p_right, p_new_left, p_new_right, p_q;
    unsigned long p_length_intervall;
    // variables for enumeration
    unsigned long p_mem_qs;
    unsigned long p_mem_lb;
    unsigned long p_mem_rb;
    unsigned long p_mem_pos;
    unsigned long p_mem_last_lb;
    unsigned long p_mem_last_rb;
    unsigned long p_mem_idx_f;
    unsigned long p_mem_idx_s;

  public:
    // constructor with arguments
    mem_iterator(cst_t& cst, std::string& p, unsigned long l, bool mums,
                 bool smems, unsigned long max_val, bool rev,
                 cst_t::node_type start, unsigned long matched_pos);

    // function to get the next possible mems (path_item)
    path_item next_path_item();

    // initializes the enumerator for the MEMs
    void init_mem_enumerator(path_item& pitem);

    // function to enumerate mems from path_items
    mem_struct next_mem_struct();
};

// computes all smems for the sequence p, that are longer than l
std::vector<path_item> compute_smems_for_seq(csa_t& bwd_csa, csa_t& fwd_csa, std::string& p, unsigned long l);

// computes all mems in an iterativ fashion.
std::vector<mem_struct> find_pattern(cst_t& cst, std::string& p, unsigned long l,
                                     bool mums, bool smems, unsigned long max_val,
                                     bool rev);

#endif
