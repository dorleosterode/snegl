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

#ifndef PROCESS_DNA_H
#define PROCESS_DNA_H

#include "sdsl_types.hpp"

struct thread_args {
  cst_t& fm_index;
  std::string& pattern;
  unsigned long l;
  FILE* out;
  std::vector<seq_info>& info;
  std::vector<std::pair<std::string, unsigned long>>& query_info;
  bitvector_t::rank_1_type& ref_rank_1;
  bitvector_t::select_1_type& ref_select_1;
  unsigned long offset;
  unsigned long nof_threads; // maybe not needed
  unsigned long thread_num;
  cst_t::node_type start_node;
  unsigned long matched_pos;
  char prev_char;
  thread_args(cst_t& fm_index,
              std::string& pattern,
              unsigned long l,
              FILE* out,
              std::vector<seq_info>& info,
              std::vector<std::pair<std::string, unsigned long>>& query_info,
              bitvector_t::rank_1_type& ref_rank_1,
              bitvector_t::select_1_type& ref_select_1,
              unsigned long offset,
              unsigned long nof_threads,
              unsigned long thread_num,
              cst_t::node_type start_node,
              unsigned long matched_pos,
              char prev_char):
    fm_index(fm_index),
    pattern(pattern),
    l(l),
    out(out),
    info(info),
    query_info(query_info),
    ref_rank_1(ref_rank_1),
    ref_select_1(ref_select_1),
    offset(offset),
    nof_threads(nof_threads),
    thread_num(thread_num),
    start_node(start_node),
    matched_pos(matched_pos),
    prev_char(prev_char) {}
};

// parallel version of process_seq
void thread_process_seq(thread_args args);

// function to compute all smems with fermi-algorithm. results are written to out_printf
void compute_smems(csa_t& bwd_csa, csa_t& fwd_csa, std::string& p, unsigned long l,
                   std::string& header, std::vector<seq_info>& info,
                   bitvector_t::rank_1_type& ref_rank_1, bitvector_t::select_1_type& ref_select_1,
                   FILE* out_printf, bool clean_up);

// function to compute mems, mums or smems. results are written to out_printf.
void process_seq(cst_t& fm_index, std::string& pattern, unsigned long l, std::string& header,
                 std::vector<seq_info>& info,
                 bitvector_t::rank_1_type& ref_rank_1, bitvector_t::select_1_type& ref_select_1,
                 bool silent, bool mums, bool smems, unsigned long max_val,
                 bool n_filter, bool rev, bool comp, FILE* out_printf, bool sort_paths, bool enum_mems,
                 bool clean_mums, bool sort, bool map_sa = false, std::string sa_name = "",
                 unsigned long nof_intervals = 0, unsigned long sof_lengths = 0);

// preprocesses the multifasta file into the right format. stores the offset information
bitvector_t preprocess_seq(std::string& old_file, std::string& new_file,
                           std::vector<seq_info>& info, bool construct_bv);

// preprocesses the multifasta file of the query into the right format. stores the offset information
unsigned long preprocess_query(std::string& old_file, std::string& new_file,
                               std::vector<std::pair<std::string, unsigned long>>& info);

// returns the last node for the lcp between p and cst. from this node
// the interval q-[i,j] can be calculated
std::pair<cst_t::node_type, unsigned long> get_starting_interval(cst_t& cst, std::string& p);

#endif
