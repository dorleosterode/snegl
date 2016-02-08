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

#ifndef UTILITY_H
#define UTILITY_H

#include <cstdlib>
#include <string>
#include <vector>
#include "sdsl_types.hpp"

struct seq_info {
  std::string header;
  unsigned long start, end;
  seq_info(std::string header, unsigned long start, unsigned long end):header(header), start(start), end(end) {}
};

struct path_item {
  long len;
  unsigned long lb, rb, p;
  path_item(long len=0, unsigned long lb=0, unsigned long rb=0, unsigned long p=0):len(len),lb(lb),rb(rb),p(p) {}
};

struct mem_struct {
  unsigned long start_ref, start_query, len;
  long abs_pos;
  unsigned long csa_pos;
  mem_struct(unsigned long start_ref=0, unsigned long start_query=0, unsigned long len=0, long abs_pos=-1, unsigned long csa_pos = 0):
    start_ref(start_ref), start_query(start_query), len(len), abs_pos(abs_pos), csa_pos(csa_pos) {}
};

// p1 comes before p2, if p1.lb < p2.lb or p1.rb < p2.rb
bool compare_path_item(path_item& p1, path_item& p2);

// m1 comes before m2, if m1.start_ref < m2.start_ref or if m1 is larger
bool compare_mem_struct_wrt_ref(std::pair<std::string, mem_struct>& m1, std::pair<std::string, mem_struct>& m2);

// m1 comes before m2, if m1.start_query < m2.start_query or if m1 is larger
bool compare_mem_struct_wrt_query(std::pair<std::string, mem_struct>& m1, std::pair<std::string, mem_struct>& m2);

unsigned long idx(unsigned long i, unsigned long len, bool rev);

char get_char(char c, bool rev, bool dna=true);

bool valid(char c, bool dna=true);

// this function corrects the given MEM (start_ref, start_query, len) w.r.t.
// info. info stores the absolute startpositions of the sequence parts of the reference.
// the start_ref is corrected to the relative startposition in the correpsponding
// sequence part. also is checked, if the MEM is crossing a sequence boundary. in that
// case the startpositions (ref and query) and the length are corrected. the correct MEM is returned
// with the additional information of the sequence-identifier.
std::vector<std::pair<std::string, mem_struct>>mem_corrected_vec(mem_struct m,
    std::vector<seq_info>& info,
    bool ref = true, unsigned long offset = 0);

// this function corrects the given MEM (start_ref, start_query, len) w.r.t.
// info. info stores the absolute startpositions of the sequence parts of the reference.
// the start_ref is corrected to the relative startposition in the correpsponding
// sequence part. also is checked, if the MEM is crossing a sequence boundary. in that
// case the startpositions (ref and query) and the length are corrected. the correct MEM is returned
// with the additional information of the sequence-identifier.
std::vector<std::pair<std::string, mem_struct>>mem_corrected(mem_struct m,
    std::vector<seq_info>& info,
    bool ref = true, unsigned long offset = 0);

// returns the header for the sequence, which contains m
std::pair<std::string, mem_struct> get_header(mem_struct m,
    std::vector<std::pair<std::string, unsigned long>>& info,
    bool ref = false, unsigned long offset = 0);

std::vector<std::pair<std::string, mem_struct>> clean_mumcands(std::vector<std::pair<std::string, mem_struct>> mums);

bool contained(mem_struct& m1, mem_struct& m2);

// iterative solution to clean mems, mums and smem from wildcard-matches.
std::vector<mem_struct> mem_without_n_it(mem_struct m, std::vector<unsigned long>& ns, bool ref);

// bitvector solution to clean mems, mums and smems from wildcard-matches.
std::vector<mem_struct> mem_without_n_bv(mem_struct& m, bitvector_t::rank_1_type& rank_1,
    bitvector_t::select_1_type& select_1, bool ref);


#endif
