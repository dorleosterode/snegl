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

#include "utility.hpp"
#include "dna_algorithms.hpp"
#include "process_dna.hpp"
#include <climits>

// returns the last node for the lcp between p and cst. from this node
// the interval q-[i,j] can be calculated
std::pair<cst_t::node_type, unsigned long> get_starting_interval(cst_t& cst, std::string& p)
{
  cst_t::node_type node = cst.root();
  cst_t::node_type prev_node = node;
  cst_t::size_type char_pos = 0;
  cst_t::size_type d = 0;
  cst_t::size_type occs = forward_search(cst, node, d, (char) tolower((int) p[0]), char_pos);
  d += 1;
  if (occs == 0)
    return std::pair<cst_t::node_type, unsigned long> (prev_node, d);
  unsigned long i = 1;
  while (occs > 0 && i < p.size()) {
    prev_node = node;
    occs = forward_search(cst, node, d, (char) tolower((int) p[i]), char_pos);
    d += 1;
    i += 1;
  }
  if (occs == 0)
    return std::pair<cst_t::node_type, unsigned long> (prev_node, d - 1);
  else
    return std::pair<cst_t::node_type, unsigned long> (node, d);
}

// all the usage of the dna_algorithms
void compute_smems(csa_t& bwd_csa, csa_t& fwd_csa, std::string& p, unsigned long l,
                   std::string& header, std::vector<seq_info>& info,
                   bitvector_t::rank_1_type& ref_rank_1, bitvector_t::select_1_type& ref_select_1,
                   FILE* out_printf, bool clean_up)
{
  bool rev = false;
  bool comp = false;
  fprintf(out_printf, "%s\n", header.c_str());
  std::vector<path_item> pot_smems = compute_smems_for_seq(bwd_csa, fwd_csa, p, l);
  std::vector<mem_struct> smems;
  std::vector<std::pair<std::string, mem_struct>> results;

  for (auto iter = pot_smems.begin(); iter != pot_smems.end(); iter ++) {
    unsigned long lb = (*iter).lb;
    unsigned long rb = (*iter).rb;
    unsigned long qs = (*iter).len;
    unsigned long pos = (*iter).p;

    unsigned long last_rb = (lb > 0)? lb - 1: 0;
    for (unsigned long k = last_rb + 1; k <= rb; k++) {
      if (pos == 0 || bwd_csa.bwt[k] != get_char(p[idx(pos - 1, p.size(), rev)], rev)) {
        mem_struct m = mem_struct(bwd_csa[k], pos, qs);
        if (clean_up) {
          // check for N's
          std::vector<mem_struct> vec_wn = mem_without_n_bv(m, ref_rank_1,
                                           ref_select_1, true);
          for (auto mem = vec_wn.begin(); mem != vec_wn.end(); mem++) {
            if ((*mem).len < l)
              continue;
            // check for borders
            std::vector<std::pair<std::string, mem_struct>> corrected_vec = mem_corrected_vec(*mem, info);
            for (auto iter_c = corrected_vec.begin(); iter_c != corrected_vec.end(); iter_c++) {
              std::pair<std::string, mem_struct> corrected = *iter_c;
              if (corrected.second.len < l)
                continue;
              results.push_back(corrected);
            }
          }
        } else {
          fprintf(out_printf, "%lu\t%lu\t%lu\n", m.start_ref + 1, m.start_query + 1, m.len);
        }
      }
    }
  }

  if (clean_up && results.size() > 0) {
    // sort the smems and discard the contained ones
    std::sort(results.begin(), results.end(), compare_mem_struct_wrt_query);
    mem_struct& prev_mem = results[0].second;
    int first_idx = 0;
    // the first smem can't be contained!
    while (prev_mem.len < l) {
      first_idx += 1;
      prev_mem = results[first_idx].second;
    }
    fprintf(out_printf, "%lu\t%lu\t%lu\n",
            prev_mem.abs_pos + 1,
            prev_mem.start_query + 1,
            prev_mem.len);
    for (int i = first_idx + 1; i < (int) results.size(); i++) {
      mem_struct& mem_wn = results[i].second;
      if (mem_wn.len < l)
        continue;
      // check the smem for containment in the smem before
      if (contained(prev_mem, mem_wn)) {
        prev_mem = mem_wn;
        continue;
      }
      prev_mem = mem_wn;
      unsigned long query_pos = mem_wn.start_query + 1;
      if (comp && rev)
        query_pos = p.size() - mem_wn.start_query;
      fprintf(out_printf, "%lu\t%lu\t%lu\n",
              mem_wn.abs_pos + 1,
              query_pos,
              mem_wn.len);
    }
  }
}

void process_seq_find_pattern(cst_t& fm_index, std::string& pattern,
                              unsigned long l, std::string& header,
                              std::vector<seq_info>& info,
                              bitvector_t::rank_1_type& ref_rank_1,
                              bitvector_t::select_1_type& ref_select_1,
                              bool silent, bool mums, bool smems, unsigned long max_val,
                              bool n_filter, bool rev, bool comp, FILE* out_printf,
                              bool clean_mums, std::vector<std::pair<std::string, mem_struct>>& mems_corrected)
{
  std::vector<mem_struct> results = find_pattern(fm_index, pattern, l, mums, smems, max_val, rev);
  for (auto iter = results.begin(); iter != results.end(); iter++) {
    mem_struct mem = *iter;
    std::vector<mem_struct> vec_wn1 = mem_without_n_bv(mem, ref_rank_1,
                                      ref_select_1, true);
    for (auto elem_wn = vec_wn1.begin(); elem_wn != vec_wn1.end(); elem_wn++) {
      mem_struct mem_wn = *elem_wn;
      if (mem_wn.len < l)
        continue;
      std::vector<std::pair<std::string, mem_struct>> corrected_p_vec = mem_corrected_vec(mem_wn, info);
      for (auto iter_c = corrected_p_vec.begin(); iter_c != corrected_p_vec.end(); iter_c++) {
        std::pair<std::string, mem_struct> corrected_p = *iter_c;
        if (corrected_p.second.len >= l) {
          if (!mums && !smems && !silent) {
            unsigned long query_pos = corrected_p.second.start_query + 1;
            if (comp && rev)
              query_pos = pattern.size() - corrected_p.second.start_query;
            fprintf(out_printf, "%s\t%lu\t%lu\t%lu\n",
                    corrected_p.first.c_str(),
                    corrected_p.second.start_ref + 1,
                    query_pos,
                    corrected_p.second.len);
          } else {
            mems_corrected.push_back(corrected_p);
          }
        }
      }
    }
  }
}

void process_seq_sort_sa(cst_t& fm_index, std::string& pattern,
                         unsigned long l, std::string& header,
                         std::vector<seq_info>& info,
                         bitvector_t::rank_1_type& ref_rank_1,
                         bitvector_t::select_1_type& ref_select_1,
                         bool silent, bool mums, bool smems, unsigned long max_val,
                         bool n_filter, bool rev, bool comp, FILE* out_printf,
                         bool enum_mems, bool clean_mums, bool sort, bool map_sa, std::string sa_name,
                         unsigned long nof_intervals, unsigned long sof_lengths,
                         std::vector<std::pair<std::string, mem_struct>>& mems_corrected)
{
  std::vector<path_item> path_items;

  mem_iterator m(fm_index, pattern, l, mums, smems, max_val, rev, fm_index.root(), 0);
  // first enumerate all path_items.
  path_item p = m.next_path_item();
  while (p.len != 0) {
    path_items.push_back(p);
    p = m.next_path_item();
  }

  if (sort)
    std::sort(path_items.begin(), path_items.end(), compare_path_item);

  if (enum_mems) {
    for (auto iter = path_items.begin(); iter != path_items.end(); iter++) {
      m.init_mem_enumerator(*iter);
      // enumerate all mems for this path item
      mem_struct mem = m.next_mem_struct();
      while (mem.len != 0) {
        if (!n_filter) {
          std::vector<mem_struct> vec_wn1 = mem_without_n_bv(mem, ref_rank_1,
                                            ref_select_1, true);
          for (auto elem_wn = vec_wn1.begin(); elem_wn != vec_wn1.end(); elem_wn++) {
            mem_struct mem_wn = *elem_wn;
            if (mem_wn.len < l)
              continue;
            std::vector<std::pair<std::string, mem_struct>> corrected_p_vec = mem_corrected_vec(mem_wn, info);
            for (auto iter_c = corrected_p_vec.begin(); iter_c != corrected_p_vec.end(); iter_c++) {
              std::pair<std::string, mem_struct> corrected_p = *iter_c;
              if (corrected_p.second.len >= l) {
                if (!mums && !smems && !silent) {
                  unsigned long query_pos = corrected_p.second.start_query + 1;
                  if (comp && rev)
                    query_pos = pattern.size() - corrected_p.second.start_query;
                  fprintf(out_printf, "%s\t%lu\t%lu\t%lu\n",
                          corrected_p.first.c_str(),
                          corrected_p.second.start_ref + 1,
                          query_pos,
                          corrected_p.second.len);
                } else {
                  mems_corrected.push_back(corrected_p);
                }
              }
            }
          }
        }
        mem = m.next_mem_struct();
      }
    }
  } else if (map_sa) {
    std::vector<path_item> tmp;
    // enumerate the mems with the sa instead of the csa. need to
    // have all sa-intervals first.
    for (auto iter = path_items.begin(); iter != path_items.end(); iter++) {
      unsigned long pos = (*iter).p;
      unsigned long lb = (*iter).lb;
      unsigned long rb = (*iter).rb;
      unsigned long qs = (*iter).len;
      unsigned long last_lb = lb;
      unsigned long last_rb = rb;

      tmp.push_back(*iter);
      nof_intervals += 1;
      sof_lengths += rb - lb + 1;
      auto parent = fm_index.parent(fm_index.node(lb, rb));
      lb = fm_index.lb(parent);
      rb = fm_index.rb(parent);
      qs = fm_index.depth(parent);

      while (qs >= l && !mums && !smems) {
        nof_intervals += 2;
        tmp.push_back(path_item(qs, lb, last_lb - 1, pos));
        tmp.push_back(path_item(qs, last_rb + 1, rb, pos));
        sof_lengths += last_lb - 1 - lb + 1;
        sof_lengths += rb - last_rb + 1 + 1;

        last_lb = lb;
        last_rb = rb;
        auto parent = fm_index.parent(fm_index.node(lb, rb));
        lb = fm_index.lb(parent);
        rb = fm_index.rb(parent);
        qs = fm_index.depth(parent);
      }
    }
    path_items = tmp;

    // iterate over the explicit SA
    int_vector_buffer<> sa(sa_name);
    std::sort(path_items.begin(), path_items.end(), compare_path_item);
    for (auto iter = path_items.begin(); iter != path_items.end(); iter++) {
      // calculate the MEM here
      for (unsigned long k = (*iter).lb; k <= (*iter).rb; k++) {
        if ((*iter).p == 1 || fm_index.csa.bwt[k] != get_char(pattern[idx((*iter).p - 2, pattern.size(), rev)], rev)) {
          mem_struct mem = mem_struct(sa[k], (*iter).p - 1, (*iter).len);
          std::vector<mem_struct> vec_wn1 = mem_without_n_bv(mem, ref_rank_1,
                                            ref_select_1, true);
          for (auto elem_wn = vec_wn1.begin(); elem_wn != vec_wn1.end(); elem_wn++) {
            mem_struct mem_wn = *elem_wn;
            if (mem_wn.len < l)
              continue;
            std::vector<std::pair<std::string, mem_struct>> corrected_p_vec = mem_corrected_vec(mem_wn, info);
            for (auto iter_c = corrected_p_vec.begin(); iter_c != corrected_p_vec.end(); iter_c++) {
              std::pair<std::string, mem_struct> corrected_p = *iter_c;
              if (corrected_p.second.len >= l) {
                if (!mums && !smems && !silent) {
                  unsigned long query_pos = corrected_p.second.start_query + 1;
                  if (comp && rev)
                    query_pos = pattern.size() - corrected_p.second.start_query;
                  fprintf(out_printf, "%s\t%lu\t%lu\t%lu\n",
                          corrected_p.first.c_str(),
                          corrected_p.second.start_ref + 1,
                          query_pos,
                          corrected_p.second.len);
                } else {
                  mems_corrected.push_back(corrected_p);
                }
              }
            }
          }
        }
      }
    }
  }
}


void process_seq_generator(cst_t& fm_index, std::string& pattern,
                           unsigned long l, std::string& header,
                           std::vector<seq_info>& info,
                           bitvector_t::rank_1_type& ref_rank_1,
                           bitvector_t::select_1_type& ref_select_1,
                           bool silent, bool mums, bool smems, unsigned long max_val,
                           bool n_filter, bool rev, bool comp, FILE* out_printf,
                           bool clean_mums, std::vector<std::pair<std::string, mem_struct>>& mems_corrected)
{
  std::vector<path_item> path_items;

  mem_iterator m(fm_index, pattern, l, mums, smems, max_val, rev, fm_index.root(), 0);
  path_item p = m.next_path_item();
  while (p.len != 0) {
    m.init_mem_enumerator(p);
    // enumerate all mems for this path item
    mem_struct mem = m.next_mem_struct();
    while (mem.len != 0) {
      if (!n_filter) {
        std::vector<mem_struct> vec_wn1 = mem_without_n_bv(mem, ref_rank_1,
                                          ref_select_1, true);
        for (auto elem_wn = vec_wn1.begin(); elem_wn != vec_wn1.end(); elem_wn++) {
          mem_struct mem_wn = *elem_wn;
          if (mem_wn.len < l)
            continue;
          std::vector<std::pair<std::string, mem_struct>> corrected_p_vec = mem_corrected_vec(mem_wn, info);
          for (auto iter_c = corrected_p_vec.begin(); iter_c != corrected_p_vec.end(); iter_c++) {
            std::pair<std::string, mem_struct> corrected_p = *iter_c;
            if (corrected_p.second.len >= l) {
              if (!mums && !smems && !silent) {
                unsigned long query_pos = corrected_p.second.start_query + 1;
                if (comp && rev)
                  query_pos = pattern.size() - corrected_p.second.start_query;
                fprintf(out_printf, "%s\t%lu\t%lu\t%lu\n",
                        corrected_p.first.c_str(),
                        corrected_p.second.start_ref + 1,
                        query_pos,
                        corrected_p.second.len);
              }

              if (mums && clean_mums) {
                unsigned long start = 0, len = 0;
                if (mem.len != corrected_p.second.len) {
                  // something has changed, while correcting and cleaning!
                  unsigned long start_diff = corrected_p.second.start_query - mem.start_query;
                  if (start_diff == 0) {
                    start = mem.start_query;
                  } else {
                    start = mem.start_query + start_diff;
                  }
                  len = corrected_p.second.len;

                  std::string s = pattern.substr(start, len);
                  mem_iterator m_iterator(fm_index, s, l, true, false, ULONG_MAX, rev, fm_index.root(), 0);
                  path_item p_item = m_iterator.next_path_item();
                  std::vector<mem_struct> results;
                  while (p_item.len != 0) {
                    m_iterator.init_mem_enumerator(p_item);
                    // enumerate all mems for this path item
                    mem_struct mem2 = m_iterator.next_mem_struct();
                    while (mem2.len != 0) {
                      results.push_back(mem2);
                      mem2 = m_iterator.next_mem_struct();
                    }
                    p_item = m_iterator.next_path_item();
                  }
                  if (results.size() != 1) {
                    continue;
                  } else
                    mems_corrected.push_back(corrected_p);
                } else
                  mems_corrected.push_back(corrected_p);
              }
              if (smems)
                mems_corrected.push_back(corrected_p);
            }
          }
        }
      } else {
        if (smems) {
          unsigned long query_pos = mem.start_query + 1;
          if (comp && rev)
            query_pos = pattern.size() - mem.start_query;
          fprintf(out_printf, "%lu\t%lu\t%lu\n", mem.start_ref + 1, query_pos, mem.len);
        } else {
          std::vector<std::pair<std::string, mem_struct>> cor_vec = mem_corrected_vec(mem, info);
          for (auto iter_c = cor_vec.begin(); iter_c != cor_vec.end(); iter_c++) {
            std::pair<std::string, mem_struct> cor = *iter_c;
            unsigned long query_pos = cor.second.start_query + 1;
            if (comp && rev)
              query_pos = pattern.size() - cor.second.start_query;
            fprintf(out_printf, "%s\t%lu\t%lu\t%lu\n",
                    cor.first.c_str(),
                    cor.second.start_ref + 1,
                    query_pos,
                    cor.second.len);
          }
        }
      }
      mem = m.next_mem_struct();
    }
    p = m.next_path_item();
  }
}

void process_seq(cst_t& fm_index, std::string& pattern, unsigned long l, std::string& header,
                 std::vector<seq_info>& info,
                 bitvector_t::rank_1_type& ref_rank_1, bitvector_t::select_1_type& ref_select_1,
                 bool silent, bool mums, bool smems, unsigned long max_val,
                 bool n_filter, bool rev, bool comp, FILE* out_printf, bool sort_paths, bool enum_mems,
                 bool clean_mums, bool sort, bool map_sa, std::string sa_name,
                 unsigned long nof_intervals, unsigned long sof_lengths)
{
  std::vector<std::pair<std::string, mem_struct>> mems_corrected;
  bool test_find_pattern = false;

  if (!silent) {
    fprintf(out_printf, "%s", header.c_str());
    if (rev)
      fprintf(out_printf, "Reverse\n");
    else
      fprintf(out_printf, "\n");
  }

  if (sort_paths)
    process_seq_sort_sa(fm_index, pattern, l, header, info, ref_rank_1,
                        ref_select_1, silent, mums, smems, max_val, n_filter,
                        rev, comp, out_printf, enum_mems, clean_mums, sort,
                        map_sa, sa_name, nof_intervals, sof_lengths, mems_corrected);
  else if (test_find_pattern)
    process_seq_find_pattern(fm_index, pattern, l, header, info, ref_rank_1,
                             ref_select_1, silent, mums, smems, max_val, n_filter,
                             rev, comp, out_printf, clean_mums,
                             mems_corrected);
  else
    process_seq_generator(fm_index, pattern, l, header, info, ref_rank_1,
                          ref_select_1, silent, mums, smems, max_val, n_filter,
                          rev, comp, out_printf, clean_mums,
                          mems_corrected);

  if (mums) {
    std::sort(mems_corrected.begin(), mems_corrected.end(), compare_mem_struct_wrt_ref);
    mems_corrected = clean_mumcands(mems_corrected);

    for (auto iter = mems_corrected.begin(); iter != mems_corrected.end(); iter++) {
      fprintf(out_printf, "%s\t%lu\t%lu\t%lu\n",
              (*iter).first.c_str(),
              (*iter).second.start_ref + 1,
              (*iter).second.start_query + 1,
              (*iter).second.len);
    }
  }

  if (smems && (mems_corrected.size() != 0)) {
    std::sort(mems_corrected.begin(), mems_corrected.end(), compare_mem_struct_wrt_query);
    mem_struct& prev_mem = mems_corrected[0].second;
    int first_idx = 0;
    // the first smem can't be contained!
    while (prev_mem.len < l) {
      first_idx += 1;
      prev_mem = mems_corrected[first_idx].second;
    }
    fprintf(out_printf, "%lu\t%lu\t%lu\n",
            prev_mem.abs_pos + 1,
            prev_mem.start_query + 1,
            prev_mem.len);
    for (int i = first_idx + 1; i < (int) mems_corrected.size(); i++) {
      mem_struct& mem_wn = mems_corrected[i].second;
      if (mem_wn.len < l)
        continue;
      // check the smem for containment in the smem before
      if (contained(prev_mem, mem_wn)) {
        prev_mem = mem_wn;
        continue;
      }
      prev_mem = mem_wn;
      unsigned long query_pos = mem_wn.start_query + 1;
      if (comp && rev)
        query_pos = pattern.size() - mem_wn.start_query;
      fprintf(out_printf, "%lu\t%lu\t%lu\n",
              mem_wn.abs_pos + 1,
              query_pos,
              mem_wn.len);
    }
  }
}

// function that returns a random character from dna
char random_dna(std::minstd_rand0& r, char dna[])
{
  std::uniform_int_distribution<> dis(0,3);
  return dna[dis(r)];
}

// preprocesses the multifasta file into the right format. stores the offset information
bitvector_t preprocess_seq(std::string& old_file, std::string& new_file,
                           std::vector<seq_info>& info, bool construct_bv)
{
  std::ifstream in(old_file);
  std::ofstream out(new_file);
  std::string line;
  std::string header;
  std::string seq;
  unsigned long seq_count = 0;
  std::string new_line;
  bool store = false;
  unsigned long resize_size = 65536;
  bit_vector bv_ns(0,0);
  if (construct_bv)
    bv_ns.resize(resize_size);

  // parameters for random_dna
  std::minstd_rand0 rd(1);
  char dna[] = { 'a', 'c', 'g', 't' };

  while (getline(in, line)) {
    if (line[0] == '>') {
      if (store) {
        // store old header with seq_len into vector
        seq_info s(header, seq_count, seq_count + seq.size());
        seq_count += seq.size();
        info.push_back(s);
        // write sequence to new_file
        out << seq;
      }
      // get new header
      std::size_t length = line.find(' ', 2);
      header = line.substr(1, length);
      store = true;
      seq.clear();
    } else {
      // transform the next line and append it to seq
      new_line.clear();
      for (int i = 0; i < (int) line.size(); i++) {
        char lower = (char) tolower((int) line[i]);
        unsigned long idx = seq_count + seq.size() + i;
        // check if bv has to be resized
        if (construct_bv && bv_ns.size() <= idx) {
          bv_ns.resize(bv_ns.capacity() + resize_size);
        }
        if (lower == 'a' || lower == 'c' || lower == 'g' || lower == 't') {
          new_line.push_back(lower);
          if (construct_bv)
            bv_ns[idx] = 0;
        } else {
          // store the information, that seq[j] was 'n' with
          // j = seq_count + seq.size() + i
          if (construct_bv)
            bv_ns[idx] = 1;
          new_line.push_back(random_dna(rd, dna));
        }
      }
      seq += new_line;
    }
  }
  // store the last header with seq_len into vector
  seq_info s(header, seq_count, seq_count + seq.size());
  seq_count += seq.size();
  info.push_back(s);
  // write sequence to new_file
  out << seq;

  // transform the bit_vector to a sparse bit_vector
  if (construct_bv)
    bv_ns.resize(seq_count);
  return bitvector_t (bv_ns);
}

// preprocesses the multifasta file of the query into the right format. stores the offset information
unsigned long preprocess_query(std::string& old_file, std::string& new_file, std::vector<std::pair<std::string, unsigned long>>& info)
{
  std::ifstream in(old_file);
  std::ofstream out(new_file);
  std::string line;
  std::string header;
  std::string seq;
  unsigned long seq_count = 0;
  bool store = false;

  while (getline(in, line)) {
    if (line[0] == '>') {
      if (store) {
        // store old header with seq_len into vector
        seq += ">";
        seq_count += seq.size();
        std::pair<std::string, unsigned long> p(header, seq_count);
        info.push_back(p);
        // write sequence to new_file
        out << seq;
      }
      // get new header
      std::size_t length = line.find(' ', 2);
      header = line.substr(1, length);
      store = true;
      seq.clear();
    } else
      seq += line;
  }
  // store the last header with seq_len into vector
  seq_count += seq.size();
  std::pair<std::string, unsigned long> p(header, seq_count);
  info.push_back(p);
  // write sequence to new_file
  out << seq;
  return seq_count;
}

std::string process_mem(mem_struct& mem, bitvector_t::rank_1_type& rank, bitvector_t::select_1_type& select,
                        unsigned long l, std::vector<seq_info>& info,
                        std::vector<std::pair<std::string, unsigned long>>& query_info, unsigned long offset,
                        FILE* out, std::string last_header)
{
  // step 1: correct wrt Ns in the reference
  std::vector<mem_struct> vec = mem_without_n_bv(mem, rank, select, true);
  // step 2: correct the position in the reference
  for (auto mem_wn = vec.begin(); mem_wn != vec.end(); mem_wn++) {
    if ((*mem_wn).len < l)
      continue;
    std::vector<std::pair<std::string, mem_struct>> cor_vec = mem_corrected_vec(*mem_wn, info);
    for (auto iter_c = cor_vec.begin(); iter_c != cor_vec.end(); iter_c++) {
      std::pair<std::string, mem_struct> cor = *iter_c;
      // step 3: correct the position in the query
      if (cor.second.len < l)
        continue;
      std::pair<std::string, mem_struct> cor2 = get_header(cor.second, query_info, false, offset);
      if (cor2.first != last_header) {
        fprintf(out, "> %s\n", cor2.first.c_str());
        last_header = cor2.first;
      }
      fprintf(out, "%s\t%lu\t%lu\t%lu\n",
              cor.first.c_str(),
              cor2.second.start_ref + 1,
              cor2.second.start_query + 1,
              cor2.second.len);
    }
  }
  return last_header;
}

void thread_process_seq(thread_args args)
{
  // step 1: find the pattern with the correct starting interval
  bool mums = false;
  bool smems = false;
  unsigned long max_val = ULONG_MAX;
  bool rev = false;

  mem_iterator m_it(args.fm_index, args.pattern, args.l, mums, smems, max_val, rev,
                    args.start_node, args.matched_pos);
  path_item p_item = m_it.next_path_item();
  std::vector<mem_struct> results;
  std::vector<path_item> path_items;
  std::string last_header = "";
  while (p_item.len != 0) {
    // found a seperator in the query
    if (p_item.len == -1) {
      // iterate over all stored path_items and output all mems!
      for (auto path = path_items.begin(); path != path_items.end(); path++) {
        m_it.init_mem_enumerator(*path);
        // enumerate all mems for this path item
        mem_struct mem = m_it.next_mem_struct();
        while (mem.len != 0) {
          if (args.thread_num != 0 && mem.start_query == 0 && mem.start_ref != 0)
            if (args.fm_index.csa.bwt[mem.csa_pos] == args.prev_char) {
              mem = m_it.next_mem_struct();
              continue;
            }

          last_header = process_mem(mem, args.ref_rank_1, args.ref_select_1,
                                    args.l, args.info, args.query_info, args.offset, args.out, last_header);
          mem = m_it.next_mem_struct();
        }
      }
      path_items.clear();
    } else {
      path_items.push_back(p_item);
    }
    p_item = m_it.next_path_item();
  }
  // iterate over all stored path_items and output all mems!
  for (auto path = path_items.begin(); path != path_items.end(); path++) {
    m_it.init_mem_enumerator(*path);
    // enumerate all mems for this path item
    mem_struct mem = m_it.next_mem_struct();
    while (mem.len != 0) {
      if (args.thread_num != 0 && mem.start_query == 0 && mem.start_ref != 0)
        if (args.fm_index.csa.bwt[mem.csa_pos] == args.prev_char) {
          mem = m_it.next_mem_struct();
          continue;
        }
      last_header = process_mem(mem, args.ref_rank_1, args.ref_select_1,
                                args.l, args.info, args.query_info, args.offset, args.out, last_header);
      mem = m_it.next_mem_struct();
    }
  }
}
