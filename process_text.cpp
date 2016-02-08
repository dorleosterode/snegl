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

#include <sdsl/suffix_arrays.hpp>
#include <sdsl/suffix_trees.hpp>
#include <fstream>
#include <cstdlib>
#include <climits>
#include "process_text.hpp"

class mem_iterator_text
{
  private:
    // arguments from the caller
    cst_huff& p_cst;
    std::string& p_pattern;
    unsigned long p_l;
    bool p_mums;
    bool p_smems;
    unsigned long p_max_val;
    bool p_rev;
    bool p_dna;
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
    mem_iterator_text(cst_huff& cst, std::string& p, unsigned long l, bool mums,
                      bool smems, unsigned long max_val, bool rev, bool dna):
      p_cst(cst), p_pattern(p), p_l(l), p_mums(mums), p_smems(smems),
      p_max_val(max_val), p_rev(rev), p_dna(dna)
    {
      // variables for the internal use
      p_plen = p_pattern.size();
      p_curr_pos = p_plen;
      p_left = 0;
      p_right = p_cst.size() - 1;
      p_new_left = 0;
      p_new_right = 0;
      p_q = 0;
      p_length_intervall = 0;
    }

    // function to get the next possible mems (path_item)
    path_item next_path_item()
    {
      while (p_curr_pos >= 1) {
        if (p_curr_pos == 0)
          return path_item();

        // check if p[curr_pos-1] is N or n. if yes, iterate until a
        // not N is found!
        while (p_curr_pos >= 1 && !valid(get_char(p_pattern[idx(p_curr_pos - 1, p_plen, p_rev)], p_rev, p_dna), p_dna)) {
          p_curr_pos -= 1;
          p_left = 0;
          p_right = p_cst.size() - 1;
          p_q = 0;
        }
        if (p_curr_pos == 0)
          return path_item();
        p_length_intervall = backward_search(p_cst.csa, p_left, p_right,
                                             get_char(p_pattern[idx(p_curr_pos - 1, p_plen, p_rev)], p_rev, p_dna),
                                             p_new_left, p_new_right);
        bool is_legal = false;

        while (p_length_intervall > 0) {
          is_legal = true;
          p_q += 1;

          // set left and right to the new values for the
          // next backward-search
          p_left = p_new_left;
          p_right = p_new_right;
          p_curr_pos -= 1;

          if (p_q >= p_l && !p_smems) {
            // if we search mums -> get the mum-candidates, if
            // interval length is exactly 1
            if (p_mums) {
              if (p_length_intervall == 1)
                return path_item(p_q, p_new_left, p_new_right, p_curr_pos + 1);
            } else {
              return path_item(p_q, p_new_left, p_new_right, p_curr_pos + 1);
            }
          }

          if (p_curr_pos >= 1) {
            // check if p[curr_pos-1] is N or n.  start with a
            // completely new intervall -> break
            if (!valid(get_char(p_pattern[idx(p_curr_pos - 1, p_plen, p_rev)], p_rev, p_dna), p_dna)) {
              if (p_smems && p_q >= p_l && is_legal) {
                p_curr_pos -= 1;
                return path_item(p_q, p_left, p_right, p_curr_pos + 2);
              }

              p_left = 0;
              p_right = p_cst.size() - 1;
              p_q = 0;
              break;
            }
            p_length_intervall = backward_search(p_cst.csa, p_left, p_right,
                                                 get_char(p_pattern[idx(p_curr_pos - 1, p_plen, p_rev)], p_rev, p_dna),
                                                 p_new_left, p_new_right);
          } else
            break;
        }

        // a parent-interval can be pushed here! check, if it is a real interval!
        if (p_smems && p_q >= p_l && is_legal)
          return path_item(p_q, p_left, p_right, p_curr_pos + 1);

        if (p_left == 0 && p_right == (p_cst.csa.size() - 1)) {
          p_curr_pos -= 1;
        } else {
          auto parent = p_cst.parent(p_cst.node(p_left, p_right));
          p_left = p_cst.lb(parent);
          p_right = p_cst.rb(parent);
          p_q = p_cst.depth(parent);
        }
      }
      return path_item();
    }

    void init_mem_enumerator(path_item& pitem)
    {
      p_mem_qs = pitem.len;
      p_mem_lb = pitem.lb;
      p_mem_rb = pitem.rb;
      p_mem_pos = pitem.p;

      p_mem_last_lb = p_mem_lb;
      p_mem_last_rb = (p_mem_lb > 0)? p_mem_lb - 1: 0;
      p_mem_idx_f = 0;
      p_mem_idx_s = 0;
    }

    // function to enumerate mems from path_items
    mem_struct next_mem_struct()
    {
      // don't return insignificant MEMs
      if ((p_mem_rb - p_mem_lb) > p_max_val)
        return mem_struct();

      while (p_mem_qs >= p_l) {
        for (unsigned long k = p_mem_lb + p_mem_idx_f; k < p_mem_last_lb; k++) {
          p_mem_idx_f += 1;
          if (p_mem_pos == 1 || p_cst.csa.bwt[k] != get_char(p_pattern[idx(p_mem_pos - 2, p_plen, p_rev)], p_rev, p_dna)) {
            return mem_struct(p_cst.csa[k], p_mem_pos - 1, p_mem_qs);
          }
        }
        // for last_rb = lb - 1, this loop is the same as lb <= k <= rb
        for (unsigned long k = p_mem_last_rb + 1 + p_mem_idx_s; k <= p_mem_rb; k++) {
          p_mem_idx_s += 1;
          if (p_mem_pos == 1 || p_cst.csa.bwt[k] != get_char(p_pattern[idx(p_mem_pos - 2, p_plen, p_rev)], p_rev, p_dna)) {
            return mem_struct(p_cst.csa[k], p_mem_pos - 1, p_mem_qs);
          }
        }
        p_mem_last_lb = p_mem_lb;
        p_mem_last_rb = p_mem_rb;

        // don't enumerate all the shorter matches, if we just
        // search for MUMs or SMEMs
        if (p_mums || p_smems) {
          break;
        }

        auto parent = p_cst.parent(p_cst.node(p_mem_lb, p_mem_rb));
        p_mem_lb = p_cst.lb(parent);
        p_mem_rb = p_cst.rb(parent);
        p_mem_qs = p_cst.depth(parent);
        p_mem_idx_f = 0;
        p_mem_idx_s = 0;
        // don't return insignificant MEMs
        if ((p_mem_rb - p_mem_lb) > p_max_val)
          break;
      }
      return mem_struct();
    }
};

void process_text(cst_huff& fm_index, std::string& pattern, unsigned long l, std::string& header,
                  std::vector<seq_info>& info,
                  bool silent, bool mums, bool smems, unsigned long max_val, bool n_filter,
                  bool rev, bool comp, bool dna, FILE* out_printf, bool clean_mums)
{
  if (!silent) {
    fprintf(out_printf, "%s", header.c_str());
    if (rev)
      fprintf(out_printf, "Reverse\n");
    else
      fprintf(out_printf, "\n");
  }

  mem_iterator_text m(fm_index, pattern, l, mums, smems, max_val, rev, dna);
  path_item p = m.next_path_item();
  std::vector<std::pair<std::string, mem_struct>> mems_corrected;
  while (p.len != 0) {
    m.init_mem_enumerator(p);
    // enumerate all mems for this path item
    mem_struct mem = m.next_mem_struct();
    while (mem.len != 0) {
      if (!n_filter) {
        std::vector<std::pair<std::string, mem_struct>> corrected_p_vec = mem_corrected_vec(mem, info);
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
                mem_iterator_text m_iterator(fm_index, s, l, true, false, ULONG_MAX, rev, dna);
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
                  mem = m.next_mem_struct();
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
      } else {
        unsigned long query_pos = mem.start_query + 1;
        if (comp && rev)
          query_pos = pattern.size() - mem.start_query;
        fprintf(out_printf, "%lu\t%lu\t%lu\n", mem.start_ref + 1, query_pos, mem.len);
      }
      mem = m.next_mem_struct();
    }
    p = m.next_path_item();
  }

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
    unsigned long query_pos = prev_mem.start_query + 1;
    if (comp && rev)
      query_pos = pattern.size() - prev_mem.start_query;
    fprintf(out_printf, "%lu\t%lu\t%lu\n",
            prev_mem.abs_pos + 1,
            query_pos,
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
      query_pos = mem_wn.start_query + 1;
      if (comp && rev)
        query_pos = pattern.size() - mem_wn.start_query;
      fprintf(out_printf, "%lu\t%lu\t%lu\n",
              mem_wn.abs_pos + 1,
              query_pos,
              mem_wn.len);
    }
  }
}

// preprocesses the multifasta file or any other file into the right
// format. stores only the offset information
int preprocess_text(std::string& old_file, std::string& new_file,
                    std::vector<seq_info>& info)
{
  std::ifstream in(old_file);
  std::ofstream out(new_file);
  std::string line;
  std::string header;
  std::string seq;
  unsigned long seq_count = 0;
  std::string new_line;

  // process the first header
  getline(in, line);
  if (line[0] != '>') {
    // this isn't the right format. don't process the sequence.
    std::cerr << "wrong format: " << old_file << "\n";
    return EXIT_FAILURE;
  }

  std::size_t length = line.find(' ', 2);
  header = line.substr(1, length);

  while (getline(in, line)) {
    if (line[0] == '>') {
      // store old header with seq_len into vector
      seq_info s(header, seq_count, seq_count + seq.size());
      seq_count += seq.size();
      info.push_back(s);
      // write sequence to new_file
      out << seq;
      // get new header
      std::size_t length = line.find(' ', 2);
      header = line.substr(1, length);
      seq.clear();
    } else {
      // transform the next line and append it to seq
      new_line.clear();
      for (int i = 0; i < (int) line.size(); i++) {
        char lower = (char) tolower((int) line[i]);
        new_line.push_back(lower);
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
  return EXIT_SUCCESS;
}
