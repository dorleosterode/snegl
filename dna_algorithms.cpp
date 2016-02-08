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

// all the algorithms to process dna sequences
#include "dna_algorithms.hpp"

// functions from class mem_iterator
// constructor with arguments
mem_iterator::mem_iterator(cst_t& cst, std::string& p, unsigned long l, bool mums,
                           bool smems, unsigned long max_val, bool rev,
                           cst_t::node_type start, unsigned long matched_pos):
  p_cst(cst), p_pattern(p), p_l(l), p_mums(mums), p_smems(smems),
  p_max_val(max_val), p_rev(rev)
{
  // variables for the internal use
  p_plen = p_pattern.size();
  p_curr_pos = p_plen;
  p_left = p_cst.lb(start);
  p_right = p_cst.rb(start);
  p_new_left = 0;
  p_new_right = 0;
  p_q = matched_pos;
  p_length_intervall = 0;
}

// function to get the next possible mems (path_item)
path_item mem_iterator::next_path_item()
{
  while (p_curr_pos >= 1) {
    if (p_curr_pos == 0)
      return path_item();

    // check if p[curr_pos-1] is N or n. if yes, iterate until a not N
    // is found!
    while (p_curr_pos >= 1 && !valid(get_char(p_pattern[idx(p_curr_pos - 1, p_plen, p_rev)], p_rev))) {
      p_curr_pos -= 1;
      p_left = 0;
      p_right = p_cst.size() - 1;
      p_q = 0;
      if (get_char(p_pattern[idx(p_curr_pos, p_plen, p_rev)], p_rev) == '>') {
        // found a seperator
        return path_item(-1);
      }
    }
    if (p_curr_pos == 0)
      return path_item();
    p_length_intervall = backward_search(p_cst.csa, p_left, p_right,
                                         get_char(p_pattern[idx(p_curr_pos - 1, p_plen, p_rev)], p_rev),
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
        if (!valid(get_char(p_pattern[idx(p_curr_pos - 1, p_plen, p_rev)], p_rev))) {
          // a parent-interval can be pushed here! check, if it is a real interval!
          if (p_smems && p_q >= p_l && is_legal) {
            unsigned long old_left = p_left;
            unsigned long old_right = p_right;
            unsigned long old_q = p_q;
            p_curr_pos -= 1;
            p_left = 0;
            p_right = p_cst.size() - 1;
            p_q = 0;

            return path_item(old_q, old_left, old_right, p_curr_pos + 2);
          }
          // found a seperator
          if (get_char(p_pattern[idx(p_curr_pos - 1, p_plen, p_rev)], p_rev) == '>') {
            p_curr_pos -= 1;
            p_left = 0;
            p_right = p_cst.size() - 1;
            p_q = 0;
            return path_item(-1);
          }

          p_left = 0;
          p_right = p_cst.size() - 1;
          p_q = 0;
          break;
        }
        p_length_intervall = backward_search(p_cst.csa, p_left, p_right,
                                             get_char(p_pattern[idx(p_curr_pos - 1, p_plen, p_rev)], p_rev),
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

void mem_iterator::init_mem_enumerator(path_item& pitem)
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
mem_struct mem_iterator::next_mem_struct()
{
  // don't return insignificant MEMs
  if ((p_mem_rb - p_mem_lb) > p_max_val)
    return mem_struct();

  while (p_mem_qs >= p_l) {
    for (unsigned long k = p_mem_lb + p_mem_idx_f; k < p_mem_last_lb; k++) {
      p_mem_idx_f += 1;
      if (p_mem_pos == 1 || p_cst.csa.bwt[k] != get_char(p_pattern[idx(p_mem_pos - 2, p_plen, p_rev)], p_rev)) {
        return mem_struct(p_cst.csa[k], p_mem_pos - 1, p_mem_qs);
      }
    }
    // for last_rb = lb - 1, this loop is the same as lb <= k <= rb
    for (unsigned long k = p_mem_last_rb + 1 + p_mem_idx_s; k <= p_mem_rb; k++) {
      p_mem_idx_s += 1;
      if (p_mem_pos == 1 || p_cst.csa.bwt[k] != get_char(p_pattern[idx(p_mem_pos - 2, p_plen, p_rev)], p_rev)) {
        return mem_struct(p_cst.csa[k], p_mem_pos - 1, p_mem_qs);
      }
    }
    p_mem_last_lb = p_mem_lb;
    p_mem_last_rb = p_mem_rb;

    // don't enumerate all the shorter matches, if we just search for
    // MUMs or SMEMs
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

// bidirectional search for csa with implicit sentinel!
unsigned long bidirectional_search_new(csa_t& csa_fwd, unsigned long l_fwd, unsigned long r_fwd,
                                       unsigned long l_bwd, unsigned long r_bwd, char c,
                                       unsigned long& l_fwd_res, unsigned long& r_fwd_res,
                                       unsigned long& l_bwd_res, unsigned long& r_bwd_res)
{
  assert(l_fwd <= r_fwd); assert(r_fwd < csa_fwd.size());
  char cc = csa_fwd.char2comp[c];
  unsigned long c_begin = csa_fwd.C[cc];
  unsigned long l_fwd_bwt = l_fwd;
  unsigned long r_fwd_bwt = r_fwd;
  if (c > 0)
    cc = csa_fwd.char2comp[c] - 1;
  l_fwd_bwt = l_fwd_bwt - (l_fwd > csa_fwd.sentinel_pos);
  r_fwd_bwt = r_fwd_bwt - (r_fwd > csa_fwd.sentinel_pos);
  auto r_s_b =  csa_fwd.wavelet_tree.lex_count(l_fwd_bwt, r_fwd_bwt+1, cc);
  unsigned long rank_l = std::get<0>(r_s_b);
  unsigned long s = std::get<1>(r_s_b), b = std::get<2>(r_s_b);
  // correct the smaller value if, the sentinel is in the interval!
  s = s + (l_fwd < csa_fwd.sentinel_pos && r_fwd >= csa_fwd.sentinel_pos);
  unsigned long rank_r = r_fwd - l_fwd - s - b + rank_l;
  l_fwd_res = c_begin + rank_l;
  r_fwd_res = c_begin + rank_r;
  assert(r_fwd_res+1 >= l_fwd_res);
  l_bwd_res = l_bwd + s;
  r_bwd_res = r_bwd - b;
  assert(r_bwd_res-l_bwd_res == r_fwd_res-l_fwd_res);
  return r_fwd_res+1-l_fwd_res;
}

// return all smems that cover the position start_pos
unsigned long smems(csa_t& fwd_csa, csa_t& bwd_csa, std::string& p, unsigned long l,
                    unsigned long start_pos, std::vector<path_item>& results)
{
  std::vector<path_item> bwd_prev, bwd_curr;
  unsigned long bwd_left, bwd_right, fwd_left, fwd_right;
  unsigned long bwd_new_left, bwd_new_right, fwd_new_left, fwd_new_right, q;
  unsigned long length_interval = 0;
  unsigned long forward_pos = 0;
  bool rev = false;
  unsigned long p_len = p.size();
  // initialize the first interval with the interval of the first character
  bwd_left = 0;
  bwd_right = bwd_csa.size() - 1;
  fwd_left = 0;
  fwd_right = fwd_csa.size() - 1;
  // fwd_cst is of the reverse of the sequence
  unsigned long old_length =  bidirectional_search_new(fwd_csa, fwd_left, fwd_right, bwd_left, bwd_right,
                              get_char(p[idx(start_pos, p_len, rev)], rev),
                              fwd_new_left, fwd_new_right, bwd_new_left, bwd_new_right);
  q = 1;
  bwd_left = bwd_new_left;
  bwd_right = bwd_new_right;
  fwd_left = fwd_new_left;
  fwd_right = fwd_new_right;
  for (unsigned long i = start_pos + 1; i < p.size(); i++) {
    if (valid(p[idx(i, p_len, rev)])) {
      // bidirectional search
      length_interval =  bidirectional_search_new(fwd_csa, fwd_left, fwd_right, bwd_left, bwd_right,
                         get_char(p[idx(i, p_len, rev)], rev),
                         fwd_new_left, fwd_new_right, bwd_new_left, bwd_new_right);
      if (length_interval != old_length) {
        bwd_prev.push_back(path_item(q, bwd_left, bwd_right, i - 1));
      }
      if (length_interval == 0) {
        break;
      }
      q += 1;
      old_length = length_interval;
      fwd_left = fwd_new_left;
      fwd_right = fwd_new_right;
      bwd_left = bwd_new_left;
      bwd_right = bwd_new_right;
    } else {
      // found an 'N'
      bwd_prev.push_back(path_item(q, bwd_left, bwd_right, i - 1));
      break;
    }
  }

  // store the last entry!
  if ((q + start_pos) == p.size()) {
    bwd_prev.push_back(path_item(q, bwd_left, bwd_right, p.size() - 1));
  }

  // this is the last position in the query sequence, that is covered
  // by all these SMEMs. This position will be returned
  forward_pos = q + start_pos;

  //the smaller intervals btw the longer prefixes should be visited
  //first.
  std::reverse(bwd_prev.begin(), bwd_prev.end());

  // search with backward-search
  int end_pos = (int) p.size();
  length_interval = 0;
  old_length = 0;
  q = 0;
  for (int i = start_pos - 1; i >= -1; i--) {
    bwd_curr.clear();
    old_length = 0;
    // iterate over all possible lcps and extend them.
    for (auto iter = bwd_prev.begin(); iter != bwd_prev.end(); iter++) {
      bwd_left = (*iter).lb;
      bwd_right = (*iter).rb;

      if (i == -1 || !valid(p[idx(i, p_len, rev)])) {
        if (bwd_curr.empty() && i < end_pos) {
          // found an smem
          end_pos = i;
          if ((*iter).len >= (long) l) {
            results.push_back(path_item((*iter).len, bwd_left, bwd_right, i + 1));
          }
        }
        continue;
      }
      length_interval = backward_search(bwd_csa, bwd_left, bwd_right,
                                        get_char(p[idx(i, p_len, rev)], rev),
                                        bwd_new_left, bwd_new_right);
      if (length_interval == 0) {
        if (bwd_curr.empty() && i < end_pos) {
          // found an smem
          end_pos = i;
          if ((*iter).len >= (long) l) {
            results.push_back(path_item((*iter).len, bwd_left, bwd_right, i + 1));
          }
        }
      }
      if (length_interval != 0 && length_interval != old_length) {
        old_length = length_interval;
        bwd_curr.push_back(path_item((*iter).len + 1, bwd_new_left, bwd_new_right, i));
      }
    }
    if (bwd_curr.empty()) {
      break;
    }
    bwd_prev = bwd_curr;
  }
  return forward_pos;
}

// computes all smems for the sequence p, that are longer than l
std::vector<path_item> compute_smems_for_seq(csa_t& bwd_csa, csa_t& fwd_csa, std::string& p, unsigned long l)
{
  unsigned long x = 0;
  std::vector<path_item> results;
  while (x < p.size()) {
    if (valid(p[x])) {
      x = smems(fwd_csa, bwd_csa, p, l, x, results);
    } else
      x++;
  }

  return results;
}

// computes all mems in an iterativ fashion.
std::vector<mem_struct> find_pattern(cst_t& cst, std::string& p, unsigned long l,
                                     bool mums, bool smems, unsigned long max_val, bool rev)
{
  unsigned long p_len = p.size();
  unsigned long curr_pos = p_len;
  unsigned long left = 0, right = 0, new_left = 0, new_right = 0, q = 0;
  unsigned long length_intervall = 0;
  std::vector<path_item> path;
  std::vector<mem_struct> mems;

  left = 0;
  right = cst.size() - 1;
  q = 0;

  while (curr_pos >= 1) {
    path.clear();

    new_left = left;
    new_right = right;

    // check if p[curr_pos-1] is N or n. if yes, iterate until a
    // not N is found!
    while (!valid(get_char(p[idx(curr_pos - 1, p_len, rev)], rev))) {
      curr_pos -= 1;
    }
    length_intervall = backward_search(cst.csa, left, right, get_char(p[idx(curr_pos - 1, p_len, rev)], rev), new_left, new_right);
    bool is_legal = false;

    while (length_intervall > 0) {
      is_legal = true;
      q += 1;
      if (q >= l && !smems) {
        // if we search mums -> get the mum-candidates, if
        // interval length is exactly 1
        if (mums) {
          if (length_intervall == 1)
            path.push_back(path_item(q, new_left, new_right, curr_pos));
        } else {
          path.push_back(path_item(q, new_left, new_right, curr_pos));
        }
      }

      left = new_left;
      right = new_right;
      curr_pos -= 1;

      if (curr_pos >= 1) {
        // check if p[curr_pos-1] is N or n.  start with a
        // completely new intervall -> break
        if (!valid(get_char(p[idx(curr_pos - 1, p_len, rev)], rev))) {
          left = 0;
          right = cst.size() - 1;
          q = 0;
          break;
        }
        length_intervall = backward_search(cst.csa, left, right, get_char(p[idx(curr_pos - 1, p_len, rev)], rev), new_left, new_right);
      } else {
        break;
      }
    }

    // a parent-interval can be pushed here! check, if it is a real interval!
    if (smems && q >= l && is_legal)
      path.push_back(path_item(q, left, right, curr_pos + 1));

    for (auto it = path.begin(); it != path.end(); ++it) {
      unsigned long last_lb, last_rb, lb, rb, pos, qs;

      qs = it->len;
      lb = it->lb;
      rb = it->rb;
      pos = it->p;

      // if the length of the interval is greater than the given
      // max_value don't enumerate mems for this interval
      if ((it->rb - it->lb) > max_val)
        continue;

      last_lb = lb;
      last_rb = lb - 1;

      while (qs >= l) {
        for (unsigned long k = lb; k < last_lb; k++) {
          if (pos == 1 || cst.csa.bwt[k] != get_char(p[idx(pos - 2, p_len, rev)], rev)) {
            mems.push_back(mem_struct(cst.csa[k], pos - 1, qs));
          }
        }
        // for last_rb = lb - 1, this loop is the same as lb <= k <= rb
        for (unsigned long k = last_rb + 1; k <= rb; k++) {
          if (pos == 1 || cst.csa.bwt[k] != get_char(p[idx(pos - 2, p_len, rev)], rev)) {
            mems.push_back(mem_struct(cst.csa[k], pos - 1, qs));
          }
        }

        last_lb = lb;
        last_rb = rb;

        // don't enumerate all the shorter matches, if we just
        // search for MUMs or SMEMs
        if (mums || smems) {
          break;
        }

        auto parent = cst.parent(cst.node(lb, rb));
        lb = cst.lb(parent);
        rb = cst.rb(parent);
        qs = cst.depth(parent);
        // if the length of the interval is greater than the given
        // max_value stop enumerating mems
        if ((it->rb - it->lb) > max_val)
          break;
      }
    }

    if (left == 0 && right == (cst.csa.size() - 1)) {
      curr_pos -= 1;
    } else {
      auto parent = cst.parent(cst.node(left, right));
      left = cst.lb(parent);
      right = cst.rb(parent);
      q = cst.depth(parent);
    }
  }

  return mems;
}
