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

#include <cassert>
#include "utility.hpp"

// m1 comes before m2, if m1.start_ref < m2.start_ref or if m1 is larger
bool compare_mem_struct_wrt_ref(std::pair<std::string, mem_struct>& m1, std::pair<std::string, mem_struct>& m2)
{
  if (m1.second.abs_pos == m2.second.abs_pos)
    return (m1.second.len > m2.second.len);
  else
    return (m1.second.abs_pos < m2.second.abs_pos);
}

// m1 comes before m2, if m1.start_query < m2.start_query or if m1 is larger
bool compare_mem_struct_wrt_query(std::pair<std::string, mem_struct>& m1, std::pair<std::string, mem_struct>& m2)
{
  if (m1.second.start_query == m2.second.start_query)
    return (m1.second.len > m2.second.len);
  else
    return (m1.second.start_query < m2.second.start_query);
}

// p1 comes before p2, if p1.lb < p2.lb or p1.rb < p2.rb
bool compare_path_item(path_item& p1, path_item& p2)
{
  return (p1.lb < p2.lb);
}

unsigned long idx(unsigned long i, unsigned long len, bool rev)
{
  return (rev? (len - i - 1): i);
}

char get_char(char c, bool rev, bool dna)
{
  if (rev && dna) {
    switch (c) {
      case 'A': case 'a':
        return 't';
      case 'C': case 'c':
        return 'g';
      case 'G': case 'g':
        return 'c';
      case 'T': case 't':
        return 'a';
      default:
        return (char) tolower((int) c);
    }
  } else
    return (char) tolower((int) c);
}

bool valid(char c, bool dna)
{
  if (!dna)
    return true;
  else
    return (c == 'A' || c == 'a' ||
            c == 'C' || c == 'c' ||
            c == 'G' || c == 'g' ||
            c == 'T' || c == 't');
}

// scan the mem m for contained borders wrt. ref or
// query. split the mem at each found border and push the mems
// into a return vector
std::vector<std::pair<std::string, mem_struct>> split_mem_at_border(mem_struct m,
    std::vector<seq_info>& info,
    int pos_in_info,
    bool ref, unsigned long offset)
{
  std::vector<std::pair<std::string, mem_struct>> ret;
  unsigned long start, other_start, end;
  if (ref) {
    start = m.start_ref + offset;
    other_start = m.start_query;
  } else {
    start = m.start_query + offset;
    other_start = m.start_ref;
  }
  end = start + m.len - 1;
  // iterate backward until m is not in the sequence
  int pos_before = pos_in_info - 1;
  while (pos_before > 1 && info[pos_before].start > start) {
    // m contains this border -> split it
    std::pair<std::string, mem_struct> p;
    p.first = info[pos_before - 1].header;
    unsigned long start_diff = 0;
    if (info[pos_before - 1].start > start) {
      start_diff = info[pos_before - 1].start - start;
    }
    unsigned long new_start = start + start_diff - info[pos_before - 1].start;
    unsigned long new_len = info[pos_before - 1].end - start + start_diff;
    if (ref)
      p.second = mem_struct(new_start, other_start + start_diff, new_len, start + start_diff);
    else
      p.second = mem_struct(other_start + start_diff, new_start, new_len, other_start + start_diff);
    ret.push_back(p);
    pos_before -= 1;
  }
  std::pair<std::string, mem_struct> p;
  if (pos_in_info == 0) {
    assert(!"reached");
    // TODO: do i need this special case? the first sequence seperator can't be contained in a MEM!!!
    // TODO: later!
    p.first = info[pos_in_info].header;
    if (ref)
      p.second = mem_struct(start, other_start, info[pos_in_info].end - start, start);
    else
      p.second = mem_struct(other_start, start, info[pos_in_info].end - start, other_start);
  } else {
    p.first = info[pos_in_info - 1].header;
    unsigned long start_diff = 0;
    if (info[pos_in_info - 1].start > start) {
      start_diff = info[pos_in_info - 1].start - start;
    }
    unsigned long new_start = start + start_diff - info[pos_in_info - 1].start;
    unsigned long new_len = info[pos_in_info - 1].end - start + start_diff;
    if (ref)
      p.second = mem_struct(new_start, other_start + start_diff, new_len, start + start_diff);
    else
      p.second = mem_struct(other_start + start_diff, new_start, new_len, other_start + start_diff);
  }
  ret.push_back(p);

  // iterate forward until m is not in the sequence
  int pos_after = pos_in_info;
  while (pos_after < (int)(info.size() - 1) && info[pos_after + 1].start <= end) {
    // m contains this border -> split it
    std::pair<std::string, mem_struct> p;
    p.first = info[pos_after].header;
    unsigned long start_diff = 0;
    if (info[pos_after].start > start)
      start_diff = info[pos_after].start - start;
    unsigned long new_start = start + start_diff - info[pos_after].start;
    unsigned long new_len = info[pos_after].end - info[pos_after].start;
    if (ref)
      p.second = mem_struct(new_start, other_start + start_diff, new_len, start + start_diff);
    else
      p.second = mem_struct(other_start + start_diff, new_start, new_len, other_start + start_diff);
    ret.push_back(p);
    pos_after += 1;
  }
  // process the last bit of the mem!
  if (pos_after < (int) info.size() && info[pos_after].start >= start) {
    std::pair<std::string, mem_struct> p;
    p.first = info[pos_after].header;
    unsigned long start_diff = 0;
    if (info[pos_after].start > start)
      start_diff = info[pos_after].start - start;
    unsigned long new_start = start + start_diff - info[pos_after].start;
    unsigned long new_len = end - info[pos_after].start + 1;
    if (ref)
      p.second = mem_struct(new_start, other_start + start_diff, new_len, start + start_diff);
    else
      p.second = mem_struct(other_start + start_diff, new_start, new_len, other_start + start_diff);
    ret.push_back(p);
  } else
    assert(!"reached");
  return ret;
}

// this function corrects the given MEM (start_ref, start_query, len) w.r.t.
// info. info stores the absolute startpositions of the sequence parts of the reference.
// the start_ref is corrected to the relative startposition in the correpsponding
// sequence part. also is checked, if the MEM is crossing a sequence boundary. in that
// case the startpositions (ref and query) and the length are corrected. the correct MEM is returned
// with the additional information of the sequence-identifier.
std::vector<std::pair<std::string, mem_struct>>mem_corrected_vec(mem_struct m,
    std::vector<seq_info>& info,
    bool ref, unsigned long offset)
{
  // search through info with a binary search
  std::vector<std::pair<std::string, mem_struct>> ret;
  std::pair<std::string, mem_struct> pr;
  long info_size = info.size();
  if (info_size == 1) {
    pr.first = "";
    if (ref)
      pr.second = mem_struct(m.start_ref + offset, m.start_query, m.len, m.start_ref + offset);
    else
      pr.second = mem_struct(m.start_ref, m.start_query + offset, m.len, m.start_ref);
    ret.push_back(pr);
    return ret;
  }
  unsigned long start, other_start;
  if (ref) {
    start = m.start_ref + offset;
    other_start = m.start_query;
  } else {
    start = m.start_query + offset;
    other_start = m.start_ref;
  }
  unsigned long end = start + m.len - 1;

  long left = 0;
  long right = info_size - 1;

  while (left <= right) {
    long middle = left + ((right - left) / 2);
    seq_info p = info[middle];

    if (start < p.start && end >= p.start) {
      return split_mem_at_border(m, info, middle, ref, offset);
    } else if (start >= p.start && start < p.end && end >= p.end) { // maybe check that middle < info.size() - 1
      return split_mem_at_border(m, info, middle + 1, ref, offset);
    } else if (middle == 0 && start < p.end && end >= p.end) {
      return split_mem_at_border(m, info, middle + 1, ref, offset);
    } else if (middle == (long)(info.size() - 1) && start < p.start) {
      return split_mem_at_border(m, info, middle, ref, offset);
    } else if (start >= p.start && end < p.end) {
      pr.first = p.header;
      unsigned long new_start = start - p.start;
      if (ref)
        pr.second = mem_struct(new_start, other_start, m.len, start);
      else
        pr.second = mem_struct(other_start, new_start, m.len, other_start);
      ret.push_back(pr);
      break;
    } else if (start < p.start) {
      right = middle - 1;
    } else if (start >= p.end) {
      left = middle + 1;
    } else {
      assert(!"reached");
      break;
    }
  }
  return ret;
}


std::pair<std::string, mem_struct> get_header(mem_struct m,
    std::vector<std::pair<std::string, unsigned long>>& info,
    bool ref, unsigned long offset)
{
  if (info.size() == 1) {
    std::pair<std::string, mem_struct> p;
    p.first = info[0].first;
    p.second = m;
    return p;
  }

  unsigned long start;
  if (ref)
    start = m.start_ref + offset;
  else
    start = m.start_query + offset;

  if (start < info[0].second) {
    std::pair<std::string, mem_struct> p;
    p.first = info[0].first;
    p.second = m;
    return p;
  }

  // binary search through info
  long left = 0;
  long right = info.size() - 1;

  while (left <= right) {
    long middle = left + ((right - left) / 2);
    std::pair<std::string, unsigned long> p = info[middle];
    if (start < p.second) {
      // check if start is in sequence[middle]
      if (middle > 0 && (start >= info[middle - 1].second)) {
        std::pair<std::string, mem_struct> ret;
        ret.first = p.first;
        if (ref)
          ret.second = mem_struct(start - info[middle - 1].second, m.start_query, m.len, m.start_ref + offset);
        else
          ret.second = mem_struct(m.start_ref, start - info[middle - 1].second, m.len, m.start_ref);
        return ret;
      } else {
        // search the lower parts
        right = middle - 1;
      }
    } else {
      if (middle < ((long) info.size() - 1) && (start < info[middle + 1].second)) {
        // start is in info[middle + 1]
        std::pair<std::string, mem_struct> ret;
        ret.first = info[middle + 1].first;
        if (ref)
          ret.second = mem_struct(start - p.second, m.start_query, m.len, m.start_ref + offset);
        else
          ret.second = mem_struct(m.start_ref, start - p.second, m.len, m.start_ref);
        return ret;
      } else {
        // search the greater parts
        left = middle + 1;
      }
    }
  }

  std::pair<std::string, mem_struct> p;
  p.first = "";
  p.second = mem_struct();
  return p;
}

// this function corrects the given MEM (start_ref, start_query, len) w.r.t.
// info. info stores the absolute startpositions of the sequence parts of the reference.
// the start_ref is corrected to the relative startposition in the correpsponding
// sequence part. also is checked, if the MEM is crossing a sequence boundary. in that
// case the startpositions (ref and query) and the length are corrected. the correct MEM is returned
// with the additional information of the sequence-identifier.
std::vector<std::pair<std::string, mem_struct>>mem_corrected(mem_struct m,
    std::vector<seq_info>& info,
    bool ref, unsigned long offset)
{

  // search through info with a binary search
  std::pair<std::string,mem_struct> ret;
  std::vector<std::pair<std::string, mem_struct>> results;
  long info_size = info.size();
  if (ref) {
    if (info_size == 1) {
      ret.first = "";
      if (ref)
        ret.second = mem_struct(m.start_ref + offset, m.start_query, m.len, m.start_ref);
      else
        ret.second = mem_struct(m.start_ref, m.start_query + offset, m.len, m.start_ref);
      results.push_back(ret);
      return results;
    }
  }

  unsigned long start, other_start;
  if (ref) {
    start = m.start_ref + offset;
    other_start = m.start_query;
  } else {
    start = m.start_query + offset;
    other_start = m.start_ref;
  }

  unsigned long end = start + m.len - 1;
  if (start < info[0].end) {
    if (end >= info[0].end) {
      // in different sequence parts
      unsigned long start_diff = info[0].end - start;
      unsigned long end_diff = end - info[0].end;
      std::pair<std::string, mem_struct> ret1;
      ret1.first = info[1].header;
      if (ref)
        ret1.second = mem_struct(start + start_diff - info[0].end, other_start + start_diff,
                                 m.len - start_diff, start + start_diff);
      else
        ret1.second = mem_struct(other_start + start_diff, start + start_diff - info[0].end,
                                 m.len - start_diff, other_start + start_diff);
      std::pair<std::string, mem_struct> ret2;
      ret2.first = info[0].header;
      if (ref)
        ret2.second = mem_struct(start, other_start, m.len - end_diff - 1, start);
      else
        ret2.second = mem_struct(other_start, start, m.len - end_diff - 1, other_start);
      results.push_back(ret1);
      results.push_back(ret2);
    } else {
      std::pair<std::string, mem_struct> ret;
      ret.first = info[0].header;
      if (ref)
        ret.second = mem_struct(start, other_start, m.len, start);
      else
        ret.second = mem_struct(other_start, start, m.len, other_start);
      results.push_back(ret);
    }
    return results;
  }

  long left = 0;
  long right = info_size - 1;

  while (left <= right) {
    long middle = left + ((right - left) / 2);
    seq_info p = info[middle];
    if (start < p.end) {
      // check if start is in sequence[middle]
      if (middle > 0 && (start >= info[middle - 1].end)) {
        if (end >= p.end) {
          // start and end are in different sequence parts
          std::pair<std::string, mem_struct> ret1;
          unsigned long start_diff = p.end - start;
          unsigned long end_diff = end - p.end;
          ret1.first = info[middle + 1].header;
          if (ref)
            ret1.second = mem_struct(start + start_diff - p.end, other_start + start_diff,
                                     m.len - start_diff, m.start_ref + offset + start_diff);
          else
            ret1.second = mem_struct(other_start + start_diff, start + start_diff - p.end,
                                     m.len - start_diff, other_start);

          std::pair<std::string, mem_struct> ret2;
          ret2.first = info[middle].header;
          if (ref)
            ret2.second = mem_struct(start - info[middle - 1].end, other_start, m.len - end_diff - 1, m.start_ref);
          else
            ret2.second = mem_struct(other_start, start - info[middle - 1].end, m.len - end_diff - 1, other_start);

          results.push_back(ret1);
          results.push_back(ret2);
        } else {
          // start and end on the same sequence part info[middle] = p
          std::pair<std::string, mem_struct> ret;
          ret.first = p.header;
          if (ref)
            ret.second = mem_struct(start - info[middle - 1].end, other_start, m.len, m.start_ref + offset);
          else
            ret.second = mem_struct(other_start, start - info[middle - 1].end, m.len, m.start_ref);
          results.push_back(ret);
        }
        break;
      } else {
        // search the lower parts
        right = middle - 1;
      }
    } else {
      if (middle < (info_size - 1) && (start < info[middle + 1].end)) {
        if (end >= info[middle + 1].end) {
          // start and end are on different sequence parts
          unsigned long end_diff = end - info[middle + 1].end;
          unsigned long start_diff = info[middle + 1].end - start;
          // shift start to the corresponding sequence info[middle + 2]
          std::pair<std::string, mem_struct> ret1;
          ret1.first = info[middle + 2].header;
          if (ref)
            ret1.second = mem_struct(start + start_diff - info[middle + 1].end, other_start + start_diff,
                                     m.len - start_diff, m.start_ref + offset + start_diff);
          else
            ret1.second = mem_struct(other_start + start_diff, start + start_diff - info[middle + 1].end,
                                     m.len - start_diff, m.start_ref + start_diff);
          std::pair<std::string, mem_struct> ret2;
          // shift end to the corresponding sequence info[middle + 1]
          ret2.first = info[middle + 1].header;
          if (ref)
            ret2.second = mem_struct(start - info[middle].end, other_start, m.len - end_diff - 1, m.start_ref + offset);
          else
            ret2.second = mem_struct(other_start, start - info[middle].end, m.len - end_diff - 1, m.start_ref);

          results.push_back(ret1);
          results.push_back(ret2);
        } else {
          // found the correspondin sequence in info[middle + 1]
          std::pair<std::string, mem_struct> ret;
          ret.first = info[middle + 1].header;
          if (ref)
            ret.second = mem_struct(start - info[middle].end, other_start, m.len, m.start_ref + offset);
          else
            ret.second = mem_struct(other_start, start - info[middle].end, m.len, m.start_ref);
          results.push_back(ret);
        }
        break;
      } else {
        // search the greater parts
        left = middle + 1;
      }
    }
  }
  return results;
}

std::vector<std::pair<std::string, mem_struct>> clean_mumcands(std::vector<std::pair<std::string, mem_struct>> mums)
{
  std::vector<std::pair<std::string, mem_struct>> tmp;
  std::pair<std::string, mem_struct> prev, curr;
  // Adapted from Stephan Kurtz's code in cleanMUMcand.c in MUMMer v3.20.
  long currentright, dbright = 0;
  bool ignorecurrent, ignoreprevious = false;
  for (auto mem_r = mums.begin(); mem_r != mums.end(); mem_r++) {
    prev = curr;
    curr = *mem_r;
    ignorecurrent = false;
    currentright = curr.second.abs_pos + curr.second.len - 1;
    if (dbright > currentright)
      ignorecurrent = true;
    else {
      if (dbright == currentright) {
        ignorecurrent = true;
        if (!ignoreprevious && prev.second.abs_pos == curr.second.abs_pos)
          ignoreprevious = true;
      } else {
        dbright = currentright;
      }
    }
    if (mem_r != mums.begin() && !ignoreprevious) {
      tmp.push_back(prev);
    }
    ignoreprevious = ignorecurrent;
  }
  if (!ignoreprevious && !mums.empty()) {
    tmp.push_back(mums.back());
  }
  return tmp;
}

bool contained(mem_struct& m1, mem_struct& m2)
{
  unsigned long s1 = m1.start_query;
  unsigned long e1 = m1.start_query + m1.len - 1;
  unsigned long s2 = m2.start_query;
  unsigned long e2 = m2.start_query + m2.len - 1;
  return (s1 <= s2 && e1 >= e2 && (s1 != s2 || e1 != e2));
}

// returns the pair (lpos, rpos) where lpos is the first position in the sequence
// of the block of Ns and rpos is the last position the block of Ns.
// returns also the interval (lborder, rborder) in the ns array.
std::pair<std::pair<long, long>,std::pair<unsigned long, unsigned long>> get_n_borders(long start, std::vector<unsigned long>& ns)
{
  unsigned long lnpos = ns[start];
  long lpos = start;
  while ((lpos > 0) && ((lnpos - 1) == ns[lpos - 1])) {
    lpos -= 1;
    lnpos -= 1;
  }
  long rpos = start;
  unsigned long rnpos = ns[start];
  while ((rpos < (long)(ns.size() - 1)) && ((rnpos + 1) == ns[rpos + 1])) {
    rpos += 1;
    rnpos += 1;
  }

  std::pair<long, long> positions(lpos, rpos);
  std::pair<unsigned long, unsigned long> interval(lnpos, rnpos);
  return std::pair<std::pair<long, long>,std::pair<unsigned long, unsigned long>> (positions, interval);
}

std::vector<mem_struct> scan_mem(mem_struct m, std::vector<unsigned long>& ns, bool ref, long pos)
{
  std::vector<mem_struct> ret;
  long last_index = ns.size() - 1;
  unsigned long pos_of_n = ns[pos];

  unsigned long start;
  if (ref)
    start = m.start_ref;
  else
    start = m.start_query;
  unsigned long end = start + m.len - 1;

  std::pair<std::pair<long, long>, std::pair<unsigned long, unsigned long>> borders = get_n_borders(pos, ns);
  // whole mem contained in N-block
  if (borders.second.first <= start && borders.second.second >= end)
    return ret;

  // scan from pos to the start
  unsigned long new_end = borders.second.first - 1;
  pos = borders.first.first - 1;
  if (pos >= 0) {
    pos_of_n = ns[pos];

    while (pos > 0 && pos_of_n >= start) {
      borders = get_n_borders(pos, ns);
      if (borders.second.second >= start) {
        // next N-block starts within mem
        unsigned long start_diff = borders.second.second - start;
        ret.push_back(mem_struct(m.start_ref + start_diff + 1,
                                 m.start_query + start_diff + 1,
                                 new_end - (start + start_diff) + 1));
      }
      pos = borders.first.first - 1;
      pos_of_n = ns[pos];
      new_end = borders.second.first - 1;
    }
  }
  if (new_end >= start) {
    ret.push_back(mem_struct(m.start_ref, m.start_query, new_end - start + 1));
  }

  // scan from pos to the end
  unsigned long new_start = borders.second.second + 1;
  pos = borders.first.second + 1;
  if (pos <= last_index) {
    pos_of_n = ns[pos];

    while (pos < last_index && pos_of_n <= end) {
      borders = get_n_borders(pos, ns);
      if (borders.second.first <= end) {
        // next N-block starts within mem
        unsigned long start_diff = new_start - start;
        ret.push_back(mem_struct(m.start_ref + start_diff,
                                 m.start_query + start_diff,
                                 borders.second.first - (start + start_diff) + 1));
      }
      pos = borders.first.second + 1;
      pos_of_n = ns[pos];
      new_start = borders.second.second + 1;
    }
  }
  if (new_start <= end) {
    unsigned long start_diff = new_start - start;
    ret.push_back(mem_struct(m.start_ref + start_diff,
                             m.start_query + start_diff,
                             end - (start + start_diff) + 1));
  }

  return ret;
}

// this function splits the mem at 'N'-symbols, if contained
std::vector<mem_struct> mem_without_n_bv(mem_struct& m, bitvector_t::rank_1_type& rank_1,
    bitvector_t::select_1_type& select_1, bool ref)
{
  unsigned long start, other_start;
  if (ref) {
    start = m.start_ref;
    other_start = m.start_query;
  } else {
    start = m.start_query;
    other_start = m.start_ref;
  }

  unsigned long end = start + m.len - 1;
  std::vector<mem_struct> res;

  // check if a 'N'-symbol is contained in m
  unsigned long before = rank_1(start);
  unsigned long after = rank_1(end + 1);

  if (before < after) {
    // process the first mem part
    unsigned long first_n, last_pos, next_pos;
    first_n = select_1(before + 1);
    if (first_n > start) {
      if (ref)
        res.push_back(mem_struct(start, other_start, first_n - start));
      else
        res.push_back(mem_struct(other_start, start, first_n - start));
      last_pos = first_n;
    } else {
      last_pos = start;
    }
    next_pos = end;
    // process all mem parts in the middle
    for (int i = 2; i <= ((int) after - (int) before); i++) {
      next_pos = select_1(before + i);

      if (next_pos > last_pos) {
        // found a possible mem [last_pos, next_pos)
        unsigned long start_diff = last_pos - start + 1;
        if (ref)
          res.push_back(mem_struct(start + start_diff, other_start + start_diff, next_pos - last_pos - 1));
        else
          res.push_back(mem_struct(other_start + start_diff, start + start_diff, next_pos - last_pos - 1));
        last_pos = next_pos;
      }
    }
    // process the last mem part
    if (last_pos < end) {
      unsigned long start_diff = last_pos - start + 1;
      if (ref)
        res.push_back(mem_struct(start + start_diff, other_start + start_diff, end - last_pos));
      else
        res.push_back(mem_struct(other_start + start_diff, start + start_diff, end - last_pos));
    }
  } else {
    res.push_back(m);
  }
  return res;
}

std::vector<mem_struct> mem_without_n_it(mem_struct m, std::vector<unsigned long>& ns, bool ref)
{
  long last_index = ns.size() - 1;
  long lborder = 0;
  long rborder = last_index;

  unsigned long new_start;
  if (ref)
    new_start = m.start_ref;
  else
    new_start = m.start_query;
  unsigned long new_end = new_start + m.len - 1;

  while (lborder <= rborder) {
    if (ns[rborder] < new_start || ns[lborder] > new_end)
      break;
    long middle = lborder + ((rborder - lborder) / 2);
    unsigned long pos_of_n = ns[middle];

    if (pos_of_n == new_end) {
      return scan_mem(m, ns, ref, middle);
    } else if (pos_of_n == new_start) {
      return scan_mem(m, ns, ref, middle);
    } else if (pos_of_n < new_end) {
      if (pos_of_n > new_start) {
        return scan_mem(m, ns, ref, middle);
      } else {
        lborder = middle + 1;
      }
    } else if (pos_of_n > new_end) {
      rborder = middle - 1;
    } else
      assert(!"reached");
  }
  std::vector<mem_struct> ret;
  ret.push_back(m);
  return ret;
}
