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

#ifndef PROCESS_TEXT_H
#define PROCESS_TEXT_H

#include "sdsl_types.hpp"
#include "utility.hpp"

using namespace sdsl;

void process_text(cst_huff& fm_index, std::string& pattern, unsigned long l, std::string& header,
                  std::vector<seq_info>& info,
                  bool silent, bool mums, bool smems, unsigned long max_val, bool n_filter,
                  bool rev, bool comp, bool dna, FILE* out_printf, bool clean_mums);

// preprocesses the multifasta file or any other file into the right
// format. stores only the offset information
int preprocess_text(std::string& old_file, std::string& new_file,
                    std::vector<seq_info>& info);
#endif
