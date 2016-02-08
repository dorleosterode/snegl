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

#ifndef SDSL_TYPES_H
#define SDSL_TYPES_H

#include <sdsl/suffix_arrays.hpp>
#include <sdsl/suffix_trees.hpp>
#include <sdsl/bit_vectors.hpp>

// typedefs for sdsl types

using namespace sdsl;

// typedef for the CST used to process normal text-data
typedef cst_sct3<csa_wt<wt_huff<>, 16, 1<<20, sa_order_sa_sampling<>>, lcp_dac<>> cst_huff;

// typedef for the used fm_index
#if 1
typedef cst_sct3<csa_wt<wt_dna<>, 16, 1<<20, sa_order_sa_sampling<>, isa_sampling<>, byte_alphabet, true>, lcp_dac<>> cst_t;
#else
typedef cst_sct3<csa_wt<wt_dna_il<7>, 16, 1<<20, sa_order_sa_sampling<>, isa_sampling<>, byte_alphabet, true>, lcp_dac<>> cst_t;
#endif

// typedef for the use of csa in the fermi-algorithm:
typedef csa_wt<wt_dna<>, 16, 1<<20, sa_order_sa_sampling<>, isa_sampling<>, byte_alphabet, true> csa_t;

//typedef for the used bitvector
typedef rrr_vector<> bitvector_t;

#endif
