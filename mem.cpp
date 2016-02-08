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

#include <fstream>
#include <stdlib.h>
#include <random>
#include <argp.h>
#include <sys/stat.h>
#include <algorithm>
#include "process_text.hpp"
#include "process_dna.hpp"
#include "sdsl_types.hpp"
#include "utility.hpp"
// for parallelisation
#include <thread>

using namespace sdsl;

template<class index_type = cst_t>
int construct_index(index_type& idx, bitvector_t& bv_ns, bitvector_t::rank_1_type& bv_ns_rank,
                    bitvector_t::select_1_type& bv_ns_select, cache_config config,
                    std::vector<seq_info>& info,
                    std::string index_name, std::string index_bv,
                    std::string text_file, std::string proc_file,
                    bool normal_text, bool map_sa, bool fermi)
{
  int err = 0;

  // try to load the index:
  std::ifstream index(index_name);
  if (index.good()) {
    idx.load(index);
    printf("loaded the index\n");

    if (!normal_text) {
      // try to load sparse bit_vector with 'N'-information
      std::ifstream sd(index_bv);
      if (sd.good()) {
        bv_ns.load(sd);
        std::ifstream rank(index_bv + "rank_1");
        if (rank.good())
          bv_ns_rank.load(rank, &bv_ns);
        else
          util::init_support(bv_ns_rank, &bv_ns);
        std::ifstream select(index_bv + "select_1");
        if (select.good())
          bv_ns_select.load(select, &bv_ns);
        else
          util::init_support(bv_ns_select, &bv_ns);
        printf("loaded the bit_vector with 'N'-information\n");
        preprocess_seq(text_file, proc_file, info, false);
        if (!fermi)
          if (std::remove(proc_file.c_str()) != 0)
            fprintf(stderr, "Error deleting file: %s\n", strerror(errno));
      } else {
        bv_ns = preprocess_seq(text_file, proc_file, info, true);
        if (!fermi)
          if (std::remove(proc_file.c_str()) != 0)
            fprintf(stderr, "Error deleting file: %s\n", strerror(errno));
        util::init_support(bv_ns_rank, &bv_ns);
        util::init_support(bv_ns_select, &bv_ns);
        store_to_file(bv_ns, index_bv);
        store_to_file(bv_ns_rank, index_bv + "_rank_1");
        store_to_file(bv_ns_select, index_bv + "_select_1");
        printf("computed the bit_vector with 'N'-information\n");
      }
    } else {
      err = preprocess_text(text_file, proc_file, info);
      if (std::remove(proc_file.c_str()) != 0)
        fprintf(stderr, "Error deleting file: %s\n", strerror(errno));
      if (err)
        return err;
    }
  } else {
    if (!normal_text) {
      bv_ns = preprocess_seq(text_file, proc_file, info, true);
      util::init_support(bv_ns_rank, &bv_ns);
      util::init_support(bv_ns_select, &bv_ns);
      printf("computed the bit_vector with 'N'-information\n");
    } else {
      err = preprocess_text(text_file, proc_file, info);
      if (err)
        return err;
    }
    printf("constructing the index\n");
    if (map_sa) {
      construct(idx, proc_file, config, 1);
      for (auto file_pair : config.file_map) {
        if (file_pair.first == conf::KEY_SA)
          continue;
        sdsl::remove(file_pair.second);
      }
    } else
      construct(idx, proc_file, 1);

    store_to_file(idx, index_name);
    if (!normal_text) {
      store_to_file(bv_ns, index_bv);
      store_to_file(bv_ns_rank, index_bv + "_rank_1");
      store_to_file(bv_ns_select, index_bv + "_select_1");
    }
    printf("constructed the index\n");
    if (!fermi)
      if (std::remove(proc_file.c_str()) != 0)
        fprintf(stderr, "Error deleting file: %s\n", strerror(errno));
  }

  return err;
}

// A description of the arguments we accept.
static char args_doc[] = "REF.fas QUERY.fas";

// the following structs and functions are used for the command-line parsing
static struct argp_option options[] = {
  {"min_length", 'l', "LENGTH", 0, "the minimal length for a MEM", 1},
  {"output", 'o', "FILE", 0, "Output to FILE instead of standard output", 1 },
  {"threads", 't', "THREADS", 0, "Number of threads that are used to process the query file", 1 },
  {"silent", 's', 0, OPTION_ARG_OPTIONAL, "produce no output", 2 },
  {"rev", 'r', 0, OPTION_ARG_OPTIONAL, "search for reverse MEMs only", 2 },
  {"both", 'b', 0, OPTION_ARG_OPTIONAL, "search for forward and reverse MEMs", 2 },
  {"comp", 'c', 0, OPTION_ARG_OPTIONAL, "report the query-position of reverse matches relative to the original query", 2 },
  {"no-filter", -1, 0, OPTION_ARG_OPTIONAL, "inhibits the filtering step", 2 },
  {"mum", -2, 0, OPTION_ARG_OPTIONAL, "search for MUMs only", 2 },
  {"max-val", -3, "NUM", 0, "maximal number of MEM in REF.fas for which MEMs should be reported", 2 },
  {"smem", -4, 0, OPTION_ARG_OPTIONAL, "search for SMEMs only", 2 },
  {"fermi", 'f', 0, OPTION_ARG_OPTIONAL, "use the fermi-algo to compute smems", 2 },
  {"normal-text", -7, 0, OPTION_ARG_OPTIONAL, "compute MEMs for other text than DNA (i.e. protein sequences or real text).", 2 },
  {"dna", -8, 0, OPTION_ARG_OPTIONAL, "only for use with --normal-text. compute MEMs for DNA without processing 'N's.", 2 },
  {"sort-paths", -9, 0, OPTION_ARG_OPTIONAL, "store the SA-intervals. Just output the MEMs if --enum-mems is set.", 2 },
  {"enum-mems", -10, 0, OPTION_ARG_OPTIONAL, "enumerate and output the MEMs for the stored SA-intervals.", 2 },
  {"clean-mums", -11, 0, OPTION_ARG_OPTIONAL, "clean the MUM-candidates and output only real MUMs.", 2 },
  {"sort", -12, 0, OPTION_ARG_OPTIONAL, "sort the SA-intervals before outputting them. Only valid with --sort-paths and --enum-mems", 2 },
  {"map-sa", -13, 0, OPTION_ARG_OPTIONAL, "the complete SA is stored and used to enumerate MEMs. Only valid with --sort-paths and --enum-mems", 2 },
  {"thread-debug", -14, 0, OPTION_ARG_OPTIONAL, "for debugging only", 2 },
  {"no-clean-mums", -15, 0, OPTION_ARG_OPTIONAL, "don't clean the MUM-candidates and output also incorrect MUMs.", 2 },
  { 0 }
};

// Used by main to communicate with parse_opt.
struct arguments {
  char* args[2];                // ref-seq & query
  unsigned long min_length;
  unsigned long max_val;
  const char* output_file;
  int silent, n_filter, mum, smems, thread_debug;
  int rev, both, comp;
  int test_smems;
  int normal_text;
  int dna;
  int path, enum_mems, clean_mums, sort;
  int map_sa, no_clean_mums;
  unsigned long threads;
};

// Parse a single option.
static error_t
parse_opt(int key, char* arg, struct argp_state* state)
{
  // Get the input argument from argp_parse, which we know is a
  // pointer to our arguments structure.
  struct arguments* arguments = (struct arguments*) state->input;

  switch (key) {
    case 'l':
      arguments->min_length = strtoul(arg, NULL, 0);
      break;
    case 's':
      arguments->silent = 1;
      break;
    case 'o':
      arguments->output_file = arg;
      break;
    case 'r':
      arguments->rev = 1;
      break;
    case 'b':
      arguments->both = 1;
      break;
    case 'c':
      arguments->comp = 1;
      break;
    case 't':
      arguments->threads = strtoul(arg, NULL, 0);
      break;
    case -1:
      arguments->n_filter = 1;
      break;
    case -2:
      arguments->mum = 1;
      break;
    case -3:
      arguments->max_val = strtoul(arg, NULL, 0);
      break;
    case -4:
      arguments->smems = 1;
      break;
    case 'f':
      arguments->test_smems = 1;
      break;
    case -7:
      arguments->normal_text = 1;
      break;
    case -8:
      arguments->dna = 1;
      break;
    case -9:
      arguments->path = 1;
      break;
    case -10:
      arguments->enum_mems = 1;
      break;
    case -11:
      arguments->clean_mums = 1;
      break;
    case -12:
      arguments->sort = 1;
      break;
    case -13:
      arguments->map_sa = 1;
      break;
    case -14:
      arguments->thread_debug = 1;
      break;
    case -15:
      arguments->no_clean_mums = 1;
      break;

    case ARGP_KEY_ARG:
      if (state->arg_num >= 2)
        // too many arguments
        argp_usage(state);

      arguments->args[state->arg_num] = arg;

      break;

    case ARGP_KEY_END:
      if (state->arg_num < 2)
        // not enough arguments
        argp_usage(state);
      break;

    default:
      return ARGP_ERR_UNKNOWN;
  }
  return 0;
}

// argp struct for the parser
static struct argp argp = { options, parse_opt, args_doc };

// checks if valid arguments are given
int check_arguments(struct arguments* arguments)
{
  int err = 0;
  if ((arguments->threads > 1 || arguments->thread_debug)
      && (arguments->mum || arguments->smems ||
          arguments->test_smems || arguments->normal_text ||
          arguments->dna || arguments->path ||
          arguments->enum_mems || arguments->sort ||
          arguments->map_sa)) {
    fprintf(stderr, "WRONG ARGUMENTS: invalid option mum, smems, fermi, normal-text, dna, sort-paths, enum-mems, sort or map-sa with t > 1\n");
    err = 1;
  } else if (arguments->mum && (arguments->smems || arguments->test_smems)) {
    fprintf(stderr, "WRONG ARGUMENTS: can't search MUMs and SMEMs simultaneously\n");
    err = 1;
  } else if (arguments->smems && arguments->test_smems) {
    fprintf(stderr, "WRONG ARGUMENTS: can't search SMEMs with two strategies simultaneously\n");
    err = 1;
  } else if (arguments->test_smems) {
    if (arguments->normal_text) {
      fprintf(stderr, "WRONG ARGUMENTS: can't use fermi algorithm with wt huff\n");
      err = 1;
    } else if (arguments->map_sa || arguments->enum_mems || arguments->sort || arguments->path) {
      fprintf(stderr, "WRONG ARGUMENTS: can't use fermi algorithm with saving and processing suffix-intervals\n");
      err = 1;
    } else if (arguments->clean_mums) {
      fprintf(stderr, "WRONG ARGUMENTS: not valid option clean-mums for computing SMEMs\n");
      err = 1;
    }
  } else if (arguments->normal_text && (arguments->path || arguments->enum_mems
                                        || arguments->sort || arguments->map_sa)) {
    fprintf(stderr, "WRONG ARGUMENTS: can't use wt huff with saving and processing suffix-intervals\n");
    err = 1;
  } else if (arguments->dna && !arguments->normal_text) {
    fprintf(stderr, "WRONG ARGUMENTS: not valid option dna without normal-text\n");
    err = 1;
  } else if (arguments->enum_mems || arguments->sort || arguments->map_sa) {
    if (!arguments->path) {
      fprintf(stderr, "WRONG ARGUMENTS: not valid option without sort-paths\n");
      err = 1;
    }
  } else if ((arguments->clean_mums || arguments->no_clean_mums) && !arguments->mum) {
    fprintf(stderr, "WRONG ARGUMENTS: not valid optein clean-mums or no-clean-mums without mum\n");
    err = 1;
  }
  return err;
}

int main(int argc, char** argv)
{
  struct arguments arguments;

  // default values
  arguments.min_length = 20;
  arguments.max_val = ULONG_MAX;
  arguments.silent = 0;
  arguments.mum = 0;
  arguments.smems = 0;
  arguments.n_filter = 0;
  arguments.output_file = "output.sdsl.out";
  arguments.rev = false;
  arguments.both = false;
  arguments.comp = false;
  arguments.test_smems = 0;
  arguments.normal_text = 0;
  arguments.dna = 0;
  arguments.path = 0;
  arguments.enum_mems = 0;
  arguments.clean_mums = 0;
  arguments.no_clean_mums = 0;
  arguments.sort = 0;
  arguments.map_sa = 0;
  arguments.threads = 1;
  arguments.thread_debug = 0;

  // parse the arguments
  argp_parse(&argp, argc, argv, 0, 0, &arguments);

  if (check_arguments(&arguments) != 0)
    return EXIT_FAILURE;

  std::string text_file = arguments.args[0];
  std::string pattern_file = arguments.args[1];
  unsigned long l = arguments.min_length;
  unsigned long max_val = arguments.max_val;
  bool silent = arguments.silent;
  bool mums = arguments.mum;
  bool smems = arguments.smems;
  bool n_filter = arguments.n_filter;
  bool rev = arguments.rev;
  bool both = arguments.both;
  bool comp = arguments.comp;
  bool test_smems = arguments.test_smems;
  bool normal_text = arguments.normal_text;
  bool dna = arguments.dna;
  bool sort_paths = arguments.path;
  bool enum_mems = arguments.enum_mems;
  bool clean_mums = arguments.clean_mums;
  if (arguments.no_clean_mums)
    clean_mums = false;
  else
    clean_mums = true; // this is the default
  bool sort = arguments.sort;
  bool map_sa = arguments.map_sa;
  unsigned long thread_num = arguments.threads;

  bool thread_query = false;
  if (thread_num > 1 || arguments.thread_debug)
    thread_query = true;

  FILE* out_printf = fopen(arguments.output_file, "w");
  std::string out_file = arguments.output_file;

  cst_huff text_index;
  cst_t fm_index;
  csa_t fwd_index;
  csa_t bwd_index;

  std::string index_name = text_file + "_index_dna";
  std::string index_name_text = text_file + "_index_text";
  std::string index_name_smems = text_file + "_csa";
  std::string index_bv = text_file + "_bv_ns";

  printf("preprocess multi-fasta\n");
  std::string proc_file = text_file + "_new";
  std::vector<seq_info> info;
  // get a sparse bit_vector for the sequence
  bitvector_t sd_bv_ns;
  bitvector_t::rank_1_type sd_bv_ns_rank_1;
  bitvector_t::select_1_type sd_bv_ns_select_1;

  std::size_t found = proc_file.find_last_of("/");
  std::string proc_file_base = proc_file.substr(found + 1);
  cache_config config = cache_config(!map_sa); // if map_sa, the files shouldn't be deleted
  std::string sa_name = cache_file_name(conf::KEY_SA, config);

  int err = 0;
  if (test_smems)
    err = construct_index(fwd_index, sd_bv_ns, sd_bv_ns_rank_1,
                          sd_bv_ns_select_1, config, info, index_name_smems,
                          index_bv, text_file, proc_file,
                          false, false, true);
  else if (normal_text)
    err = construct_index(text_index, sd_bv_ns, sd_bv_ns_rank_1,
                          sd_bv_ns_select_1, config, info, index_name_text,
                          index_bv, text_file, proc_file,
                          true, false, false);
  else
    err = construct_index(fm_index, sd_bv_ns, sd_bv_ns_rank_1,
                          sd_bv_ns_select_1, config, info, index_name,
                          index_bv, text_file, proc_file,
                          false, map_sa, false);
  if (err) {
    fprintf(stderr, "error occured while constructing or loading index\n");
    return err;
  }

  // check if sa_name exists, if i want to map the explicit sa
  if (map_sa) {
    struct stat buffer;
    if (stat(sa_name.c_str(), &buffer) != 0) {
      fprintf(stderr, "sa file can't be opened");
      return EXIT_FAILURE;
    }
  }

  // construct also the reverse index!
  if (test_smems) {
    std::string rev_file = proc_file + "_rev";
    std::string rev_file_index = index_name_smems + "_rev";

    std::ifstream index(rev_file_index);
    if (index.good()) {
      bwd_index.load(index);
      printf("loaded the reverse index\n");
      if (std::remove(proc_file.c_str()) != 0)
        fprintf(stderr, "Error deleting file: %s\n", strerror(errno));
      if (std::remove(rev_file.c_str()) != 0)
        fprintf(stderr, "Error deleting file: %s\n", strerror(errno));
    } else {
      std::ifstream to_reverse(proc_file);
      if (to_reverse) {
        to_reverse.seekg(0, to_reverse.end);
        int length = to_reverse.tellg();
        to_reverse.seekg(0, to_reverse.beg);

        char* buffer = new char [length + 1];
        to_reverse.read(buffer, length);
        to_reverse.close();

        // revert the buffer content here and write it to the new file
        for (int i = 0; i < length / 2; ++i) {
          char tmp = buffer[i];
          buffer[i] = buffer[length - 1 - i];
          buffer[length - 1 - i] = tmp;
        }
        buffer[length] = '\0';
        std::ofstream reversed(rev_file);
        reversed << buffer;

        delete[] buffer;
      }

      printf("constructing the reverse index\n");
      construct(bwd_index, rev_file, 1);

      store_to_file(bwd_index, rev_file_index);
      printf("constructed the reverse index\n");
      if (std::remove(proc_file.c_str()) != 0)
        fprintf(stderr, "Error deleting file: %s\n", strerror(errno));
      if (std::remove(rev_file.c_str()) != 0)
        fprintf(stderr, "Error deleting file: %s\n", strerror(errno));
    }
  }

  // write space requierment to disk:
  std::ofstream json_cst("cst.json");
  std::ofstream json_csa("csa.json");
  std::ofstream json_csa_rev("csarev.json");
  std::ofstream json_text("text.json");
  write_structure<JSON_FORMAT>(fm_index, json_cst);
  write_structure<JSON_FORMAT>(text_index, json_text);
  write_structure<JSON_FORMAT>(fwd_index, json_csa);
  write_structure<JSON_FORMAT>(bwd_index, json_csa_rev);

  // variables for some stats, that are important for map_sa
  unsigned long nof_intervals = 0;
  unsigned long sof_lengths = 0;

  if (thread_query) {
    // read in multifasta-file and preprocess it. store the offset-information for each entry
    std::cout << "preprocess query multi-fasta\n";
    std::string pat_file = pattern_file + "_new";
    std::vector<std::pair<std::string, unsigned long>> query_info;
    unsigned long query_size = preprocess_query(pattern_file, pat_file, query_info);
    unsigned long chunk_size = query_size/thread_num;
    // seperate the single chunks in different files
    std::cout << "chunk the query file. chunk size: " << chunk_size << "\n";
    std::ifstream queries(pat_file);
    std::vector<std::string> pats(thread_num);
    for (unsigned long i = 0; i < (thread_num - 1); i++) {
      char* buf = new char [chunk_size + 1];
      queries.read(buf, chunk_size);
      buf[chunk_size] = '\0';
      pats[i] = std::string(buf);
      delete[] buf;
    }
    unsigned long rest_size = query_size - (chunk_size * (thread_num - 1));
    char* buf = new char [rest_size + 1];
    queries.read(buf, rest_size);
    buf[rest_size] = '\0';
    pats[thread_num - 1] = std::string(buf);
    delete[] buf;
    queries.close();

    // step 1: calculate the starting interval
    std::cout << "calculate the starting intervals\n";
    std::vector<std::pair<cst_t::node_type, unsigned long>> nodes(thread_num);
    for (int i = 0; i < (int)(thread_num - 1); i++) {
      nodes[i] = get_starting_interval(fm_index, pats[i + 1]);
    }
    nodes[thread_num - 1] = std::pair<cst_t::node_type, unsigned long> (fm_index.root(), 0);

    // step 2: call process_seq with all the needed information
    // for thread_num threads!
    std::vector<std::thread> threads(thread_num);
    std::vector<FILE*> fileps(thread_num);
    for (int i = 0; i < (int) thread_num; i++) {
      char prev_char;
      if (i > 0)
        prev_char = (char) tolower((int) pats[i-1].back());
      else
        prev_char = 'z';

      std::string thread_out = out_file + "_" + std::to_string(i);
      FILE* thread_out_filep = fopen(thread_out.c_str(), "w");
      fileps[i] = thread_out_filep;
      thread_args args(fm_index, pats[i], l, thread_out_filep,
                       info, query_info, sd_bv_ns_rank_1, sd_bv_ns_select_1,
                       chunk_size * i, thread_num, i, nodes[i].first,
                       nodes[i].second, prev_char);
      threads[i] = std::thread(thread_process_seq, args);
    }

    // step 3: join all the started threads
    for (int i = 0; i < (int) thread_num; i++) {
      threads[i].join();
    }

    // close all file-pointer
    for (int i = 0; i < (int) thread_num; i++) {
      fclose(fileps[i]);
    }

    // TODO: forgot something?
  } else {
    // iterate over the query-file and search each query
    std::ifstream in(pattern_file);
    std::string line;
    std::string pattern;
    std::string header;
    bool get_seq = false;
    unsigned long seq_counter = 0;

    if (smems || test_smems) {
      printf("processing the reads\n");
    }
    while (std::getline(in, line)) {
      if (line[0] == '>') {
        if (get_seq) {
          seq_counter += 1;
          if (!smems && !test_smems) {
            printf("process sequencenr. %lu\n", seq_counter);
          }
          // process the sequences before processing new header
          if (!normal_text) {
            if (both && !test_smems) {
              process_seq(fm_index, pattern, l, header, info, sd_bv_ns_rank_1, sd_bv_ns_select_1,
                          silent, mums, smems, max_val, n_filter, false, comp, out_printf,
                          sort_paths, enum_mems, clean_mums, sort);
              process_seq(fm_index, pattern, l, header, info, sd_bv_ns_rank_1, sd_bv_ns_select_1,
                          silent, mums, smems, max_val, n_filter, true, comp, out_printf,
                          sort_paths, enum_mems, clean_mums, sort);
            } else if (test_smems) {
              compute_smems(fwd_index, bwd_index, pattern, l,
                            header, info, sd_bv_ns_rank_1,
                            sd_bv_ns_select_1,
                            out_printf, !n_filter);
            } else {
              process_seq(fm_index, pattern, l, header, info, sd_bv_ns_rank_1, sd_bv_ns_select_1,
                          silent, mums, smems, max_val, n_filter, rev, comp, out_printf,
                          sort_paths, enum_mems, clean_mums, sort, map_sa, sa_name);
            }
          } else {
            if (both) {
              process_text(text_index, pattern, l, header, info, silent, mums,
                           smems, max_val, n_filter, false, comp, dna, out_printf, clean_mums);
              process_text(text_index, pattern, l, header, info, silent, mums,
                           smems, max_val, n_filter, true, comp, dna, out_printf, clean_mums);
            } else
              process_text(text_index, pattern, l, header, info, silent, mums,
                           smems, max_val, n_filter, rev, comp, dna, out_printf, clean_mums);
          }
        }
        header = line;
        get_seq = true;
        pattern.clear();
      } else if (get_seq) {
        pattern += line;
      }
    }
    seq_counter += 1;
    if (!smems && !test_smems) {
      printf("process sequencenr. %lu\n", seq_counter);
    }

    // process the last sequence
    if (!normal_text) {
      if (both && !test_smems) {
        process_seq(fm_index, pattern, l, header, info, sd_bv_ns_rank_1, sd_bv_ns_select_1,
                    silent, mums, smems, max_val, n_filter, false, comp, out_printf,
                    sort_paths, enum_mems, clean_mums, sort);
        process_seq(fm_index, pattern, l, header, info, sd_bv_ns_rank_1, sd_bv_ns_select_1,
                    silent, mums, smems, max_val, n_filter, true, comp, out_printf,
                    sort_paths, enum_mems, clean_mums, sort);
      } else if (test_smems) {
        compute_smems(fwd_index, bwd_index, pattern, l,
                      header, info, sd_bv_ns_rank_1,
                      sd_bv_ns_select_1,
                      out_printf, !n_filter);
      } else {
        process_seq(fm_index, pattern, l, header, info, sd_bv_ns_rank_1, sd_bv_ns_select_1,
                    silent, mums, smems, max_val, n_filter, rev, comp, out_printf,
                    sort_paths, enum_mems, clean_mums, sort, map_sa, sa_name);
      }
    } else {
      if (both) {
        process_text(text_index, pattern, l, header, info, silent, mums,
                     smems, max_val, n_filter, false, comp, dna, out_printf, clean_mums);
        process_text(text_index, pattern, l, header, info, silent, mums,
                     smems, max_val, n_filter, true, comp, dna, out_printf, clean_mums);
      } else
        process_text(text_index, pattern, l, header, info, silent, mums,
                     smems, max_val, n_filter, rev, comp, dna, out_printf, clean_mums);
    }

    if (map_sa)
      fprintf(stdout, "nof_intervals: %lu\nsof_lengths: %lu\n", nof_intervals, sof_lengths);
  }
  return 0;
}
