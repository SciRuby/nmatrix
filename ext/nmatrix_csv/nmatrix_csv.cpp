/////////////////////////////////////////////////////////////////////
// = NMatrix
//
// A linear algebra library for scientific computation in Ruby.
// NMatrix is part of SciRuby.
//
// NMatrix was originally inspired by and derived from NArray, by
// Masahiro Tanaka: http://narray.rubyforge.org
//
// == Copyright Information
//
// SciRuby is Copyright (c) 2010 - 2014, Ruby Science Foundation
// NMatrix is Copyright (c) 2012 - 2014, John Woods and the Ruby Science Foundation
//
// Please see LICENSE.txt for additional copyright notices.
//
// == Contributing
//
// By contributing source code to SciRuby, you agree to be bound by
// our Contributor Agreement:
//
// * https://github.com/SciRuby/sciruby/wiki/Contributor-Agreement
//
// == nmatrix_csv.cpp
//
// Main file for nmatrix_csv extension
//
#include <ruby.h>
#include "nmatrix.h"
#include <ruby/io.h>
#include <cstdio>
#include <vector>
#include <algorithm> // std::max

#include "csv.h"

using namespace std;

static size_t dim2 = 0;

const size_t BUFFERSIZE = 1024;
char* buffer = NULL;

vector< vector<string> > unfilled_data; // raw data read from csv_file

// push field read into vector
static void field_callback(void* s, size_t len, void* data) {
  unfilled_data.back().push_back(string((char*)s, len));
}

// push a new vector into unfilled_data when a new row starts
static void row_callback(int c, void* data) {
  size_t last_sz = unfilled_data.back().size();
  dim2 = std::max(last_sz, dim2);
  unfilled_data.push_back(vector<std::string>());
}

// only called inside cpp
// free data 
void free_data(FILE* p_file, csv_parser* p_parser)
{
  fclose(p_file);
  csv_free(p_parser);
  delete [] buffer;
  unfilled_data.clear();
}

// only called inside cpp
// actually call csv_parse, csv_fini to do the parse job
// return 0 on success, -1 otherwise
static int actually_parse(FILE* p_file, csv_parser* p_parser) {

  long bytes_read, parser_processed;
  int result;

  while(!feof(p_file))
  {
    bytes_read = fread((void*)buffer, sizeof(char), BUFFERSIZE, p_file);
    if(bytes_read != BUFFERSIZE && ferror(p_file)) {
      rb_warning("fread error\n");
      return -1;
    }
    parser_processed = csv_parse(p_parser, buffer, bytes_read,
              field_callback,
              row_callback,
              NULL);
    if(parser_processed != bytes_read) { // error handling
      rb_warning("error encountered executing csv_parse: %s\n",
                csv_strerror(csv_error(p_parser)));
      return -1;
    }
  }

  result = csv_fini(p_parser, field_callback, row_callback, NULL);

  if(result != 0) {
    rb_warning("error encountered executing csv_fini: %s\n",
              csv_strerror(csv_error(p_parser)));
    return -1;
  }

  if(dim2 == 0 || unfilled_data.size() - 1 == 0) { // no element read
    return -1;
  }

  return 0;
}

// only called inside this cpp
static VALUE build_return_array()
{
  VALUE ret_arr = rb_ary_new();
  VALUE ret_dimension = rb_ary_new2(2);
  rb_ary_push(ret_dimension, INT2FIX(unfilled_data.size() - 1));
  rb_ary_push(ret_dimension, INT2FIX(dim2));
  VALUE ret_data = rb_ary_new();
  for(size_t i = 0; i < unfilled_data.size() - 1; i++)
  {
    for(size_t j = 0; j < unfilled_data[i].size(); j++) {
      rb_ary_push(ret_data, 
        rb_str_new(unfilled_data[i][j].c_str(), unfilled_data[i][j].size()));
    }
    for(size_t j = 0; j < dim2 - unfilled_data[i].size(); j++) {
      rb_ary_push(ret_data, Qnil);
    }
  }
  // TODO: return actual array
  rb_ary_push(ret_arr, ret_dimension);
  rb_ary_push(ret_arr, ret_data);
  return ret_arr;
}

/*
 * call-seq Csv.parse_file(filename, delim) -> Array
 *
 * Parse a csv file into a nmatrix using strict mode in libcsv,
 * any malformed structure will be viewed as an error.
 * If parse successfully, an array of [shape, data] will be returned,
 * shape is the dimension of the matrix, data will be the matrix data
 * as 1D - array. Any unaligned row will be filled with nil to pad with
 * longest row.
 *
 * *Arguments* :
 * - +filename+ -> string file path for given csv file
 * - +delim+ -> delimiter of csv file, comma by default
 *
 */
static VALUE parse_file(VALUE self, VALUE str, VALUE delim)
{
  int result;
  csv_parser parser;
  dim2 = 0;

  if((result = csv_init(&parser, CSV_STRICT | CSV_STRICT_FINI)) != 0) {
    rb_raise(rb_eRuntimeError, 
      "error encountered when initializing csv_parser\n");
    return Qnil;
  }

  csv_set_delim(&parser, RSTRING_PTR(delim)[0]);

  const char* filename = RSTRING_PTR(str);
  FILE* p_file;
  p_file = fopen(filename, "rb");
  if(!p_file) { 
    csv_free(&parser);
    rb_raise(rb_eRuntimeError, "open file error\n");
    return Qnil;
  }

  buffer = new char[BUFFERSIZE];

  unfilled_data.push_back(vector<std::string>());

  result = actually_parse(p_file, &parser);
  if(result == -1) {
    free_data(p_file, &parser);
    return Qnil;
  }

  VALUE ret_arr = build_return_array();
  // Cleaning
  free_data(p_file, &parser);

  return ret_arr;
}

extern "C"  void Init_nmatrix_csv() {
  static VALUE mCsv;
  mCsv = rb_define_module("Csv");
  rb_define_module_function(mCsv, "parse_file", (METHOD)parse_file, 2);
}