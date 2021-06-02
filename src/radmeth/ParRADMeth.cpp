/*    Copyright (C) 2013 University of Southern California and
 *                       Egor Dolzhenko
 *                       Andrew D Smith
 *
 *    Authors: Andrew D. Smith and Egor Dolzhenko
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <thread>

// GSL headers
#include <gsl/gsl_cdf.h>

// smithlab headers
#include "OptionParser.hpp"
#include "smithlab_os.hpp"
#include "smithlab_utils.hpp"

// Local headers.
#include "regression.hpp"
#include "combine_pvals.hpp"
#include "merge.hpp"

// MPI Header
#include <mpi.h>

//OpenMP Header
#include <omp.h>

//memcpy Header
#include <string.h>

//uncomment this line to measure performance
//#define MEASURE 1

//uncomment this line to measure performance by phases
//#define DEBUG 1
#define simple_t unsigned int


using std::string;
using std::vector;
using std::istringstream;
using std::ostringstream;
using std::cerr;
using std::endl;
using std::istream;
using std::ostream;

static bool
lt_locus_pval(const PvalLocus &r1, const PvalLocus &r2) {
  return r1.combined_pval < r2.combined_pval;
}

static bool
ls_locus_position(const PvalLocus &r1, const PvalLocus &r2) {
  return r1.pos < r2.pos;
}

void
fdr(vector<PvalLocus> &loci) {

      std::sort(loci.begin(), loci.end(), lt_locus_pval);

      for (size_t ind = 0; ind < loci.size(); ++ind) {
        const double current_score = loci[ind].combined_pval;

        //Assign a new one.
        const double corrected_pval = loci.size()*current_score/(ind + 1);
        loci[ind].corrected_pval = corrected_pval;
      }

      for (vector<PvalLocus>::reverse_iterator
            it = loci.rbegin() + 1; it != loci.rend(); ++it) {

        const PvalLocus &prev_locus = *(it - 1);
        PvalLocus &cur_locus = *(it);

        cur_locus.corrected_pval =
              std::min(prev_locus.corrected_pval, cur_locus.corrected_pval);
      }

      for (vector<PvalLocus>::iterator it = loci.begin();
            it != loci.end(); ++it) {
        PvalLocus &cur_locus = *(it);
        if (cur_locus.corrected_pval > 1.0)
          cur_locus.corrected_pval = 1.0;
      }

      // Restore original order
      std::sort(loci.begin(), loci.end(), ls_locus_position);
}

// Splits a string using white-space characters as delimeters.
static vector<string>
split(string input) {
  istringstream iss(input);
  string token;
  vector<string> tokens;

  while (iss >> token)
    tokens.push_back(token);

  return tokens;
}

// Given the maximum likelihood estimates of the full and reduced models, the
// function outputs the p-value of the log-likelihood ratio. *Note* that it is
// assumed that the reduced model has one fewer factor than the reduced model.
double
loglikratio_test(double null_loglik, double full_loglik) {

  // The log-likelihood ratio statistic.
  const double log_lik_stat = -2*(null_loglik - full_loglik);

  // It is assumed that null model has one fewer factor than the full model.
  // Hence the number of degrees of freedom is 1.
  const size_t degrees_of_freedom = 1;

  // Log-likelihood ratio statistic has a chi-sqare distribution.
  double chisq_p = gsl_cdf_chisq_P(log_lik_stat, degrees_of_freedom);
  const double pval = 1.0 - chisq_p;

  return pval;
}

bool
has_low_coverage(const Regression &reg, size_t test_factor) {

  bool is_covered_in_test_factor_samples = false;
  bool is_covered_in_other_samples = false;

  for (size_t sample = 0; sample < reg.design.sample_names.size(); ++sample) {
    if (reg.design.matrix[sample][test_factor] == 1) {
      if (reg.props.total[sample] != 0)
        is_covered_in_test_factor_samples = true;
    } else {
      if (reg.props.total[sample] != 0)
        is_covered_in_other_samples = true;
    }
  }

  return !is_covered_in_test_factor_samples || !is_covered_in_other_samples;
}

bool
has_extreme_counts(const Regression &reg) {

  bool is_maximally_methylated = true;
  bool is_unmethylated = true;

  for (size_t sample = 0; sample < reg.design.sample_names.size(); ++sample) {
    if (reg.props.total[sample] != reg.props.meth[sample])
      is_maximally_methylated = false;

    if (reg.props.meth[sample] != 0)
      is_unmethylated = false;
  }



  return is_maximally_methylated || is_unmethylated;
}

void read_props_line(std::string this_line,  std::vector<std::string> &chrom,
                      std::vector<size_t> &position, std::vector<std::string> &strand,
                      std::vector<std::string> &context, std::vector<simple_t> &total,
                      std::vector<simple_t> &meth, size_t props_expected) {

  // Skip lines contining only the newline character (e.g. the last line of the
  // proportion table).
  if(this_line.empty())
    return;

  size_t token_end;
  size_t token_start = 0;
  size_t token_colon;

  // Get the row name (which must be specified like this: "chr:position") and
  // parse it.
  token_end = this_line.find("\t", token_start);
  if (token_end == string::npos)
    throw SMITHLABException("Each row in the count table must start with "
                            "a line chromosome:position:strand:context "
                            "and then a tab character but no tab found." );
  string row_name_encoding = this_line.substr(token_start, token_end - token_start);

  // Every row must start an identifier consisiting of genomic loci of the
  // corresponding site. Here we check this identifier has the correct number
  // of colons.
  const size_t num_colon =
            std::count(row_name_encoding.begin(), row_name_encoding.end(), ':');

  if (num_colon != 3)
    throw SMITHLABException("Each row in the count table must start with "
                            "a line chromosome:position:strand:context."
                            "Got \"" + row_name_encoding + "\" instead." );

  // First parse the row identifier.
  token_colon = row_name_encoding.find(":", token_start);
  string chrom_str = row_name_encoding.substr(token_start, token_colon - token_start);
  if (chrom_str.empty())
    throw SMITHLABException("Error parsing " + row_name_encoding +
                            ": chromosome name is missing.");
  chrom.push_back(chrom_str);
  token_start = token_colon + 1;

  token_colon = row_name_encoding.find(":", token_start);
  string position_encoding = row_name_encoding.substr(token_start, token_colon - token_start);
  try{
    position.push_back((size_t)stol(position_encoding));
  }catch(...){
    throw SMITHLABException("The token \"" +position_encoding + "\" "
                            "does not encode a natural number");
  }
  token_start = token_colon + 1;

  token_colon = row_name_encoding.find(":", token_start);
  string strand_str = row_name_encoding.substr(token_start, token_colon - token_start);
  strand.push_back(strand_str);
  token_start = token_colon + 1;

  string context_str = row_name_encoding.substr(token_start);
  context.push_back(context_str);

  // After parsing the row identifier, parse count proportions.
  size_t total_count, meth_count;
  token_start = token_end + 1;
  string props_encoding = this_line.substr(token_start);
  const char *props_chars = props_encoding.c_str();
  token_start = 0;
  token_end = 0;
  size_t props_count = 0;
  try {
    while (props_chars[token_end] != '\0') {
      if (props_chars[token_end] == '\t'){
        string prop_encoding = props_encoding.substr(token_start, token_end - token_start);
        if (props_count % 2){
          //props_count is odd then push_back to meth
          meth_count = stoi(prop_encoding);
          meth.push_back( (simple_t) meth_count);
        } else {
          //props_count is even then push_back to total
          total_count = stoi(prop_encoding);
          total.push_back( (simple_t) total_count);
        }
        props_count++;
        token_start = token_end + 1;
      }

      token_end++;
    }
  } catch(...){
    throw SMITHLABException("This row's proportions are not"
                            "natural numbers:\n" + this_line);
  }


  if (props_count != props_expected)
    throw SMITHLABException("This row does not encode proportions"
                            "correctly:\n" + this_line);
  return;
}


int
main(int argc, char *argv[]) {

  try {

    #ifdef MEASURE
    double local_measure, global_measure;
    std::chrono::time_point<std::chrono::system_clock> measure_start, measure_end;
    std::chrono::duration<double> measure_seconds;
    measure_start = std::chrono::system_clock::now();
    #endif

    #ifdef DEBUG
    double local_time, global_time;
    std::chrono::time_point<std::chrono::system_clock> start, seg_start, end;
    std::chrono::duration<double> elapsed_seconds, elapsed_segment;
    start = std::chrono::system_clock::now();
    seg_start = std::chrono::system_clock::now();
    #endif

    MPI_Init(&argc, &argv);

    int rank, numprocs;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    const string prog_name = strip_path(argv[0]);

    const string main_help_message =
      "Uage: " + prog_name + " [COMMAND] [PARAMETERS]\n\n"
      "Available commands: \n"
      "  regression  Calculates multi-factor differential methylation scores.\n"
      "  adjust      Adjusts the p-value of each site based on the p-value of "
                     "its neighbors.\n"
      "  merge       Combines significantly differentially methylated CpGs into"
                     " DMRs.\n";

    if (argc == 1) {
      cerr << "Analysis of differential methylation in multi-factor bisulfite "
              "sequencing experiments.\n"
           << main_help_message;
      MPI_Finalize();
      return EXIT_SUCCESS;
    }

    // The first command-line argument is the name of the command to run.
    const string command_name = argv[1];

    // Run beta-binoimial regression using the specified table with proportions
    // and design matrix.
    if (command_name == "regression") {
      string outfile;
      string test_factor_name;
      bool VERBOSE = false;
      std::ifstream table_file;

      OptionParser opt_parse(prog_name + "\t" + command_name, "Calculates "
                             "multi-factor differential methylation scores.",
                             "<design-matrix> <data-matrix>");

      opt_parse.add_opt("out", 'o', "output file (default: stdout)",
                        false, outfile);

      opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);

      opt_parse.add_opt("factor", 'f', "a factor to test",
                        true, test_factor_name);

      vector<string> leftover_args;
      opt_parse.parse(argc - 1, (const char**)(argv + 1), leftover_args);

      if(!rank){
        if (argc == 2 || opt_parse.help_requested()) {
          cerr << opt_parse.help_message() << endl
               << opt_parse.about_message() << endl;
          MPI_Abort(MPI_COMM_WORLD, EXIT_SUCCESS);
        }
        if (opt_parse.about_requested()) {
          cerr << opt_parse.about_message() << endl;
          MPI_Abort(MPI_COMM_WORLD, EXIT_SUCCESS);
        }
        if (opt_parse.option_missing()) {
          cerr << opt_parse.option_missing_message() << endl;
          MPI_Abort(MPI_COMM_WORLD, EXIT_SUCCESS);
        }
        if (leftover_args.size() != 2) {
          cerr << opt_parse.help_message() << endl;
          MPI_Abort(MPI_COMM_WORLD, EXIT_SUCCESS);
        }
      }

      const string design_filename(leftover_args.front());
      const string table_filename(leftover_args.back());

      std::ifstream design_file(design_filename.c_str());
      if (!design_file)
        throw SMITHLABException("could not open file: " + design_filename);

      //PROBLEMA AQUÍ, STREAMS DEFINIDAS EN UN SCOPE MUY PEQUEÑO
      //ADEMAS SI ELIJO ESTE APROACH NO NECESITO EL STRING DE SALIDA, USO OTRO DE MPI
      if(!rank){
        table_file = std::ifstream(table_filename.c_str());
        if (!table_file)
          throw SMITHLABException("could not open file: " + table_filename);
      }

      Regression full_regression;            // Initialize the full design
      design_file >> full_regression.design; // matrix from file.

      // Check that the provided test factor name exists and find its index.
      // Here we identify test factors with their indexes to simplify naming.
      vector<string>::const_iterator test_factor_it =
        std::find(full_regression.design.factor_names.begin(),
                  full_regression.design.factor_names.end(), test_factor_name);

      if(!rank){
        if (test_factor_it == full_regression.design.factor_names.end())
          throw SMITHLABException("Error: " + test_factor_name +
                                  " is not a part of the design specification.");
      }

      size_t test_factor = test_factor_it -
                                  full_regression.design.factor_names.begin();

      Regression null_regression;
      null_regression.design = full_regression.design;
      remove_factor(null_regression.design, test_factor);

      // Make sure that the first line of the proportion table file contains
      // names of the samples. Throw an exception if the names or their order
      // in the proportion table does not match those in the full design matrix.
      if(!rank){
        string sample_names_encoding;
        getline(table_file, sample_names_encoding);

        if (full_regression.design.sample_names != split(sample_names_encoding))
          throw SMITHLABException(sample_names_encoding + " does not match factor "
                                  "names or their order in the design matrix. "
                                  "Please verify that the design matrix and the "
                                  "proportion table are correctly formatted.");
      }

      //Hasta aqui todo es igual en todas las versiones, ahora la cosa cambia
      //En la versión secuencial se mete al loop de computacion, leyendo linea a linea
      //En mi version paralela se lee todo el archivo en P0 y lo distribuye antes del loop computacional
      //Vamos a hacer un nuevo aproach mixto

      //PRIMERO leo el archivo paralelamente con MPI
      //Cada proceso tiene su chunk de lineas, se lee todo el archivo pero no se necesita distribuir nada
      //Utilizas la mecanica del overlap del ejemplo para que cada chunk contenga lineas completas
      MPI_File in;
      const char *table_parallel = table_filename.c_str();
      if ( MPI_File_open(MPI_COMM_WORLD, table_parallel, MPI_MODE_RDONLY, MPI_INFO_NULL, &in) ){
        throw SMITHLABException("could not open file: " + table_filename);
      }

      MPI_Offset globalstart;
      long mysize;
      char *chunk;
      MPI_Offset globalend;
      MPI_Offset filesize;
      const int overlap = 50 + ( full_regression.design.sample_names.size() * 30 );

      /* figure out who reads what */
      MPI_File_get_size(in, &filesize);
      filesize--;  /* get rid of text file eof */
      mysize = filesize/numprocs;
      globalstart = rank * mysize;
      globalend   = globalstart + mysize - 1;
      if (rank == numprocs-1) globalend = filesize-1;

      /* add overlap to the end of everyone's chunk except last proc... */
      if (rank != numprocs-1)
          globalend += overlap;

      mysize =  globalend - globalstart + 1;

      /* allocate memory */
      chunk = (char*) malloc( (mysize + 1)*sizeof(char));

      /* everyone reads in their part */
      MPI_File_read_at_all(in, globalstart, chunk, mysize, MPI_CHAR, MPI_STATUS_IGNORE);
      chunk[mysize] = '\0';

      MPI_File_close(&in);

      int locstart=0, locend=mysize-1;
      while(chunk[locstart] != '\n') locstart++;
      locstart++;
      if (rank != numprocs-1) {
          locend-=overlap;
          locend++;
          while(chunk[locend] != '\n') locend++;
      }
      mysize = locend-locstart+1;
      chunk[locend + 1] = '\0';

      //En este punto tienes un array de chars en cada proceso,asegurate de que acabe con '\0' que tu ejemplo no lo hace
      //SEGUNDO Lo transformas en string
      //string lines_string( &(chunk[locstart]) );

      //TERCERO Transformas la string en istringstream
      //istringstream lines_stream(lines_string);

      //La salida del loop no la vas escribiendo en un archivo
      //La escribes en un ostringstream
      ostringstream output_lines;

      #ifdef DEBUG
      end = std::chrono::system_clock::now();
      elapsed_segment = end - seg_start;
      local_time = elapsed_segment.count();
      MPI_Reduce(&local_time, &global_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
      if (!rank){
        std::cout << "Tiempo en la fase 1(Preparacion con MPI I/O): " << global_time << " segundos\n";
      }
      MPI_Barrier(MPI_COMM_WORLD);
      seg_start = std::chrono::system_clock::now();
      #endif

      std::vector<std::vector<std::string>> chrom_v;
      std::vector<std::vector<size_t>> position_v;
      std::vector<std::vector<std::string>> strand_v;
      std::vector<std::vector<std::string>> context_v;
      std::vector<std::vector<simple_t>> array_total_v;
      std::vector<std::vector<simple_t>> array_meth_v;

      //one vector per thread
      #pragma omp parallel shared(chrom_v, position_v, strand_v, context_v, array_total_v, array_meth_v)
      {

        if ( !( omp_get_thread_num() ) ){
          chrom_v.resize(omp_get_num_threads());
          position_v.resize(omp_get_num_threads());
          strand_v.resize(omp_get_num_threads());
          context_v.resize(omp_get_num_threads());
          array_total_v.resize(omp_get_num_threads());
          array_meth_v.resize(omp_get_num_threads());
        }

      }

      #pragma omp parallel shared(chrom_v, position_v, strand_v, context_v, array_total_v, array_meth_v, chunk, locstart, locend, mysize) firstprivate(full_regression)
      {
        int num_threads = omp_get_num_threads();
        int thread_num = omp_get_thread_num();
        long th_size = mysize / num_threads;
        long th_start = locstart + ( th_size * thread_num );
        long th_end = locstart + ( th_size * ( thread_num + 1 ) );
        if ( thread_num == (num_threads - 1) ){
          th_end = locend;
        } else {
          th_end += overlap;
        }
        long th_chunk_size = th_end + 1 - th_start;
        char* th_chunk =  (char*) malloc( (th_chunk_size + 1) * sizeof(char));
        memcpy( th_chunk, &( chunk[th_start] ), th_chunk_size );
        th_start = 0;
        th_end = th_chunk_size - 1;
        if(thread_num){
          while(th_chunk[th_start] != '\n') th_start++;
          th_start++;
        }
        if( thread_num != (num_threads - 1) ){
          th_end -= overlap;
          while(th_chunk[th_end] != '\n') th_end++;
        }
        th_chunk[th_end + 1] = '\0';


        string lines_string( &(th_chunk[th_start]) );

        std::vector<std::string> th_chrom;
        std::vector<size_t> th_position;
        std::vector<std::string> th_strand;
        std::vector<std::string> th_context;
        std::vector<simple_t> th_array_total;
        std::vector<simple_t> th_array_meth;

        size_t str_start = 0;
        size_t str_end;

        while( 1 ){
          string this_line;
          if ((str_end = lines_string.find("\n", str_start)) == string::npos) {
            if (!(this_line = lines_string.substr(str_start)).empty()) {
              read_props_line(this_line, th_chrom, th_position, th_strand, th_context, th_array_total, th_array_meth, full_regression.design.sample_names.size() * 2);
            }
            break;
          }

          this_line = lines_string.substr(str_start, str_end - str_start);
          read_props_line(this_line, th_chrom, th_position, th_strand, th_context, th_array_total, th_array_meth, full_regression.design.sample_names.size() * 2);
          str_start = str_end + 1;
        }

        //places its partial result on its main vector slot
        chrom_v[omp_get_thread_num()] = th_chrom;
        position_v[omp_get_thread_num()] = th_position;
        strand_v[omp_get_thread_num()] = th_strand;
        context_v[omp_get_thread_num()] = th_context;
        array_total_v[omp_get_thread_num()] = th_array_total;
        array_meth_v[omp_get_thread_num()] = th_array_meth;

      }

      //concat all parts
      std::vector<std::string> chrom;
      for(long unsigned int i=0; i<chrom_v.size(); i++){
        chrom.insert(chrom.end(), chrom_v[i].begin(), chrom_v[i].end());
        chrom_v[i].clear();
      }

      std::vector<size_t> position;
      for(long unsigned int i=0; i<position_v.size(); i++){
        position.insert(position.end(), position_v[i].begin(), position_v[i].end());
        position_v[i].clear();
      }

      std::vector<std::string> strand;
      for(long unsigned int i=0; i<strand_v.size(); i++){
        strand.insert(strand.end(), strand_v[i].begin(), strand_v[i].end());
        strand_v[i].clear();
      }

      std::vector<std::string> context;
      for(long unsigned int i=0; i<context_v.size(); i++){
        context.insert(context.end(), context_v[i].begin(), context_v[i].end());
        context_v[i].clear();
      }

      std::vector<simple_t> array_total;
      for(long unsigned int i=0; i<array_total_v.size(); i++){
        array_total.insert(array_total.end(), array_total_v[i].begin(), array_total_v[i].end());
        array_total_v[i].clear();
      }

      std::vector<simple_t> array_meth;
      for(long unsigned int i=0; i<array_meth_v.size(); i++){
        array_meth.insert(array_meth.end(), array_meth_v[i].begin(), array_meth_v[i].end());
        array_meth_v[i].clear();
      }

      //get number of lines on each proc
      size_t lines = chrom.size();



      #ifdef DEBUG
      end = std::chrono::system_clock::now();
      elapsed_segment = end - seg_start;
      local_time = elapsed_segment.count();
      MPI_Reduce(&local_time, &global_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
      if (!rank){
        std::cout << "Tiempo en la fase 2(Lectura): " << global_time << " segundos\n";
      }
      MPI_Barrier(MPI_COMM_WORLD);
      seg_start = std::chrono::system_clock::now();
      #endif


      double *fit_array = (double*) malloc(sizeof(double) * lines);

      if(fit_array == NULL)
        throw SMITHLABException("Too much memory required(fit_array)");

      // Performing the log-likelihood ratio test on proportions from each row
      // of the proportion table.
      #pragma omp parallel for firstprivate(full_regression, null_regression) shared(array_total, array_meth, fit_array, test_factor) schedule(dynamic)
      for (size_t i = 0; i < (lines); ++i) {

        full_regression.props.total.clear();
        full_regression.props.meth.clear();
        for(int k = 0; k < (int) full_regression.design.sample_names.size(); ++k){
          full_regression.props.total.push_back((size_t) array_total[(i * full_regression.design.sample_names.size()) + k]);
          full_regression.props.meth.push_back((size_t) array_meth[(i * full_regression.design.sample_names.size()) + k]);
        }



        // Do not perform the test if there's no coverage in either all case or
        // all control samples. Also do not test if the site is completely
        // methylated or completely unmethylated across all samples.
        if (has_low_coverage(full_regression, test_factor)) {
          fit_array[i] = -1;
        }
        else if (has_extreme_counts(full_regression)) {
          fit_array[i] = -1;
        }
        else {
          fit(full_regression);
          null_regression.props = full_regression.props;
          fit(null_regression);
          const double pval = loglikratio_test(null_regression.max_loglik,
                                         full_regression.max_loglik);

          // If error occured in the fitting algorithm (i.e. p-val is nan or
          // -nan).
          fit_array[i] = ( (pval != pval) ? -1 : pval);
        }
      }

      #ifdef DEBUG
      end = std::chrono::system_clock::now();
      elapsed_segment = end - seg_start;
      local_time = elapsed_segment.count();
      MPI_Reduce(&local_time, &global_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
      if (!rank){
        std::cout << "Tiempo en fase 3(loop computacional): " << global_time << " segundos\n";
      }
      MPI_Barrier(MPI_COMM_WORLD);
      seg_start = std::chrono::system_clock::now();
      #endif

      //escritura
      for (size_t i = 0; i < (lines); ++i) {

        size_t coverage_factor = 0, coverage_rest = 0,
               meth_factor = 0, meth_rest = 0;

        for(size_t s = 0; s < full_regression.design.sample_names.size(); ++s) {
          if(full_regression.design.matrix[s][test_factor] != 0) {
            coverage_factor += array_total[s + (i * full_regression.design.sample_names.size())];
            meth_factor += array_meth[s + (i * full_regression.design.sample_names.size())];
          } else {
            coverage_rest += array_total[s + (i * full_regression.design.sample_names.size())];
            meth_rest += array_meth[s + (i * full_regression.design.sample_names.size())];
          }
        }

        output_lines << chrom[i] << "\t"
            << position[i] << "\t"
            << strand[i] << "\t"
            << context[i] << "\t";

        output_lines << fit_array[i];

        output_lines << "\t" << coverage_factor << "\t" << meth_factor
            << "\t" << coverage_rest << "\t" << meth_rest << endl;

      }


      #ifdef DEBUG
      end = std::chrono::system_clock::now();
      elapsed_segment = end - seg_start;
      local_time = elapsed_segment.count();
      MPI_Reduce(&local_time, &global_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
      if (!rank){
        std::cout << "Tiempo en fase 4(Escritura a string): " << global_time << " segundos\n";
      }
      MPI_Barrier(MPI_COMM_WORLD);
      seg_start = std::chrono::system_clock::now();
      #endif

      //Asi que despues del loop cojes en cada proceso la string que te dio y calculas la lenght
      string output_string = output_lines.str();
      size_t output_length = output_string.length();

      //Con estas length haces Allgather y todos los procesos saben donde hacer la escritura paralela en el archivo de salida
      size_t *procs_length = (size_t *) malloc( numprocs * sizeof( size_t ) );
      MPI_Allgather( &output_length, 1, MPI_AINT, procs_length, 1, MPI_AINT, MPI_COMM_WORLD);

      #ifdef DEBUG
      end = std::chrono::system_clock::now();
      elapsed_segment = end - seg_start;
      if (!rank){
        std::cout << "Tiempo en fase 5(MPI_Allgather): " << elapsed_segment.count() << " segundos\n";
      }
      MPI_Barrier(MPI_COMM_WORLD);
      seg_start = std::chrono::system_clock::now();
      #endif

      //Calculas donde empieza a escribir cada proceso
      size_t write_offset = 0;
      for(int i = 0; i < rank; i++)
        write_offset+=procs_length[i];

      //abres el archivo de salida paralelamente
      MPI_File out;
      if ( MPI_File_open(MPI_COMM_WORLD, outfile.c_str(), MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &out) ){
        throw SMITHLABException("could not open file: " + outfile);
      }

      //cada proceso tiene toda la informacion para escribir su parte de la salida donde le toca
      MPI_File_write_at_all(out, (MPI_Offset)write_offset, output_string.c_str(), output_length, MPI_CHAR, MPI_STATUS_IGNORE);

      //Despues de escribir cierras el archivo
      MPI_File_close(&out);

      #ifdef DEBUG
      end = std::chrono::system_clock::now();
      elapsed_segment = end - seg_start;
      local_time = elapsed_segment.count();
      MPI_Reduce(&local_time, &global_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
      if (!rank){
        std::cout << "Tiempo en fase 6(Escritura con MPI I/O): " << global_time << " segundos\n";
      }
      seg_start = std::chrono::system_clock::now();
      elapsed_seconds = seg_start - start;
      if (!rank){
        std::cout << "Tiempo total: " << elapsed_seconds.count() << " segundos\n";
      }
      #endif

      #ifdef MEASURE
      measure_end = std::chrono::system_clock::now();
      measure_seconds = measure_end - measure_start;
      local_measure = measure_seconds.count();
      MPI_Reduce(&local_measure, &global_measure, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
      if (!rank){
        std::cout << "Tiempo total: " << global_measure << " segundos\n";
      }
      #endif

    // Combine p-values using the Z test.
    } else if (command_name == "adjust" && !rank) {
      string outfile;
      string bin_spec = "1:200:1";

      /****************** GET COMMAND LINE ARGUMENTS ***************************/
      OptionParser opt_parse(prog_name + "\t" + command_name, "computes "
                             "adjusted p-values using autocorrelation",
                             "<regression-output>");
      opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)",
            false , outfile);
      opt_parse.add_opt("bins", 'b', "corrlation bin specification",
            false , bin_spec);
      vector<string> leftover_args;
      opt_parse.parse(argc - 1, (const char**)(argv + 1), leftover_args);
      if (argc == 2 || opt_parse.help_requested()) {
        cerr << opt_parse.help_message() << endl;
        return EXIT_SUCCESS;
      }
      if (opt_parse.about_requested()) {
        cerr << opt_parse.about_message() << endl;
        return EXIT_SUCCESS;
      }
      if (opt_parse.option_missing()) {
        cerr << opt_parse.option_missing_message() << endl;
        return EXIT_SUCCESS;
      }
      if (leftover_args.size() != 1) {
        cerr << opt_parse.help_message() << endl;
        return EXIT_SUCCESS;
      }
      const string bed_filename = leftover_args.front();
      /*************************************************************************/

      BinForDistance bin_for_dist(bin_spec);

      std::ifstream bed_file(bed_filename.c_str());

      if (!bed_file)
        throw "could not open file: " + bed_filename;

      cerr << "Loading input file." << endl;

      // Read in all p-value loci. The loci that are not correspond to valid
      // p-values (i.e. values in [0, 1]) are skipped.
      vector<PvalLocus> pvals;
      std::string input_line, prev_chrom;

      size_t chrom_offset = 0;

      while(getline(bed_file, input_line)) {

        try {
          std::istringstream iss(input_line);
          iss.exceptions(std::ios::failbit);
          std::string chrom, sign, name;
          size_t position;
          double pval;
          iss >> chrom >> position >> sign >> name >> pval;

          // Skip loci that do not correspond to valid p-values.
          if (0 <= pval && pval <= 1) {
            // locus is on new chrom.
            if (!prev_chrom.empty() && prev_chrom != chrom)
              chrom_offset += pvals.back().pos;

            PvalLocus plocus;
            plocus.raw_pval = pval;
            plocus.pos = chrom_offset +
            bin_for_dist.max_dist() + 1 + position;

            pvals.push_back(plocus);
            prev_chrom = chrom;
          }
        } catch (std::exception const & err) {
          std::cerr << err.what() << std::endl << "Couldn't parse the line \""
                    << input_line << "\"." << std::endl;
          std::terminate();
        }
      }

      cerr << "[done]" << endl;

      cerr << "Combining p-values." << endl;
      combine_pvals(pvals, bin_for_dist);
      cerr << "[done]" << endl;

      cerr << "Running multiple test adjustment." << endl;
      fdr(pvals);
      cerr << "[done]" << endl;

      std::ofstream of;
      if (!outfile.empty()) of.open(outfile.c_str());
        std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

      std::ifstream original_bed_file(bed_filename.c_str());

      update_pval_loci(original_bed_file, pvals, out);

      //TODO: Check that the regions do not overlap & sorted

    } else if (command_name == "merge" && !rank) {
      /* FILES */
      string outfile;
      string bin_spec = "1:200:25";
      double cutoff = 0.01;

      /****************** GET COMMAND LINE ARGUMENTS ***************************/
      OptionParser opt_parse("dmrs", "a program to merge significantly "
                             "differentially methylated CpGs into DMRs",
                             "<bed-file-in-radmeth-format>");
      opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)",
            false , outfile);
      opt_parse.add_opt("cutoff", 'p', "P-value cutoff (default: 0.01)",
            false , cutoff);
      vector<string> leftover_args;
      opt_parse.parse(argc - 1, (const char**)(argv + 1), leftover_args);
      if (argc == 1 || opt_parse.help_requested()) {
        cerr << opt_parse.help_message() << endl;
        return EXIT_SUCCESS;
      }
      if (opt_parse.about_requested()) {
        cerr << opt_parse.about_message() << endl;
        return EXIT_SUCCESS;
      }
      if (opt_parse.option_missing()) {
        cerr << opt_parse.option_missing_message() << endl;
        return EXIT_SUCCESS;
      }
      if (leftover_args.size() != 1) {
        cerr << opt_parse.help_message() << endl;
        return EXIT_SUCCESS;
      }
      const string bed_filename = leftover_args.front();
      /************************************************************************/

      std::ofstream of;
      if (!outfile.empty()) of.open(outfile.c_str());
      std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

      std::ifstream bed_file(bed_filename.c_str());

      if (!bed_file)
        throw "could not open file: " + bed_filename;

      merge(bed_file, out, cutoff);
    } else {
      cerr << "ERROR: \"" << command_name << "\" is not a valid command.\n"
           << main_help_message;
      return EXIT_SUCCESS;
    }
    /****************** END COMMAND LINE OPTIONS *****************/

  }
  catch (const SMITHLABException &e) {
    cerr << e.what() << endl;
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    return EXIT_FAILURE;
  }
  MPI_Finalize();
  return EXIT_SUCCESS;
}
