#ifndef UTIL_H
#define UTIL_H

#include <string>
#include <vector>
#include <iostream>

namespace annz_util {
  
  // Parse a line of input data (e.g. magnitudes/errors and optional redshift).
  bool parse_input_data(const std::string data, const int nIn, const int nOut, 
                      bool & truesPresent, std::vector<double> & inputs, 
                      std::vector<double> & errors, std::vector<double> & trues, std::vector<double> & sigma2s);
  
  // Skip past comments/blank lines in the given input file stream.
  void skip_comments(std::ifstream& fstr);
  
  // Locate the next non-comment and non-blank line from the given file stream.
  // The line is returned in the given string.
  bool get_next_dataline(std::ifstream& file, std::string& line);
  
  // Extract the flag and data from the given line.
  int parse_param_line(const std::string& buffer, std::string& flag, std::vector<std::string>& data);
  
  // String-numeric conversions.
  int string_to_int(const std::string s);
  double string_to_double(const std::string s);

  //------------------------------------------------------------
  // Output error message on bad file format.
  inline void file_error(std::string f, std::string l) {
    std::cerr << "ERROR: Bad file format in " << f << ": " << l << std::endl;
  }

}

#endif // UTIL_H
