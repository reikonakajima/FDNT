#include <sstream>
#include <fstream>
#include <iostream>
#include "util.h"

bool annz_util::parse_input_data(const std::string data, const int nIn, const int nOut, 
		     bool& truesPresent, std::vector<double>& inputs, 
		     std::vector<double>& errors, std::vector<double>& trues) {

  // Clear data vectors.
  inputs.clear();
  errors.clear();
  trues.clear();

  // Attach a stringstream to the data string.
  std::stringstream ss(data, std::stringstream::in);

  // Temporary data storage.
  double temp;
  
  // Read inputs (expect n_inputs values).
  for (int ii = 0; ii < nIn; ++ii) {
    if (ss >> temp) {
      inputs.push_back(temp);
    } else {
      std::cerr << "ERROR: Bad input data (insufficient fields)" << std::endl;
      return false;
    }
  }
  
  // Read errors (again, expect n_inputs values).
  for (int ii = 0; ii < nIn; ++ii) {
    if (ss >> temp) {
      errors.push_back(temp);
    } else {
      std::cerr << "ERROR: Bad input data (insufficient fields)" << std::endl;
      return false;
    }
  }
  
  // Is there any more?
  if (ss >> temp) {
    // Assume it's evaluation data.
    trues.push_back(temp);
    
    // Look for the rest of the outputs.
    for (int ii = 1; ii < nOut; ++ii) {
      if (ss >> temp) {
	trues.push_back(temp);
      } else {
	std::cerr << "ERROR: Bad input data (insufficient fields)" << std::endl;
	return false;
      }
    }

    // Ensure we are now at the end of the buffer.
    ss >> std::ws;
    if (!ss.eof()) {
      std::cerr << "ERROR: Bad input data (too many fields)" << std::endl;
      return false;
    }

    truesPresent = true;

  } else {
    truesPresent = false;
  }
  
  return true;
  
}

//-------------------------------------------------
// Skip past comments/ blank lines in a file.
void annz_util::skip_comments(std::ifstream& fstr)
{
  int place = fstr.tellg();
  char a = fstr.get();
  std::string word1;
  fstr >> word1;
  std::cout << "got character " << a << " and word " << word1 << std::endl;

  if(word1.empty() || word1[0]=='#' || a=='#')
  {
    std::string word2;
    getline(fstr, word2);
    if(fstr.eof()) return;
    std::cout << "skipping line " << a << " // " << word1 << " // " << word2 << std::endl;
    skip_comments(fstr);
    return;
  }
  if(a=='\n') place++;
  fstr.seekg(place);
}

// Locate the next non-comment and non-blank line from the given file stream.
bool annz_util::get_next_dataline(std::ifstream& file, std::string& line) {
  
  skip_comments(file);

  // Return fail if end of file.
  if (file.eof()) return false;
  
  // Put the next line in 'line'
  getline(file, line);
  return true;
}

// Extract the flag and data from the given line.
int annz_util::parse_param_line(const std::string& buffer, std::string& flag, std::vector<std::string>& data) {
  
  std::stringstream ss(buffer, std::stringstream::in);

  // First item is flag.
  ss >> flag;
  
  ss >> std::ws;
  std::string temp;
  data.clear();
  while (!ss.eof()) {
    ss >> temp;
    data.push_back(temp);
    ss >> std::ws;
  }
}

// Read int from string.
int annz_util::string_to_int(const std::string s) {
  int i = 0;
  std::stringstream ss(s, std::stringstream::in);
  ss >> i;
  return i;
}

// Read double from string.
double annz_util::string_to_double(const std::string s) {
  double d = 0.0;
  std::stringstream ss(s, std::stringstream::in);
  ss >> d;
  return d;
}
