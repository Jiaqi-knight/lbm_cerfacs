#include "ParamFile.hpp"
#include "Output.hpp"
#include "GlobalVars.hpp"

using namespace std;

void readParamFile(ifstream& in_run_input_file)
{ 
  string text_line;
  string keyword, value;

  // Now read in parameters
  while (getline (in_run_input_file, text_line)) {
    if (readLine(text_line, keyword, value)){
      analyzeLine(keyword, value.c_str());
    }
  }
  printInfo();
}

void analyzeLine(string & keyword, const char * value){

  if (!keyword.compare("lbm model")){
    string tmp_value = string(value);
    if (!tmp_value.compare("d2q9") and !tmp_value.compare("d2q5")){
      cout << "LBM model not in (D2Q4, D2Q5, D2Q9)" << endl;
      exit (EXIT_FAILURE);
    }
    StringToUpperCase(tmp_value);
    LBMmodel_ = string(tmp_value.c_str());
    if (!LBMmodel_.compare("D2Q4")){
      Q_ = 5;   
    }
    else if (!LBMmodel_.compare("D2Q5")){
      Q_ = 5;   
    }
    else if (!LBMmodel_.compare("D2Q9")){
      Q_ = 9;
    }
  }
  if (!keyword.compare("regularisation")){
    string tmp_value = string(value);
    if (!tmp_value.compare("false") and !tmp_value.compare("true")){
      cout << "Regularisation not in (True, False)" << endl;
      exit (EXIT_FAILURE);
    }
    reg_ = (tmp_value == "false") ? false : true;
  }
  else if (!keyword.compare("reference temperature (k)")) sscanf(value,"%lf",&Tref_);
  else if (!keyword.compare("reference pressure (pa)"))   sscanf(value,"%lf",&Pref_);
  else if (!keyword.compare("kinematic viscosity (m^2/s)")) sscanf(value,"%lf",&Nu_);
  else if (!keyword.compare("number of iterations"))        sscanf(value,"%d",&nite_);
  else if (!keyword.compare("simulation time (s)"))        sscanf(value,"%lf",&simTime_);
  else if (!keyword.compare("periodicity in output file"))  sscanf(value,"%d",&freqOut_);
  else if (!keyword.compare("periodicity in integrated quantities"))  sscanf(value,"%d",&freqIntQuant_);
  else if (!keyword.compare("configuration file")) config_file = string(value);
  else if (!keyword.compare("rgas"))        sscanf(value,"%lf",&Rgas_);
  else if (!keyword.compare("gamma"))       sscanf(value,"%lf",&Gamma_);
  else if (!keyword.compare("aero distribution development order"))   sscanf(value,"%d",&aeroOrder_);
  else if (!keyword.compare("cryogenic ratio"))       sscanf(value,"%lf",&cryoRatio_);
  else if (!keyword.compare("number of threads"))       sscanf(value,"%d",&nThreads_);
}


bool readLine(string & str, string & keyword, string & value) {

  const string delimiters(" \t\n\r");
  string::size_type pos, last_pos;

  if(str.empty()) return false; // empty str 

  replaceTabsAndReturns(str);
  pos = str.find_first_of("#");
  if (pos == 0)  return false;  // a full comment line
  if (pos != string::npos) str.erase(pos); // remove comment at end if necessary

  pos = str.find("//");
  if (pos != string::npos) str.erase(pos); // remove comment at end if necessary
    
  // find the : sign and split string
  string name_part, value_part;
  pos = str.find(":");
  if (pos == string::npos) return false;

  name_part = str.substr(0, pos);
  value_part = str.substr(pos+1, string::npos);

  // remove left white space or tabs
  last_pos = name_part.find_first_not_of(delimiters, 0);
  pos = name_part.find_first_of(delimiters, last_pos);
  keyword = name_part.substr(last_pos, name_part.length() - last_pos);
  keyword = trim(keyword);
  StringToLowerCase(keyword);

  // remove left white space or tabs
  last_pos = value_part.find_first_not_of(delimiters, 0);
  pos = value_part.find_first_of(delimiters, last_pos);
  value = value_part.substr(last_pos, value_part.length() - last_pos);
  value = trim(value);
  if (keyword.compare("configuration file")) StringToLowerCase(value);

  return true;
}

string trim(string & str)
{
  if(str.empty())
      return str;
  size_t firstScan = str.find_first_not_of(' ');
  size_t first     = firstScan == string::npos ? str.length() : firstScan;
  size_t last      = str.find_last_not_of(' ');
  return str.substr(first, last-first+1);
}

string replaceTabsAndReturns(string & str)
{
  size_t pos;
  if(str.empty())
      return str;

  while (1){
    pos = str.find_first_of("\t");
    if (pos == string::npos) break;
    str.replace(pos,1," ");      
  }

  while (1){
    pos = str.find_first_of("\r");
    if (pos == string::npos) break;
    str.replace(pos,1," ");      
  }

  return str;
}