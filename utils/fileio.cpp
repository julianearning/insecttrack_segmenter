#include "fileio.hpp"
#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <ios>
#include <cmath>
#include <cstring>
#include <unordered_map>


const char * noise_color = "#ff0000";

const char * colors [] = {  "#696969", "#7f0000", "#006400", "#808000", "#483d8b", "#008b8b", "#4682b4", "#000080", "#d2691e", "#7f007f", "#8fbc8f", "#ff4500", "#ffa500", "#00ff00", "#00ff7f", "#e9967a", "#dc143c", "#00ffff", "#0000ff", "#a020f0", "#adff2f", "#d8bfd8", "#1e90ff", "#db7093", "#f0e68c", "#ffff54",  "#90ee90", "#ff1493", "#7b68ee", "#ee82ee" };
 

void hex_to_rgb(const char * hex, int * r, int * g, int * b) {
  /* char* number = new char[2];
  number[0] = *(hex+1);
  number[1] = *(hex+2);
  *r = (int)std::stoul(number, nullptr, 16);
  number[0] = *(hex+3);
  number[1] = *(hex+4);
  *g = std::stoul(number, nullptr, 16);
  number[0] = *(hex+5);
  number[1] = *(hex+6);
  *b = std::stoul(number, nullptr, 16); */
  // Check if the input hex string is valid
  if(hex[0] == '#' && strlen(hex) == 7) {
      // Convert each pair of hex characters to integers
      *r = std::stoi(std::string(hex + 1, 2), nullptr, 16);
      *g = std::stoi(std::string(hex + 3, 2), nullptr, 16);
      *b = std::stoi(std::string(hex + 5, 2), nullptr, 16);
  } else {
      // If the input is invalid, set RGB values to 0
      *r = *g = *b = 0;
  }
}



//-------------------------------------------------------------------
// Finds min and max Point of *cloud*
// the Points are returned in *min* and *max*.
//-------------------------------------------------------------------
void find_min_max_point(const PointCloud &cloud, Point &min, Point &max) {
  min = max = cloud[0];
  for (std::vector<Point>::const_iterator p = cloud.begin(); p != cloud.end(); ++p) {
    for(size_t i = 0; i<cloud[0].data.size(); i++) {
      if (min.data.at(i) > p->data.at(i)) {
        min.data.at(i) = p->data.at(i);
      } else if (max.data.at(i) < p->data.at(i)) {
        max.data.at(i) = p->data.at(i);
      }
    }
    if (min.t > p->t) {
      min.t = p->t;
    } else if (max.t < p->t) {
      max.t = p->t;
    }
  }
}


// Split string *input* into substrings by *delimiter*. The result is
// returned in *result*
void split(const std::string &input, std::vector<std::string> &result,
           const char delimiter) {
  std::stringstream ss(input);
  std::string element;

  while (std::getline(ss, element, delimiter)) {
    result.push_back(element);
  }
}

bool isNullOrWhitespace(const std::string& str) {
  return str.empty()
      || std::all_of(str.begin(), str.end(), [](char c) {
      return std::isspace(static_cast<unsigned char>(c));
  });
}

int FileParser::parse_file(std::string filename, PointCloud & pc, double scale_factor_t, bool is4d) {
    const char delimiter = ',';
    std::ifstream infile(filename);
    std::string line;
    std::vector<std::string> items;
    if (infile.fail()) throw std::exception();
    while (!infile.eof()) {
        std::getline(infile, line, '\n');
        if(!isNullOrWhitespace(line)) {
          Point point;
          items.clear();
          split(line, items, delimiter);
          point.data.push_back((double)std::stoi(items[0].c_str()));
          point.data.push_back((double)std::stoi(items[1].c_str()));
          if(is4d) point.data.push_back((double)std::stoi(items[2].c_str()));
          //point.t = std::stod(items[2].c_str()) / scale_factor_t;
          // since when t is scaled you might get something like xxxxx.000000000000001. I want to compare for equality so its a problem. So I round... which means that the scaling might be off a little bit so I need to store the original data
          if(is4d) {
            original_time_information.push_back(std::stod(items[3].c_str()));
            point.t = std::round((std::stod(items[3].c_str()) * scale_factor_t) * 10.0 ) / 10.0;
          } else {
            original_time_information.push_back(std::stod(items[2].c_str()));
            point.t = std::round((std::stod(items[2].c_str()) * scale_factor_t) * 10.0 ) / 10.0;
          } 
          pc.push_back(point);
        }
    }
    return 0;
}



int FileParser::write_csv(std::string filename, PointCloud & pc, std::vector<int> & labels, bool append, std::vector<size_t>* noise_points) {
  int r=0,g=0,b=0;
  // set ostream f to either file or stdout
  std::streambuf* buf;
  std::ofstream fhelp;
  if(!filename.empty()) {
    if(append==true) {
      fhelp.open(filename.c_str(), std::ios::out | std::ios::app);
    } else {
      fhelp.open(filename.c_str());
    }
    if (fhelp.fail())
      return 1;
    buf = fhelp.rdbuf();
  } else {
    buf = std::cout.rdbuf();
  }
  std::ostream f(buf);

  int i=0;
  int noise_points_idx=0;
  for (PointCloud::const_iterator it = pc.begin(); it != pc.end(); ++it) {
    if(labels.at(i)==-1) hex_to_rgb(noise_color, &r, &g, &b);
    else hex_to_rgb(colors[labels.at(i)%30], &r, &g, &b);
    for(size_t k = 0; k<pc.at(0).data.size(); k++) {
      f << it->data.at(k)<<",";
      //f << it->x << "," << it->y << "," << original_time_information.at(i) << ",";)
    }
    f << original_time_information.at(i) << ",";
    if( (noise_points_idx<noise_points->size()) && (noise_points->at(noise_points_idx)==i) ) { 
      f << labels.at(i)<<", "<<r<<", "<<g<<", "<<b<<", "<<1<<"\n";
      noise_points_idx++;
    } else {
      f << labels.at(i)<<", "<<r<<", "<<g<<", "<<b<<", "<<0<<"\n";
    }
    i++;
  }
  return 0;
}



int FileParser::write_csv(std::string filename, PointCloud & pc, std::vector<int> & labels, bool append) {
  int r=0,g=0,b=0;
  // set ostream f to either file or stdout
  std::streambuf* buf;
  std::ofstream fhelp;
  if(!filename.empty()) {
    if(append==true) {
      fhelp.open(filename.c_str(), std::ios::out | std::ios::app);
    } else {
      fhelp.open(filename.c_str());
    }
    if (fhelp.fail())
      return 1;
    buf = fhelp.rdbuf();
  } else {
    buf = std::cout.rdbuf();
  }
  std::ostream f(buf);

  int i=0;
  for (PointCloud::const_iterator it = pc.begin(); it != pc.end(); ++it) {
    if(labels.at(i)==-1) hex_to_rgb(noise_color, &r, &g, &b);
    else hex_to_rgb(colors[labels.at(i)%30], &r, &g, &b);
    //else hex_to_rgb(colors[0], &r, &g, &b);
    //f << it->x << "," << it->y << "," << it->t << ",";
    for(size_t k = 0; k<pc[0].data.size(); k++) {
      f << it->data.at(k)<<",";
      //f << it->x << "," << it->y << "," << original_time_information.at(i) << ",";
    }
    f << original_time_information.at(i) << ",";
    //f << it->t<<",";
    //f << it->x << "," << it->y << "," << original_time_information.at(i) << ",";
    f << labels.at(i)<<", "<<r<<", "<<g<<", "<<b<<"\n";
    i++;
  }
  return 0;
}

int FileParser::write_loess(std::string filename, std::vector<PointCloud> * data, std::vector<int> * labels, bool append, double scale_factor_t) {
  // set ostream f to either file or stdout
  std::streambuf* buf;
  std::ofstream fhelp;
  if(!filename.empty()) {
    if(append==true) {
      fhelp.open(filename.c_str(), std::ios::out | std::ios::app);
    } else {
      fhelp.open(filename.c_str());
    }
    if (fhelp.fail())
      return 1;
    buf = fhelp.rdbuf();
  } else {
    buf = std::cout.rdbuf();
  }
  std::ostream f(buf);

  for(size_t i = 0; i<data->size(); i++) {
    for(size_t j = 0; j<data->at(i).size(); j++) {
      for(size_t k = 0; k<data->at(i).at(j).data.size(); k++) {
        f << data->at(i).at(j).data.at(k)<<",";
      }
      f << data->at(i).at(j).t/scale_factor_t<< "," <<labels->at(i)<<"\n";
      //f << original_time_information.at(i) << "," <<labels->at(i)<<"\n";
      //f << data->at(i).at(j).x << ", " << data->at(i).at(j).y << ", " << data->at(i).at(j).t/scale_factor_t << ", " <<labels->at(i)<<"\n";
    }
  }

  return 0;
}


int FileParser::write_csv_no_label(std::string filename, PointCloud & pc) {
  // set ostream f to either file or stdout
  std::streambuf* buf;
  std::ofstream fhelp;
  if(!filename.empty()) {
    fhelp.open(filename.c_str());
    if (fhelp.fail())
      return 1;
    buf = fhelp.rdbuf();
  } else {
    buf = std::cout.rdbuf();
  }
  std::ostream f(buf);

  int i=0;
  for (PointCloud::const_iterator it = pc.begin(); it != pc.end(); ++it) {
    //f << it->x << "," << it->y << "," << it->t << "\n";
    for(size_t k = 0; k<pc[0].data.size(); k++) {
      f << it->data.at(k)<<",";
      //f << it->x << "," << it->y << "," << original_time_information.at(i) << ",";
    }
    f << original_time_information.at(i) << "\n";
    i++;
  }
  return 0;
}


size_t FileParser::find_pos_start(std::ifstream & is, Point & p, bool is4d) {
  size_t ndims = is4d ? 4 : 3;
  std::vector<std::string> items;
  const char delimiter = ',';
  double goal_t=original_time_information.at(0); 
  std::string line;

  is.seekg(-1, std::ios_base::end);  // go to end of file and then one back to avoid EOF
  size_t num_characters = is.tellg();                    
  size_t curr_search_position = (size_t)(num_characters/2);   // binary search so start looking in the middle
  is.seekg(curr_search_position, std::ios::beg); // move cursor to that position
  while (is.peek() != '\n') {   // find the first newline by walking backwards
    is.seekg(-1, std::ios::cur) ;
  }
  int position_newline = is.tellg();    // store the position of that newline

  is.seekg(1, std::ios::cur) ;        // move cursor forward to avoid the newline
  std::getline(is, line, '\n');       // read line
  items.clear();
  split(line, items, delimiter);      // split by the , 
  double curr_t;
  if(is4d) curr_t = std::stod(items[3].c_str());   // read the time information
  else curr_t = std::stod(items[2].c_str());
  
  size_t lo=0; 
  size_t hi=num_characters;

  while(((curr_t)!=(goal_t)) && (is.tellg()>0) && (is.tellg()<num_characters)) {
    if(curr_t > goal_t) { // move left
      if(hi<position_newline) break;
      //curr_search_position = (int)(position_newline/2);
      hi = position_newline;
      curr_search_position = lo+(int)((hi-lo)/2);
    }
    else if(curr_t < goal_t) { // move right 
      if(lo>position_newline) break;
      lo = position_newline;
      curr_search_position = position_newline+(int)((num_characters-position_newline)/2);
    }

    is.seekg(curr_search_position, std::ios::beg); // move cursor to that position
    while ((is.peek() != '\n') && (is.tellg()!=0)) {   // find the first newline by walking backwards
      is.seekg(-1, std::ios::cur) ;
    }
    position_newline = is.tellg();    // store the position of that newline
    //std::cout<<position_newline<<"\n";
  
    is.seekg(1, std::ios::cur) ;        // move cursor forward to avoid the newline
    std::getline(is, line, '\n');       // read line
    items.clear();
    split(line, items, delimiter);      // split by the , 
    if(is4d) curr_t = std::stod(items[3].c_str());   // read the time information
    else curr_t = std::stod(items[2].c_str());
  }

  // there could be several points with the same timestamp so we actually have to find the matching point, but all have to be in the overlapping part so go to the beginning

  if(((double)std::stoi(items[0].c_str())==p.data.at(0)) && ((double)std::stoi(items[1].c_str())==p.data.at(1))) {
    if(is4d && ((double)std::stoi(items[2].c_str())==p.data.at(2))) return position_newline+1;
    else if(!is4d) return position_newline;
  }

  // if we get here we have a problem

  // move back until there is no more timestamp with the same value
  // std::cout<<is.tellg()<<"\n";
  //std::cout<<position_newline-7<<"\n";
  is.seekg(position_newline-7, std::ios::beg);
  //std::cout<<is.tellg()<<"\n";
  int old_position_newline=position_newline;
  while((curr_t==goal_t)&& (is.tellg()>0)) {
    old_position_newline=position_newline;
    while ((is.peek() != '\n') && (is.tellg()>0)) {   // find the first newline by walking backwards
      is.seekg(-1, std::ios::cur) ;
    }
    position_newline = is.tellg();    // store the position of that newline
    is.seekg(1, std::ios::cur) ;        // move cursor forward to avoid the newline
    std::getline(is, line, '\n');       // read line
    items.clear();
    split(line, items, delimiter);      // split by the , 
    if(is4d) curr_t = std::stod(items[3].c_str());   // read the time information
    else curr_t = std::stod(items[2].c_str());
  }

  is.seekg(0, std::ios::beg);

  return old_position_newline+1;

}



int FileParser::merge_csv(std::string filename, PointCloud & pc, std::vector<int> & labels, double scale_factor_t, bool is4d) {
  std::ifstream is;
  is.open(filename.c_str());

  if(is.peek() == std::ifstream::traits_type::eof()) {   // first timeslot can be written normally

    ssize_t max_label = -1;

    for(size_t i = 0; i<labels.size(); i++) {
      if(max_label < labels.at(i)) {
        max_label = labels.at(i);
      }
    }

    std::ofstream os_max_label;
    os_max_label.open(".max_label",std::ios::trunc);
    os_max_label << max_label;
    os_max_label.close();

    this->write_csv(filename, pc, labels, false);

  } else {

    std::ifstream is_max_label;
    is_max_label.open(".max_label");
    int max_label;
    is_max_label >> max_label;
    is_max_label.close();


    size_t start_pos = find_pos_start(is, pc.at(0),  is4d); // finds start position of the overlap with binary search, probably has bugs
    is.clear();
    is.seekg(start_pos, std::ios::beg);

    std::vector<std::string> items;
    const char delimiter = ',';
    std::string line;
    std::vector<int> labels_outfile;
    PointCloud old_pc;

    // read the overelapping section from the output file
    while (!is.eof()) {  
      std::getline(is, line, '\n');
      //std::cout<<line<<"\n";
      if (line[0] == '#' || line.empty() || line.find_first_not_of("\n\r\t ") == std::string::npos) continue;
      if(!isNullOrWhitespace(line)) {
        Point point;
        items.clear();
        split(line, items, delimiter);
        point.data.push_back((double)std::stoi(items[0].c_str()));
        point.data.push_back((double)std::stoi(items[1].c_str()));
        if(is4d) point.data.push_back((double)std::stoi(items[2].c_str()));
        //point.t = std::stod(items[2].c_str()) / scale_factor_t;
        // since when t is scaled you might get something like xxxxx.000000000000001. I want to compare for equality so its a problem. So I round... which means that the scaling might be off a little bit so I need to store the original data
        if(is4d) {
          point.t = std::round((std::stod(items[3].c_str()) * scale_factor_t) * 10.0 ) / 10.0;
          labels_outfile.push_back(std::stoi(items[4].c_str()));
        } else {
          point.t = std::round((std::stod(items[2].c_str()) * scale_factor_t) * 10.0 ) / 10.0;
          labels_outfile.push_back(std::stoi(items[3].c_str()));
        } 
        old_pc.push_back(point);
      }
    }

    std::unordered_map<int, int> cluster_merged;
    ssize_t curr_max_label=-1;

    for(size_t i = 0; i<old_pc.size(); i++) {
      if(labels.at(i)!=-1) {
        cluster_merged[labels.at(i)] = labels_outfile.at(i);
      }
    }
    for(size_t i = 0; i<pc.size(); i++) {
      if(labels.at(i)!=-1) {
        if(cluster_merged.find(labels.at(i)) != cluster_merged.end()) {
          labels.at(i) = cluster_merged[labels.at(i)];
        } else {
          labels.at(i) = labels.at(i)+max_label+1; 
        }
      }

      if(curr_max_label < labels.at(i)) {
        curr_max_label = labels.at(i);
      }
    }

    is.close();

    curr_max_label = (curr_max_label>max_label) ? curr_max_label : max_label;

    std::ofstream os_max_label;
    os_max_label.open(".max_label",std::ios::trunc);
    os_max_label << std::to_string(curr_max_label);
    os_max_label.close();

    this->write_csv(filename, pc, labels, true);
  }

 return 1;
}


/* int FileParser::merge_csv(std::string filename, PointCloud & pc, std::vector<int> & labels, double scale_factor_t, bool is4d) {
  std::ifstream in;
  in.open(filename.c_str());

  if(in.peek() == std::ifstream::traits_type::eof()) {   // first timeslot can be written normally

    in.close();
    this->write_csv(filename, pc, labels, false);

  } else {  // merging

    int label_first=labels.at(0);
    Point point_first(pc.at(0));
    int label=0;
    int max_label=0;
    const char delimiter = ',';
    std::ifstream infile(filename);
    std::string line;
    std::vector<std::string> items;
    if (infile.fail()) throw std::exception();
    while (!infile.eof()) {  // get data until the first Point was identified
      std::getline(infile, line, '\n');
      if (line[0] == '#' || line.empty() || line.find_first_not_of("\n\r\t ") == std::string::npos) continue;
      if(!isNullOrWhitespace(line)) {
        Point point;
        items.clear();
        split(line, items, delimiter);

        point.data.push_back((double)std::stoi(items[0].c_str()));
        point.data.push_back((double)std::stoi(items[1].c_str()));
        if(is4d) point.data.push_back((double)std::stoi(items[2].c_str()));
        //point.t = std::stod(items[2].c_str()) / scale_factor_t;
        // since when t is scaled you might get something like xxxxx.000000000000001. I want to compare for equality so its a problem. So I round... which means that the scaling might be off a little bit so I need to store the original data
        if(is4d) {
          original_time_information.push_back(std::stod(items[3].c_str()));
          point.t = std::round((std::stod(items[3].c_str()) * scale_factor_t) * 10.0 ) / 10.0;
          label = std::stoi(items[4].c_str());
        } else {
          original_time_information.push_back(std::stod(items[2].c_str()));
          point.t = std::round((std::stod(items[2].c_str()) * scale_factor_t) * 10.0 ) / 10.0;
          label = std::stoi(items[3].c_str());
        } 

        if(label > max_label) {
          max_label=label;
        }
        if(point==point_first) {
          break;
        }
      }
    }
    std::vector<int> merged_labels;  
    std::vector<bool> changed;
    int count=0;
    for(size_t i = 0; i<labels.size(); i++) {  // merge labels
      changed.push_back(false);
      if(labels.at(i) == label_first) {
        labels.at(i) = label;
        changed.at(i) = true;
        count++;
      }
    }
    merged_labels.push_back(label);
    changed.at(0) = true;

    int index = 1;

    while (!infile.eof()) {
      std::getline(infile, line, '\n');
      if(!isNullOrWhitespace(line)) {
        Point point;
        items.clear();
        split(line, items, delimiter);

        point.data.push_back((double)std::stoi(items[0].c_str()));
        point.data.push_back((double)std::stoi(items[1].c_str()));
        if(is4d) point.data.push_back((double)std::stoi(items[2].c_str()));
        //point.t = std::stod(items[2].c_str()) / scale_factor_t;
        // since when t is scaled you might get something like xxxxx.000000000000001. I want to compare for equality so its a problem. So I round... which means that the scaling might be off a little bit so I need to store the original data
        if(is4d) {
          original_time_information.push_back(std::stod(items[3].c_str()));
          point.t = std::round((std::stod(items[3].c_str()) * scale_factor_t) * 10.0 ) / 10.0;
          label = std::stoi(items[4].c_str());
        } else {
          original_time_information.push_back(std::stod(items[2].c_str()));
          point.t = std::round((std::stod(items[2].c_str()) * scale_factor_t) * 10.0 ) / 10.0;
          label = std::stoi(items[3].c_str());
        }

        if(label > max_label) {
          max_label=label;
        }
        if(std::find(merged_labels.begin(), merged_labels.end(),  label) == merged_labels.end()) {
          label_first = labels.at(index);
          for(size_t i = 0; i<labels.size(); i++) {
            if(labels.at(i) == label_first) {
              labels.at(i) = label;
              changed.at(i) = true;
            }
          }
          merged_labels.push_back(label_first);
        }
        index++;
      }
    }


    for(size_t i = index; i<labels.size(); i++) {
      if(!changed.at(i)) {
        labels.at(i) = labels.at(i)+1;
      }
    }    

    this->write_csv(filename, pc, labels, true);

  }
  return 0;
} */

//-------------------------------------------------------------------
// saves smoothen cloud as gnuplot script.
// The script is saved as file 'debug_smoothed.gnuplot' in the current
// working directory. It contains the original PointCloud *cloud* in
// black and the smoothed PointCloud *cloud_smooth* in red.
// *fname* is the name of the new file. The function returns false if
// an error occurred.
//-------------------------------------------------------------------
bool FileParser::debug_gnuplot(PointCloud &cloud, PointCloud &cloud_downsampled,
                   const char *fname, bool is4d) {
  // saves cloud in gnuplot file
  std::ofstream of(fname);
  of << std::fixed;  // set float style
  if (!of.is_open()) {
    std::cerr << "[Error] could not save under '" << fname << std::endl;
    return false;
  }
  // find ranges for each axis
  Point min, max;
  find_min_max_point(cloud, min, max);

  // Write header
 
  // when max and min are the same, the script can't be ploted, so the range
  // must be changed
  if (max.data.at(0) > min.data.at(0)) {
    of << "set xrange [" << min.data.at(0) << ":" << max.data.at(0) << "]\n";
  } else {
    of << "set xrange [" << (min.data.at(0) - 1.0) << ":" << (max.data.at(0) + 1.0) << "]\n";
  }
  if (max.data.at(1) > min.data.at(1)) {
    of << "set yrange [" << min.data.at(1) << ":" << max.data.at(1) << "]\n";
  } else {
    of << "set yrange [" << (min.data.at(1) - 1.0) << ":" << (max.data.at(1) + 1.0) << "]\n";
  }
  if(!is4d) {
    if (max.t > min.t) {
      of << "set zrange [" << min.t << ":" << max.t << "]\n ";
    } else {
      of << "set zrange [" << (min.t - 1.0) << ":" << (max.t + 1.0) << "]\n ";
    }
  } else {
    if (max.data.at(2) > min.data.at(2)) {
      of << "set zrange [" << min.data.at(2) << ":" << max.data.at(2) << "]\n ";
    } else {
      of << "set zrange [" << (min.data.at(2) - 1.0) << ":" << (max.data.at(2) + 1.0) << "]\n ";
    }
  }
  of << "splot ";
  
  of << "'-' with points lc 'black' title 'original', '-' with points lc "
        "'red' title 'smoothed'\n";
  for (PointCloud::const_iterator it = cloud.begin(); it != cloud.end(); ++it) {
    /* for(size_t k = 0; k<cloud[0].data.size(); k++) {
      //of << it->x << " " << it->y<<" "<<it->t<< std::endl;
      of << it->data.at(k)<<" ";
    }
    of << it->t<<std::endl; */
    of << *it << "\n";
    //std::cout<<*it<<"\n";
  }
  of << "e\n";

  size_t idx=0;
  for (PointCloud::const_iterator it = cloud_downsampled.begin();
       it != cloud_downsampled.end(); ++it) {
    /* for(size_t k = 0; k<cloud[0].data.size(); k++) {
      //of << it->x << " " << it->y<<" "<<it->t<< std::endl;
      of << it->data.at(k)<<" ";
    } */
    of << *it<<"\n";
  }
  of << "e\npause mouse keypress\n";

  of.close();
  return true;
}


