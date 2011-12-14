#ifndef FILEPREFIX_H
#define FILEPREFIX_H

#include <string>

struct FilePrefix
{
  FilePrefix(const std::string& filePrefix)
  {
    // Strip any existing extension
    std::string filename(filePrefix);
    std::string::size_type idx;

    idx = filename.rfind('.');

    if(idx != std::string::npos)
    {
	//std::string extension = filename.substr(idx+1);
	prefix = filename.substr(0, idx);
    }
    else
    {
	// No extension found
	prefix = filePrefix;
    }
    
    std::cout << "Changed prefix from " << filePrefix << " to " << prefix << std::endl; 
  }
  
  std::string prefix;
};

#endif
