/*    nuAnalysis class header file */


class nuAnalysis {

  bool               Debug;    // Global debug flag
  bool            FileFlag;    // Global flag to say read file
  bool           ChainFlag;    // Global flag to say chain files in dir
  std::string ROOTfilename;    // ROOT input file name
  std::string CHAINdirname;    // Director for chain of ROOT files

public:
  nuAnalysis(bool Dbg = true,
	     std::string Rfn = "Initaliation dummy",
	     std::string Cdn = "Initaliation dummy");

  ~nuAnalysis() {  }

private:
  void createChains() ;
  
};
