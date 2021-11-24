/*      -------- nuAnalysis:  Run control class header file --------      */
// RunControl call is a singleton.  Its role is to set flags and to allow
// various file names to be defined from the command line.
// The various flags can be printed with teh option "-h".  The flags are
// defined as follows:
//              -h : Generates "help" printout to show use of flags etc.
//              -d : Sets debug flag: RunControl::Debug = true
//   -f <filename> : Sets ROOT <filename> to be read (single file).  If -f
//                   is specified, <filename> must exist or execution is
//                   terminated
//   -c <dir name> : Sets directory containing ROOT files to be chained. If -c
//                   is specified, <dir name> must exist or execution is
//                   terminated
// Note that -f and -c can be used together, all files will be read.
class RunControl {

private:
  static RunControl* instance;

  bool               Debug;    // Global debug flag
  bool            FileFlag;    // Global flag to say read file
  bool           ChainFlag;    // Global flag to say chain files in dir
  std::string ROOTfilename;    // ROOT input file name
  std::string CHAINdirname;    // Director for chain of ROOT files
  
  RunControl(bool Dbg = true,
	     std::string Rfn = "Initaliation dummy",
	     std::string Cdn = "Initaliation dummy");

  ~RunControl() { };
  
public:
  static RunControl* getInstance();

  bool               getDebug(){ return Debug; };
  std::string getROOTfilename(){ return ROOTfilename; }
  std::string getCHAINdirname(){ return CHAINdirname; }
  
  void print();

  void ParseArgs(int nArgs, char *ArgV[]);
  
};
