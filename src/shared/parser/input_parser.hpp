
#ifndef SHAXIPP_INPUT_PARSER_HPP_
#define SHAXIPP_INPUT_PARSER_HPP_

template <typename ...parts>
struct InputParser : parts...
{
  InputParser(int argc, char *argv[])
  {
    if (argc != 2){
      throw std::runtime_error
	("Wrong # of cmd line args, should be: ./exe <path-to-inputfile>");
    }

    const std::string inputFile = argv[1];
    std::cout << "*** Parsing input file ***\n";
    this->parseGeneral(inputFile);
    this->parseIo(inputFile);
    this->parseMaterial(inputFile);
    this->parseForcing(inputFile);
    this->parseRom(inputFile);
    this->parseSampling(inputFile);
  }
};

#endif
