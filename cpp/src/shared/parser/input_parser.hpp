
#ifndef SHAXIPP_INPUT_PARSER_HPP_
#define SHAXIPP_INPUT_PARSER_HPP_

template <typename ...Mixins>
struct InputParser : Mixins...
{
  InputParser(int argc, char *argv[])
  {
    if (argc != 2){
      throw std::runtime_error("Wrong # of cmd line args: ./exe <input_file>");
    }

    const std::string inputFile = argv[1];
    std::cout << "*** Parsing input ***\n";
    this->parseGeneral(inputFile);
    this->parseIo(inputFile);
    this->parseMaterial(inputFile);
    this->parseForcing(inputFile);
    this->parseRom(inputFile);
    this->parseSampling(inputFile);
  }
};

#endif
