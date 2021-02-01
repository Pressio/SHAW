
#include "./fom/fom_problem_rank_one_forcing.hpp"

int main(int argc, char *argv[])
{
  Kokkos::initialize (argc, argv);
  {
    using parser_t = kokkosapp::commonTypes::parser_t;
    parser_t parser(argc, argv);

    std::vector<bool> vb;

    // general section
    vb.push_back(parser.getMeshDir() == "fullMeshString");
    vb.push_back(parser.checkDispersion() == false);
    vb.push_back(parser.checkCfl() == false);
    vb.push_back(parser.getTimeStepSize() == 1.5);
    vb.push_back(parser.getNumSteps() == 100);

    // i/o section
    vb.push_back(parser.enableSnapshotMatrix() == true);
    vb.push_back(parser.writeSnapshotsBinary() == false);
    vb.push_back(parser.getSnapshotFreq(dofId::vp) == 3);
    vb.push_back(parser.getSnapshotFreq(dofId::sp) == 4);
    vb.push_back(parser.getSnapshotFileName(dofId::vp) == "snaps_vp");
    vb.push_back(parser.getSnapshotFileName(dofId::sp) == "snaps_sp");

    vb.push_back(parser.enableSeismogram() == true);
    vb.push_back(parser.writeSeismogramBinary() == true);
    vb.push_back(parser.getSeismogramFileName() == "myfileName");
    vb.push_back(parser.getSeismoFreq() == 12 );

    using ang_t = typename parser_t::receivers_loc_t;
    ang_t goldAngles = {25, 50, 120, 160};
    vb.push_back(parser.getSeismoReceiversAnglesDeg() == goldAngles );

    // forcing
    vb.push_back(parser.getSourceSignalKind() == signalKind::sinusoid);
    vb.push_back(parser.getSourceProperty("depth") == 1111.);
    vb.push_back(parser.getSourceProperty("angle") == 88.);
    vb.push_back(parser.getSourceProperty("period") == 43.);
    vb.push_back(parser.getSourceProperty("delay") == 12.);


    if (std::none_of(vb.begin(), vb.end(), std::logical_not<bool>()))
      std::puts("PASS");
    else
      std::puts("FALSE");
  }
  Kokkos::finalize();

  return 0;
}
