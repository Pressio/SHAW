
#include "./shared/all.hpp"
#include "./kokkos/types.hpp"

int main(int argc, char *argv[])
{
  Kokkos::initialize (argc, argv);
  {
    using parser_t = kokkosapp::commonTypes::parser_type;
    parser_t parser(argc, argv);

    std::vector<bool> vb;
    vb.push_back(parser.getMeshDir() == "fullMeshString");
    vb.push_back(parser.checkDispersion() == false);
    vb.push_back(parser.checkCfl() == false);
    vb.push_back(parser.getTimeStepSize() == 1.5);
    vb.push_back(parser.getNumSteps() == 100);

    // these are the default settings for i/o
    vb.push_back(parser.enableSnapshotMatrix() == false);
    vb.push_back(parser.writeSnapshotsBinary() == true);
    vb.push_back(parser.getSnapshotFreq(dofId::vp) == 0);
    vb.push_back(parser.getSnapshotFreq(dofId::sp) == 0);
    vb.push_back(parser.getSnapshotFileName(dofId::vp) == "snaps_vp");
    vb.push_back(parser.getSnapshotFileName(dofId::sp) == "snaps_sp");

    vb.push_back(parser.enableSeismogram() == false);
    vb.push_back(parser.writeSeismogramBinary() == true);
    vb.push_back(parser.getSeismogramFileName() == "seismogram");
    vb.push_back(parser.getSeismoFreq() == 0);

    using ang_t = typename parser_t::receivers_loc_t;
    ang_t goldAngles = {5,30,60,90,120,150,175};
    vb.push_back(parser.getSeismoReceiversAnglesDeg() == goldAngles );

    // forcing
    vb.push_back(parser.getSourceSignalKind() == signalKind::ricker);
    vb.push_back(parser.viewDepths()[0]== 1100.);
    vb.push_back(parser.viewAngles()[0]== 88.);
    vb.push_back(parser.viewPeriods()[0] == 40.);
    vb.push_back(parser.viewDelays()[0]== 10.);

    // material
    vb.push_back(parser.getMaterialModelKind() == materialModelKind::prem);

    if (std::none_of(vb.begin(), vb.end(), std::logical_not<bool>()))
      std::puts("PASS");
    else
      std::puts("FAIL");
  }
  Kokkos::finalize();

  return 0;
}
