
#include "./shared/all.hpp"
#include "./kokkos/types.hpp"

int main(int argc, char *argv[])
{
  Kokkos::initialize (argc, argv);
  {
    using parser_t = kokkosapp::commonTypes::parser_type;
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
    vb.push_back(parser.viewDepths()[0]== 1111.);
    vb.push_back(parser.viewAngles()[0]== 88.);
    vb.push_back(parser.viewPeriods()[0] == 43.);
    vb.push_back(parser.viewDelays()[0]== 12.);

    // // material
    vb.push_back(parser.getMaterialModelKind() == materialModelKind::unilayer);

    const auto a1 = parser.viewDiscontinuityDepthsKm();
    vb.push_back(a1.size()==1);
    vb.push_back(a1[0]==0.);

    const auto a2 = parser.viewDensityParametrization();
    vb.push_back(a2.size()==1);
    vb.push_back(a2.front().size()==3);
    vb.push_back(a2[0][0]==2000.);
    vb.push_back(a2[0][1]==0.);
    vb.push_back(a2[0][2]==0.);

    const auto a3 = parser.viewVelocityParametrization();
    vb.push_back(a3.size()==1);
    vb.push_back(a3.front().size()==3);
    vb.push_back(a3[0][0]==5000.);
    vb.push_back(a3[0][1]==0.);
    vb.push_back(a3[0][2]==0.);

    if (std::none_of(vb.begin(), vb.end(), std::logical_not<bool>()))
      std::puts("PASS");
    else
      std::puts("FAIL");
  }
  Kokkos::finalize();

  return 0;
}
