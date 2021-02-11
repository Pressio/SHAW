
#ifndef ALL_HPP_
#define ALL_HPP_

#include <fstream>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <vector>

#include "./constants.hpp"

#include "./enums/dof_id_enum.hpp"
#include "./enums/supported_signal_enums.hpp"
#include "./enums/supported_material_model_enums.hpp"
#include "./enums/supported_samplable_params_enums.hpp"

#include "./complexity.hpp"
#include "./various/print_perf.hpp"

#ifdef SHW_HAVE_TPL_EIGEN
#include "./meta/meta_eigen.hpp"
#endif
#include "./meta/meta_kokkos.hpp"
#include "./various/equality.hpp"
#include "./various/angular_helpers.hpp"

#include "./parser/parser_general_section.hpp"
#include "./parser/parser_io_section.hpp"
#include "./parser/parser_material_model.hpp"
#include "./parser/parser_forcing_section.hpp"
#include "./parser/input_parser.hpp"

#include "./io/matrix_write.hpp"
#include "./io/matrix_read.hpp"
#include "./io/read_basis.hpp"
#include "./io/read_reference_state.hpp"
#include "./io/vector_write.hpp"
#include "./io/vector_read.hpp"

#include "./checkers/check_dispersion_criterion.hpp"
#include "./checkers/check_cfl.hpp"

#include "./mesh_helpers/read_graph_file.hpp"
#include "./mesh_helpers/mesh_info.hpp"
#include "./mesh_helpers/read_vpcoeff_file.hpp"

#include "./material_models/material_model_unilayer.hpp"
#include "./material_models/material_model_bilayer.hpp"
#include "./material_models/material_model_prem.hpp"
#include "./material_models/material_model_create.hpp"

#include "./forcing/signal.hpp"
#include "./forcing/rank_one_forcing.hpp"
#include "./forcing/rank_two_forcing.hpp"

#include "./various/map_nominal_location_to_velocity_grid_point.hpp"
#include "./collectors/state_observer.hpp"
#include "./collectors/seismogram.hpp"

#endif
