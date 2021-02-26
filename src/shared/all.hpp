/*
//@HEADER
// ************************************************************************
//
// all.hpp
//                     		Pressio/SHAW
//                         Copyright 2019
// National Technology & Engineering Solutions of Sandia, LLC (NTESS)
//
// Under the terms of Contract DE-NA0003525 with NTESS, the
// U.S. Government retains certain rights in this software.
//
// Pressio is licensed under BSD-3-Clause terms of use:
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its
// contributors may be used to endorse or promote products derived
// from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Francesco Rizzi (fnrizzi@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

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

#include "./various/signal.hpp"
#include "./nominal_to_grid_mappers/map_point_source_to_velocity_grid_point.hpp"
#include "./nominal_to_grid_mappers/map_nominal_location_to_velocity_grid_point.hpp"

#endif
