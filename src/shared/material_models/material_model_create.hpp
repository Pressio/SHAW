
#ifndef MATERIAL_MODEL_CREATE_HPP_
#define MATERIAL_MODEL_CREATE_HPP_

template<typename scalar_t, typename parser_t, typename mesh_info_t>
std::shared_ptr<MaterialModelBase<scalar_t>>
createMaterialModel(const parser_t & parser, const mesh_info_t & meshInfo)
{
  auto kind = parser.getMaterialModelKind();
  if (kind == materialModelKind::unilayer){
    using ret_t = UnilayerMaterialModel<scalar_t, parser_t>;
    return std::make_shared<ret_t>(parser, meshInfo);
  }
  else if (kind == materialModelKind::bilayer){
    using ret_t = BilayerMaterialModel<scalar_t, parser_t>;
    return std::make_shared<ret_t>(parser, meshInfo);
  }
  else if (kind == materialModelKind::prem){
    using ret_t = PremMaterialModel<scalar_t, parser_t>;
    return std::make_shared<ret_t>(parser, meshInfo);
  }
  else
    throw std::runtime_error("Cannot create material object, invalid material model kind");
}

#endif
