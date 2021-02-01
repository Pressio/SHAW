
#ifndef MATERIAL_MODEL_CREATE_HPP_
#define MATERIAL_MODEL_CREATE_HPP_

template<typename scalar_t, typename parser_t>
std::shared_ptr<MaterialModelBase<scalar_t>> createMaterialModel(const parser_t & parser)
{
  auto kind = parser.getMaterialModelKind();
  if (kind == materialModelKind::unilayer){
    using ret_t = UnilayerMaterialModel<scalar_t, parser_t>;
    return std::make_shared<ret_t>(parser);
  }
  else if (kind == materialModelKind::bilayer){
    using ret_t = BilayerMaterialModel<scalar_t, parser_t>;
    return std::make_shared<ret_t>(parser);
  }
  else if (kind == materialModelKind::prem){
    using ret_t = PremMaterialModel<scalar_t, parser_t>;
    return std::make_shared<ret_t>(parser);
  }
  else if (kind == materialModelKind::ak135f){
    using ret_t = Ak135fMaterialModel<scalar_t, parser_t>;
    return std::make_shared<ret_t>(parser);
  }
  else
    throw std::runtime_error("Cannot create material object, invalid material model kind");
}

#endif
