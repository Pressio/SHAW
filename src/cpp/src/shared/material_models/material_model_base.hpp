
#ifndef MATERIAL_MODEL_BASE_HPP_
#define MATERIAL_MODEL_BASE_HPP_

template <typename scalar_t>
class MaterialModelBase
{
public:
  virtual void computeAt(const scalar_t & radiusFromEarthCenterMeters,
			 const scalar_t & angleRadians,
			 scalar_t & density,
			 scalar_t & vs) const = 0;
};

#endif
