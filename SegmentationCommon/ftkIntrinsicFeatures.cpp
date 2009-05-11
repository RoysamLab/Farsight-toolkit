#include "ftkIntrinsicFeatures.h"

namespace ftk
{
//CONSTANT STATIC VALUES FOR EACH FEATURE:
FeatureInfoType IntrinsicFeatures::Info[N] = { 
	{ "volume", "units", "description" },
	{ "integrated_intensity", "units", "description" },
	{ "eccentricity", "units", "description" },
	{ "elongation", "units", "description" },
	{ "orientation", "units", "description" },
	{ "bounding_box_volume", "units", "description" },
	{ "sum", "units", "description" },
	{ "mean", "units", "description" },
	{ "median", "units", "description" },
	{ "minimum", "units", "description" },
	{ "maximum", "units", "description" },
	{ "sigma", "units", "description" },
	{ "variance", "units", "description" },
	{ "surface_gradient", "units", "description" },
	{ "interior_gradient", "units", "description" },
	{ "surface_intensity", "units", "description" },
	{ "interior_intensity", "units", "description" },
	{ "intensity_ratio", "units", "description" },
	{ "radius_variation", "units", "description" },
	{ "surface_area", "units", "description" },
	{ "shape", "units", "description" },
	{ "shared_boundary", "units", "description" },
	{ "skew", "units", "description" },
	{ "energy", "units", "description" },
	{ "entropy", "units", "description" },
	{ "t_energy", "units", "description" },
	{ "t_entropy", "units", "description" },
	{ "inverse_diff_moment", "units", "description" },
	{ "inertia", "units", "description" },
	{ "cluster_shade", "units", "description" },
	{ "cluster_prominence", "units", "description" }
};

void IntrinsicFeatures::Print()
{

}

} //end namespace ftk 

