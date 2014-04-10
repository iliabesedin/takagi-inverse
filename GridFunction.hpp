#ifndef _GRID_FUNCTION_HPP_
#define _GRID_FUNCTION_HPP_

#include "NDVectorArray.hpp"
#include <map>

template <size_t dim> using GridIndex = VectorArrayIndex<dim>;
template <size_t dim> using GridCoordinate = std::array<float, dim>;
template <size_t dim> using GridLinearCombination = std::map<GridIndex<dim>, float>;

template <typename T, size_t dim> class GridFunction : public NDVector<T, dim> {
public:
    GridCoordinate<dim> minbounds;
    GridCoordinate<dim> stepsize;
private:
    std::vector<T> data;
    template <size_t d> void NLinearCoefficientsPartial(GridLinearCombination<dim>& result, const GridCoordinate<d>& coordinates, const GridIndex<dim-d> PartialIndeces, const float& PartialWeight) const {
	float position = (coordinates[d-1]-minbounds[d-1])/stepsize[d-1];
	//TODO:usable exceptions
	if (position < 0) throw("Out of bounds.");
	if (position > this->dsize(d-1)) throw("Out of bounds.");
	size_t li = floor(position);
	size_t ri = ceil(position);
	float rw = (position - li)*PartialWeight;
	float lw = (ri - position)*PartialWeight;
	GridIndex<dim-d+1> li_full, ri_full;
	li_full[0] = li; ri_full[0] = ri;
	std::copy(PartialIndeces.begin(), PartialIndeces.end(), std::next(li_full.begin()));
	std::copy(PartialIndeces.begin(), PartialIndeces.end(), std::next(ri_full.begin()));
	
	GridCoordinate<d-1> child_coordinates;
	std::copy(coordinates.begin(), std::prev(coordinates.end()), child_coordinates.begin());
	NLinearCoefficientsPartial(result, child_coordinates, li_full, lw);
	NLinearCoefficientsPartial(result, child_coordinates, ri_full, rw);
	return;
    };
    void NLinearCoefficientsPartial(GridLinearCombination<dim>& result, const GridCoordinate<1>& coordinates, const GridIndex<dim-1> PartialIndeces, const float& PartialWeight) const {
    	float position = (coordinates[0]-minbounds[0])/stepsize[0];
	//TODO:usable exceptions
	if (position < 0) throw("Out of bounds.");
	if (position > this->dsize(0)) throw("Out of bounds.");
	size_t li = floor(position);
	size_t ri = ceil(position);
	float rw = (position - li)*PartialWeight;
	float lw = (ri - position)*PartialWeight;
	GridIndex<dim> li_full, ri_full;
	li_full[0] = li; ri_full[0] = ri;
	std::copy(PartialIndeces.begin(), PartialIndeces.end(), std::next(li_full.begin()));
	std::copy(PartialIndeces.begin(), PartialIndeces.end(), std::next(ri_full.begin()));
	
	result[li_full] = lw;
	result[ri_full] = rw;
	return;
    };
public:
    GridLinearCombination<dim>  NLinearCoefficients (const GridCoordinate<dim>& coordinates) const {
	GridLinearCombination<dim> result;
	NLinearCoefficientsPartial(result, coordinates, GridIndex<0>(), 1);
	return result;
    };
    T NLinearValue (const GridCoordinate<dim>& coordinates) const {
	GridLinearCombination<dim> points = NLinearCoefficients (coordinates);
	T value = 0;
	for (typename GridLinearCombination<dim>::const_iterator point = points.begin(); point != points.end(); point++)
	    value += this->operator() (point->first) * point->second;
	return value;
    };
};

#endif
