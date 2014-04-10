#ifndef _GRID_FUNCTION_HPP_
#define _GRID_FUNCTION_HPP_

#include <vector>
#include <map>
#include <array>
#include <algorithm>

template <size_t dim> using GridIndex = std::array<size_t, dim>;
template <size_t dim> using GridCoordinate = std::array<float, dim>;
template <size_t dim> using GridLinearCombination = std::map<GridIndex<dim>, float>;

template <typename T, size_t dim> class GridFunction {
public:
    GridCoordinate<dim> minbounds;
    GridCoordinate<dim> stepsize;
    const GridCoordinate<dim>& size() const { return stepno; };
    void set_size(const GridIndex<dim>& new_size) { 
	size_t scalar_size = 1;
	for (size_t d=0; d< dim; d++)
	    scalar_size *= new_size[d];
	data.resize(scalar_size);
	stepno = new_size;
    };
private:
    std::vector<T> data;
    GridIndex<dim> stepno;
    template <size_t d> size_t ND2Scalar(const GridIndex<d>& index) const {
	GridIndex<d-1> child_index;
	std::copy(std::next(index.begin()), index.end(), child_index.begin());
	return stepno[0] * this->ND2Scalar(child_index) + index[0];
    };
    size_t ND2Scalar(const GridIndex<1>& index) const { return index[0]; };
    template <size_t d> void NLinearCoefficientsPartial(GridLinearCombination<dim>& result, const GridCoordinate<d>& coordinates, const GridIndex<dim-d> PartialIndeces, const float& PartialWeight) const {
	float position = (coordinates[d-1]-minbounds[d-1])/stepsize[d-1];
	//TODO:usable exceptions
	if (position < 0) throw("Out of bounds.");
	if (position > stepno[d-1]) throw("Out of bounds.");
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
	if (position > stepno[0]) throw("Out of bounds.");
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
    const T& operator() (const GridIndex<dim>& indeces) const { return data[ND2Scalar(indeces)]; };
    T& operator() (const GridIndex<dim>& indeces) { return data[ND2Scalar(indeces)]; };
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
//    GridFunction (const GridCoordinate& minbounds, const GridCoordinate& stepsize, const )
};

#endif
