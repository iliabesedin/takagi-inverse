#ifndef _NDVECTORARRAY_HPP_
#define _NDVECTORARRAY_HPP_

#include <array>
#include <vector>
#include <algorithm>

template <size_t dim> using VectorArrayIndex = std::array<size_t, dim>;

template <size_t dim> class NDVectorArray {
public:
    virtual const VectorArrayIndex<dim>& size() const;
    virtual const size_t& dsize(const size_t& d) const;
    template <size_t d> size_t ND2Scalar(const VectorArrayIndex<d>& index) const {
	GridIndex<d-1> child_index;
	std::copy(std::next(index.begin()), index.end(), child_index.begin());
	return dsize(dim-d) * this->ND2Scalar(child_index) + index[0];
    };
    size_t ND2Scalar(const VectorArrayIndex<1>& index) const { return index[0]; };
};

template <typename T, size_t dim> class NDVector : public NDVectorArray<dim> {
private:
    std::vector<T> data;
    VectorArrayIndex<dim> elemno;
public:
    const VectorArrayIndex<dim>& size() const { return elemno; };
    const size_t& dsize(const size_t& d) const { return elemno[d]; };
    void set_size(const VectorArrayIndex<dim>& new_size) { { 
	size_t scalar_size = 1;
	for (size_t d=0; d< dim; d++)
	    scalar_size *= new_size[d];
	data.resize(scalar_size);
	elemno = new_size;
    };
    const T& operator() (const VectorArrayIndex<dim>& indeces) const { return data[ND2Scalar(indeces)]; };
    T& operator() (const VectorArrayIndex<dim>& indeces) { return data[ND2Scalar(indeces)]; };
};

#endif
