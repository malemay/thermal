#include <cmath>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
//
arma::Col<arma::uword> cellFromXY (const arma::mat & coords, int rast_rows, int rast_cols, arma::vec extent);

// [[Rcpp::export]]
double assess_transform_cpp(const arma::vec & xval, const arma::vec & yval, const arma::mat & coords,
		double theta, double htrans, double vtrans, int nrows, int ncols, arma::vec extent,
		int min_overlap) {

	// Creating the transformation matrix from the supplied transformation parameters
	arma::mat tmat = {{cos(theta),  sin(theta),     0.0},
			  {-sin(theta), cos(theta),     0.0},
			  {htrans,      vtrans,         1.0}};

	// Performing the matrix product
	arma::mat new_coords = coords * tmat;

	// Creating a vector that will hold the indices of the coordinates that are within the raster
	arma::Col<arma::uword> coords_indices(new_coords.n_rows);
	int j = 0;

	// Identifying the rows that contain coordinates within the bounds of the target raster
	for(arma::uword i = 0; i < new_coords.n_rows; i++) {
		if(new_coords(i, 0) > extent[0] && new_coords(i, 0) < extent[1] && new_coords(i, 1) > extent[2] && new_coords(i, 1) < extent[3])
			coords_indices[j++] = i;
	}

	// We return a correlation of 0 if there is no overlap
	if(j < min_overlap) return 0.0;

	// Keeping only the the relevant values
	coords_indices = coords_indices.head_rows(j);

	// Keeping only the coordinates that are within the bounds of the raster
	new_coords = new_coords.rows(coords_indices);
	
	// Computing the indices of the target raster that correspond to the coordinates
	arma::Col<arma::uword> x_indices = cellFromXY(new_coords, nrows, ncols, extent);

	// Computing the correlation from the relevant raster values
	arma::mat output = cor(xval.rows(x_indices), yval.rows(coords_indices));

	// Returning the result
	return output[0];
}

arma::Col<arma::uword> cellFromXY (const arma::mat & coords, int rast_rows, int rast_cols, arma::vec extent) {
// size of x and y should be the same

	size_t size = coords.n_rows;
	arma::Col<arma::uword> cells(size);

	double xmin = extent[0];
	double xmax = extent[1];
	double ymin = extent[2];
	double ymax = extent[3];

	double yr_inv = rast_rows / (ymax - ymin);
	double xr_inv = rast_cols / (xmax - xmin);

	for (size_t i = 0; i < size; i++) {
		// cannot use trunc here because trunc(-0.1) == 0
		long row = std::floor((ymax - coords(i, 1)) * yr_inv);
		// points in between rows go to the row below
		// except for the last row, when they must go up
		if (coords(i, 1) == ymin) {
			row = rast_rows - 1 ;
		}

		long col = std::floor((coords(i, 0) - xmin) * xr_inv);
		// as for rows above. Go right, except for last column
		if (coords(i, 0) == xmax) {
			col = rast_cols - 1 ;
		}

		cells[i] = row * rast_cols + col;
	}

	return cells;
}

