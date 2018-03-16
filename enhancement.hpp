/*
Fingerprint enhancer
Copyright (C) 2013-2015  Raja Lehtihet

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#ifndef ENHANCEMENT_HPP
#define ENHANCEMENT_HPP

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <vector>

#ifndef M_PI
#define M_PI 3.1415926535f
#endif

#ifdef min
#undef min
#endif

#ifdef max
#undef max
#endif


typedef float				real32;
typedef double				real64;

namespace fprint
{

struct vec2
{
	float x, y;

	vec2() : x(0.0f), y(0.0f) {}
	vec2(const float v) : x(v), y(v) {}
	vec2(const float ix, const float iy) : x(ix), y(iy) {}

	operator float *(){ return &x; }
	operator const float *() const { return &x; }

	void operator += (const vec2 &v);
	void operator -= (const vec2 &v);
	void operator *= (const float s);
	void operator *= (const vec2 &v);
	void operator /= (const float s);
	void operator /= (const vec2 &v);
};

vec2 operator + (const vec2 &u, const vec2 &v);
vec2 operator + (const vec2 &v, const float s);
vec2 operator - (const vec2 &u, const vec2 &v);
vec2 operator - (const vec2 &v, const float s);
vec2 operator - (const vec2 &v);
vec2 operator * (const vec2 &u, const vec2 &v);
vec2 operator * (const float s, const vec2 &v);
vec2 operator * (const vec2 &v, const float s);
vec2 operator / (const vec2 &u, const vec2 &v);
vec2 operator / (const vec2 &v, const float s);

bool operator == (const vec2 &u, const vec2 &v);
bool operator != (const vec2 &u, const vec2 &v);

float length(const vec2 &v);

///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////

typedef std::vector<vec2> vec2_vec;

class edge2d
{
protected:
	uint32_t	s[2];
	int32_t		f[2];	// facets
public:
	edge2d() { s[0] = s[1] = 0; f[0] = -1; f[1] = -1; }
	edge2d(uint32_t s0, uint32_t s1, int32_t f0 = -1, int32_t f1 = -1)
	{
		s[0] = std::min(s0, s1); s[1] = std::max(s0, s1);
		f[0] = f0; f[1] = f1;
	}

	edge2d(const edge2d &e) { s[0] = e.s[0]; s[1] = e.s[1]; f[0] = e.f[0]; f[1] = e.f[1]; }

	void		get_sites(uint32_t &s0, uint32_t &s1) const	{ s0 = s[0]; s1 = s[1]; }
	void		set_sites(uint32_t s0, uint32_t s1)			{  s[0] = std::min(s0, s1); s[1] = std::max(s0, s1); }
	void		get_facets(int32_t &f0, int32_t &f1)		{ f0 = f[0]; f1 = f[1]; }
	void		set_facets(int32_t f0, int32_t f1)			{ f[0] = f0; f[1] = f1; }
};

typedef std::vector<edge2d> edge2d_vec;

class facet2d
{
protected:
	std::vector<uint32_t>		sites;
public:
	facet2d() {}
	facet2d(const facet2d &f) { sites	= f.sites; }

	void		get_sites(std::vector<uint32_t> &sites) const { sites = this->sites; }
	void		set_sites(const std::vector<uint32_t> &sites) { this->sites = sites; }
	void		add_site(uint32_t s) { sites.push_back(s); }
};

typedef std::vector<facet2d> facet2d_vec;

bool build_del2d(const vec2_vec &in_sites, edge2d_vec &edges, facet2d_vec &facets);

///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////

///
/// image 2D matrix object
///
struct matrix
{
public:
	matrix();
	~matrix();

	matrix(uint32_t w, uint32_t h);
	matrix(uint32_t w, uint32_t h, real32 *d);

	uint32_t	get_width() { return width; }
	uint32_t	get_height() { return height; }

	void		set(uint32_t x, uint32_t y, real32 r) { data[y * width + x] = r; }
	real32		get(uint32_t x, uint32_t y) { return data[y * width + x]; }

	void		set(uint32_t i, real32 r) { data[i] = r; }
	real32		get(uint32_t i) { return data[i]; }
	real32*		data;
protected:
	uint32_t	width, height;

};

matrix* load_bmp(const char *file_name);
matrix* load_tiff(const char *file_name);
matrix* load_image(const char *file_name);

matrix* normalize(matrix* m);
matrix* normalize2(matrix* m);

matrix* gaussian(uint32_t s, real32 sigma);
matrix* gabor(uint32_t s, real32 lambda, real32 theta);

matrix* convolute(matrix *image, matrix *filter, bool mul_coef);

void get_mins(matrix* mat, std::vector<vec2> &mins, int win);


///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////

///
/// fingerprint related objects
///

class enhancer
{
public:
	enhancer();
	~enhancer();

	void			set(matrix* m, uint32_t bs, real32 thres);

	uint32_t		get_columns() { return columns; }
	uint32_t		get_rows() { return rows; }
	uint32_t		get_block_size() { return block_size; }

	uint32_t		reselect(real32 thres);
	bool			is_garbage(uint32_t c, uint32_t r); // x, y

	real32			variance(uint32_t c, uint32_t r); // x, y
	real32			orientation(uint32_t c, uint32_t r)
	{
		return blocks[c + r * columns].angle;
	}

	real32			ortho_orientation(uint32_t c, uint32_t r)
	{
		return blocks[c + r * columns].ortho_angle;
	}

	real32			lambda(uint32_t c, uint32_t r)
	{
		return blocks[c + r * columns].lambda;
	}

	matrix*			create_variance_image();
	matrix*			get_orig() { return orig; }
	matrix*			get_gradient_X() { return gradX; }
	matrix*			get_gradient_Y() { return gradY; }
	matrix*			get_gaussian() { return gauss; }
	matrix*			get_level1() { return level1; }
	matrix*			get_level2() { return level2; }

	const std::vector<vec2>&	get_minima() const { return minima; }
	const edge2d_vec&			get_edges() const { return edges; }
	const facet2d_vec&			get_faces() const { return faces; }

protected:

	struct block
	{
		uint32_t	x, y, w, h;		// block bounding box (start/end):  x <= i < x + w!
		bool		garbage;		// set if this block is garbage
		real32		mean;			// block mean
		real32		variance;		// block variance
		real32		angle;			// block orientation
		real32		ortho_angle;	// block orthogonal orientation
		real32		lambda;			// distance between ridges

		std::vector<size_t>	minima;	// minima points position
		matrix*		orig;
		matrix*		gabor;
		matrix*		level1;
		matrix*		level2;

		block() : orig(0), gabor(0), level1(0), level2(0) {}
		~block()
		{
			if( orig )
				delete orig;
			if( gabor )
				delete gabor;
			if( level1 )
				delete level1;
			if( level2 )
				delete level2;
		}
	};

	real32			compute_block_mean(const block &b);
	real32			compute_block_variance(const block &b);
	void			compute_orientations();
	void			elimnate_facets( const facet2d_vec &facets, uint32_t j, edge2d_vec &redges );
	void			compute_lambda_by_block(std::vector<real32> &lengths, block &b);
	void			compute_lambda();
	void			compute_gabor();
	void			compute_gabor_level1();
	void			smooth_orientation();
	void			smooth_lambda();

	std::vector<block>	blocks;
	matrix*			gaussian_mask;
	matrix*			orig;
	matrix*			gauss;
	matrix*			gradX;
	matrix*			gradY;
	matrix*			level1;
	matrix*			level2;
	matrix*			counts;
	uint32_t		columns;
	uint32_t		rows;
	uint32_t		block_size;
	real32			threshold;

	real32			min_variance;		// minimum variance
	real32			max_variance;		// maximum variance

	std::vector<vec2>	minima;

	edge2d_vec	edges;					// delaunay edges
	facet2d_vec	faces;					// delaunay facets
};
}
#endif // ENHANCEMENT_HPP
