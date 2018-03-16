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

#include "enhancement.hpp"

#include <cstdio>
#include <fstream>
#include <map>
#include <iostream>
#include <string>
#include <sstream>
#include <cfloat>
#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <memory>

namespace fprint {

void vec2::operator += (const vec2 &v)
{
	x += v.x;
	y += v.y;
}

void vec2::operator -= (const vec2 &v)
{
	x -= v.x;
	y -= v.y;
}

void vec2::operator *= (const float s)
{
	x *= s;
	y *= s;
}

void vec2::operator *= (const vec2 &v)
{
	x *= v.x;
	y *= v.y;
}

void vec2::operator /= (const float s)
{
	x /= s;
	y /= s;
}

void vec2::operator /= (const vec2 &v)
{
	x /= v.x;
	y /= v.y;
}

vec2 operator + (const vec2 &u, const vec2 &v)
{
	return vec2(u.x + v.x, u.y + v.y);
}

vec2 operator + (const vec2 &v, const float s)
{
	return vec2(v.x + s, v.y + s);
}

vec2 operator - (const vec2 &u, const vec2 &v)
{
	return vec2(u.x - v.x, u.y - v.y);
}

vec2 operator - (const vec2 &v, const float s)
{
	return vec2(v.x - s, v.y - s);
}

vec2 operator - (const vec2 &v)
{
	return vec2(-v.x, -v.y);
}

vec2 operator * (const vec2 &u, const vec2 &v)
{
	return vec2(u.x * v.x, u.y * v.y);
}

vec2 operator * (const float s, const vec2 &v)
{
	return vec2(v.x * s, v.y * s);
}

vec2 operator * (const vec2 &v, const float s)
{
	return vec2(v.x * s, v.y * s);
}

vec2 operator / (const vec2 &u, const vec2 &v)
{
	return vec2(u.x / v.x, u.y / v.y);
}

vec2 operator / (const vec2 &v, const float s)
{
	return vec2(v.x / s, v.y / s);
}

bool operator == (const vec2 &u, const vec2 &v)
{
	return (u.x == v.x && u.y == v.y);
}

bool operator != (const vec2 &u, const vec2 &v)
{
	return (u.x != v.x || u.y != v.y);
}

float length(const vec2 &v)
{
	return sqrtf(v.x * v.x + v.y * v.y);
}

////////////////////////////////////////////////////////////////////////////////

bool operator < (const edge2d &el, const edge2d &er)
{
	uint32_t sl[2], sr[2];
	el.get_sites(sl[0], sl[1]);
	er.get_sites(sr[0], sr[1]);
	if( sl[0] < sr[0] )
		return true;
	if( sl[0] > sr[0] )
		return false;
	if( sl[1] < sr[1] )
		return true;
	return false;
}


struct edge2d_cmp
{
	bool operator() (const edge2d &el, const edge2d &er)
	{
		return el < er;
	}
};

typedef std::map<edge2d, edge2d, edge2d_cmp> edge2d_map;
typedef std::vector<edge2d> edge2d_vec;

///////////////////////////////////////////////////////////////////////////////

bool write_qh2d_input(char *file_name, const vec2_vec &sites)
{
	std::fstream ofs;
	ofs.open(file_name, std::ios_base::out);
	if( ofs.fail() )
	{
		std::cout << "unable to open " << file_name << " for output" << std::endl;
		return false;
	}

	ofs << "2" << std::endl;
	ofs << sites.size() << std::endl;

	for(size_t i = 0; i < sites.size(); i++)
		ofs << sites[i].x << " " << sites[i].y << std::endl;

	ofs.close();
	return true;
}

uint32_t string_to_int(const std::string &str)
{
	std::istringstream iss(str);
	uint32_t num;
	iss >> num;
	return num;
}

bool read_qh2d_output(char *file_name, edge2d_map &edges, facet2d_vec &facets)
{
	edges.clear();

	std::fstream ifs;
	ifs.open(file_name, std::ios_base::in);
	if( ifs.fail() )
	{
		std::cout << "unable to open " << file_name << " for input" << std::endl;
		return false;
	}

	std::string str;
	// read the dimension (discard)
	std::getline(ifs, str);

	// read the number of points
	std::getline(ifs, str);
	uint32_t num_points, num_facets;

	std::istringstream iss(str);

	iss >> num_points >> num_facets;

	for( size_t i = 0; i < num_points; i++ )
		std::getline(ifs, str);

	for( size_t i = 0; i < num_facets; i++ )
	{
		std::getline(ifs, str);
		std::istringstream iss(str);
		int32_t np;
		iss >> np;
		std::vector<int32_t> points;

		facet2d	f;

		for( int p = 0; p < np; p++ )
		{
			int32_t id;
			iss >> id;
			points.push_back(id);
			f.add_site(id);
		}

		facets.push_back(f);

		for( int p = 0; p < np; p++ )
		{
			edge2d e(points[p], points[(p + 1) % points.size()]);
			if( edges.find(e) == edges.end() )
			{
				edges[e] = e;
			}

			int32_t f0, f1;
			edges[e].get_facets(f0, f1);
			if( f0 == -1 && f1 != facets.size())
				edges[e].set_facets(facets.size(), f1);
			else if( f1 == -1 && f0 != facets.size())
				edges[e].set_facets(f0, facets.size());
		}
	}
	ifs.close();

	return true;
}

bool build_del2d(const vec2_vec &in_sites, edge2d_vec &edges, facet2d_vec &facets)
{
	char infile_name[TMP_MAX];
	char outfile_name[TMP_MAX];

	edges.resize(0);

	if( in_sites.size() > 3 )
	{
		tmpnam(infile_name);
		tmpnam(outfile_name);

		write_qh2d_input(outfile_name, in_sites);

		std::string command("qhull d o TO ");
		command += infile_name;
		command += " < ";
		command += outfile_name;

		::system(command.c_str());

		edge2d_map edges_m;
		read_qh2d_output(infile_name, edges_m, facets);

		for( edge2d_map::iterator e = edges_m.begin(); e != edges_m.end(); e++ )
		{
			edges.push_back(e->second);
		}

		command = std::string("rm ");
		system((command + infile_name).c_str());
		system((command + outfile_name).c_str());
	}
	return true;
}

////////////////////////////////////////////////////////////////////////////////

matrix::matrix() : width(0), height(0), data(NULL) {}
matrix::~matrix() { delete[] data; }

matrix::matrix(uint32_t w, uint32_t h) : width(w), height(h)
{
	data	= new real32[w * h];
	memset(data, 0, sizeof(float) * w * h);
}

matrix::matrix(uint32_t w, uint32_t h, real32 *d) : width(w), height(h)
{
	data	= new real32[w * h];
	memcpy(data, d, w * h * sizeof(real32));
}


matrix* normalize(matrix* m)
{
	// find min/max
	real32	min_val	= FLT_MAX;
	real32	max_val	= -FLT_MAX;

	for( uint32_t i = 0; i < m->get_height() * m->get_width(); i++ )
	{
		real32	val = m->get(i);
		min_val	= std::min(min_val, val);
		max_val	= std::max(max_val, val);
	}

	// normalize
	real32	mx	= 1.0f - max_val;
	mx	= ( mx > 0 ) ? max_val : 1.0f;

	real32	n	= mx / (max_val - min_val);
	for( uint32_t i = 0; i < m->get_height() * m->get_width(); i++ )
	{
		real32	val = (m->get(i) - min_val) * n;
		m->set(i, val);
	}

	return m;
}

matrix* normalize2(matrix* m)
{
	// find min/max
	real32	min_val	= FLT_MAX;
	real32	max_val	= -FLT_MAX;

	for( uint32_t i = 0; i < m->get_height() * m->get_width(); i++ )
	{
		real32	val = m->get(i);
		min_val	= std::min(min_val, val);
		max_val	= std::max(max_val, val);
	}

	// normalize
	real32	mx	= 1.0f;

	real32	n	= mx / (max_val - min_val);
	for( uint32_t i = 0; i < m->get_height() * m->get_width(); i++ )
	{
		real32	val = (m->get(i) - min_val) * n;
		m->set(i, val);
	}

	return m;
}

matrix* gaussian(uint32_t s, real32 sigma)
{
	if( !(s & 1) )
		s++;

	real32	data[s * s];

	real32	invSigSq	= (fabs(sigma) < 0.001f) ? 0.0f : (1.0f / (sigma * sigma));
	real32	h_fs	= real32(s) * 0.5f;

	for( int32_t y = 0; y < s; y++ )
	{
		real32	fy = real32(y) - h_fs;
		real32	fySq	= fy * fy;
		for( int32_t x = 0; x < s; x++ )
		{
			real32 fx = real32(x) - h_fs;
			real32 fxSq	= fx * fx;
			data[x + s * y]	= exp( -(fxSq + fySq) * invSigSq );
		}
	}
	return new matrix(s, s, data);
}

matrix* gabor(uint32_t s, real32 lambda, real32 theta)
{
	real32 w2 = real32(s) * 0.5f;
	real32 h2 = real32(s) * 0.5f;

	real32 sigma_x = 4.0f;
	real32 sigma_y = 4.0f;

	matrix	*m	= new matrix(s, s);

	for( uint32_t y = 0; y < s; y++ )
	{
		for( uint32_t x = 0; x < s; x++ )
		{
			real32 xp = real32(x) - w2;
			real32 yp = real32(y) - h2;
			real32 x0 = xp * cos(theta) + yp * sin(theta);
			real32 y0 = -xp * sin(theta) + yp * cos(theta);
			real32 r = ((x0 * x0) / (sigma_x * sigma_x)) + ((y0 * y0) / (sigma_y * sigma_y));
			real32 g = exp(-0.5f * r) * cos((2 * M_PI * x0) / lambda);

			m->set(x, y, g);
		}
	}

	return m;
}

matrix* convolute(matrix *image, matrix *filter, bool mul_coef)
{
	matrix* res	= new matrix(image->get_width(), image->get_height());
	int32_t	fc	= filter->get_width() >> 1;
	int32_t	iw	= image->get_width();
	int32_t	ih	= image->get_height();

	for( uint32_t y = 0; y < image->get_height(); y++ )
	{
		for( uint32_t x = 0; x < image->get_width(); x++ )
		{
			real32	val		= real32(0);
			real32	coef	= real32(0);

			int32_t	y_min	= std::max<int32_t>(0, y - fc);
			int32_t	y_max	= std::min<int32_t>(image->get_height() - 1, y + fc);
			int32_t	x_min	= std::max<int32_t>(0, x - fc);
			int32_t	x_max	= std::min<int32_t>(image->get_width() - 1, x + fc);

			for( int32_t fy = y_min; fy <= y_max; fy++ )
			{
				for( int32_t fx = x_min; fx <= x_max; fx++ )
				{

					real32 fv = filter->get(fx - x_min, fy - y_min);
					coef += fv;
					val += image->get(fx, fy) * fv;
				}
			}

			if( mul_coef && coef != 0.0f )
				val /= coef;
			res->set(x, y, val);
		}
	}

	return res;
}

void get_mins(matrix* mat, std::vector<vec2> &mins, int win)
{
	int32_t	hw = win >> 1;

	for( int32_t y = 0; y < mat->get_height(); y++ )
	{
		for( int32_t x = 0; x < mat->get_width(); x++ )
		{
			real32	pv = mat->get(x, y);
			real32	min_val = FLT_MAX;

			int32_t	y_min	= std::max<int32_t>(0, y - hw);
			int32_t	y_max	= std::min<int32_t>(mat->get_height() - 1, y + hw);
			int32_t	x_min	= std::max<int32_t>(0, x - hw);
			int32_t	x_max	= std::min<int32_t>(mat->get_width() - 1, x + hw);

			for( int32_t fy = y_min; fy <= y_max; fy++ )
			{
				for( int32_t fx = x_min; fx <= x_max; fx++ )
				{

					if( fx == x && fy == y )
						continue;

					float v = mat->get(fx, fy);

					if( v < min_val )
						min_val = v;
				}
			}

			if( pv < min_val )
			{
				mins.push_back(vec2(real32(x),real32(y)));
			}
		}
	}
}
////////////////////////////////////////////////////////////////////////////////

enhancer::enhancer() :
	gaussian_mask(0), orig(0), gauss(0), gradX(0), gradY(0), level1(0), level2(0), counts(0)
{

}

enhancer::~enhancer()
{
	if( gaussian_mask )
		delete gaussian_mask;
	if( orig )
		delete orig;
	if( gauss )
		delete gauss;
	if( gradX )
		delete gradX;
	if( gradY )
		delete gradY;
	if( level1 )
		delete level1;
	if( level2 )
		delete level2;
	if( counts )
		delete counts;
}

real32 enhancer::compute_block_mean(const block &b)
{
	uint32_t	x, y;
	real32	mean = 0;

	for( y = b.y; y < (b.y + b.h); y++ )
	{
		for( x = b.x; x < (b.x + b.w); x++ )
		{
			mean += gauss->get(x, y);
		}
	}

	mean /= real32(b.w * b.h);

	return mean;
}

real32 enhancer::compute_block_variance(const block& b)
{
	uint32_t	x, y;

	real32	mean		= b.mean;
	real32	variance	= 0;

	for( y = b.y; y < (b.y + b.h); y++ )
	{
		for( x = b.x; x < (b.x + b.w); x++ )
		{
			real32	i = gauss->get(x, y);
			variance += (i - mean) * (i - mean);
		}
	}

	variance	/= (b.w * b.h);

	return variance;
}


void enhancer::set(matrix* m, uint32_t block_size, real32 thres)
{
	matrix noise(m->get_width(), m->get_height());

	for( uint32_t y = 0; y < m->get_height(); y++ )
	{
		for( uint32_t x = 0; x < m->get_width(); x++ )
		{
			real32 val	= (rand() & 0x1F) / 256.0f;
			noise.set(x, y, val);
		}
	}

	orig	= m;

	// first compute the gaussian
	fprint::matrix*	gauss_kernel	= fprint::gaussian(3, 1.8f);

	this->block_size	= block_size;

	// add the noise
	matrix* noise2	= convolute(&noise, gauss_kernel, true);
	for( uint32_t y = 0; y < m->get_height(); y++ )
	{
		for( uint32_t x = 0; x < m->get_width(); x++ )
		{
			real32 n	= noise2->get(x, y);
			real32 v	= orig->get(x, y);
			orig->set(x, y, v + n);

		}
	}

	fprint::normalize(orig);

	delete noise2;

	gauss	= convolute(orig, gauss_kernel, true);

	// compute blocks and variances
	uint32_t		j = 0;
	uint32_t		i = 0;

	blocks.clear();

	min_variance	= FLT_MAX;
	max_variance	= -FLT_MAX;

	for( uint32_t y = 0; y < m->get_height(); y += block_size, j++ )
	{
		i	= 0;
		for( uint32_t x = 0; x < m->get_width(); x += block_size, i++ )
		{
			block	b;
			b.x	= x;
			b.y	= y;
			b.w = std::min(m->get_width() - x, block_size);
			b.h = std::min(m->get_height() - y, block_size);

			b.garbage	= false;

			b.mean	= compute_block_mean(b);
			// variance computation should always follow mean computation!!!
			b.variance	= compute_block_variance(b);

			min_variance	= std::min(b.variance, min_variance);
			max_variance	= std::max(b.variance, max_variance);

			blocks.push_back(b);
		}
	}

	columns	= i;
	rows	= j;

	// select the blocks basing on threshold
	reselect(thres);

	// compute orientations
	compute_orientations();

	// compute minima
	vec2_vec mins;

	fprint::get_mins(gauss, mins, 5);

	for( size_t i = 0; i < mins.size(); i++ )
	{
		uint32_t	r	= (mins[i].y / block_size);
		uint32_t	c	= (mins[i].x / block_size);
		uint32_t	bidx	= r * columns + c;
		if( !blocks[bidx].garbage )
		{
			blocks[bidx].minima.push_back(minima.size());
			minima.push_back(mins[i]);
		}
	}

	// build delaunay
	build_del2d(minima, edges, faces);

	// compute lambda
	compute_lambda();

	// compute Gabor
	if( level1 )
		delete level1;
	if( level2 )
		delete level2;
	if( counts )
		delete counts;

	level1	= new matrix(columns * block_size /* * 3*/, rows * block_size /* * 3 */);
	level2	= new matrix(columns * block_size, rows * block_size);
	counts	= new matrix(columns * block_size, rows * block_size);
	compute_gabor();
}

uint32_t enhancer::reselect(real32 thres)
{
	uint32_t	num	= blocks.size();

	for( size_t b = 0; b < blocks.size(); b++ )
	{
		real32	val	= (blocks[b].variance - min_variance) / (max_variance - min_variance);
		if( val < thres )
		{
			blocks[b].garbage	= true;
			num--;
		} else blocks[b].garbage = false;
	}

	return num;
}

bool enhancer::is_garbage(uint32_t c, uint32_t r)
{
	return blocks[r * columns + c].garbage;
}

real32 enhancer::variance(uint32_t c, uint32_t r)
{
	return blocks[r * columns + c].variance;
}

matrix* enhancer::create_variance_image()
{
	matrix*	img	= new matrix(orig->get_width(), orig->get_height());

	for( size_t i = 0; i < blocks.size(); i++ )
	{
		real32	val	= (blocks[i].variance - min_variance) / (max_variance - min_variance);
		for( uint32_t y = blocks[i].y; y < blocks[i].y + blocks[i].h; y++ )
		{
			for( uint32_t x = blocks[i].x; x < blocks[i].x + blocks[i].w; x++ )
			{
				img->set(x, y, val);
			}
		}
	}

	return img;
}

void enhancer::compute_orientations()
{
	fprint::matrix sobel_X(3, 3);
	sobel_X.set(0, 0, -1);
	sobel_X.set(1, 0, -2);
	sobel_X.set(2, 0, -1);

	sobel_X.set(0, 1, 0);
	sobel_X.set(1, 1, 0);
	sobel_X.set(2, 1, 0);

	sobel_X.set(0, 2, 1);
	sobel_X.set(1, 2, 2);
	sobel_X.set(2, 2, 1);

	fprint::matrix sobel_Y(3, 3);
	sobel_Y.set(0, 0, -1);
	sobel_Y.set(0, 1, -2);
	sobel_Y.set(0, 2, -1);

	sobel_Y.set(1, 0, 0);
	sobel_Y.set(1, 1, 0);
	sobel_Y.set(1, 2, 0);

	sobel_Y.set(2, 0, 1);
	sobel_Y.set(2, 1, 2);
	sobel_Y.set(2, 2, 1);

	gradX = fprint::convolute(gauss, &sobel_X, true);
	gradY = fprint::convolute(gauss, &sobel_Y, true);

	for( size_t i = 0; i < blocks.size(); i++ )
	{
		real32	vx	= 0;
		real32	vy	= 0;

		if( !blocks[i].garbage )
		{
			int32_t start_y	= std::max<int32_t>(0, blocks[i].y - block_size);
			int32_t start_x	= std::max<int32_t>(0, blocks[i].x - block_size);
			int32_t	end_y	= std::min<int32_t>(gradY->get_height(), blocks[i].y + blocks[i].h + block_size);
			int32_t	end_x	= std::min<int32_t>(gradX->get_width(), blocks[i].x + blocks[i].w + block_size);

			for( uint32_t y = start_y; y < end_y; y++ )
			{
				for( uint32_t x = start_x; x < end_x; x++ )
				{

					real32 gx = gradX->get(x, y);
					real32 gy = gradY->get(x, y);

					vx	+= 2 * gx * gy;
					vy	+= gx * gx - gy * gy;
				}
			}

			blocks[i].angle	= atan(tan(M_PI * 0.5 + 0.5 * atan2( vx, vy )));
			blocks[i].ortho_angle	= atan(tan(0.5 * atan2( vx, vy )));
		}
	}

	//smooth_orientation();
}

void enhancer::smooth_orientation()
{
	real32	o[rows][columns];

	for( uint32_t r = 0; r < rows; r++ )
	{
		for( uint32_t c = 0; c < columns; c++ )
		{
			real32	oval	= 0;
			real32	factor	= 0;

			o[r][c]	= -FLT_MAX;

			for( int32_t sy = -1; sy <= 1; sy++ )
			{
				for( int32_t sx = -1; sx <= 1; sx++ )
				{
					int32_t ssx	= c + sx;
					int32_t ssy	= r + sy;

					if( ssx >= 0 && ssy >= 0 && ssx < columns && ssy < rows )
					{
						if( !blocks[ssy * columns + ssx].garbage )
						{
							oval	+= orientation(ssx, ssy);
							factor++;
						}
					}
				}
			}

			if( factor >= 4.0f )
				o[r][c]	= oval / factor;
		}
	}

	for( uint32_t r = 0; r < rows; r++ )
	{
		for( uint32_t c = 0; c < columns; c++ )
		{
			if( o[r][c] != -FLT_MAX )
			{
				blocks[r * columns + c].angle	= o[r][c];
				blocks[r * columns + c].ortho_angle	= M_PI * 0.5 + o[r][c];
			}
		}
	}
}

//
// eliminate the facets around a vertex keeping only the edges
//
void enhancer::elimnate_facets( const facet2d_vec &facets, uint32_t j, edge2d_vec &redges )
{

	for( size_t i=0; i < facets.size(); i++ )
	{
		std::vector<uint32_t> sl;
		facets[i].get_sites(sl);

		if ( j == sl[0] )
		{
			redges.push_back(edge2d(sl[1],sl[2]));
		}
		if ( j == sl[1] )
		{
			redges.push_back(edge2d(sl[2],sl[0]));
		}
		if ( j == sl[2] )
		{
			redges.push_back(edge2d(sl[1],sl[0]));
		}
	}
}

void intersect_sites_fn(const vec2 &pt0,const vec2 &pt1, const vec2 &pt3, const vec2 &pt2, vec2 &intersect_s)
{

	real32 pt0_sign	= - (pt0.x - pt2.x) * (pt3.y - pt2.y) + (pt3.x - pt2.x) * (pt0.y - pt2.y);
	real32 denom0	= (pt1.x - pt0.x) * (pt3.y - pt2.y) - (pt3.x - pt2.x) * (pt1.y - pt0.y);
	real32 ua		= pt0_sign / denom0;
	if( fabs(denom0) < 0.01f )
		std::cout << "*************" << std::endl << "Denom ERROR!!!!" << std::endl << "***************" << std::endl;
	intersect_s.x	= pt0.x + ua * (pt1.x - pt0.x);
	intersect_s.y	= pt0.y + ua * (pt1.y - pt0.y);

}

void is_intersect_segments(const vec2 &pt0, const float &teta, const edge2d_vec &redges,const vec2_vec &sitesM, vec2_vec &intersections)
{
	vec2 pt1, pt4;
	pt1.x= pt0.x - ( cos((M_PI/2) +teta) * 100 );
	pt1.y= pt0.y - ( sin(teta + (M_PI/2)) * 100 );
	pt4.x= pt0.x + ( cos((M_PI/2) +teta) * 100);
	pt4.y= pt0.y + ( sin(teta + (M_PI/2)) * 100);

	uint32_t sll[2];
	vec2 intersect_s;
	int k=0; int m=0;

	for( size_t j = 0; j <  redges.size(); j++ )
	{
		redges[j].get_sites(sll[0],sll[1]);
		vec2 pt2 = sitesM[sll[0]];
		vec2 pt3 = sitesM[sll[1]];
		real32 pt2_sign=(pt2.x-pt0.x)*(pt1.y-pt0.y)-(pt1.x-pt0.x)*(pt2.y-pt0.y);
		real32 pt3_sign=(pt3.x-pt0.x)*(pt1.y-pt0.y)-(pt1.x-pt0.x)*(pt3.y-pt0.y);
		real32 pt0_sign=(pt0.x-pt2.x)*(pt3.y-pt2.y)-(pt3.x-pt2.x)*(pt0.y-pt2.y);
		real32 pt1_sign=(pt1.x-pt2.x)*(pt3.y-pt2.y)-(pt3.x-pt2.x)*(pt1.y-pt2.y);

		real32 pt24_sign=(pt2.x-pt0.x)*(pt4.y-pt0.y)-(pt4.x-pt0.x)*(pt2.y-pt0.y);
		real32 pt34_sign=(pt3.x-pt0.x)*(pt4.y-pt0.y)-(pt4.x-pt0.x)*(pt3.y-pt0.y);
		//float pt0_sign=(pt0.x-pt2.x)*(pt3.y-pt2.y)-(pt3.x-pt2.x)*(pt0.y-pt2.y);
		real32 pt4_sign=(pt4.x-pt2.x)*(pt3.y-pt2.y)-(pt3.x-pt2.x)*(pt4.y-pt2.y);

		//compute denom to get intersection point coordinates

		if(((pt2_sign*pt3_sign)<0.0f)&&((pt1_sign*pt0_sign) < 0.0f))
		{
			k++;m++;//continue;
			//edges1_rem.push_back(edge2d(sll[0],sll[1]));

			intersect_sites_fn(pt0, pt1, pt3, pt2, intersect_s);
			intersections.push_back(intersect_s);

			//std::cout << "\t(" << pt2.x << ", " << pt2.y << ")-("  << pt3.x << ", " << pt3.y << ") with (" << pt0.x << ", " << pt0.y << ")-("  << pt1.x << ", " << pt1.y << ")" << std::endl;
		}

		if(((pt24_sign*pt34_sign)<0.0f)&&((pt4_sign*pt0_sign) < 0.0f))
		{
			k++;m++;//continue;
			//edges1_rem.push_back(edge2d(sll[0],sll[1]));

			intersect_sites_fn(pt0, pt4, pt3, pt2, intersect_s);
			intersections.push_back(intersect_s);

			//std::cout << "\t(" << pt2.x << ", " << pt2.y << ")-("  << pt3.x << ", " << pt3.y << ") with (" << pt0.x << ", " << pt0.y << ")-("  << pt1.x << ", " << pt1.y << ")" << std::endl;
		}
	}
}

void enhancer::compute_lambda_by_block(std::vector<real32> &lengths, block &b )
{
	b.lambda = 0;
	std::sort(lengths.begin(), lengths.end());
	if( lengths.size()!=0 )
	{
		b.lambda = lengths[lengths.size() >> 1];
		if( std::abs(b.lambda) < 0.01)
			b.lambda = 0;
	}
}

void enhancer::compute_lambda()
{
	vec2_vec	intersect_sites;	// intersection points temporary variable
	edge2d_vec	redges;				// remaining edges
	std::vector<real32>	lengths;

	for( size_t i = 0; i < blocks.size(); i++ )
	{
		if( !(blocks[i].garbage) )
		{
			lengths.resize(0);

			real32 theta			= blocks[i].angle;

			if( blocks[i].orig )
				delete blocks[i].orig;

			// create and fill the image with black
			blocks[i].orig			= new matrix(block_size * 3, block_size * 3);
			assert(blocks[i].orig);

			for( uint32_t y = 0; y < block_size * 3; y++ )
				for( uint32_t x = 0; x < block_size * 3; x++ )
					blocks[i].orig->set(x, y, 0);

			// copy 3x3 blocks from original image
			int32_t		start_cpy_x		= std::max<int32_t>(0, blocks[i].x - block_size);
			int32_t		start_cpy_y		= std::max<int32_t>(0, blocks[i].y - block_size);
			int32_t		end_cpy_y		= std::min<int32_t>(orig->get_height(), blocks[i].y + 2 * block_size);
			int32_t		end_cpy_x		= std::min<int32_t>(orig->get_width(), blocks[i].x + 2 * block_size);
			int32_t		offset_x		= (blocks[i].x - block_size < 0) ? block_size : 0;//start_cpy_x - blocks[i].x + block_size;
			int32_t		offset_y		= (blocks[i].y - block_size < 0) ? block_size : 0;//start_cpy_y - blocks[i].y + block_size;

			for( int32_t y = start_cpy_y, o_y = offset_y; y < end_cpy_y; y++, o_y++ )
			{
				for( int32_t x = start_cpy_x, o_x = offset_x; x < end_cpy_x; x++, o_x++ )
				{
					real32 val = orig->get(x, y);
					blocks[i].orig->set(o_x, o_y, val);
				}
			}

			// normalize the block
			fprint::normalize(blocks[i].orig);

			// compute lambda
			for( size_t j = 0; j < blocks[i].minima.size(); j++ )
			{
				redges.resize(0);
				intersect_sites.resize(0);

				elimnate_facets(faces, blocks[i].minima[j], redges);

				is_intersect_segments(minima[blocks[i].minima[j]], theta, redges, minima, intersect_sites);

				for( size_t l = 0; l < intersect_sites.size(); l++ )
				{
					vec2	pt0 = minima[blocks[i].minima[j]];
					vec2	pt1	= intersect_sites[l];
					real32	len = length(pt0 - pt1);
					lengths.push_back(len);
				}

				compute_lambda_by_block(lengths, blocks[i]);
			}
		}
		else
		{
			blocks[i].lambda	= -1;
		}
	}

	smooth_lambda();
}

void enhancer::smooth_lambda()
{
	real32	l[rows][columns];

	for( uint32_t r = 0; r < rows; r++ )
	{
		for( uint32_t c = 0; c < columns; c++ )
		{
			real32	lval	= 0;
			real32	factor	= 0;

			l[r][c]	= -FLT_MAX;

			for( int32_t sy = -1; sy <= 1; sy++ )
			{
				for( int32_t sx = -1; sx <= 1; sx++ )
				{
					int32_t ssx	= c + sx;
					int32_t ssy	= r + sy;

					if( ssx >= 0 && ssy >= 0 && ssx < columns && ssy < rows )
					{
						if( !blocks[ssy * columns + ssx].garbage && lambda(ssx, ssy) > 0 )
						{
							lval	+= lambda(ssx, ssy);
							factor++;
						}
					}
				}
			}

			if( factor >= 4.0f )
				l[r][c]	= lval / factor;
		}
	}

	for( uint32_t r = 0; r < rows; r++ )
	{
		for( uint32_t c = 0; c < columns; c++ )
		{
			uint32_t	idx	= r * columns + c;
			if( lambda(c, r) > 0 && l[r][c] != -FLT_MAX )
			{
				blocks[idx].lambda	= l[r][c];
			} //else blocks[idx].garbage	= true;
		}
	}
}

void enhancer::compute_gabor()
{
	// generate gaussian mask
	gaussian_mask	= fprint::gaussian(block_size * 3, block_size * 0.8f);
	fprint::normalize(gaussian_mask);

	//fprint::normalize(orig);
	for( size_t i = 0; i < blocks.size(); i++ )
	{
		if( !(blocks[i].garbage) )
		{
			// compute the gabor filter
			real32 theta		= blocks[i].angle;
			real32 lambda		= blocks[i].lambda;

			if( blocks[i].gabor )
				delete blocks[i].gabor;

			blocks[i].gabor		= /*gaussian(15, 2.8);*/gabor(9, lambda, M_PI - theta);

			blocks[i].level1	= convolute(blocks[i].orig, blocks[i].gabor, false);

			for( uint32_t h = 0; h < block_size * 3; h++ )
			{
				for( uint32_t w = 0; w < block_size * 3; w++ )
				{
					real32	g	= 1.0f;//gaussian_mask->get(w, h);
					real32	l	= blocks[i].level1->get(w, h);
					blocks[i].level1->set(w, h, l * g);
				}
			}

			fprint::normalize(blocks[i].level1);
		}
	}

	for( size_t i = 0; i < blocks.size(); i++ )
	{
		if( !(blocks[i].garbage) )
		{
			// copy to level 1 image
			for( uint32_t y = 0; y < /*block_size * 3*/blocks[i].h; y++ )
			{
				for( uint32_t x = 0; x < /*block_size * 3*/ blocks[i].w; x++ )
				{
					real32	res	= blocks[i].level1->get(x + block_size, y + block_size);
					//res	/= count;
					level1->set(blocks[i].x /* * 3*/ + x, blocks[i].y /* * 3*/ + y, res);

					int32_t	ix, iy;
					ix	= blocks[i].x + x - block_size;
					iy	= blocks[i].y + y - block_size;
					if( ix >= 0 && iy >= 0 && ix < block_size * columns && iy < block_size * rows )
					{
						real32 valCol = level2->get(ix, iy);
						real32 valCount	= counts->get(ix, iy);

						level2->set(ix, iy, valCol + res);
						valCount++;
						counts->set(ix, iy, valCount);
					}

				}
			}
		}
	}

	for( size_t i = 0; i < blocks.size(); i++ )
	{
		if( !(blocks[i].garbage) )
		{
			uint32_t h	= blocks[i].h + blocks[i].y;
			uint32_t w	= blocks[i].w + blocks[i].x;

			matrix	m(block_size, block_size);

			for( uint32_t y = blocks[i].y; y < h; y++ )
			{
				for( uint32_t x = blocks[i].x; x < w; x++ )
				{
					real32	res	= level2->get(x, y);
					real32	count = counts->get(x, y);

					if( count > 0 )
					{
						//level2->set(x, y, res/count);
						m.set(x - blocks[i].x, y - blocks[i].y, res / count);
					} else {
						level2->set(x, y, 0.0f);
					}
				}
			}

			fprint::normalize2(&m);
			for( uint32_t y = 0; y < blocks[i].h; y++ )
			{
				for( uint32_t x = 0; x < blocks[i].w; x++ )
				{
					real32	res	= m.get(x, y);
					level2->set(blocks[i].x + x, blocks[i].y + y, res);
				}
			}
		}
	}

	for( size_t i = 0; i < blocks.size(); i++ )
	{
		if( blocks[i].garbage )
		{
			uint32_t h	= blocks[i].h + blocks[i].y;
			uint32_t w	= blocks[i].w + blocks[i].x;

			for( uint32_t y = blocks[i].y; y < h; y++ )
			{
				for( uint32_t x = blocks[i].x; x < w; x++ )
				{
					level2->set(x, y, 0);
				}
			}
		}
	}

	fprint::normalize(level2);

	for( size_t i = 0; i < blocks.size(); i++ )
	{
		if( !blocks[i].garbage )
		{
			uint32_t h	= blocks[i].h + blocks[i].y;
			uint32_t w	= blocks[i].w + blocks[i].x;

			for( uint32_t y = blocks[i].y; y < h; y++ )
			{
				for( uint32_t x = blocks[i].x; x < w; x++ )
				{
					real32 bin = level1->get(x, y);
					//if( bin > 0.35f )
					//level2->set(x, y, 1.0f);
					//else level2->set(x, y, 0.0f);
					level2->set(x, y, bin);
				}
			}
		}
	}

}

}	// namespace fprint
