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

#include <iostream>
#include <cstdio>
#include <cstring>
#include <cctype>
#include <cfloat>

#include <cmath>

#include <malloc.h>

#include <tiffio.h>

#include "enhancement.hpp"

extern "C" unsigned char* read_bmp (FILE* fd, int* w, int* h, int* bps, int* spp,
									int* xres, int* yres,
									unsigned char** color_table, int* color_table_size,
									int* color_table_elements);

namespace fprint {
matrix* load_bmp(const char *file_name)
{
	matrix*	img;
	FILE*	f = fopen(file_name, "rb");
	if( f == NULL )
		return NULL;

	unsigned char* clr_tbl = 0;
	int clr_tbl_size = 0, clr_tbl_elems = 0;

	int32_t	w, h, bps, spp, xres, yres;

	uint8_t *cdata = read_bmp(f, &w, &h, &bps, &spp,
							  &xres, &yres, &clr_tbl,
							  &clr_tbl_size, &clr_tbl_elems);
	fclose(f);

	if( !cdata )
		return NULL;

	if( clr_tbl)
		free(clr_tbl);

	if( bps != 8 || spp != 1 )
	{
		free(cdata);
		return NULL;
	}

	std::cout << "image loaded. width: " << w << ", height: " << h
			  << ". bps: " << bps << ", spp: " << spp << std::endl;


	// convert to float
	real32 *data	= new real32[w * h];

	for( uint32 i = 0; i < w * h; i++ )
	{
		real32	r	= cdata[i];
		r		*= 1.0f / 255.0f;
		data[i]	= r;
	}

	img	= new matrix(w, h, data);

	delete [] data;
	// no color table anyway or RGB* ?
	//if (!clr_tbl || spp >= 3)
	return img;

}

matrix* load_tiff(const char *file_name)
{
	TIFF* tif = TIFFOpen(file_name, "r");
	if (tif)
	{
		uint32	w = 0, h = 0, bps = 0, spp = 0;
		size_t	npixels;
		uint32*	raster;
		matrix*	img	 = NULL;
		uint16	orient = 0;

		TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &w);
		TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h);
		TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bps);
		TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &spp);
		TIFFGetField(tif, TIFFTAG_ORIENTATION, &orient);

		std::cout << "bps: " << bps << ", spp: " << spp << ", orient: " << orient << std::endl;

		if( bps != 8 )
		{
			TIFFClose(tif);
			return img;
		}

		npixels = w * h;
		raster = (uint32*) _TIFFmalloc(npixels * sizeof (uint32));
		if (raster != NULL)
		{
			if (TIFFReadRGBAImage(tif, w, h, raster, 0))
			{
				// convert to float
				real32 *data	= new real32[w * h];

				for( uint32 p = 0; p < h; p++ )
				{
					uint32	r_ph	= (h - 1 - p) * w;
					uint32	_ph		= p * w;

					for( uint32 i = 0; i < w; i++ )
					{
						real32	r	= (raster[i + _ph] & 0xFF);
						r		*= 1.0f / 255.0f;
						data[r_ph + i]	= r;
					}
				}

				img	= new matrix(w, h, data);

				delete [] data;
			}
			_TIFFfree(raster);
		}
		TIFFClose(tif);

		return img;
	}

	return NULL;

}


matrix* load_image(const char *file_name)
{
	std::string	str;
	size_t		len = strlen(file_name);

	for( size_t i = 0; i < len; i++ )
		str += std::tolower(file_name[i]);

	// find the last dot
	size_t i = str.size() - 1;
	for( ; i >= 0; i-- )
	{
		if(str[i] == '.')
			break;
	}

	std::string ext;
	for( ; i < str.size(); ++i )
		ext += str.c_str()[i];

	if( ext == ".tiff" || ext == ".tif" )
		return load_tiff(file_name);
	else if( ext == ".bmp" )
		return load_bmp(file_name);

	return NULL;
}
}	// namespace fprint
