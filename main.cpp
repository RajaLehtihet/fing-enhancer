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
#include <string>
#include <getopt.h>
#include <unistd.h>
#include <tiffio.h>
#include <enhancement.hpp>

using namespace std;
using namespace fprint;

void enhancement(const std::string& filename, int option);

int main(int argc, char* argv[])
{
	static const char *opt_string = "f:";

	int opt = getopt(argc, argv, opt_string);

	while ( opt != -1)
	{
		switch(opt)
		{
		case 'f':
		{
			enhancement(optarg, opt);
			break;
		}

		default :
			return 1;
		}

		opt = getopt(argc, argv, opt_string);
	}

	return 0;
}

void enhancement(const std::string& filename, int option)
{
	matrix *image = NULL;
	enhancer e;

	image = fprint::load_image(filename.c_str());

	//Check the image loading
	if( image == NULL )
	{
		std::cout << "invalid format" << std::endl;
	}
	cout << "the image file name is : " << filename << endl;

	//Enhancement following the algorithm steps
	e.set(image, 15, 0.05f);
	int w = e.get_level2()->get_width();
	int h = e.get_level2()->get_height();

	//Convert image data to float
	uint8_t *data	= new uint8_t[w * h];
	for( uint32 i = 0; i < w * h; i++ )
	{
		real32	r	= e.get_level2()->data[i];
		r		*= 255.0f;
		data[i]	= static_cast<uint8_t>(r);
	}

	//Prepare the enhanced image for writing
	TIFF* output_image = TIFFOpen("enhanced_image.tiff", "w");

	TIFFSetField(output_image, TIFFTAG_IMAGEWIDTH, w);
	TIFFSetField(output_image, TIFFTAG_IMAGELENGTH, h);
	TIFFSetField(output_image, TIFFTAG_BITSPERSAMPLE, 8);
	TIFFSetField(output_image, TIFFTAG_SAMPLESPERPIXEL, 1);
	TIFFSetField(output_image, TIFFTAG_ORIENTATION, (int)ORIENTATION_TOPLEFT);
	TIFFSetField(output_image, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
	TIFFSetField(output_image, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
	TIFFSetField(output_image, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);

	tsize_t image_s;
	if( (image_s = TIFFWriteEncodedStrip(output_image, 0, data, w * h * sizeof(uint8_t))) == -1)
	{
		std::cerr << "Unable to write tiff file: " << "enhanced_image.tiff" << std::endl;
	}
	else
	{
		std::cout << "Image is saved! size is : " << image_s << std::endl;
	}

	TIFFClose(output_image);

	delete [] data;
}
