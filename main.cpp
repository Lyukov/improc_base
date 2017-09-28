#include "imageformats.hpp"
#include "imageio.hpp"

float ToGray(ColorFloatPixel pixel)
{
	return pixel.b * 0.114f + pixel.g * 0.587f + pixel.r * 0.299f;
}

void TestFunc(char *inputfilename, char *outputfilename)
{
	ColorFloatImage image = ImageIO::FileToColorFloatImage(inputfilename);
	GrayscaleFloatImage res(image.Width(), image.Height());
	
	for (int j = 0; j < res.Height(); j++)
		for (int i = 0; i < res.Width(); i++)
			res(i, j) = ToGray(image(res.Width() - 1 - i, j));

	ImageIO::ImageToFile(image, outputfilename);
}

int main(int argc, char **argv)
{
	if (argc < 3)
        return 0;
    
   // ColorFloatImage res = ImageIO::FileToColorFloatImage(argv[1]);
   // ImageIO::ImageToFile(res, argv[2]);
	TestFunc(argv[1], argv[2]);
	return 0;
}
