#include "../Task1/imageformats.hpp"
#include "../Task1/imageio.cpp"
#include "../Task1/imageio.hpp"
#include "../Task1/pixelformats.hpp"

#include <string.h>
#include <cmath>
#include <iostream>

namespace Metrics
{
double MSE( const ColorFloatImage& image1, const ColorFloatImage& image2 )
{
    assert( image1.Height() == image2.Height() && image1.Width() == image2.Width() );
    double sum = 0.0;
    for( int j = 0; j < image1.Height(); ++j )
    {
        for( int i = 0; i < image1.Width(); ++i )
        {
            ColorFloatPixel diff = image1( i, j ) + -1.f * image2( i, j );
            sum += diff.r * diff.r + diff.g * diff.g + diff.b * diff.b;
        }
    }
    sum /= image1.Height() * image1.Width();
    return sum;
}

double PSNR( const ColorFloatImage& image1, const ColorFloatImage& image2 )
{
    assert( image1.Height() == image2.Height() && image1.Width() == image2.Width() );
    double mse = MSE( image1, image2 );
    double S = 255.0;
    return 10.0 * log( S * S / mse );
}
}

int main( int argc, char** argv )
{
    if( argc < 2 )
    {
        std::cout << "Bad command" << std::endl;
    }
    if( !strcmp( argv[1], "mse" ) )
    {
        if( argc < 4 )
        {
            std::cout << "Bad command" << std::endl;
        }
        ColorFloatImage image1 = ImageIO::FileToColorFloatImage( argv[2] );
        ColorFloatImage image2 = ImageIO::FileToColorFloatImage( argv[3] );
        double mse = Metrics::MSE( image1, image2 );
        std::cout << "MSE = " << mse << std::endl;
    } else if( !strcmp( argv[1], "psnr" ) )
    {
        if( argc < 4 )
        {
            std::cout << "Bad command" << std::endl;
        }
        ColorFloatImage image1 = ImageIO::FileToColorFloatImage( argv[2] );
        ColorFloatImage image2 = ImageIO::FileToColorFloatImage( argv[3] );
        double psnr = Metrics::PSNR( image1, image2 );
        std::cout << "PSNR = " << psnr << std::endl;
    }
    return 0;
}
