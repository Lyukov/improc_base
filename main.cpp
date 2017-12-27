#include "../Task1/imageformats.hpp"
#include "../Task1/imageio.cpp"
#include "../Task1/imageio.hpp"
#include "../Task1/pixelformats.hpp"

#include <string.h>
#include <cmath>
#include <iostream>

namespace Metrics
{
//

double MSE( const ColorFloatImage& image1, const ColorFloatImage& image2 )
{
    assert( image1.Height() == image2.Height() && image1.Width() == image2.Width() );
    double sum = 0.0;
    for( int j = 0; j < image1.Height(); ++j )
    {
        for( int i = 0; i < image1.Width(); ++i )
        {
            ColorFloatPixel diff = image1( i, j ) - image2( i, j );
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

// Mean value of a part of the image
ColorFloatPixel Mean( const ColorFloatImage& image,
                      int x1 = 0,
                      int x2 = -1,
                      int y1 = 0,
                      int y2 = -1 )
{
    if( x2 == -1 )
    {
        x2 = image.Width();
        y2 = image.Height();
    }
    ColorFloatPixel mean = ColorFloatPixel( 0.0 );
    for( int j = y1; j < y2; ++j )
    {
        for( int i = x1; i < x2; ++i )
        {
            mean += image( i, j );
        }
    }
    const int N = ( x2 - x1 ) * ( y2 - y1 );
    mean = mean * ( 1.0 / N );
    return mean;
}

double StandardDerivation( const ColorFloatImage& image,
                           ColorFloatPixel mean = ColorFloatPixel( -1.0f ),
                           int x1 = 0,
                           int x2 = -1,
                           int y1 = 0,
                           int y2 = -1 )
{
    if( x2 == -1 )
    {
        x2 = image.Width();
        y2 = image.Height();
    }
    if( mean.r == -1.0f )
    {
        mean = Mean( image, x1, x2, y1, y2 );
    }
    double sigma = 0.0;
    for( int j = y1; j < y2; ++j )
    {
        for( int i = x1; i < x2; ++i )
        {
            ColorFloatPixel diff = image( i, j ) - mean;
            sigma += diff.r * diff.r + diff.g * diff.g + diff.b * diff.b;
        }
    }
    const int N = ( x2 - x1 ) * ( y2 - y1 );
    if( N == 1 )
    {
        return 0;
    }
    return sqrt( sigma / ( N - 1 ) );
}

double Covariance( const ColorFloatImage& image1,
                   const ColorFloatImage& image2,
                   ColorFloatPixel mean1 = ColorFloatPixel( -1.0f ),
                   ColorFloatPixel mean2 = ColorFloatPixel( -1.0f ),
                   int x1 = 0,
                   int x2 = -1,
                   int y1 = 0,
                   int y2 = -1 )
{
    if( x2 == -1 )
    {
        x2 = image1.Width();
        y2 = image1.Height();
    }
    if( mean1.r == -1.0f )
    {
        mean1 = Mean( image1, x1, x2, y1, y2 );
    }
    if( mean2.r == -1.0f )
    {
        mean2 = Mean( image2, x1, x2, y1, y2 );
    }
    double covariance = 0.0;
    for( int j = y1; j < y2; ++j )
    {
        for( int i = x1; i < x2; ++i )
        {
            ColorFloatPixel diff1 = image1( i, j ) - mean1;
            ColorFloatPixel diff2 = image2( i, j ) - mean2;
            covariance += diff1.r * diff2.r + diff1.g * diff2.g + diff1.b * diff2.b;
        }
    }
    const int N = ( x2 - x1 ) * ( y2 - y1 );
    if( N == 1 )
    {
        return 0;
    }
    covariance /= ( N - 1 );
    return covariance;
}

double SSIM( const ColorFloatImage& image1,
             const ColorFloatImage& image2,
             int x1 = 0,
             int x2 = -1,
             int y1 = 0,
             int y2 = -1 )
{
    if( x2 == -1 )
    {
        x2 = image1.Width();
        y2 = image1.Height();
    }

    ColorFloatPixel mean1 = Mean( image1, x1, x2, y1, y2 );
    double sigma1 = StandardDerivation( image1, mean1, x1, x2, y1, y2 );
    ColorFloatPixel mean2 = Mean( image2, x1, x2, y1, y2 );
    double sigma2 = StandardDerivation( image2, mean2, x1, x2, y1, y2 );
    double covariance = Covariance( image1, image2, mean1, mean2, x1, x2, y1, y2 );

    const double c1 = (0.01 * 255) * (0.01 * 255);
    const double c2 = (0.03 * 255) * (0.03 * 255);

    double m1_m2 = mean1.r * mean2.r + mean1.g * mean2.g + mean1.b * mean2.b;
    double m1_sqr = mean1.r * mean1.r + mean1.g * mean1.g + mean1.b * mean1.b;
    double m2_sqr = mean2.r * mean2.r + mean2.g * mean2.g + mean2.b * mean2.b;
    double ssim = ( ( 2 * m1_m2 + c1 ) * ( 2 * covariance + c2 ) ) /
                  ( ( m1_sqr + m2_sqr + c1 ) * ( sigma1 * sigma1 + sigma2 * sigma2 + c2 ) );
    return ssim;
}

double MSSIM( const ColorFloatImage& image1, const ColorFloatImage& image2 )
{
    assert( image1.Height() == image2.Height() && image1.Width() == image2.Width() );
    double mssim = 0.0;
    for( int j = 0; j < image1.Height(); j += 8 )
    {
        int j2 = j + 8 > image1.Height() ? image1.Height() : j + 8;
        for( int i = 0; i < image1.Width(); i += 8 )
        {
            int i2 = i + 8 > image1.Width() ? image1.Width() : i + 8;
            mssim += SSIM( image1, image2, i, i2, j, j2 );
        }
    }
    int P = ( ( image1.Width() + 7 ) / 8 ) * ( ( image1.Height() + 7 ) / 8 );
    mssim /= P;
    return mssim;
}

// end of namespace
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
        std::cout << "MSE = " << Metrics::MSE( image1, image2 ) << std::endl;
    }
    else if( !strcmp( argv[1], "psnr" ) )
    {
        if( argc < 4 )
        {
            std::cout << "Bad command" << std::endl;
        }
        ColorFloatImage image1 = ImageIO::FileToColorFloatImage( argv[2] );
        ColorFloatImage image2 = ImageIO::FileToColorFloatImage( argv[3] );
        std::cout << "PSNR = " << Metrics::PSNR( image1, image2 ) << std::endl;
    }
    else if( !strcmp( argv[1], "ssim" ) )
    {
        if( argc < 4 )
        {
            std::cout << "Bad command" << std::endl;
        }
        ColorFloatImage image1 = ImageIO::FileToColorFloatImage( argv[2] );
        ColorFloatImage image2 = ImageIO::FileToColorFloatImage( argv[3] );
        std::cout << "SSIM = " << Metrics::SSIM( image1, image2 ) << std::endl;
    }
    else if( !strcmp( argv[1], "mssim" ) )
    {
        if( argc < 4 )
        {
            std::cout << "Bad command" << std::endl;
        }
        ColorFloatImage image1 = ImageIO::FileToColorFloatImage( argv[2] );
        ColorFloatImage image2 = ImageIO::FileToColorFloatImage( argv[3] );
        std::cout << "MSSIM = " << Metrics::MSSIM( image1, image2 ) << std::endl;
    }
    else if( !strcmp( argv[1], "metrics" ) )
    {
        if( argc < 4 )
        {
            std::cout << "Bad command" << std::endl;
        }
        ColorFloatImage image1 = ImageIO::FileToColorFloatImage( argv[2] );
        ColorFloatImage image2 = ImageIO::FileToColorFloatImage( argv[3] );
        std::cout << "MSE   = " << Metrics::MSE( image1, image2 ) << std::endl;
        std::cout << "PSNR  = " << Metrics::PSNR( image1, image2 ) << std::endl;
        std::cout << "SSIM  = " << Metrics::SSIM( image1, image2 ) << std::endl;
        std::cout << "MSSIM = " << Metrics::MSSIM( image1, image2 ) << std::endl;
    }
    return 0;
}
