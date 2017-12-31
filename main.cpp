#include "imageformats.hpp"
#include "imageio.hpp"
#include "task1.hpp"

#include <string.h>
#include <cmath>
#include <iostream>
#include <queue>

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
            sum += sqr( diff.r ) + sqr( diff.g ) + sqr( diff.b );
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
    return 10.0 * log10( S * S / mse );
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
            sigma += sqr( diff.r ) + sqr( diff.g ) + sqr( diff.b );
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

    const double c1 = ( 0.01 * 255 ) * ( 0.01 * 255 );
    const double c2 = ( 0.03 * 255 ) * ( 0.03 * 255 );

    double m1_m2 = mean1.r * mean2.r + mean1.g * mean2.g + mean1.b * mean2.b;
    double m1_sqr = sqr( mean1.r ) + sqr( mean1.g ) + sqr( mean1.b );
    double m2_sqr = sqr( mean2.r ) + sqr( mean2.g ) + sqr( mean2.b );
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

GrayscaleByteImage Canny( const GrayscaleFloatImage& image,
                          float sigma,
                          float thr_high,
                          float thr_low )
{
    int win_size = ceilf( sigma * 3.f );
    GrayscaleFloatImage kernelDx( win_size * 2 + 1, win_size * 2 + 1 );
    for( int j = 0; j < kernelDx.Height(); ++j )
    {
        for( int i = 0; i < kernelDx.Width(); ++i )
        {
            kernelDx( i, j ) = GaussDx( i - win_size, j - win_size, sigma );
        }
    }
    GrayscaleFloatImage kernelDy = Mirror( Rotate( kernelDx, 90, true ), 'x' );  // Transpose
    GrayscaleFloatImage image_dx = Convolution( image, kernelDx );
    GrayscaleFloatImage image_dy = Convolution( image, kernelDy );
    GrayscaleFloatImage gradient( image.Width(), image.Height() );
    GrayscaleFloatImage gradient_direction( image.Width(), image.Height() );
    float max_gradient = 0.0f;
    for( int j = 0; j < image.Height(); ++j )
    {
        for( int i = 0; i < image.Width(); ++i )
        {
            gradient( i, j ) = hypot( image_dy( i, j ), image_dx( i, j ) );
            if( gradient( i, j ) > max_gradient )
            {
                max_gradient = gradient( i, j );
            }
            gradient_direction( i, j ) = atan2( image_dy( i, j ), image_dx( i, j ) );
        }
    }

    // Non-maximum suppression
    GrayscaleFloatImage preserved( image.Width(), image.Height() );
    for( int j = 0; j < image.Height(); ++j )
    {
        for( int i = 0; i < image.Width(); ++i )
        {
            float direction = gradient_direction( i, j ) * 180.0 / M_PI;
            if( direction < 0 )
            {
                direction += 180.0;
            }
            float grad1, grad2;
            if( direction <= 22.5 || direction >= 157.5 )
            {
                grad1 = gradient( i - 1, j );
                grad2 = gradient( i + 1, j );
            }
            else if( direction <= 67.5 )
            {
                grad1 = gradient( i + 1, j + 1 );
                grad2 = gradient( i - 1, j - 1 );
            }
            else if( direction <= 112.5 )
            {
                grad1 = gradient( i, j + 1 );
                grad2 = gradient( i, j - 1 );
            }
            else
            {
                grad1 = gradient( i + 1, j - 1 );
                grad2 = gradient( i - 1, j + 1 );
            }
            float grad = gradient( i, j );
            if( grad > grad1 && grad > grad2 )
            {
                preserved( i, j ) = grad;
            }
            else
            {
                preserved( i, j ) = 0.0;
            }
        }
    }

    // Thresholding
    GrayscaleByteImage thresholded( image.Width(), image.Height() );
    std::queue<std::pair<int, int> > Q;
    const unsigned char WEAK = 127;
    const unsigned char STRONG = 255;
    for( int j = 0; j < preserved.Height(); ++j )
    {
        for( int i = 0; i < preserved.Width(); ++i )
        {
            if( preserved( i, j ) < thr_low * max_gradient )
            {
                thresholded( i, j ) = 0;
            }
            else if( preserved( i, j ) < thr_high * max_gradient )
            {
                thresholded( i, j ) = WEAK;
            }
            else
            {
                thresholded( i, j ) = STRONG;
                Q.push( std::pair<int, int>( i, j ) );
            }
        }
    }
    // Hysteresis
    int movedir[8][2] = {{-1, -1}, {-1, 0}, {-1, 1}, {0, -1}, {0, 1}, {1, -1}, {1, 0}, {1, 1}};
    while( !Q.empty() )
    {
        int x = Q.front().first;
        int y = Q.front().second;
        Q.pop();
        for( int k = 0; k < 8; ++k )
        {
            int i = x + movedir[k][0];
            int j = y + movedir[k][1];
            if( thresholded( i, j ) == WEAK )
            {
                thresholded( i, j ) = STRONG;
                Q.push( std::pair<int, int>( i, j ) );
            }
        }
    }
   for( int j = 0; j < thresholded.Height(); ++j )
    {
        for( int i = 0; i < thresholded.Width(); ++i )
        {
            if( thresholded( i, j ) == WEAK )
            {
                thresholded( i, j ) = 0;
            }
        }
    }
    return thresholded;
}

float GaborFunction( float x,
                     float y,
                     float sigma,
                     float gamma,
                     float theta,
                     float lambda,
                     float psi )
{
    float x1 = x * cos( theta ) + y * sin( theta );
    float y1 = -x * sin( theta ) + y * cos( theta );
    return exp( -( sqr( x1 ) + sqr( gamma * y1 ) ) / ( 2.0 * sqr( sigma ) ) ) *
           cos( 2.0 * M_PI * x1 / lambda + psi );
}

template <typename PixelType>
ImageBase<PixelType> GaborFilter( const ImageBase<PixelType>& image,
                                  float sigma,
                                  float gamma,
                                  float theta,
                                  float lambda,
                                  float psi )
{
    int win_size = ceilf( sigma * 5.f );
    GrayscaleFloatImage kernel( win_size * 2 + 1, win_size * 2 + 1 );
    for( int j = 0; j < kernel.Height(); ++j )
    {
        for( int i = 0; i < kernel.Width(); ++i )
        {
            kernel( i, j ) =
                GaborFunction( i - win_size, j - win_size, sigma, gamma, theta, lambda, psi );
        }
    }
    ImageBase<PixelType> result = Convolution( image, kernel );
    //  result.for_each_pixel( []( PixelType p ) { return p * 0.5f + PixelType( 128 ); } );
    return result;
}

int main( int argc, char** argv )
{
    if( argc < 2 )
    {
        std::cout << "Bad command" << std::endl;
        return -1;
    }
    if( !strcmp( argv[1], "mse" ) )
    {
        if( argc < 4 )
        {
            std::cout << "Bad command" << std::endl;
            return -1;
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
            return -1;
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
            return -1;
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
            return -1;
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
            return -1;
        }
        ColorFloatImage image1 = ImageIO::FileToColorFloatImage( argv[2] );
        ColorFloatImage image2 = ImageIO::FileToColorFloatImage( argv[3] );
        std::cout << "MSE   = " << Metrics::MSE( image1, image2 ) << std::endl;
        std::cout << "PSNR  = " << Metrics::PSNR( image1, image2 ) << std::endl;
        std::cout << "SSIM  = " << Metrics::SSIM( image1, image2 ) << std::endl;
        std::cout << "MSSIM = " << Metrics::MSSIM( image1, image2 ) << std::endl;
    }
    else if( !strcmp( argv[1], "canny" ) )
    {
        if( argc < 7 )
        {
            std::cout << "Bad command" << std::endl;
            return -1;
        }
        float sigma = atof( argv[2] );
        float thr_high = atof( argv[3] );
        float thr_low = atof( argv[4] );
        GrayscaleFloatImage image = ImageIO::FileToGrayscaleFloatImage( argv[5] );
        GrayscaleByteImage canny = Canny( image, sigma, thr_high, thr_low );
        ImageIO::ImageToFile( canny, argv[6] );
    }
    else if( !strcmp( argv[1], "gabor" ) )
    {
        if( argc < 9 )
        {
            std::cout << "Bad command" << std::endl;
            return -1;
        }
        float sigma = atof( argv[2] );
        float gamma = atof( argv[3] );
        float theta = atof( argv[4] );
        theta *= M_PI / 180.f;
        float lambda = atof( argv[5] );
        float psi = atof( argv[6] );
        GrayscaleFloatImage image = ImageIO::FileToGrayscaleFloatImage( argv[7] );
        GrayscaleFloatImage gabor = GaborFilter( image, sigma, gamma, theta, lambda, psi );
        ImageIO::ImageToFile( gabor, argv[8] );
    }
    return 0;
}
