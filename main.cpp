#define _USE_MATH_DEFINES

#include <math.h>
#include <string.h>
#include <algorithm>
#include <vector>
#include "imageformats.hpp"
#include "imageio.hpp"

template <typename T>
inline T sqr( const T &x )
{
    return x * x;
}

float Gauss( float x, float y, float sigma )
{
    return exp( -( sqr( x ) + sqr( y ) ) / ( 2.f * sqr( sigma ) ) ) / ( 2.f * M_PI * sqr( sigma ) );
}

inline float GaussDx( float x, float y, float sigma )
{
    return Gauss( x, y, sigma ) * ( -x / sqr( sigma ) );
}

GrayscaleFloatImage ToGrayscale( const ColorFloatImage &image )
{
    GrayscaleFloatImage result( image.Width(), image.Height() );
    for( int j = 0; j < result.Height(); ++j )
    {
        for( int i = 0; i < result.Width(); ++i )
        {
            result( i, j ) = image( i, j ).toGray();
        }
    }
    return result;
}

template <typename PixelType>
ImageBase<PixelType> Mirror( const ImageBase<PixelType> &image, char axis = 'x' )
{
    assert( axis == 'x' || axis == 'y' );
    ImageBase<PixelType> result = ImageBase<PixelType>( image.Width(), image.Height() );
    for( int j = 0; j < image.Height(); ++j )
    {
        for( int i = 0; i < image.Width(); ++i )
        {
            result( i, j ) = ( axis == 'x' ) ? image( image.Width() - 1 - i, j ) :
                                               image( i, image.Height() - 1 - j );
        }
    }
    return result;
}

template <typename PixelType>
ImageBase<PixelType> Rotate( const ImageBase<PixelType> &image,
                             int angle = 90,
                             bool direction = true )  // true = clockwise, false = counterclockwise
{
    if( direction == false )
    {
        angle = 360 - angle;
    }
    angle %= 360;
    if( angle == 0 )
    {
        return image.Copy();
    }
    int width = ( angle == 0 || angle == 180 ) ? image.Width() : image.Height();
    int height = ( angle == 0 || angle == 180 ) ? image.Height() : image.Width();
    ImageBase<PixelType> result = ImageBase<PixelType>( width, height );
    for( int j = 0; j < height; ++j )
    {
        for( int i = 0; i < width; ++i )
        {
            PixelType pixel;
            switch( angle )
            {
            case 90:
                pixel = image( j, width - i - 1 );
                break;
            case 180:
                pixel = image( width - i - 1, height - j - 1 );
                break;
            case 270:
                pixel = image( height - j - 1, i );
                break;
            }
            result( i, j ) = pixel;
        }
    }
    return result;
}

template <typename PixelType>
ImageBase<PixelType> MedianFilter( const ImageBase<PixelType> &image, int radius = 1 )
{
    assert( radius > 0 );
    ImageBase<PixelType> result = ImageBase<PixelType>( image.Width(), image.Height() );
    std::vector<PixelType> window( sqr( 2 * radius + 1 ) );
    for( int j = 0; j < result.Height(); ++j )
    {
        for( int i = 0; i < result.Width(); ++i )
        {
            window.clear();
            for( int k = -radius; k < radius + 1; ++k )
            {
                for( int l = -radius; l < radius + 1; ++l )
                {
                    window.push_back( image( i + l, j + k ) );
                }
            }
            PixelType pixel;
            std::sort( window.begin(), window.end(),
                       []( const PixelType &p1, const PixelType &p2 ) { return p1.r < p2.r; } );
            pixel.r = window.at( window.size() / 2 + 1 ).r;
            std::sort( window.begin(), window.end(),
                       []( const PixelType &p1, const PixelType &p2 ) { return p1.g < p2.g; } );
            pixel.g = window.at( window.size() / 2 + 1 ).g;
            std::sort( window.begin(), window.end(),
                       []( const PixelType &p1, const PixelType &p2 ) { return p1.b < p2.b; } );
            pixel.b = window.at( window.size() / 2 + 1 ).b;
            result( i, j ) = pixel;
        }
    }
    return result;
}

template <typename PixelType, typename KernelType>
ImageBase<PixelType> Convolution( const ImageBase<PixelType> &image,
                                  const ImageBase<KernelType> &kernel )
{
    ImageBase<PixelType> result( image.Width(), image.Height() );
    int kerw = kernel.Width();
    int kerh = kernel.Height();
    for( int j = 0; j < result.Height(); ++j )
    {
        for( int i = 0; i < result.Width(); ++i )
        {
            PixelType res( 0 );
            for( int l = -kerh / 2; l <= kerh / 2; ++l )
            {
                for( int k = -kerw / 2; k <= kerw / 2; ++k )
                {
                    res += image( i - k, j - l ) * kernel( k + kerw / 2, l + kerw / 2 );
                }
            }
            result( i, j ) = res;
        }
    }
    return result;
}

template <typename PixelType>
ImageBase<PixelType> SobelFilter( const ImageBase<PixelType> &image, char axis = 'x' )
{
    assert( axis == 'x' || axis == 'y' );
    GrayscaleFloatImage kernel( 3, 3 );
    if( axis == 'x' )
    {
        float kerdata[9] = {-1.f, 0.f, 1.f,  //
                            -2.f, 0.f, 2.f,  //
                            -1.f, 0.f, 1.f};
        kernel = kerdata;
    }
    else if( axis == 'y' )
    {
        float kerdata[9] = {-1.f, -2.f, -1.f,  //
                            0.f,  0.f,  0.f,   //
                            1.f,  2.f,  1.f};
        kernel = kerdata;
    }
    ImageBase<PixelType> result = Convolution( image, kernel );
    result.for_each_pixel( []( PixelType p ) { return p + PixelType( 128 ); } );
    return result;
}

template <typename PixelType>
ImageBase<PixelType> GammaCorrection( const ImageBase<PixelType> &image, float gamma = 1.f )
{
    ImageBase<PixelType> result( image.Width(), image.Height() );
    for( int j = 0; j < result.Height(); ++j )
    {
        for( int i = 0; i < result.Width(); ++i )
        {
            PixelType gamma_px( gamma );
            result( i, j ) = apply(
                []( float x, float gamma ) -> float { return pow( ( x / 255.f ), gamma ) * 255.f; },
                image( i, j ), gamma_px );
        }
    }
    return result;
}

template <typename PixelType>
ImageBase<PixelType> GaussianFilter( const ImageBase<PixelType> &image,
                                     float sigma = 1.f,
                                     float gamma = 1.f )
{
    int win_size = ceil( sigma * 4.f );
    GrayscaleFloatImage kernel( win_size * 2 + 1, win_size * 2 + 1 );
    for( int j = 0; j < kernel.Height(); ++j )
    {
        for( int i = 0; i < kernel.Width(); ++i )
        {
            kernel( i, j ) = Gauss( i - win_size, j - win_size, sigma );
        }
    }
    if( gamma == 1.f )
    {
        return Convolution( image, kernel );
    }
    else
    {
        return GammaCorrection( Convolution( image, kernel ), gamma );
    }
}

template <typename PixelType>
ImageBase<PixelType> Gradient( const ImageBase<PixelType> &image, float sigma )
{
    int win_size = ceil( sigma * 4.f );
    GrayscaleFloatImage kernelDx( win_size * 2 + 1, win_size * 2 + 1 );
    for( int j = 0; j < kernelDx.Height(); ++j )
    {
        for( int i = 0; i < kernelDx.Width(); ++i )
        {
            kernelDx( i, j ) = GaussDx( i - win_size, j - win_size, sigma );
        }
    }
    GrayscaleFloatImage kernelDy = Mirror( Rotate( kernelDx, 90, true ), 'x' );  // Transpose
    ImageBase<PixelType> image_dx = Convolution( image, kernelDx );
    ImageBase<PixelType> image_dy = Convolution( image, kernelDy );
    ImageBase<PixelType> result( image.Width(), image.Height() );
    for( int j = 0; j < result.Height(); ++j )
    {
        for( int i = 0; i < result.Width(); ++i )
        {
            result( i, j ) = apply( []( float p1, float p2 ) -> float { return hypot( p1, p2 ); },
                                    image_dx( i, j ), image_dy( i, j ) );
        }
    }
    return result;
}

int main( int argc, char **argv )
{
    ColorFloatImage image = ImageIO::FileToColorFloatImage( argv[1] );
    ColorFloatImage image5 = GaussianFilter( image, 0.5f, 2.f );
    ImageIO::ImageToFile( image5, argv[2] );
    return 0;
}
