#pragma once

#define _USE_MATH_DEFINES

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <algorithm>
#include <iostream>
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
    int win_size = ceilf( sigma * 3.f );
    GrayscaleFloatImage kernel( win_size * 2 + 1, win_size * 2 + 1 );
    float sum = 0.f;
    for( int j = 0; j < kernel.Height(); ++j )
    {
        for( int i = 0; i < kernel.Width(); ++i )
        {
            kernel( i, j ) = Gauss( i - win_size, j - win_size, sigma );
            sum += kernel(i, j);
        }
    }
    for( int j = 0; j < kernel.Height(); ++j )
    {
        for( int i = 0; i < kernel.Width(); ++i )
        {
            kernel( i, j ) /= sum;
        }
    }
    if( sigma == 0.f ) {
		kernel( 0, 0 ) = 1.f;
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

inline float abs( float x ) { return ( x > 0 ) ? x : -x; }
template <typename PixelType>
PixelType BilinearInterpolate( const ImageBase<PixelType> &image, float a, float b )
{
    float i = floorf( a );
    float j = floorf( b );
    return image( i, j ) * ( i + 1 - a ) * ( j + 1 - b ) +
           image( i + 1, j ) * ( a - i ) * ( j + 1 - b ) +
           image( i, j + 1 ) * ( i + 1 - a ) * ( b - j ) +
           image( i + 1, j + 1 ) * ( a - i ) * ( b - j );
}

template <typename PixelType>
ImageBase<PixelType> Rotate( const ImageBase<PixelType> &image, float angle, bool direction = true )
{
    if( direction == false )
    {
        angle = -angle;
    }
    angle *= M_PI / 180.f;
    float sina = sin( angle );
    float cosa = cos( angle );
    int h = image.Height();
    int w = image.Width();
    int H = ceil( (float)h * abs( cosa ) + (float)w * abs( sina ) );
    int W = ceil( (float)h * abs( sina ) + (float)w * abs( cosa ) );
    ImageBase<PixelType> result( W, H );
    struct point
    {
        float x;
        float y;
    };
    std::vector<point> v( 4 );
    v[0].x = -cosa * w / 2 + sina * h / 2;
    v[1].x = cosa * w / 2 + sina * h / 2;
    v[2].x = -cosa * w / 2 - sina * h / 2;
    v[3].x = cosa * w / 2 - sina * h / 2;
    v[0].y = -sina * w / 2 - cosa * h / 2;
    v[1].y = sina * w / 2 - cosa * h / 2;
    v[2].y = -sina * w / 2 + cosa * h / 2;
    v[3].y = sina * w / 2 + cosa * h / 2;
    std::sort( v.begin(), v.end(), []( const point &p1, const point &p2 ) { return p1.x < p2.x; } );
    std::sort( v.begin() + 1, v.end() - 1,
               []( const point &p1, const point &p2 ) { return p1.y < p2.y; } );
    float c[4];
    if( abs( sina ) < 1e-6f || abs( cosa ) < 1e-6f )
    {
        c[0] = c[1] = c[2] = c[3] = 0.f;
    }
    else
    {
        c[0] = ( v[1].x - v[0].x ) / ( v[1].y - v[0].y );
        c[1] = ( v[2].x - v[0].x ) / ( v[2].y - v[0].y );
        c[2] = ( v[1].x - v[3].x ) / ( v[1].y - v[3].y );
        c[3] = ( v[2].x - v[3].x ) / ( v[2].y - v[3].y );
    }
    for( int y = -H / 2; y < H / 2; ++y )
    {
        int l, r;
        if( y <= v[0].y )
            l = v[0].x + c[0] * ( y - v[0].y );
        else
            l = v[0].x + c[1] * ( y - v[0].y );
        if( y <= v[3].y )
            r = v[3].x + c[2] * ( y - v[3].y );
        else
            r = v[3].x + c[3] * ( y - v[3].y );
        for( int x = l; x <= r; ++x )
        {
            float a = cosa * x + sina * y;
            float b = -sina * x + cosa * y;
            //result( x + W / 2, y + H / 2 ) = image( a + w / 2, b + h / 2 ); // Without bilinear interpolation
            result( x + W / 2, y + H / 2 ) = BilinearInterpolate( image, a + w / 2, b + h / 2 );
        }
    }
    return result;
}

template <typename PixelType>
ImageBase<PixelType> Scale( const ImageBase<PixelType> &image, float x, float y = -1.f )
{
    if( y == -1.f )
        y = x;
    assert( x > 0 && y > 0 );
    int W = floorf( image.Width() * x );
    int H = floorf( image.Height() * y );
    ImageBase<PixelType> result( W, H );
    for( int j = 0; j < H; ++j )
    {
        for( int i = 0; i < W; ++i )
        {
            float a = (float)i / x;
            float b = (float)j / y;
            // result( i, j ) = image( a, b ); // Without bilinear interpolation
            result( i, j ) = BilinearInterpolate( image, a, b );
        }
    }
    return result;
}