#include <string.h>
#include <algorithm>
#include <vector>
#include "imageformats.hpp"
#include "imageio.hpp"

float ToGray( ColorFloatPixel pixel )
{
    return pixel.b * 0.114f + pixel.g * 0.587f + pixel.r * 0.299f;
}

template <typename T>
inline T sqr( const T &x )
{
    return x * x;
}

template <typename PixelType>
ImageBase<PixelType> Mirror( const ImageBase<PixelType> &image, const char axis = 'x' )
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

int main( int argc, char **argv )
{
    /* if( !strcmp( argv[1], "mirror" ) )
     {
     }*/
    ColorByteImage image = ImageIO::FileToColorByteImage( argv[1] );
    ColorByteImage image1 = MedianFilter(image, 1);
    ImageIO::ImageToFile( image1, "lena1.bmp" );
    image1 = MedianFilter(image, 2);
    ImageIO::ImageToFile( image1, "lena2.bmp" );
    image1 = MedianFilter(image, 3);
    ImageIO::ImageToFile( image1, "lena3.bmp" );
    return 0;
}
