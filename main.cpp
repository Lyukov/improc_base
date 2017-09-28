#include <string.h>
#include "imageformats.hpp"
#include "imageio.hpp"

float ToGray( ColorFloatPixel pixel )
{
    return pixel.b * 0.114f + pixel.g * 0.587f + pixel.r * 0.299f;
}

template <typename PixelType>
ImageBase<PixelType> Mirror( const ImageBase<PixelType> &image, const char axis = 'x' )
{
    assert( axis == 'x' || axis == 'y' );
    ImageBase<PixelType> result = ImageBase<PixelType>( image.Width(), image.Height() );
    if( axis == 'x' )
    {
        for( int j = 0; j < image.Height(); ++j )
        {
            for( int i = 0; i < image.Width(); ++i )
            {
                result( i, j ) = image( image.Width() - 1 - i, j );
            }
        }
    }
    else if( axis == 'y' )
    {
        for( int j = 0; j < image.Height(); ++j )
        {
            for( int i = 0; i < image.Width(); ++i )
            {
                result( i, j ) = image( i, image.Height() - 1 - j );
            }
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
    ColorByteImage image1 = Mirror( image, 'x' );
    ColorByteImage image2 = Mirror( image, 'y' );
    ImageIO::ImageToFile( image1, argv[2] );
    ImageIO::ImageToFile( image2, argv[3] );
    return 0;
}
