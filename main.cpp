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

int main( int argc, char **argv )
{
    /* if( !strcmp( argv[1], "mirror" ) )
     {
     }*/
    ColorByteImage image = ImageIO::FileToColorByteImage( argv[1] );
    ColorByteImage image1 = Mirror(image, 'y');
    ImageIO::ImageToFile(image1, "qweer.bmp");
    return 0;
}
