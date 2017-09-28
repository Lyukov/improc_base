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
    else if( axis == 'y' )  // Здесь между дублированием кода и проверкой условия на каждом пикселе
                            // я выбрал дублирование
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
    switch( angle )
    {
    case 0:
        return image.Copy();
    case 90:
    {
        ImageBase<PixelType> result = ImageBase<PixelType>( image.Height(), image.Width() );
        for( int j = 0; j < image.Height(); ++j )
        {
            for( int i = 0; i < image.Width(); ++i )
            {
                result( image.Height() - j - 1, i ) = image( i, j );
            }
        }
        return result;
    }
    case 180:
    {
        ImageBase<PixelType> result = ImageBase<PixelType>( image.Width(), image.Height() );
        for( int j = 0; j < image.Height(); ++j )
        {
            for( int i = 0; i < image.Width(); ++i )
            {
                result( image.Width() - i - 1, image.Height() - j - 1 ) = image( i, j );
            }
        }
        return result;
    }
    case 270:
    {
        ImageBase<PixelType> result = ImageBase<PixelType>( image.Height(), image.Width() );
        for( int j = 0; j < image.Height(); ++j )
        {
            for( int i = 0; i < image.Width(); ++i )
            {
                result( j, image.Width() - i - 1 ) = image( i, j );
            }
        }
        return result;
    }
    }
    return image.Copy();
}

int main( int argc, char **argv )
{
    
    return 0;
}
