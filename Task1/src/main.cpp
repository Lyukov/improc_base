#include <random>
#include "task1.hpp"
int main(int argc, char **argv)
{
    std::mt19937 generator(42);
    int N = 1000;
    ColorByteImage noise(N, N);
    for(int i = 0; i < N; ++i)
    {
        for(int j = 0; j < N; ++j)
        {
            noise(i, j) = ColorBytePixel(generator() % 256, generator() % 256, generator() % 256);
        }
    }
    ImageIO::ImageToFile(noise, "noise1000.bmp");
    return 0;

    if(!strcmp(argv[1], "mirror"))
    {
        ColorByteImage image = ImageIO::FileToColorByteImage(argv[3]);
        ImageIO::ImageToFile(Mirror(image, argv[2][0]), argv[4]);
    }
    else if(!strcmp(argv[1], "rotate"))
    {
        bool direction = true;
        if(!strcmp(argv[2], "cw"))
        {
            direction = true;
        }
        else if(!strcmp(argv[2], "ccw"))
        {
            direction = false;
        }
        float angle = atof(argv[3]);
        if(angle == 90.f || angle == 180.f || angle == 270.f)
        {
            ColorByteImage image = ImageIO::FileToColorByteImage(argv[4]);
            image = Rotate(image, int(angle), direction);
            ImageIO::ImageToFile(image, argv[5]);
        }
        else
        {
            ColorFloatImage image = ImageIO::FileToColorFloatImage(argv[4]);
            image = Rotate(image, angle, direction);
            ImageIO::ImageToFile(image, argv[5]);
        }
    }
    else if(!strcmp(argv[1], "sobel"))
    {
        ColorFloatImage image = ImageIO::FileToColorFloatImage(argv[3]);
        ImageIO::ImageToFile(SobelFilter(image, argv[2][0]), argv[4]);
    }
    else if(!strcmp(argv[1], "median"))
    {
        int radius = atoi(argv[2]);
        ColorByteImage image = ImageIO::FileToColorByteImage(argv[3]);
        ImageIO::ImageToFile(MedianFilter(image, radius), argv[4]);
    }
    else if(!strcmp(argv[1], "gauss"))
    {
        float sigma = atof(argv[2]);
        float gamma = atof(argv[3]);
        ColorFloatImage image = ImageIO::FileToColorFloatImage(argv[4]);
        ImageIO::ImageToFile(GaussianFilter(image, sigma, gamma), argv[5]);
    }
    else if(!strcmp(argv[1], "gauss_edges"))
    {
        float sigma = atof(argv[2]);
        ColorFloatImage image = ImageIO::FileToColorFloatImage(argv[3]);
        ColorFloatImage gauss = GaussianFilter(image, sigma);
        image = GaussianFilter(image, sigma / 2.f);
        ColorFloatImage res(image.Width(), image.Height());
        for(int j = 0; j < image.Height(); ++j)
        {
            for(int i = 0; i < image.Width(); ++i)
            {
                res(i, j) = image(i, j) + gauss(i, j) * -1.f;
            }
        }
        GrayscaleFloatImage result = ToGrayscale(res);
        if(argc > 5)
        {
            float x = atof(argv[5]);
            for(int j = 0; j < result.Height(); ++j)
            {
                for(int i = 0; i < result.Width(); ++i)
                {
                    result(i, j) = abs(result(i, j)) * x;
                }
            }
        }
        ImageIO::ImageToFile(result, argv[4]);
    }
    else if(!strcmp(argv[1], "gradient"))
    {
        float sigma = atof(argv[2]);

        ColorFloatImage image = ImageIO::FileToColorFloatImage(argv[3]);
        GrayscaleFloatImage result = ToGrayscale(Gradient(image, sigma));
        if(argc > 5)
        {
            float x = atof(argv[5]);
            for(int j = 0; j < result.Height(); ++j)
            {
                for(int i = 0; i < result.Width(); ++i)
                {
                    result(i, j) *= x;
                }
            }
        }
        ImageIO::ImageToFile(result, argv[4]);
    }
    else if(!strcmp(argv[1], "gamma"))
    {
        float gamma = atof(argv[2]);
        ColorFloatImage image = ImageIO::FileToColorFloatImage(argv[3]);
        ImageIO::ImageToFile(GammaCorrection(image, gamma), argv[4]);
    }
    else if(!strcmp(argv[1], "scale"))
    {
        float x = atof(argv[2]);
        if(argc > 5)
        {
            float y = atof(argv[3]);
            ColorFloatImage image = ImageIO::FileToColorFloatImage(argv[4]);
            ImageIO::ImageToFile(Scale(image, x, y), argv[5]);
        }
        else
        {
            ColorFloatImage image = ImageIO::FileToColorFloatImage(argv[3]);
            ImageIO::ImageToFile(Scale(image, x), argv[4]);
        }
    }
    else if(!strcmp(argv[1], "grayscale"))
    {
        GrayscaleByteImage image = ImageIO::FileToGrayscaleByteImage(argv[2]);
        ImageIO::ImageToFile(image, argv[3]);
    }

    return 0;
}
