#include "TXLib.h"
#include <math.h>

RGBQUAD* Video_memory = NULL;

const bool BLINK_COURSOR = 0;

const int NUMBER_OF_GRAPHS = 1;

const int WINDOW_WIDTH = 600;
const int WINDOW_HEIGHT = 600;

const double GRAPH_HEIGHT = 0.8;
const double BUTTON_DISTANÑE = 0.1;

const double RADIUS_OF_SPHERE = 153;

const COLORREF WINDOW_BACKGROUND = RGB(0, 0, 0);

const COLORREF SHPERE_COLOR = RGB(128, 255, 0);
const COLORREF SQUARE_COLOR = RGB(74, 74, 74);
const COLORREF COORD_LINE_COLOR = RGB(144, 144, 144);
const COLORREF GRID_LINE_COLOR = RGB(144, 144, 144);

const COLORREF VECTOR_COLOR = RGB(255, 255, 255);

const COLORREF TEXT_COLOR = RGB(215, 215, 215);
const COLORREF ARROW_COLOR = RGB(144, 144, 144);
const COLORREF DOT_COLOR = RGB(255, 0, 0);

const double R_COLOR = 0.501960784;
const double G_COLOR = 0;
const double B_COLOR = 0.501960784;

const double R_SOURCE_LIGHT_COLOR = 1;
const double G_SOURCE_LIGHT_COLOR = 0;
const double B_SOURCE_LIGHT_COLOR = 0;

const int    COORD_L_THICKNESS = 3;
const double ARROW_THICKNESS = 0.25;

const int GRID_L_THICKNESS = 1;
const int SPLIT_DENCITY = 25;

#pragma once

class Vector
{
public:
    double x, y, z;

    Vector(double x_coord, double y_coord, double z_coord)
    {
        x = x_coord;
        y = y_coord;
        z = z_coord;
    }

    static Vector Vector_Normal(Vector* vector)
    {
        return *vector / sqrt(vector->x * vector->x + vector->y * vector->y + vector->z * vector->z);
    }

    Vector operator + (Vector b)
    {
        Vector result_vector(x + b.x, y + b.y, z + b.z);
        return result_vector;
    }
    Vector operator - (Vector b)
    {
        Vector result_vector(x - b.x, y - b.y, z - b.z);
        return result_vector;
    }
    Vector operator * (double lambda)
    {
        Vector result_vector(lambda * x, lambda * y, lambda * z);
        return result_vector;
    }
    Vector operator / (double lambda)
    {
        Vector result_vector(x / lambda, y / lambda, z / lambda);
        return result_vector;
    }
    double operator ^ (Vector b)
    {
        return (x * b.x + y * b.y + z * b.z);
    }
};

class Source
{
public:
    double fong, lambert;

    Source(double intensity, double flare_component)
    {
        lambert = intensity;
        fong = flare_component;
    }
};

class Sphere_Parameters
{
public:
    static double Intensity(Vector* radius, Vector* source)
    {
        return Vector::Vector_Normal(radius) ^ Vector::Vector_Normal(source);
    }

    static double Depth(double x, double y, double radius)
    {
        return sqrt(radius * radius - x * x - y * y);
    }

    static double Scalar_Flare(Vector* view, Vector* rv, Vector* source)
    {
        return pow((Vector::Vector_Normal(&((*rv) * 2 - *source)) ^ Vector::Vector_Normal(view)), 35);
    }
};

class Normalize
{
public:
    static Vector Normalize_Check(Vector* color)
    {
        color->x = Top_Limit(&(color->x));
        color->y = Top_Limit(&(color->y));
        color->z = Top_Limit(&(color->z));
        return *color;
    }

private:
    static double Top_Limit(double* number)
    {
        if (*number > 1)
        {
            return 1;
        }
        else
        {
            return abs(*number);
        }
    }
};

class Window
{
public:
    static void Window_Init(double size_X, double size_Y, COLORREF background, bool coursor)
    {
        txCreateWindow(size_X, size_Y, FALSE);
        txSetFillColor(background);
        txClear();
        txTextCursor(coursor);
    }
};

class Figure
{
public:
    static void Figure_Main()
    {
        Window::Window_Init(WINDOW_WIDTH, WINDOW_HEIGHT, WINDOW_BACKGROUND, BLINK_COURSOR);

        double x1 = WINDOW_WIDTH / 2;
        double y1 = WINDOW_HEIGHT / 2;

        double COORD_X1 = -500;
        double COORD_Y1 = 0;
        double COORD_Z1 = 0;

        double COORD_X2 = 80;
        double COORD_Y2 = -200;

        double z;

        for (double angle = 0; !txGetAsyncKeyState(VK_SPACE); angle += 0.0174)
        {


            for (int x = -RADIUS_OF_SPHERE; x <= RADIUS_OF_SPHERE; x++)
            {
                for (int y = -RADIUS_OF_SPHERE; y <= RADIUS_OF_SPHERE; y++)
                {
                    if (x * x + y * y <= RADIUS_OF_SPHERE * RADIUS_OF_SPHERE)
                    {
                        z = Sphere_Parameters::Depth(x, y, RADIUS_OF_SPHERE);

                        Vector radius_of_sphere(x, y, z);
                        Vector view_point(-x, -y, z);
                        Vector vector_source1(COORD_X1, COORD_Y1, z);
                        Vector vector_source2(COORD_X2, COORD_Y2, z);

                        Source source1(Sphere_Parameters::Intensity(&radius_of_sphere, &vector_source1),
                            Sphere_Parameters::Scalar_Flare(&view_point, &radius_of_sphere, &vector_source1));

                        Source source2(Sphere_Parameters::Intensity(&radius_of_sphere, &vector_source2),
                            Sphere_Parameters::Scalar_Flare(&view_point, &radius_of_sphere, &vector_source2));

                        Vector color_sum((source1.fong + source1.lambert + source2.fong + source2.lambert) * (R_SOURCE_LIGHT_COLOR + R_COLOR),
                            (source1.fong + source1.lambert + source2.fong + source2.lambert) * (G_SOURCE_LIGHT_COLOR + G_COLOR),
                            (source1.fong + source1.lambert + source2.fong + source2.lambert) * (B_SOURCE_LIGHT_COLOR + B_COLOR));

                        if (source1.lambert + source2.lambert > 0)
                        {
                            color_sum = Normalize::Normalize_Check(&color_sum);
                            txSetPixel(x + x1, y + y1, RGB(z * color_sum.x / 0.6, z * color_sum.y / 0.6, z * color_sum.z / 0.6));
                        }
                        else
                        {
                            txSetPixel(x + x1, y + y1, RGB(0, 0, 0));
                        }
                    }
                }
            }
            COORD_X1 = -500 * cos(angle * 5);
        }
    }
};

int main()
{
    Video_memory = txVideoMemory();
    txBegin();
    Figure::Figure_Main();
    txEnd();
    return 0;
}