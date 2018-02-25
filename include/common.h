
#ifndef SCARLET_STUDY_COMMON_H
#define SCARLET_STUDY_COMMON_H

# include <stdlib.h>
# include <stdio.h>
# include <string.h>
# include <math.h>
# include <time.h>



#define RANDOM_SEED 42

/*
 * helpers:
 */
double fRand(double fMin, double fMax) {
    double f = (double) rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

void fillArrayWithRandomValues(double *array, long length) {
    for (long i = 0; i < length; i++) {
        array[i] = fRand(-1000.0, +1000.0);
    }
}


/*
 * internal methods:
 */




int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
double r8_abs ( double x );
double r8_max ( double x, double y );
double r8_min ( double x, double y );
void r8vec_copy ( int n, double a1[], double a2[] );
void timestamp ( );
double wtime ( );







/******************************************************************************/

int i4_max ( int i1, int i2 )

/******************************************************************************/
/*
  Purpose:

    I4_MAX returns the maximum of two I4's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 August 2006

  Author:

    John Burkardt

  Parameters:

    Input, int I1, I2, are two integers to be compared.

    Output, int I4_MAX, the larger of I1 and I2.
*/
{
    int value;

    if ( i2 < i1 )
    {
        value = i1;
    }
    else
    {
        value = i2;
    }
    return value;
}
/******************************************************************************/

int i4_min ( int i1, int i2 )

/******************************************************************************/
/*
  Purpose:

    I4_MIN returns the smaller of two I4's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 August 2006

  Author:

    John Burkardt

  Parameters:

    Input, int I1, I2, two integers to be compared.

    Output, int I4_MIN, the smaller of I1 and I2.
*/
{
    int value;

    if ( i1 < i2 )
    {
        value = i1;
    }
    else
    {
        value = i2;
    }
    return value;
}
/******************************************************************************/

double r8_abs ( double x )

/******************************************************************************/
/*
  Purpose:

    R8_ABS returns the absolute value of an R8.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, double X, the quantity whose absolute value is desired.

    Output, double R8_ABS, the absolute value of X.
*/
{
    double value;

    if ( 0.0 <= x )
    {
        value = + x;
    }
    else
    {
        value = - x;
    }
    return value;
}
/******************************************************************************/

double r8_max ( double x, double y )

/******************************************************************************/
/*
  Purpose:

    R8_MAX returns the maximum of two R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, double X, Y, the quantities to compare.

    Output, double R8_MAX, the maximum of X and Y.
*/
{
    double value;

    if ( y < x )
    {
        value = x;
    }
    else
    {
        value = y;
    }
    return value;
}
/******************************************************************************/

double r8_min ( double x, double y )

/******************************************************************************/
/*
  Purpose:

    R8_MIN returns the minimum of two R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, double X, Y, the quantities to compare.

    Output, double R8_MIN, the minimum of X and Y.
*/
{
    double value;

    if ( y < x )
    {
        value = y;
    }
    else
    {
        value = x;
    }
    return value;
}
/******************************************************************************/

void r8vec_copy ( int n, double a1[], double a2[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_COPY copies an R8VEC.

  Discussion:

    An R8VEC is a vector of R8"s.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 July 2005

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vectors.

    Input, double A1[N], the vector to be copied.

    Input, double A2[N], the copy of A1.
*/
{
    int i;

    for ( i = 0; i < n; i++ )
    {
        a2[i] = a1[i];
    }
    return;
}
/******************************************************************************/

void timestamp ( void )

/******************************************************************************/
/*
  Purpose:

    TIMESTAMP prints the current YMDHMS date as a time stamp.

  Example:

    31 May 2001 09:45:54 AM

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 September 2003

  Author:

    John Burkardt

  Parameters:

    None
*/
{
# define TIME_SIZE 40

    static char time_buffer[TIME_SIZE];
    const struct tm *tm;
    size_t len;
    time_t now;

    now = time ( NULL );
    tm = localtime ( &now );

    len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

    fprintf ( stdout, "%s\n", time_buffer );

    return;
# undef TIME_SIZE
}
/******************************************************************************/

double wtime ( )

/******************************************************************************/
/*
  Purpose:

    WTIME reports the elapsed wallclock time.

  Discussion:

    The reliability of this function depends in part on the value of
    CLOCKS_PER_SECOND.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    27 April 2009

  Author:

    John Burkardt

  Parameters:

    Output, double WTIME, the a reading of the wall clock timer,
    in seconds.
*/
{
    double value;

    value = ( double ) clock ( )
            / ( double ) CLOCKS_PER_SEC;

    return value;
}


#endif //SCARLET_STUDY_COMMON_H
