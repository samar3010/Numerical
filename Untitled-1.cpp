#include <iostream>
#include <cmath>
using namespace std;
const double eps = 0.2;
double f(double x)
{
return -13 - (20 * x) + (19 * pow(x, 2)) - (3 * pow(x, 3));
}
double falsePosition(double xl, double xu)
{
double iter = 0;
double xr = 0;
double xrOld = 0;
double error = 0;
do
{
xrOld = xr;
xr = xu - (f(xu) * (xl - xu)) / (f(xl) - f(xu));
error = abs((xr - xrOld) / xr) * 100;
cout << "iteration=" << iter << " | xl=" << xl << " | f(xl)=" << f(xl)
<< " | xu=" << xu << " | f(xu)"
<< f(xu) << " | xr=" << xr << " | f(xr)=" << f(xr) << " Error%="
<< error << endl << endl;
if (f(xl) * f(xr) > 0)
{
xl = xr;
}
else
{
xu = xr;
}
iter++;
} while (error > eps);
return xr;
}
int main(){
float xl, xu, root;
cout << "Enter lower value = ";
cin >> xl;
cout << "Enter upper value = ";
cin >> xu;
if (f(xl) * f(xu) < 0)
{
root = falsePosition(xl, xu);
cout << "root= " << root;
}
else
{
cout << "Not correct xl and xu" << endl;
}
return 0;
}