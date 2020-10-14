#include<iostream>
#include<iomanip>
#include<fstream>
#include<cmath>
#include<ctime>
#include<boost/array.hpp>
#include<boost/numeric/odeint.hpp>

#define COST_FUNCTION_PRECISION 7

using namespace std;
using namespace boost::numeric::odeint;
typedef boost::array<double, 8> stateType;

vector<double> rFunc; //orbit r(t) for each pair of initital conditions

//gravitational constant, masses and mass parameter
double G,M1,M2,m;
//oblate spheroid axes and moments of inertia
double a1,b1,c1, I1x,I1y,I1z;
//ellipsoid axes and moments of inertia
double a2,b2,c2, I2x,I2y,I2z;

//load the inputs into a vector<double>
vector<double> load_inputs(const char *path)
{
    ifstream fpInputs(path);
    if (!fpInputs.is_open())
    {
        cout << "Inputs file not found. Exiting...\n";
        exit(EXIT_FAILURE);
    }
    double temp;
    vector<double> vec;
    while (fpInputs >> temp)
        vec.push_back(temp);
    fpInputs.close();
    return vec;
}

////////////////////////////////////////////////////////////////////////////////

//potential's r-derivative
double dV_dr(const stateType &x)
{
    double trace1 = (I1x + I1y + I1z)/M1;
    double trace2 = (I2x + I2y + I2z)/M2;
    return (G*M1*M2/(x[0]*x[0]))*(1.0 + (3.0/(2.0*x[0]*x[0]))*( trace1 + trace2 -
           (3.0/2.0)*( (I1x + I1y - cos(2.0*x[2])*(I1y - I1x))/M1 +
                       (I2x + I2y - cos(2.0*x[3])*(I2y - I2x))/M2 ) ) );
}

//potential's phi1-derivative
double dV_dphi1(const stateType &x)
{
    return 3.0*G*M2*sin(2.0*x[2])*(I1y - I1x)/(2.0*x[0]*x[0]*x[0]);
}

//potential's phi2-derivate
double dV_dphi2(const stateType &x)
{
    return 3.0*G*M1*sin(2.0*x[3])*(I2y - I2x)/(2.0*x[0]*x[0]*x[0]);
}

//8 odes in total
void odes(const stateType &x, stateType &dxdt, double t)
{
    //calculate the potential partial derivatives first
    double coeff1 = dV_dr(x);
    double coeff2 = dV_dphi1(x);
    double coeff3 = dV_dphi2(x);

    dxdt[0] = x[4];
    dxdt[1] = x[5];
    dxdt[2] = x[6];
    dxdt[3] = x[7];
    dxdt[4] = x[0]*x[5]*x[5] - coeff1/m;
    dxdt[5] = (coeff2 + coeff3)/(m*x[0]*x[0]) - 2.0*x[4]*x[5]/x[0];
    dxdt[6] = -(1.0 + m*x[0]*x[0]/I1z)*coeff2/(m*x[0]*x[0]) - coeff3/(m*x[0]*x[0]) + 2.0*x[4]*x[5]/x[0];
    dxdt[7] = -(1.0 + m*x[0]*x[0]/I2z)*coeff3/(m*x[0]*x[0]) - coeff2/(m*x[0]*x[0]) + 2.0*x[4]*x[5]/x[0];
}

//observe function (used to communicate with the odes function at each time step)
void observe(const stateType &x, const double t)
{
    //approximate collision detection for the 2 bodies, assuming a circumscribed
    //sphere for each body, the radius of which is equal to the longest semi axis
    //a (a1 for body 1 and a2 for body 2)
    if (x[0] <= a1 + a2)
    {
        cout << " Collision detection. Redefine the initial condition grid (r,vthita)\n";
        exit(EXIT_FAILURE);
    }
    rFunc.push_back(x[0]);
}

//Simpson's rule for computing the integral of the function (r(t) - r0)/r0
double simpson(const vector<double> &rFunc, double dt, double r0)
{
    double I0 = (rFunc[0] - r0)/r0;
    double I1 = (rFunc[rFunc.size() - 1] - r0)/r0;
    double sum = 0.0;
    for (int i = 1; i < rFunc.size() - 1; ++i)
    {
        if (i%2 != 0)
            sum += 4.0*(rFunc[i] - r0)/r0;
        else
            sum += 2.0*(rFunc[i] - r0)/r0;
    }
    return (dt/3.0)*(I0 + sum + I1); //integral value
}

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
    //inputs vector
    vector<double> inputs = load_inputs("../resources/sys_inputs.txt");

    //store the inputs into the corresponding variables
    G = inputs[0];
    M1 = inputs[1];
    M2 = inputs[2];
    m = M1*M2/(M1 + M2);

    a1 = inputs[3];
    b1 = inputs[4];
    c1 = inputs[5];
    I1x = M1*(b1*b1 + c1*c1)/5.0;
    I1y = M1*(a1*a1 + c1*c1)/5.0;
    I1z = M1*(a1*a1 + b1*b1)/5.0;

    a2 = inputs[6];
    b2 = inputs[7];
    c2 = inputs[8];
    I2x = M2*(b2*b2 + c2*c2)/5.0;
    I2y = M2*(a2*a2 + c2*c2)/5.0;
    I2z = M2*(a2*a2 + b2*b2)/5.0;

    //grid's r - bounds and step
    double r0 = inputs[22];
    double rmax = inputs[23];
    double dr = inputs[24];

    //grid's vthita - bounds and step
    double vthita0 = inputs[25];
    double vthitamax = inputs[26];
    double dvthita = inputs[27];

    //time parameters
    double t0 = inputs[19];
    double tmax = inputs[20];
    double dt = inputs[21];
    
    ofstream fpJ("../resources/cost_function_data.txt"); //output file
    fpJ << fixed << setprecision(COST_FUNCTION_PRECISION);
    
    cout << "Traversing the grid (r,vthita)...\n";
    for (double vthita = vthita0; vthita <= vthitamax + dvthita; vthita += dvthita)
    {
        for (double r = r0; r <= rmax + dr; r += dr)
        {
            printf("\r%.3lf / %.3lf, %.8lf / %.8lf",r,rmax,vthita,vthitamax);
            fflush(stdout);

            //r gets value from the for loop
            double thita = inputs[12];
            double phi1 = inputs[13];
            double phi2 = inputs[14];

            double vr = inputs[15];
            //vthita gets value from the for loop
            double vphi1 = inputs[17];
            double vphi2 = inputs[18];

            double rStart = r; //r(0)
            double vthitaStart = vthita; //vthita(0)

            stateType initConditions = {rStart,thita,phi1,phi2, vr,vthitaStart,vphi1,vphi2};
            runge_kutta_dopri5<stateType> stepper;
            integrate_const(stepper, odes, initConditions, t0,tmax,dt, observe);
            double J = fabs(simpson(rFunc,dt,r)); //cost function calculation
            rFunc.clear();
            fpJ << r << " " << vthita << " " << J << endl;
        }
    }
    cout << "\nCost function evaluated successfully.\n";
    fpJ.close();
    return 0;
}
