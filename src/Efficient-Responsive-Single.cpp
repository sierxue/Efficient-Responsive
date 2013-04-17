//============================================================================
// Name        :    Efficient-Responsive-Single.cpp
// Author      :    Tong WANG
// Email       :    tong.wang@nus.edu.sg
// Version     :    v2.0 (originall written on 2009-05-05, re-organized on 2013-04-17)
// Copyright   :    ...
// Description :    code for the Efficient-Responsive paper
//                  For Responsiveness paper Case 1: single order at either epoch 1 (efficient) or epoch 2 (responsive)
//                  --binary choice between E and R
//                  --for variable cost differential \delta and/or fixed cost differential K
//============================================================================

#include <iostream>
#include <fstream>
#include <iomanip> //required by setprecision()
#include <cmath>
#include <vector>

using namespace std;

//---------------------------------------------------------------------------

#define  STEP 100                   //discrete each standard deviation into STEP=100 segments
#define  N 4                        //only consider +/- 4 standard deviations for normal distributions

//---------------------------------------------------------------------------

double  b;                          //substitutability
double  c_E, c_R, delta, K;         //unit cost for E and R, variable cost differential, fixed cost of responsiveness;
double  u, var, stdev;              //mean demand, variance, stdev of signal epsilon
vector<double>  phi(2*N*STEP+1);    //standard Normal p.d.f.

ofstream file, file2;               //output files


//============================================================================

//first order condition of the Efficient firm for the E-R case
double foc1_ER(double qq)
{
    double out = (2-b)*u + b*c_R - 2*c_E - 2*(2-b*b)*qq;
    double intg = 0;

    for (int j=0;j<=2*N*STEP;j++)
    {
        double ee = ((double)j/STEP-N)*stdev;
        if (ee + u - c_R - b*qq<0)
            intg += (ee + u - c_R - 2*b*qq) * phi[j] ;
    }

    return (out + b*intg/STEP)/2;
}

//============================================================================


int main()
{

    //Open output file
    //file.open("Efficient-Responsive-Single.txt");

    //if (! file)
    //{
        //if fail to open the file
    //    cerr << "can't open output file Efficient-Responsive-Single.txt!" << endl;
    //    exit(EXIT_FAILURE);
    //}

    file2.open("Efficient-Responsive-Single-Nash.txt");

    if (! file2)
    {
        //if fail to open the file
        cerr << "can't open output file Efficient-Responsive-Single-Nash.txt!" << endl;
        exit(EXIT_FAILURE);
    }

    
    cout << setprecision(10);
    //file << setprecision(10);
    file2 << setprecision(5);



    cout << "u\tva\tb\tc_E\tc_R\tK\tq_EE\tPi_EE\tCS_EE\tTS_EE\tq_RR\tPi_RR\tCS_RR\tTS_RR\tqE_ER\tqR_ER\tPiE_ER\tPiR_ER\tCS_ER\tTS_ER\tNashEq\tCS*\tTS*" << endl;
    //file << "u\tvar\tb\tc_E\tc_R\tK\tq_EE\tPi_EE\tCS_EE\tTS_EE\tq_RR\tPi_RR\tCS_RR\tTS_RR\tqE_ER\tqR_ER\tPiE_ER\tPiR_ER\tCS_ER\tTS_ER\tNashEq\tCS*\tTS*" << endl;

    file2 << "K\tdelta\tb\tNashE\tNashR\tCS01\tCS12\tTS01\tTS12\tCoop01\tCoop12" << endl;

    //initialize array for standard normal pdf
    for (int i=0; i<=2*N*STEP; i++)
        phi[i] = (1/sqrt(2.0*M_PI))*exp(-pow(((double)i/STEP-N),2.0)/2);

    //initialize default values of the parameters
    u  = 10;
    var = 16;

    delta = 1;
    c_E = 1;
    c_R = c_E*delta;
    K = 0;

    b = 1;


    //initialize the vectors to save threshold values
    const int II=1000;
    const int JJ=100;
    const int KK=5;

    //thresholds are such that if the uncertainty $var$ exceeds the threshold, the corresponding optimal strategy switches to R from E
    vector< vector<double> > NashE(KK+1, vector<double>(JJ+1));     //threshold for the optimal Nash strategy if the other player chooses E
    vector< vector<double> > NashR(KK+1, vector<double>(JJ+1));     //threshold for the optimal Nash strategy if the other player chooses R
    vector< vector<double> > CS01(KK+1, vector<double>(JJ+1));
    vector< vector<double> > CS12(KK+1, vector<double>(JJ+1));
    vector< vector<double> > TS01(KK+1, vector<double>(JJ+1));
    vector< vector<double> > TS12(KK+1, vector<double>(JJ+1));
    vector< vector<double> > Coop01(KK+1, vector<double>(JJ+1));
    vector< vector<double> > Coop12(KK+1, vector<double>(JJ+1));


    for (int kk=0; kk<=KK; kk++)            //for-loop for delta or K
    {
        //K = kk;                 //iterate over multiple values for K
        
        delta = 1+0.1*kk;     //iterate over multiple values for delta
        c_R = c_E*delta;

        for (int jj=0; jj<=JJ; jj++)        //for b
        {
            b = jj*0.01;

            for (int ii=0; ii<=II; ii++)    //for var
            {
                var = ii*0.02;


                //file << u << "\t" << var << "\t"  << b << "\t"  << c_E << "\t"  << c_R << "\t"  << K << "\t" ;
                cout << u << "\t" << var << "\t"  << b << "\t"  << c_E << "\t"  << c_R << "\t"  << K << "\t" ;

                //stdev of the signal (info available to the Responsive firm)
                stdev = sqrt(var);


                //for case E-E
                double	q_EE = (u-c_E)/(2+b);
                double	Pi_EE = pow((u-c_E)/(2+b), 2.0);
                double	CS_EE = (1+b)*Pi_EE;
                double	TS_EE = (3+b)*Pi_EE;

                //file << q_EE << "\t"  << Pi_EE << "\t"  << CS_EE << "\t"  << TS_EE << "\t";
                cout << q_EE << "\t"  << Pi_EE << "\t"  << CS_EE << "\t"  << TS_EE << "\t";


                //for case R-R
                double	q_RR=0, Pi_RR=0, CS_RR, TS_RR;

                for (int i=0;i<=2*N*STEP;i++)
                {
                    double e = ((double)i/STEP-N)*stdev;
                    if (u + e - c_R>0)
                    {
                        q_RR += (u+e-c_R)/(2+b) * phi[i];
                        Pi_RR += pow((u+e-c_R)/(2+b), 2.0) * phi[i];
                    }
                }

                q_RR /= STEP;
                Pi_RR /= STEP;
                CS_RR = (1+b)*Pi_RR;
                Pi_RR -= K;       //account for fixed cost
                TS_RR = CS_RR + 2*Pi_RR;

                //file << q_RR << "\t"  << Pi_RR << "\t"  << CS_RR << "\t"  << TS_RR << "\t";
                cout << q_RR << "\t"  << Pi_RR << "\t"  << CS_RR << "\t"  << TS_RR << "\t";



                //for case E-R (firm 1 E, firm 2 R)
                double	qR_ER=0, qE_ER=0, PiR_ER=0, PiE_ER=0, CS_ER=0, TS_ER=0;

                //search for optimal q1 until foc1=0
                double q1=0, q2=0;
                double lb=0, ub=(2-b)*(u-c_E)/(2*(2-b*b));      //lower and upper bounds used in the searching algorithm

                if ((ub>lb)&&(foc1_ER(q1)>0))
                {
                    q1=(lb+ub)/2;
                    double temp=foc1_ER(q1);
                    for (;ub-lb>0.001;)
                    {
                        if (temp>0)
                            lb = q1;
                        else
                            ub = q1;

                        q1 = (lb+ub)/2;
                        temp = foc1_ER(q1);
                    }
                }


                //for given q1, calculating expected q2, Pi1, Pi2 with regard to the second signal e
                for (int j=0; j<=2*N*STEP; j++)
                {
                    double e = ((double)j/STEP-N)*stdev;        //realization of the second signal e

                    if (e + u - c_R - b*q1<0)
                        q2 = 0;
                    else
                        q2 = (e + u - c_R - b*q1)/2;

                    qR_ER += q2 * phi[j];
                    PiE_ER += (e + u - q1 - b*q2 - c_E)*q1 * phi[j] ;
                    PiR_ER += (e + u - q2 - b*q1 - c_R)*q2 * phi[j];
                    CS_ER += 0.5*( pow(q1,2.0) + pow(q2,2.0) + 2*b*q1*q2 ) * phi[j];
                }


                qE_ER = q1;
                qR_ER /= STEP;
                PiE_ER /= STEP;
                PiR_ER /= STEP;
                CS_ER /= STEP;
                PiR_ER -= K;       //account for fixed cost
                TS_ER = CS_ER + PiR_ER + PiE_ER;

                //file << qE_ER << "\t"  << qR_ER << "\t"  << PiE_ER << "\t"  << PiR_ER << "\t"  << CS_ER << "\t"  << TS_ER << "\t";
                cout << qE_ER << "\t"  << qR_ER << "\t"  << PiE_ER << "\t"  << PiR_ER << "\t"  << CS_ER << "\t"  << TS_ER << "\t";



                //Nash equilibrium analysis
                if ((Pi_EE>=PiR_ER)&&(PiE_ER>=Pi_RR))
                {
                    cout << 0 << "\t";          //E-E
                    //file << 0 << "\t";          //E-E
                }
                else if ((Pi_EE<PiR_ER)&&(PiE_ER<Pi_RR))
                {
                    cout << 2 << "\t";          //R-R
                    //file << 2 << "\t";          //R-R
                }
                else if ((Pi_EE>=PiR_ER)&&(PiE_ER<Pi_RR))
                {
                    cout << 0.5 << "\t";        //E-E / R-R
                    //file << 0.5 << "\t";        //E-E / R-R
                }
                else if ((Pi_EE<PiR_ER)&&(PiE_ER>=Pi_RR))
                {
                    cout << 1 << "\t";          //E-R
                    //file << 1 << "\t";          //E-R
                }
                else
                {
                    cout << "???\t";            //error
                    file << "???\t";            //error
                }

                //best strategy for maximizing Consumer Surplus
                if ((CS_EE>=CS_ER)&&(CS_EE>=CS_RR))
                {
                    cout << 0 << "\t";          //E-E
                    //file << 0 << "\t";          //E-E
                }
                else if ((CS_ER>=CS_EE)&&(CS_ER>=CS_RR))
                {
                    cout << 1 << "\t";          //E-R
                    //file << 1 << "\t";          //E-R
                }
                else
                {
                    cout << 2 << "\t";          //R-R
                    //file << 2 << "\t";          //R-R
                }

                //best strategy for maximizing Total Surplus
                if ((TS_EE>=TS_ER)&&(TS_EE>=TS_RR))
                {
                    cout << 0 << endl;          //E-E
                    //file << 0 << endl;          //E-E
                }
                else if ((TS_ER>=TS_EE)&&(TS_ER>=TS_RR))
                {
                    cout << 1 << endl;          //E-R
                    //file << 1 << endl;          //E-R
                }
                else
                {
                    cout << 2 << endl;          //R-R
                    //file << 2 << endl;          //R-R
                }




                ////////////////////////////////////////////////////////////////////////
                //keep track of the thresholds on $var$
                
                if (Pi_EE>=PiR_ER)
                    NashE[kk][jj] = var;

                if (PiE_ER>=Pi_RR)
                    NashR[kk][jj] = var;

                if ((CS_EE>=CS_ER)&&(CS_EE>=CS_RR))
                    CS01[kk][jj] = var;

                if ((CS_EE>=CS_RR)||(CS_ER>=CS_RR))
                    CS12[kk][jj] = var;

                if ((TS_EE>=TS_ER)&&(TS_EE>=TS_RR))
                    TS01[kk][jj] = var;

                if ((TS_EE>=TS_RR)||(TS_ER>=TS_RR))
                    TS12[kk][jj] = var;

                if ((Pi_EE*2>=PiE_ER+PiR_ER)&&(Pi_EE>=Pi_RR))
                    Coop01[kk][jj] = var;

                if ((Pi_EE>=Pi_RR)||(PiE_ER+PiR_ER>=2*Pi_RR))
                    Coop12[kk][jj] = var;


            }   //end-for var

        }   //end-for b

        //output the thresholds to file2
        for (int jj=0; jj<=JJ; jj++)
            file2 << K << "\t" << delta << "\t" << jj*0.01 << "\t" << NashE[kk][jj]  << "\t" << NashR[kk][jj]  << "\t" << CS01[kk][jj]  << "\t" << CS12[kk][jj]  << "\t" << TS01[kk][jj]  << "\t" << TS12[kk][jj]  << "\t" << Coop01[kk][jj]  << "\t" << Coop12[kk][jj] << endl;

    }   //end big parameter loop



    //file.close();
    file2.close();
    
    return 0;

}
