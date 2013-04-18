//============================================================================
// Name        :    Efficient-Responsive-Multi.cpp
// Author      :    Tong WANG
// Email       :    tong.wang@nus.edu.sg
// Version     :    v2.0 (originall written on 2009-05-05, re-organized on 2013-04-17)
// Copyright   :    ...
// Description :    code for the Efficient-Responsive paper
//                  For Responsiveness paper Case 2: possible multi order at both epoch 1 and epoch 2 for the Responsive firm
//                  --binary choice between E and R
//                  --variable cost differential \delta and fixed cost differential K are NOT yet implemented
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
double  c_E, c_R1, c_R2;            //unit cost for E and R;
double  delta, K;                   //variable cost differential, fixed cost of responsiveness;
double  u, var, stdev;              //mean demand, variance, stdev of signal epsilon
vector<double>  phi(2*N*STEP+1);    //standard Normal p.d.f.

ofstream file, file2;               //output files


//============================================================================

//first order condition of the Efficient firm for the E-R case
double foc1_ER(double qq)
{
    double out = (2-b)*u + b*c_R2 - 2*c_E - 2*(2-b*b)*qq;
    double intg = 0;

    for (int j=0;j<=2*N*STEP;j++)
    {
        double ee = ((double)j/STEP-N)*stdev;
        if (ee + u - c_R2 - b*qq<0)
            intg += (ee + u - c_R2 - 2*b*qq) * phi[j] ;
    }

    return (out + b*intg/STEP)/2;
}


double Pi1_RR(double qq1, double qq2)
{
    double intg = 0;

    if (qq1>qq2)
    {
        for (int j=0;j<=2*N*STEP;j++)
        {
            double ee = ((double)j/STEP-N)*stdev;
            if (ee + u - c_R2 <= b*qq1 + 2*qq2)
                intg += (ee + u - qq1 - b*qq2)*qq1 * phi[j] ;
            else if ((ee + u - c_R2 > b*qq1 + 2*qq2)&&(ee + u - c_R2 <= (2+b)*qq1))
                intg += (0.5*(2-b)*(ee+u) + 0.5*b*c_R2 - (1-b*b/2)*qq1)*qq1 * phi[j] ;
            else
                intg += ( pow((ee+u-c_R2)/(2+b),2.0) + c_R2*qq1 ) * phi[j] ;
        }
    }
    else
    {
        for (int j=0;j<=2*N*STEP;j++)
        {
            double ee = ((double)j/STEP-N)*stdev;
            if (ee + u - c_R2 <= b*qq2 + 2*qq1)
                intg += (ee + u - qq1 - b*qq2)*qq1 * phi[j] ;
            else if ((ee + u - c_R2 > b*qq2 + 2*qq1)&&(ee + u - c_R2 <= (2+b)*qq2))
                intg += (0.25*pow(ee+u-c_R2-b*qq2,2.0) + c_R2*qq1) * phi[j] ;
            else
                intg += ( pow((ee+u-c_R2)/(2+b),2.0) + c_R2*qq1 ) * phi[j] ;
        }
        
    }

    return intg/STEP - c_R1*qq1;
}


//============================================================================


int main()
{

    /*********************************************************
     //DO NOT OUTPUT DETAILED NUMBER IN THE PRODUCTION VERSION
     //Open output file
     file.open("Efficient-Responsive-Multi.txt");
     
     if (! file)
     {
     //if fail to open the file
     cerr << "can't open output file Efficient-Responsive-Multi.txt!" << endl;
     exit(EXIT_FAILURE);
     }
     *********************************************************/

    file2.open("Efficient-Responsive-Multi-Nash.txt");

    if (! file2)
    {
        //if fail to open the file
        cerr << "can't open output file Efficient-Responsive-Multi-Nash.txt!" << endl;
        exit(EXIT_FAILURE);
    }

    
    cout << setprecision(5);
    //file << setprecision(10);
    file2 << setprecision(5);



    //cout << "u\tvar\tb\tc_E\tc_R1\tc_R2\tK\tq_EE\tPi_EE\tCS_EE\tTS_EE\tq_RR\tPi_RR\tCS_RR\tTS_RR\tqE_ER\tqR_ER\tPiE_ER\tPiR_ER\tCS_ER\tTS_ER\tNashEq\tq1_RRm\tQ1_RRm\tQ2_RRm\tPi1_RRm\tPi2_RRm\tPi_RRm\tCS_RRm\tTS_RRm\tNashEqm\tCSm*\tTSm*" << endl;
    //file << "u\tvar\tb\tc_E\tc_R1\tc_R2\tK\tq_EE\tPi_EE\tCS_EE\tTS_EE\tq_RR\tPi_RR\tCS_RR\tTS_RR\tqE_ER\tqR_ER\tPiE_ER\tPiR_ER\tCS_ER\tTS_ER\tNashEq\tq1_RRm\tQ1_RRm\tQ2_RRm\tPi1_RRm\tPi2_RRm\tPi_RRm\tCS_RRm\tTS_RRm\tNashEqm\tCSm*\tTSm*" << endl;

    file2 << "K\tdelta\tb\tNashE\tNashR\tNashEm\tNashRm\tCS01m\tCS12m\tTS01m\tTS12m\tCoop01m\tCoop12m" << endl;

    //initialize array for standard normal pdf
    for (int i=0; i<=2*N*STEP; i++)
        phi[i] = (1/sqrt(2.0*M_PI))*exp(-pow(((double)i/STEP-N),2.0)/2);

    //initialize default values of the parameters
    u  = 10;
    var = 16;

    c_E = 1;
    c_R1 = 1;
    c_R2 = 1;
    
    delta = 1;  //THE CODE IS NOT READY FOR VARIABLE/FIXED COST DIFFERENTIAL
    K = 0;      //DO NOT CHANGE THESE VALUES!

    b = 1;


    //initialize the vectors to save threshold values
    const int II=1000;
    const int JJ=100;
    const int KK=0;

    //thresholds are such that if the uncertainty $var$ exceeds the threshold, the corresponding optimal strategy switches to R from E
    vector< vector<double> > NashE(KK+1, vector<double>(JJ+1));     //threshold for the optimal Nash strategy if the other player chooses E
    vector< vector<double> > NashR(KK+1, vector<double>(JJ+1));     //threshold for the optimal Nash strategy if the other player chooses R
    vector< vector<double> > NashEm(KK+1, vector<double>(JJ+1));
    vector< vector<double> > NashRm(KK+1, vector<double>(JJ+1));
    vector< vector<double> > CS01m(KK+1, vector<double>(JJ+1));
    vector< vector<double> > CS12m(KK+1, vector<double>(JJ+1));
    vector< vector<double> > TS01m(KK+1, vector<double>(JJ+1));
    vector< vector<double> > TS12m(KK+1, vector<double>(JJ+1));
    vector< vector<double> > Coop01m(KK+1, vector<double>(JJ+1));
    vector< vector<double> > Coop12m(KK+1, vector<double>(JJ+1));


    for (int kk=0; kk<=KK; kk++)            //empty for-loop left for future introduction of K or delta 
    {
        for (int jj=0; jj<=JJ; jj++)        //for b
        {
            b = jj*0.01;

            cout << "Calculating K=" << K << ", delta=" << delta << ", b=" << b << " ... " << endl;

            for (int ii=0; ii<=II; ii++)    //for var
            {
                var = ii*0.02;


                //file << u << "\t" << var << "\t"  << b << "\t"  << c_E << "\t"  << c_R1 << "\t"  << c_R2 << "\t"  << K << "\t" ;
                //cout << u << "\t" << var << "\t"  << b << "\t"  << c_E << "\t"  << c_R1 << "\t"  << c_R2 << "\t"  << K << "\t";

                //stdev of the signal (info available to the Responsive firm)
                stdev = sqrt(var);


                //for case E-E
                double	q_EE = (u-c_E)/(2+b);
                double	Pi_EE = pow((u-c_E)/(2+b), 2.0);
                double	CS_EE = (1+b)*Pi_EE;
                double	TS_EE = (3+b)*Pi_EE;

                //file << q_EE << "\t"  << Pi_EE << "\t"  << CS_EE << "\t"  << TS_EE << "\t";
                //cout << q_EE << "\t"  << Pi_EE << "\t"  << CS_EE << "\t"  << TS_EE << "\t";


                //for case R-R
                double	q_RR=0, Pi_RR=0, CS_RR, TS_RR;

                for (int i=0;i<=2*N*STEP;i++)
                {
                    double e = ((double)i/STEP-N)*stdev;
                    if (u + e - c_R2>0)
                    {
                        q_RR += (u+e-c_R2)/(2+b) * phi[i];
                        Pi_RR += pow((u+e-c_R2)/(2+b), 2.0) * phi[i];
                    }
                }

                q_RR /= STEP;
                Pi_RR /= STEP;
                CS_RR = (1+b)*Pi_RR;
                Pi_RR -= K;       //account for fixed cost
                TS_RR = CS_RR + 2*Pi_RR;

                //file << q_RR << "\t"  << Pi_RR << "\t"  << CS_RR << "\t"  << TS_RR << "\t";
                //cout << q_RR << "\t"  << Pi_RR << "\t"  << CS_RR << "\t"  << TS_RR << "\t";



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

                    if (e + u - c_R2 - b*q1<0)
                        q2 = 0;
                    else
                        q2 = (e + u - c_R2 - b*q1)/2;

                    qR_ER += q2 * phi[j];
                    PiE_ER += (e + u - q1 - b*q2 - c_E)*q1 * phi[j] ;
                    PiR_ER += (e + u - q2 - b*q1 - c_R2)*q2 * phi[j];
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
                //cout << qE_ER << "\t"  << qR_ER << "\t"  << PiE_ER << "\t"  << PiR_ER << "\t"  << CS_ER << "\t"  << TS_ER << "\t";


                /*********************************************************
                //DO NOT OUTPUT DETAILED NUMBER IN THE PRODUCTION VERSION

                //Nash equilibrium analysis
                if ((Pi_EE>=PiR_ER)&&(PiE_ER>=Pi_RR))
                {
                    cout << 0 << "\t";          //E-E
                    file << 0 << "\t";          //E-E
                }
                else if ((Pi_EE<PiR_ER)&&(PiE_ER<Pi_RR))
                {
                    cout << 2 << "\t";          //R-R
                    file << 2 << "\t";          //R-R
                }
                else if ((Pi_EE>=PiR_ER)&&(PiE_ER<Pi_RR))
                {
                    cout << 0.5 << "\t";        //E-E / R-R
                    file << 0.5 << "\t";        //E-E / R-R
                }
                else if ((Pi_EE<PiR_ER)&&(PiE_ER>=Pi_RR))
                {
                    cout << 1 << "\t";          //E-R
                    file << 1 << "\t";          //E-R
                }
                else
                {
                    cout << "???\t";            //error
                    file << "???\t";            //error
                }
                *********************************************************/


                
                
                
                
                //for case R-Rm
                double q1_RRm=0, Q1_RRm=0, Q2_RRm=0, Pi1_RRm=0, Pi2_RRm=0, Pi_RRm=0, CS_RRm=0, TS_RRm=0;
                {
                    
                    //search for q1 until foc1=0
                    double q1=0, q2=0;
                    double lb=0, ub=(2-b)*(u-c_E)/(2*(2-b*b)); //qE_ER  //lower and upper bounds used in the searching algorithm
                    
                    //line search q1
                    double qq1, max=-1.0e20;
                    for (qq1=lb; qq1<=ub; qq1+=0.001)
                    {
                        double temp = Pi1_RR(qq1,q2);
                        
                        if (temp>max)
                        {
                            q1 = qq1;
                            max = temp;
                        }
                    }
                    
                    
                    //for given q1, calculating expected Q1, Q2, Pi1, Pi2 with regard to the second signal e
                    for (int j=0; j<=2*N*STEP; j++)
                    {
                        double e = ((double)j/STEP-N)*stdev; //realization of the second signal e
                        double Q1,Q2;
                        
                        if (e + u - c_R2 <= b*q1 + 2*q2)
                        {
                            Q1 = q1;
                            Q2 = q2;
                        }
                        else if ((e + u - c_R2 > b*q1 + 2*q2)&&(e + u - c_R2 <= (2+b)*q1))
                        {
                            Q1 = q1;
                            Q2 = (e + u - c_R2 - b*q1)/2;
                        }
                        else
                        {
                            Q1 = (e + u - c_R2)/(2+b);
                            Q2 = (e + u - c_R2)/(2+b);
                        }
                        
                        
                        //take expectation again with regard to e
                        Q1_RRm += Q1 * phi[j];
                        Q2_RRm += Q2 * phi[j];
                        Pi1_RRm += ((e + u - Q1 - b*Q2)*Q1 - c_R1*q1 - c_R2*(Q1-q1))  * phi[j];
                        Pi2_RRm += ((e + u - Q2 - b*Q1)*Q2 - c_R1*q2 - c_R2*(Q2-q2))  * phi[j] ;
                        CS_RRm += 0.5*( pow(Q1,2.0) + pow(Q2,2.0) + 2*b*Q1*Q2 )  * phi[j];
                        
                    }
                    
                    
                    q1_RRm = q1;
                    Q1_RRm /= STEP;
                    Q2_RRm /= STEP;
                    Pi1_RRm /= STEP;
                    Pi2_RRm /= STEP;
                    CS_RRm /= STEP;
                    TS_RRm = CS_RRm + Pi1_RRm + Pi2_RRm;
                    
                    Pi_RRm = (Pi1_RRm+Pi2_RRm)/2;
                    
                    //file << q1_RRm << "\t" << Q1_RRm << "\t" << Q2_RRm << "\t" << Pi1_RRm << "\t" << Pi2_RRm << "\t" << Pi_RRm << "\t" << CS_RRm << "\t" << TS_RRm << "\t";
                    //cout << q1_RRm << "\t" << Q1_RRm << "\t" << Q2_RRm << "\t" << Pi1_RRm << "\t" << Pi2_RRm << "\t" << Pi_RRm << "\t" << CS_RRm << "\t" << TS_RRm << "\t";
                    
                }
                
                /*********************************************************
                //DO NOT OUTPUT DETAILED NUMBER IN THE PRODUCTION VERSION
                //Nash equilibrium analysis for multiple ordering
                if ((Pi_EE>=PiR_ER)&&(PiE_ER>=Pi_RRm))
                {
                    cout << 0 << "\t";          //E-E
                    file << 0 << "\t";          //E-E
                }
                else if ((Pi_EE<PiR_ER)&&(PiE_ER<Pi_RRm))
                {
                    cout << 2 << "\t";          //R-R
                    file << 2 << "\t";          //R-R
                }
                else if ((Pi_EE>=PiR_ER)&&(PiE_ER<Pi_RRm))
                {
                    cout << 0.5 << "\t";        //E-E / R-R
                    file << 0.5 << "\t";        //E-E / R-R
                }
                else if ((Pi_EE<PiR_ER)&&(PiE_ER>=Pi_RRm))
                {
                    cout << 1 << "\t";          //E-R
                    file << 1 << "\t";          //E-R
                }
                else
                {
                    cout << "???\t";            //error
                    file << "???\t";            //error
                }

                
                
                
                //best strategy for maximizing Consumer Surplus
                if ((CS_EE>=CS_ER)&&(CS_EE>=CS_RRm))
                {
                    cout << 0 << "\t";          //E-E
                    file << 0 << "\t";          //E-E
                }
                else if ((CS_ER>=CS_EE)&&(CS_ER>=CS_RRm))
                {
                    cout << 1 << "\t";          //E-R
                    file << 1 << "\t";          //E-R
                }
                else
                {
                    cout << 2 << "\t";          //R-R
                    file << 2 << "\t";          //R-R
                }

                //best strategy for maximizing Total Surplus
                if ((TS_EE>=TS_ER)&&(TS_EE>=TS_RRm))
                {
                    cout << 0 << endl;          //E-E
                    file << 0 << endl;          //E-E
                }
                else if ((TS_ER>=TS_EE)&&(TS_ER>=TS_RRm))
                {
                    cout << 1 << endl;          //E-R
                    file << 1 << endl;          //E-R
                }
                else
                {
                    cout << 2 << endl;          //R-R
                    file << 2 << endl;          //R-R
                }
                ***********************************************************/



                ////////////////////////////////////////////////////////////////////////
                //keep track of the thresholds on $var$
                
                if (Pi_EE>=PiR_ER)
                    NashE[kk][jj] = var;
                
                if (PiE_ER>=Pi_RR)
                    NashR[kk][jj] = var;

                if (Pi_EE>=PiR_ER)
                    NashEm[kk][jj] = var;

                if (PiE_ER>=Pi_RRm)
                    NashRm[kk][jj] = var;

                if ((CS_EE>=CS_ER)&&(CS_EE>=CS_RRm))
                    CS01m[kk][jj] = var;

                if ((CS_EE>=CS_RRm)||(CS_ER>=CS_RRm))
                    CS12m[kk][jj] = var;

                if ((TS_EE>=TS_ER)&&(TS_EE>=TS_RRm))
                    TS01m[kk][jj] = var;

                if ((TS_EE>=TS_RRm)||(TS_ER>=TS_RRm))
                    TS12m[kk][jj] = var;

                if ((Pi_EE*2>=PiE_ER+PiR_ER)&&(Pi_EE>=Pi_RRm))
                    Coop01m[kk][jj] = var;

                if ((Pi_EE>=Pi_RRm)||(PiE_ER+PiR_ER>=2*Pi_RRm))
                    Coop12m[kk][jj] = var;


            }   //end-for var

        }   //end-for b

        //output the thresholds to file2
        for (int jj=0; jj<=JJ; jj++)
            file2 << K << "\t" << delta << "\t" << jj*0.01 << "\t" << NashE[kk][jj]  << "\t" << NashR[kk][jj]  << "\t"  << NashEm[kk][jj]  << "\t" << NashRm[kk][jj]  << "\t" << CS01m[kk][jj]  << "\t" << CS12m[kk][jj]  << "\t" << TS01m[kk][jj]  << "\t" << TS12m[kk][jj]  << "\t" << Coop01m[kk][jj]  << "\t" << Coop12m[kk][jj] << endl;

    }   //end big parameter loop



    //file.close();
    file2.close();
    
    return 0;

}
