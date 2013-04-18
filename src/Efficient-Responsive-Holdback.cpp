//============================================================================
// Name        :    Efficient-Responsive-Holdback.cpp
// Author      :    Tong WANG
// Email       :    tong.wang@nus.edu.sg
// Version     :    v2.0 (originall written on 2009-05-05, re-organized on 2013-04-17)
// Copyright   :    ...
// Description :    code for the Efficient-Responsive paper
//                  For Responsiveness paper Case 1: single order at either epoch 1 (efficient) or epoch 2 (responsive)
//                  --with Hold-back (efficient firm can hold back and salvage some committed quantity)
//                  --variable cost differential delta and fixed cost differential K are NOT yet implemented
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
double	s;                          //salvage value
double  u, var, stdev;              //mean demand, variance, stdev of signal epsilon
vector<double>  phi(2*N*STEP+1);    //standard Normal p.d.f.

ofstream file, file2;               //output files


//============================================================================

//first order condition of the firm 1 for the E-Eh case
double foc1_EEh(double qq)
{
	double intg = 0;
    
	for (int j=0;j<=2*N*STEP;j++)
    {
        double ee = ((double)j/STEP-N)*stdev;
        if (ee + u - s - (2+b)*qq>=0)
            intg += (ee + u - s - (2+b)*qq) * phi[j] ;
    }
    
	return s-c_E + intg/STEP;
}

//first order condition of the Efficient firm for the R-Eh case
//Case I:  Q1 <= (c2-s)/(2-b)
double foc1_ERh_I(double qq1)
{
	double intg = 0;
    
	for (int j=0;j<=2*N*STEP;j++)
	{
		double ee = ((double)j/STEP-N)*stdev;
		if ( (ee+u>=2*qq1+s)&&(ee+u<b*qq1+c_R))
			intg += (ee + u - s - 2*qq1) * phi[j] ;
		else if ( ee+u >= b*qq1+c_R)
			intg += (0.5*(2-b)*(ee + u) - s + 0.5*b*c_R - (2-b*b)*qq1) * phi[j] ;
	}
    
	return s-c_E + intg/STEP;
}

//Case II:  Q1 > (c2-s)/(2-b)
double foc1_ERh_II(double qq1)
{
	double intg = 0;
    
	for (int j=0;j<=2*N*STEP;j++)
	{
		double ee = ((double)j/STEP-N)*stdev;
		if (ee + u >= (2+b)*qq1 + (2*s-b*c_R)/(2-b))
			intg += (0.5*(2-b)*(ee + u) - s + 0.5*b*c_R - (2-b*b)*qq1) * phi[j] ;
	}
    
	return s-c_E + intg/STEP;
}

//============================================================================


int main()
{
    /*********************************************************
    //DO NOT OUTPUT DETAILED NUMBER IN THE PRODUCTION VERSION
    //Open output file
    file.open("Efficient-Responsive-Holdback.txt");

    if (! file)
    {
        //if fail to open the file
        cerr << "can't open output file Efficient-Responsive-Holdback.txt!" << endl;
        exit(EXIT_FAILURE);
    }
    *********************************************************/

    file2.open("Efficient-Responsive-Holdback-Nash.txt");

    if (! file2)
    {
        //if fail to open the file
        cerr << "can't open output file Efficient-Responsive-Holdback-Nash.txt!" << endl;
        exit(EXIT_FAILURE);
    }

    
    cout << setprecision(5);
    //file << setprecision(10);
    file2 << setprecision(5);



    //cout << "u\tvar\tb\tc_E\tc_R\tK\ts\tQ_EE\tq_EE\tPi_EE\tCS_EE\tTS_EE\tq_RR\tPi_RR\tCS_RR\tTS_RR\tQ1_ER\tqE_ER\tqR_ER\tPiE_ER\tPiR_ER\tCS_ER\tTS_ER\tNashEq\tCS*\tTS*" << endl;
    //file << "u\tvar\tb\tc_E\tc_R\tK\ts\tQ_EE\tq_EE\tPi_EE\tCS_EE\tTS_EE\tq_RR\tPi_RR\tCS_RR\tTS_RR\tQ1_ER\tqE_ER\tqR_ER\tPiE_ER\tPiR_ER\tCS_ER\tTS_ER\tNashEq\tCS*\tTS*" << endl;

    file2 << "K\tdelta\ts\tb\tNashEh\tNashRh\tCS01h\tCS12h\tTS01h\tTS12h\tCoop01h\tCoop12h" << endl;

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
    s = -10;

    b = 1;

    
    vector<double> s_array = {0.9, 0.5, 0, -0.5, -1, -1.5, -2, -2.5, -3, -3.5, -4, -5, -10, -20} ;

    //initialize the vectors to save threshold values
    const int II = 1000;
    const int JJ = 100;
    const int KK = s_array.size();

    //thresholds are such that if the uncertainty $var$ exceeds the threshold, the corresponding optimal strategy switches to R from E
    vector< vector<double> > NashEh(KK+1, vector<double>(JJ+1));     //threshold for the optimal Nash strategy if the other player chooses E
    vector< vector<double> > NashRh(KK+1, vector<double>(JJ+1));     //threshold for the optimal Nash strategy if the other player chooses R
    vector< vector<double> > CS01h(KK+1, vector<double>(JJ+1));
    vector< vector<double> > CS12h(KK+1, vector<double>(JJ+1));
    vector< vector<double> > TS01h(KK+1, vector<double>(JJ+1));
    vector< vector<double> > TS12h(KK+1, vector<double>(JJ+1));
    vector< vector<double> > Coop01h(KK+1, vector<double>(JJ+1));
    vector< vector<double> > Coop12h(KK+1, vector<double>(JJ+1));



    for (int kk=0; kk<KK; kk++)            //for-loop for delta or K
    {
        s = s_array[kk];
        
        //K = kk;                 //iterate over multiple values for K
        
        //delta = 1+0.1*kk;     //iterate over multiple values for delta
        //c_R = c_E*delta;
        
        for (int jj=0; jj<=JJ; jj++)        //for b
        {
            b = jj*0.01;

            cout << "Calculating K=" << K << ", delta=" << delta << ", s=" << s << ", b=" << b << " ... " << endl;

            for (int ii=0; ii<=II; ii++)    //for var
            {
                var = ii*0.02;


                //file << u << "\t" << var << "\t"  << b << "\t"  << c_E << "\t"  << c_R << "\t"  << K << "\t"  << s << "\t";
                //cout << u << "\t" << var << "\t"  << b << "\t"  << c_E << "\t"  << c_R << "\t"  << K << "\t"  << s << "\t";

                //stdev of the signal (info available to the Responsive firm)
                stdev = sqrt(var);


				////////////////////////////////////////////////////////////////////////
				//for case E-Eh
                
				double q_EE=0, Q_EE, Pi_EE=0, CS_EE=0, TS_EE;

                //search for Q
                double lb=0, ub=u;	//lower and upper bounds used in the searching algorithm
                double Q=lb,q;
                
                if ((ub>lb)&&(foc1_EEh(Q)>0))
                {
                    Q=(lb+ub)/2;
                    double temp=foc1_EEh(Q);
                    
                    for (;ub-lb>0.001;)
                    {
                        if (temp>0)
                            lb = Q;
                        else
                            ub = Q;
                        
                        Q = (lb+ub)/2;
                        temp = foc1_EEh(Q);
                    }
                }
                
                
                //for given Q, calculating expected q1, q2, Pi1, Pi2 with regard to the second signal e
                for (int j=0;j<=2*N*STEP;j++)
                {
                    double e = ((double)j/STEP-N)*stdev; //realization of the second signal e
                    
                    if (e + u < s)
                        q = 0;
                    else if (e + u < (2+b)*Q+s )
                        q = (e+u-s)/(2+b);
                    else
                        q = Q;
                    
                    
                    q_EE += q * phi[j];
                    Pi_EE += ((e + u - q - b*q - s)*q + (s-c_E)*Q )* phi[j];
                    CS_EE += (1+b)*q*q * phi[j];
                }
                
                
                Q_EE=Q;
                q_EE /= STEP;
                Pi_EE /= STEP;
                CS_EE /= STEP;
                TS_EE = 2*Pi_EE+CS_EE;
                
                //file << Q_EE << "\t"  << q_EE << "\t"  << Pi_EE << "\t"  << CS_EE << "\t"  << TS_EE << "\t";
                //cout << Q_EE << "\t"  << q_EE << "\t"  << Pi_EE << "\t"  << CS_EE << "\t"  << TS_EE << "\t";
 


                ///////////////////////////////////////////////////////////////////
                //for case R-Rh (the same as R-R)
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
                TS_RR = (3+b)*Pi_RR;
                
                //file << q_RR << "\t"  << Pi_RR << "\t"  << CS_RR << "\t"  << TS_RR << "\t";
                //cout << q_RR << "\t"  << Pi_RR << "\t"  << CS_RR << "\t"  << TS_RR << "\t";


                
                ///////////////////////////////////////////////////////////////////
                //for case E-Rh (firm 1 E, firm 2 R)
                double	Q1_ER,qR_ER,qE_ER,PiR_ER,PiE_ER,CS_ER,TS_ER;
                
                //Case I:  Q1 <= (c2-s)/(2-b)
                double Q1L=0,q1L=0,q2L=0;
                double qR_L=0,qE_L=0,PiR_L=0,PiE_L=0,CS_L=0;

                //search for Q1L
                double lb2=0, ub2=(c_R-s)/(2-b);	//lower and upper bounds used in the searching algorithm
                Q1L=lb2;
                
                if ((ub2>lb2)&&(foc1_ERh_I(Q1L)>0))
                {
                    Q1L=(lb2+ub2)/2;
                    double temp=foc1_ERh_I(Q1L);
                    
                    for (;ub2-lb2>0.001;)
                    {
                        if (temp>0)
                            lb2 = Q1L;
                        else
                            ub2 = Q1L;
                        
                        Q1L = (lb2+ub2)/2;
                        temp = foc1_ERh_I(Q1L);
                    }
                }
                
                
                //for given Q1, calculating expected q1, q2, Pi1, Pi2 with regard to the second signal e
                for (int j=0;j<=2*N*STEP;j++)
                {
                    double e = ((double)j/STEP-N)*stdev; //realization of the second signal e
                    
                    if (e + u < s)
                    {
                        q1L = 0;
                        q2L = 0;
                    }
                    else if (e + u < 2*Q1L+s )
                    {
                        q1L = (e+u-s)/2;
                        q2L = 0;
                    }
                    else if (e + u < b*Q1L+c_R )
                    {
                        q1L = Q1L;
                        q2L = 0;
                    }
                    else
                    {
                        q1L = Q1L;
                        q2L = (e+u-c_R-b*Q1L)/2;
                    }
                    
                    
                    
                    qE_L += q1L * phi[j];
                    qR_L += q2L * phi[j];
                    PiE_L += ((e + u - q1L - b*q2L - s)*q1L + (s-c_E)*Q1L )* phi[j];
                    PiR_L += (e + u - q2L - b*q1L - c_R)*q2L * phi[j] ;
                    CS_L += 0.5*( pow(q1L,2.0) + pow(q2L,2.0) + 2*b*q1L*q2L ) * phi[j];
                }
                
                qE_L /= STEP;
                qR_L /= STEP;
                PiE_L /= STEP;
                PiR_L /= STEP;
                CS_L /= STEP;

                
                    
                
                //Case II:  Q1 > (c2-s)/(2-b)
                double Q1H=0,q1H=0,q2H=0;
                double qR_H=0,qE_H=0,PiR_H=0,PiE_H=0,CS_H=0;

                //search for Q1H
                double lb3=(c_R-s)/(2-b), ub3=u;	//lower and upper bounds used in the searching algorithm
                Q1H=lb3;
                
                if ((ub3>lb3)&&(foc1_ERh_II(Q1H)>0))
                {
                    Q1H=(lb3+ub3)/2;
                    double temp=foc1_ERh_II(Q1H);
                    
                    for (;ub3-lb3>0.001;)
                    {
                        if (temp>0)
                            lb3 = Q1H;
                        else
                            ub3 = Q1H;
                        
                        Q1H = (lb3+ub3)/2;
                        temp = foc1_ERh_II(Q1H);
                    }
                }
                
                //for given Q1H, calculating expected q1, q2, Pi1, Pi2 with regard to the second signal e
                for (int j=0;j<=2*N*STEP;j++)
                {
                    double e = ((double)j/STEP-N)*stdev; //realization of the second signal e
                    
                    if (e + u < s)
                    {
                        q1H = 0;
                        q2H = 0;
                    }
                    else if (e + u < (2*c_R-b*s)/(2-b))
                    {
                        q1H = (e+u-s)/2;
                        q2H = 0;
                    }
                    else if (e + u < (2+b)*Q1H + (2*s-b*c_R)/(2-b))
                    {
                        q1H = ( 2*(e+u-s) - b*(e+u-c_R)) / (4-b*b);
                        q2H = ( 2*(e+u-c_R) - b*(e+u-s)) / (4-b*b);
                    }
                    else
                    {
                        q1H = Q1H;
                        q2H = (e+u-c_R-b*Q1H)/2;
                    }
                    
                    
                    qE_H += q1H * phi[j];
                    qR_H += q2H * phi[j];
                    PiE_H += ((e + u - q1H - b*q2H - s)*q1H + (s-c_E)*Q1H )* phi[j];
                    PiR_H += (e + u - q2H - b*q1H - c_R)*q2H * phi[j] ;
                    CS_H += 0.5*( pow(q1H,2.0) + pow(q2H,2.0) + 2*b*q1H*q2H ) * phi[j];
                }
                
                
                qE_H /= STEP;
                qR_H /= STEP;
                PiE_H /= STEP;
                PiR_H /= STEP;
                CS_H /= STEP;
                
                
                //comparison of case I and II
                if (PiE_H > PiE_L)
                {
                    Q1_ER=Q1H;
                    qE_ER=qE_H;
                    qR_ER=qR_H;
                    PiE_ER=PiE_H;
                    PiR_ER=PiR_H;
                    CS_ER=CS_H;
                }
                else
                {
                    Q1_ER=Q1L;
                    qE_ER=qE_L;
                    qR_ER=qR_L;
                    PiE_ER=PiE_L;
                    PiR_ER=PiR_L;
                    CS_ER=CS_L;
                }
                
                TS_ER = CS_ER + PiR_ER + PiE_ER;
                
                
                //file << Q1_ER << "\t" << qE_ER << "\t"  << qR_ER << "\t"  << PiE_ER << "\t"  << PiR_ER << "\t"  << CS_ER << "\t"  << TS_ER << "\t";
                //cout << Q1_ER << "\t" << qE_ER << "\t"  << qR_ER << "\t"  << PiE_ER << "\t"  << PiR_ER << "\t"  << CS_ER << "\t"  << TS_ER << "\t";

                
                

                /*********************************************************
                //DO NOT OUTPUT DETAILED NUMBER IN THE PRODUCTION VERSION

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

                *********************************************************/



                ////////////////////////////////////////////////////////////////////////
                //keep track of the thresholds on $var$
                
                if (Pi_EE>=PiR_ER)
                    NashEh[kk][jj] = var;

                if (PiE_ER>=Pi_RR)
                    NashRh[kk][jj] = var;

                if ((CS_EE>=CS_ER)&&(CS_EE>=CS_RR))
                    CS01h[kk][jj] = var;

                if ((CS_EE>=CS_RR)||(CS_ER>=CS_RR))
                    CS12h[kk][jj] = var;

                if ((TS_EE>=TS_ER)&&(TS_EE>=TS_RR))
                    TS01h[kk][jj] = var;

                if ((TS_EE>=TS_RR)||(TS_ER>=TS_RR))
                    TS12h[kk][jj] = var;

                if ((Pi_EE*2>=PiE_ER+PiR_ER)&&(Pi_EE>=Pi_RR))
                    Coop01h[kk][jj] = var;

                if ((Pi_EE>=Pi_RR)||(PiE_ER+PiR_ER>=2*Pi_RR))
                    Coop12h[kk][jj] = var;


            }   //end-for var

        }   //end-for b

        //output the thresholds to file2
        for (int jj=0; jj<=JJ; jj++)
            file2 << K << "\t" << delta << "\t" << s << "\t" << jj*0.01 << "\t" << NashEh[kk][jj]  << "\t" << NashRh[kk][jj]  << "\t" << CS01h[kk][jj]  << "\t" << CS12h[kk][jj]  << "\t" << TS01h[kk][jj]  << "\t" << TS12h[kk][jj]  << "\t" << Coop01h[kk][jj]  << "\t" << Coop12h[kk][jj] << endl;

    }   //end big parameter loop



    //file.close();
    file2.close();
    
    return 0;

}
