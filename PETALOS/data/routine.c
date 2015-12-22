#include <stdio.h>
#include <stdlib.h>
#include <math.h>


int main (void){
  double *l; //List of the last Nperiod noise samples
  double noise, h0, h, n;
  int i, i_period, i_period2, i_period3, j, counter1, counter2;
  
  int step; // the step of the signal, in ps

  int T0; // the updating time for the noise in ps

  int T; // T0/step 

  int Nperiod; //How many noise sample we  want to convolute

  int N; //The time length of the convolution width (in ps)

  int TotalTime; //Total Simulated Time, in unit of step

  double noise_std_dev; // Noise std deviation
  
  //DEBUG FILE *out;
  //DEBUG FILE *out2;

  T0=500;
  step=5;
  T=(int)(T0/step);
  Nperiod = 16;
  N=T*Nperiod;
  TotalTime=64000;
  noise_std_dev=1.466E-3;
  noise =0.0;

  if ((l=(double *)malloc(sizeof(double)*Nperiod))==NULL){
    printf("memory allocation error on l\n");
    return 1;
  };


  /*Now I guess you'll have a for routine where you convolute the 1pe signal with the events distribution
   * I Add a first routine to initialize the noise signal (Nperiod steps).
   */
  //DEBUG out=fopen("noise_data.dat", "w");
  //DEBUG out2=fopen("h0_data.dat", "w");
  srand(1);
  for (i=0; i<Nperiod; i++){
    l[i]=((double)(rand()%1024-512))/512.0*noise_std_dev*(double)T*sqrt(3); //Noise sample
  };
  counter1=0;
  counter2=0;
  i_period=0;
  i_period2=0;
  i_period3=0;
  for (i =0; i<TotalTime; i++){
    counter1++;
    if (T==counter1){
      counter1=0;
      i_period++;
      counter2++;
      if (Nperiod==counter2){
        counter2=0;
        i_period3++;
      };
    };

    if (i_period*T==i){
      i_period2++;
      if(Nperiod==i_period2) i_period2=0;
      l[i_period2]=((double)(rand()%1024-512))/512.0*noise_std_dev*(double)T*sqrt(3); //Noise sample
    };
    noise=0.0;
    for(j=0; j<Nperiod; j++){
      n=((double)(i-(i_period3*Nperiod+j)*T));
      if(n==0) h0=2.0*1.0/(2.0*(double)T);
      else h0=sin(2.0*M_PI*n/(2.0*(double)T))/(M_PI*n);
      noise+=l[(i_period2+j+Nperiod/2) % Nperiod]*h0;
    };
  //DEBUG fprintf(out, "%lf\n", noise); //Here I write on a file, you can just add the value of "noise" to the signal sample.
  };
  //DEBUG fclose(out); 
  //DEBUG fclose(out2); 
 return 0;
}

  

