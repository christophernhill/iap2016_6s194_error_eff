#include <stdio.h>
#include <stdlib.h>
#include <math.h>
// CPP Macros to control whether we truncate
// - initial conditions
// - time varying part
#define TRUNC_IC
#define TRUNC_VARYING
#define TRUNC_ON
// #undef  TRUNC_ON

/* Helper macros */
#define HEX__(n) 0x##n##LU
#define B8__(x) ((x&0x0000000FLU)?1:0) \
+((x&0x000000F0LU)?2:0) \
+((x&0x00000F00LU)?4:0) \
+((x&0x0000F000LU)?8:0) \
+((x&0x000F0000LU)?16:0) \
+((x&0x00F00000LU)?32:0) \
+((x&0x0F000000LU)?64:0) \
+((x&0xF0000000LU)?128:0)

/* User macros */
#define B8(d) ((unsigned char)B8__(HEX__(d)))
#define B16(dmsb,dlsb) (((unsigned short)B8(dmsb)<<8) \
+ B8(dlsb))
#define B32(dmsb,db2,db3,dlsb) (((unsigned long)B8(dmsb)<<24) \
+ ((unsigned long)B8(db2)<<16) \
+ ((unsigned long)B8(db3)<<8) \
+ B8(dlsb))
#define B64(dmsb,db2,db3,db4,db5,db6,db7,dlsb) (  \
  ((unsigned long)B8(dmsb)<<56)                   \
+ ((unsigned long)B8(db2)<<48)                    \
+ ((unsigned long)B8(db3)<<40)                    \
+ ((unsigned long)B8(db4)<<32)                    \
+ ((unsigned long)B8(db5)<<24)                    \
+ ((unsigned long)B8(db6)<<16)                    \
+ ((unsigned long)B8(db7)<<8)                     \
+ B8(dlsb)                                        \
  )


typedef union {
  double f;
  char   fc[8];
} double_cast;

double bitTrunc(double var,int truncOff){
 double_cast tmpVal;
 tmpVal.f=var;
    tmpVal.fc[0]=0;
    tmpVal.fc[1]=0;
    tmpVal.fc[2]=0;
    tmpVal.fc[3]=0;
    tmpVal.fc[4]=0;
    tmpVal.fc[5]=0;
 // tmpVal.fc[7]=0;
 return tmpVal.f;
}

int main(int argc, char *argv[]){
 int truncOff;   // How far from end to truncate
 int nread;

 // Where to truncate from
 truncOff = 0;
 if ( argc > 1 ) {
  nread=sscanf(argv[1],"%d",&truncOff);
  if ( nread != 1 ) {
   fprintf(stderr,"Usage: ./a.out truncOff\n");
   exit(-1);
  }
 }

 // Variables
 // Problem size
 int Nl=100;
 int Nt=100;
 int NtCycl=10;  // Assumed to be integer when checking reference

 // Index variables initial values
 int jNM1orig;
 int jNM1=0;
 int jN=1;
 int jNP1=2;

 //
 long unsigned const number = B64(00000000,
                                  00000000,
                                  00000000,
                                  00000000,
                                  00000000,
                                  00000000,
                                  00000001,
                                  00000101);
 printf("number %lu\n",number);

 // Loop counters
 int i, j;

 // Scalars
 double_cast dt2=(double_cast)pow((double)1.0/(double)Nt,(double)2);
 double_cast dx2=(double_cast)pow((double)1.0/(double)Nl,(double)2);
 double_cast c2=(double_cast)((double)1.0);

 // Time steppig uVals array
 double_cast uVals[3][Nl];

 j=0;
 double vTmp;
 for (i=0;i<Nl;++i) {
  vTmp=sin(2.*M_PI/((double)Nl)*(double)((double)i+0.5))*cos(2.*M_PI/(double)Nt*(double)(j));
#ifdef TRUNC_ON
  vTmp=bitTrunc(vTmp,truncOff);
#endif
  uVals[jNM1][i]=(double_cast)(vTmp);
 }
 j=1;
 for (i=0;i<Nl;++i) {
  vTmp=sin(2.*M_PI/((double)Nl)*(double)((double)i+0.5))*cos(2.*M_PI/(double)Nt*(double)(j));
#ifdef TRUNC_ON
  vTmp=bitTrunc(vTmp,truncOff);
#endif
  uVals[jN][i]=(double_cast)(vTmp);
 }

 printf("%le\n",M_PI);
 for (i=0;i<Nl;++i) {
  printf("%le\n",uVals[jN][i].f);
 }
 
 // Calculate reference solution (at end time)
 // time runs fron 0 to dt*(Nt-1)
 double uRef[Nl];
 j=Nt-1;
 for (i=0;i<Nl;++i) {
  vTmp=sin(2.*M_PI/((double)Nl)*(double)((double)i+0.5))*cos(2.*M_PI/(double)Nt*(double)(j));
  // vTmp=bitTrunc(vTmp,truncOff);
  uRef[i]=(double)(vTmp);
 }

 // Calculate numerical solution
 double_cast v1, v2;
 int iM1, iP1;
 for (j=2;j<Nt*NtCycl;++j) {
  for (i=0;i<Nl;++i) {
   iM1=i-1; if ( iM1 <   0 ) { iM1 = Nl-1; }
   iP1=i+1; if ( iP1 == Nl ) { iP1 = 0;    }
   v1=(double_cast)(2.*uVals[jN][i].f-uVals[jNM1][i].f);
   v2=(double_cast)(
       dt2.f*c2.f*(
        uVals[jN][iP1].f +
        uVals[jN][iM1].f -
        2.*uVals[jN][i].f
       )/dx2.f
    );
    vTmp=v1.f+v2.f;
#ifdef TRUNC_ON
    vTmp=bitTrunc(vTmp,truncOff);
#endif
    uVals[jNP1][i].f=vTmp;
  }
  jNM1orig=jNM1;
  jNM1=jN;
  jN=jNP1;
  jNP1=jNM1orig;
 }

//calculate l1,l2 and infinity norm (Jenny)
 double diff[Nl];
 double l1, l2, infinity_norm;

 for (i=0;i<Nl;++i) {
  printf("%le %le %le\n",uVals[jN][i].f,uRef[i],uRef[i]-uVals[jN][i].f);
  diff[i] = uRef[i]-uVals[jN][i].f;
 }

 l1 = 0.0;
 l2 = 0.0;
 infinity_norm = diff[0];
 
 for (i=0;i<Nl;i++){
  l1 = l1+fabs(diff[i]);
  l2 = l2+pow(diff[i], 2);
  if (fabs(diff[i])>infinity_norm){
    infinity_norm = fabs(diff[i]);
  }
 }
 l2 = pow(l2, 0.5);

 printf("l1 =%le l2 =%le infinity_norm =%le\n", l1, l2, infinity_norm);
 
}
