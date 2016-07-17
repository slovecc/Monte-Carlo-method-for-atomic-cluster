/*transizione di fase Lennard Jones per un cluster Argon-Xenon
Per la compilazione:
gcc -lm -o LJ_Ar-Xe LJ_Ar-Xe.c ran2.c                 
 */ 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define SEED -1L
#define sigma_Xe 1.1490
#define epsi_Xe 1.9484
#define sigma_Ar 1.0
#define epsi_Ar 1.0
#define sigma_XA 1.0745
#define epsi_XA 1.3958
//dichiarazione delle funzioni
void writeCoord(double * rx_Xe, double * ry_Xe, double * rz_Xe,
                double * rx_Ar, double * ry_Ar, double * rz_Ar) ;

void ask(int *N,double *L,double *T,double *dr,double *rc2,int *nCycles,int *nEq,int *ntsteps);

double total_e ( double * rx_Xe, double * ry_Xe, double * rz_Xe,
                 double * rx_Ar, double * ry_Ar, double * rz_Ar, 
                 int N, double L,double rc2, 
                 double ecut, 
                 double *vir_Xe,double * vir_Ar, double * vir_XA );

void init ( double * rx_Xe, double * ry_Xe, double * rz_Xe,
            double * rx_Ar, double * ry_Ar, double * rz_Ar, double L);

//variabili globali
FILE *file;
int N, nCycles, nEq, ntsteps;
double L, T, dr, rc2;

/*_________funzione main___________*/

void main (  ) {

  extern float ran2();
//variabili Xenon  
  double * rx_Xe, * ry_Xe, * rz_Xe;
  double vir_Xe;
  double dx_Xe,dy_Xe,dz_Xe, rxold_Xe,ryold_Xe,rzold_Xe; 
  double  rr3_Xe, ecut_Xe;
  
//variabili Argon
  double * rx_Ar, * ry_Ar, * rz_Ar;
  double vir_Ar;
  double dx_Ar,dy_Ar,dz_Ar, rxold_Ar,ryold_Ar,rzold_Ar; 
  double   rr3_Ar, ecut_Ar;
 
//variabili interazione Xenon-Argon
  double vir_XA;
   double  rr3_XA, ecut_XA;
  

  double E_new, E_old, vir_old;
  double ecut, esum, vir_sum, Press;
  double spost;
  int c,a;
  double rho, V;
  int i,j,k;
  int  nSamp, nAcc;
  double tfactor=0.999;//fattore di riduzione della temperatura
  double Tinit, RatioAcc; 
  long idum;
  double change_x,change_y, change_z;//variabili per lo scambio

   Tinit=T;

  file=fopen("LJ_Ar-Xe.data","w");

// inizializzazione generatore
  idum=SEED;
  ran2(&idum);

/*_________chiama ask_________*/
ask(&N,&L,&T,&dr,&rc2,&nCycles,&nEq,&ntsteps);
  rho=((double)N)/(2*pow(L,3)); //N/2 perchè considero metà delle particelle alla volta


//Calcolo del potenziale di taglio
//Xenon
  rr3_Xe = (pow(sigma_Xe,3))/(rc2*rc2*rc2);
  ecut_Xe = 4*epsi_Xe*(rr3_Xe*rr3_Xe*rr3_Xe*rr3_Xe-rr3_Xe*rr3_Xe);
   
//Argon
  rr3_Ar = (pow(sigma_Ar,3))/(rc2*rc2*rc2);
  ecut_Ar = 4*epsi_Ar*(rr3_Ar*rr3_Ar*rr3_Ar*rr3_Ar-rr3_Ar*rr3_Ar);
  
//Xenon-Argon
  rr3_XA = (pow(sigma_XA,3))/(rc2*rc2*rc2);
  ecut_XA = 4*epsi_XA*(rr3_XA*rr3_XA*rr3_XA*rr3_XA-rr3_XA*rr3_XA);

//ecut totale  
  ecut=ecut_Xe+ecut_Ar+ecut_XA; 

// Quadrato del raggio di troncamento
  rc2*=rc2;

//Per efficienza computazionale uso 1/T
  T = 1.0/T;

//calcolo del volume
  V = L*L*L;
  
// Informazioni stampate sullo schermo
      fprintf(stdout,"NVT Metropolis Monte Carlo Simulation"
	    " of the Lennard-Jones fluid.\n"
	    "---------------------------------------------\n");
  fprintf(stdout,"# L = %.5lf; rho = %.5lf; N = %i; rc = %.5lf\n",
	  L,rho,N,sqrt(rc2));
  fprintf(stdout,"# nCycles %i, nEq %i, dR %.5lf\n\n",
	  nCycles,nEq,dr);
  fprintf(stdout,"\nResults:\n");

 fprintf(stdout," AcceptionRatio   Energy      virial    pressure   Temperature\n");
 
//Numero totale di cicli=numero di cicli di "equilibrazione"+numero di cicli di produzione
  nCycles+=nEq;

//allocozione per i vettori delle posizioni per Xe e Ar
//Xenon
  rx_Xe = (double*)malloc(0.5*N*sizeof(double));
  ry_Xe = (double*)malloc(0.5*N*sizeof(double));
  rz_Xe = (double*)malloc(0.5*N*sizeof(double));
//Argon
  rx_Ar = (double*)malloc(0.5*N*sizeof(double));
  ry_Ar = (double*)malloc(0.5*N*sizeof(double));
  rz_Ar = (double*)malloc(0.5*N*sizeof(double)); 


/*________chiamata alla funzione init_________*/
  init(rx_Xe,ry_Xe,rz_Xe,rx_Ar,ry_Ar,rz_Ar,L);
//loop sulle temperature
for(k=0;k<=ntsteps;k++){ 
/*_________chiamata alla funzione total_e_____*/
  E_old = total_e(rx_Xe, ry_Xe, rz_Xe, rx_Ar, ry_Ar, rz_Ar,
                      N,L,rc2,ecut,
                      &vir_Xe, &vir_Ar, &vir_XA);
  nAcc = 0;
  esum = 0.0;
  nSamp = 0;
  vir_sum = 0.0;
  for (c=0;c<nCycles;c++) {
     spost=ran2(&idum); //tramite spost decido se effettuare uno spostamento o scambio di posto
      if(spost<=0.45){//sposto lo Xe
    //scelgo casualmente un atomo i
      i=(int)(ran2(&idum)*N*0.5);
    //calcolo lo spostamento
      dx_Xe = dr*(0.5-ran2(&idum));
      dy_Xe = dr*(0.5-ran2(&idum));
      dz_Xe = dr*(0.5-ran2(&idum));

    //salvo la posizione corrente
     rxold_Xe=rx_Xe[i];
     ryold_Xe=ry_Xe[i];
     rzold_Xe=rz_Xe[i];

    //sposto l'atomo i
     rx_Xe[i]+=dx_Xe;
     ry_Xe[i]+=dy_Xe;
     rz_Xe[i]+=dz_Xe;
   

    //applico le condizioni al contorno periodiche
    if (rx_Xe[i]<0.0) rx_Xe[i]+=L;
    if (rx_Xe[i]>L)   rx_Xe[i]-=L;
    if (ry_Xe[i]<0.0) ry_Xe[i]+=L;
    if (ry_Xe[i]>L)   ry_Xe[i]-=L;
    if (rz_Xe[i]<0.0) rz_Xe[i]+=L;
    if (rz_Xe[i]>L)   rz_Xe[i]-=L;

    //calcolo la nuova energia
    E_new = total_e(rx_Xe, ry_Xe, rz_Xe, rx_Ar, ry_Ar, rz_Ar,
                      N,L,rc2,ecut,
                      &vir_Xe, &vir_Ar, &vir_XA);
      
    //accetto?
    if (ran2(&idum) < exp(-T*(E_new-E_old))) {
      E_old=E_new;
      vir_old=vir_Xe+vir_Ar+vir_XA;
      nAcc++;
    }
    //se non accetto riassegno la posizione precedente
    else {
      rx_Xe[i]=rxold_Xe;
      ry_Xe[i]=ryold_Xe;
      rz_Xe[i]=rzold_Xe;
    }

    //accumulo i risultati sia in caso di accetazione che nel caso opposto
    if (c>nEq) {
      esum+=E_old;
      vir_sum+=vir_old;
      nSamp++;
    }
   }//chiusura dello spostamento dello Xe
 
// sposto solo l' Ar

   else if((0.45<spost)&&(spost<=0.9)){
   //scielgo casualmente un atomo i
    i=(int)(ran2(&idum)*N*0.5);
    //calcolo lo spostamento
    dx_Ar = dr*(0.5-ran2(&idum));
    dy_Ar = dr*(0.5-ran2(&idum));
    dz_Ar = dr*(0.5-ran2(&idum));

    //salvo la posizione corrente
    rxold_Ar=rx_Ar[i];
    ryold_Ar=ry_Ar[i];
    rzold_Ar=rz_Ar[i];

    //sposto l'atomo i
    rx_Ar[i]+=dx_Ar;
    ry_Ar[i]+=dy_Ar;
    rz_Ar[i]+=dz_Ar;
  

    //applico le condizioni al contorno periodiche
    if (rx_Ar[i]<0.0) rx_Ar[i]+=L;
    if (rx_Ar[i]>L)   rx_Ar[i]-=L;
    if (ry_Ar[i]<0.0) ry_Ar[i]+=L;
    if (ry_Ar[i]>L)   ry_Ar[i]-=L;
    if (rz_Ar[i]<0.0) rz_Ar[i]+=L;
    if (rz_Ar[i]>L)   rz_Ar[i]-=L;

    //calcolo la nuova energia
    E_new = total_e(rx_Xe, ry_Xe, rz_Xe, rx_Ar, ry_Ar, rz_Ar,
                      N,L,rc2,ecut,
                      &vir_Xe, &vir_Ar, &vir_XA);
       
    //accetto?
    if (ran2(&idum) < exp(-T*(E_new-E_old))) {
      E_old=E_new;
      vir_old=vir_Xe+vir_Ar+vir_XA;
      nAcc++;
    }
     //se non accetto riassegno la posizione precedente
    else {
      rx_Ar[i]=rxold_Ar;
      ry_Ar[i]=ryold_Ar;
      rz_Ar[i]=rzold_Ar;
    }

    //accumulo i risultati sia in caso di accetazione che nel caso opposto
    if (c>nEq) {
      esum+=E_old;
      vir_sum+=vir_old;
      nSamp++;
    }
   }//chiusura dello spostamento dello Ar
    //effettuo lo scambio di due atomi di specie diverse
   else if(spost>0.9){
   //scelgo due atomi a caso
    i=(int)(ran2(&idum)*N*0.5);
    j=(int)(ran2(&idum)*N*0.5);

    //salvo la posizione corrente
    rxold_Ar=rx_Ar[i];
    ryold_Ar=ry_Ar[i];
    rzold_Ar=rz_Ar[i];
    rxold_Xe=rx_Xe[j];
    ryold_Xe=ry_Xe[j];
    rzold_Xe=rz_Xe[j];

    /* scambio di posto */
   
    change_x=rx_Xe[j];
    change_y=ry_Xe[j];
    change_z=rz_Xe[j];
    rx_Xe[j]=rx_Ar[i];
    ry_Xe[j]=ry_Ar[i];
    rz_Xe[j]=rz_Ar[i];
    rx_Ar[i]=change_x;
    ry_Ar[i]=change_y;
    rz_Ar[i]=change_z;

    //calcolo la nuova energia
    E_new = total_e(rx_Xe, ry_Xe, rz_Xe, rx_Ar, ry_Ar, rz_Ar,
                      N,L,rc2,ecut,
                      &vir_Xe, &vir_Ar, &vir_XA);
      
   //accetto?
    if (ran2(&idum) < exp(-T*(E_new-E_old))) {
      E_old=E_new;
      vir_old=vir_Xe+vir_Ar+vir_XA;
      nAcc++;
    }
   //se non accetto riassegno la posizione precedente
    else {
      rx_Ar[i]=rxold_Ar;
      ry_Ar[i]=ryold_Ar;
      rz_Ar[i]=rzold_Ar;
      rx_Xe[j]=rxold_Xe;
      ry_Xe[j]=ryold_Xe;
      rz_Xe[j]=rzold_Xe;
    }

    //accumulo i risultati sia in caso di accetazione che nel caso opposto
    if (c>nEq) {
      esum+=E_old;
      vir_sum+=vir_old;
      nSamp++;
    }
   }
  }//chiusura degli ncycles!

  //Calcolo il rapporto di accettazione e le quantità medie
      
      RatioAcc=((double)nAcc)/(nCycles);
      esum/=nSamp;
      vir_sum/=nSamp;
      Press=rho*(T-vir_sum);
  
      fprintf(stdout," %14.8f   %7.5f   %7.5f   %7.5f    %7.5f\n", RatioAcc,esum,vir_sum,Press,1./T);

fprintf(file, "%d %12.4f %12.4f  %12.4f \n", j, 1.0/T, esum/nSamp, Press);

if(k%50==0)//salva ogni 100 
   { 
     writeCoord(rx_Xe,ry_Xe,rz_Xe,rx_Ar,ry_Ar,rz_Ar);
   }
T=1.0/T;
T*=tfactor;//annealing schedule
T=1.0/T;
rc2=2.8717; //ripristina prima di iniziare un altro loop

 }//chiusura loop sulle Temperature
}//chiusura del main


/*______subroutine writeCoord()_______*/
//in questa subroutine mi trascrivo le posizioni degli atomi alla fine del ciclo fissata la T
//conterrà gli snapshots in ogni file numerato diversamente
void writeCoord(double * rx_Xe, double * ry_Xe, double * rz_Xe,
                double * rx_Ar, double * ry_Ar, double * rz_Ar) 
  {
    int i;
    char filename[100];
    static int num=0;
    sprintf(filename, "coords_%d.dat",num);
    FILE* output=fopen (filename, "w");
    for(i=0;i<N/2;i++){         
	fprintf(output, "%lf  %lf  %lf  %lf  %lf  %lf   \n", rx_Xe[i],ry_Xe[i],rz_Xe[i], rx_Ar[i],ry_Ar[i],rz_Ar[i]);
		}
    num++;
    fclose(output);
  }


/*_________subroutine ask_________*/
void ask(int *N,double *L,double *T,double *dr,double *rc2,int *nCycles,int *nEq,int *ntsteps){
    fprintf(stderr, "Number of atoms (suggested:144) N=");
    scanf("%i",N);
    fprintf(stderr, "side cube (suggested 40->11.7813): L=");
    scanf("%lf",L);
    fprintf(stderr, "Requested temperature  (suggested 3.4249) T=");
    scanf("%lf",T);
    fprintf(stderr, "step dr (suggested:0.1) dr=");
    scanf("%lf",dr);
    fprintf(stderr, "r cut (suggested:9.75->2.8717) rc=");
    scanf("%lf",rc2);
    fprintf(stderr, "Number of cycles (suggested:1000) Ncycles=");
    scanf("%i",nCycles);
    fprintf(stderr, "Number of equilibritation steps (suggested:1000) Neq=");
    scanf("%i",nEq);
    
    fprintf(stderr, "Number of steps of temperature (suggested:1000) NTsteps=");
    scanf("%i",ntsteps);
    }


/*__ Souroutine per il calcolo dell'energia e del viriale ___*/
double total_e ( double * rx_Xe, double * ry_Xe, double * rz_Xe,
                 double * rx_Ar, double * ry_Ar, double * rz_Ar, 
                 int N, double L,double rc2, 
                 double ecut, 
                 double *vir_Xe,double * vir_Ar, double * vir_XA ){
   int i,j;
   double dx_Xe, dy_Xe, dz_Xe, r2_Xe, r6i_Xe;
   double dx_Ar, dy_Ar, dz_Ar, r2_Ar, r6i_Ar;
   double dx_XA, dy_XA, dz_XA, r2_XA, r6i_XA;
   double e_Xe = 0.0, e_XA = 0.0,e_Ar = 0.0,hL=L/2.0;
   double e_tot;

   *vir_Xe=0.0;
   *vir_XA=0.0;
   *vir_Ar=0.0;

//calcolo energia per Xenon
   for (i=0;i<=(N-2)/2;i++) {
     for (j=i+1;j<N/2;j++) {
	dx_Xe  = (rx_Xe[i]-rx_Xe[j]);
	dy_Xe  = (ry_Xe[i]-ry_Xe[j]);
	dz_Xe  = (rz_Xe[i]-rz_Xe[j]);
	//Convenzione dell'immagine minima
	if (dx_Xe>hL)       dx_Xe-=L;
	else if (dx_Xe<-hL) dx_Xe+=L;
	if (dy_Xe>hL)       dy_Xe-=L;
	else if (dy_Xe<-hL) dy_Xe+=L;
	if (dz_Xe>hL)       dz_Xe-=L;
	else if (dz_Xe<-hL) dz_Xe+=L;
	r2_Xe = dx_Xe*dx_Xe + dy_Xe*dy_Xe + dz_Xe*dz_Xe;
	if (r2_Xe<rc2) {
	  r6i_Xe   = (pow(sigma_Xe,3))/(r2_Xe*r2_Xe*r2_Xe);
	  e_Xe   += 4*(epsi_Xe)*(r6i_Xe*r6i_Xe - r6i_Xe);
	  *vir_Xe += 48*epsi_Xe*(r6i_Xe*r6i_Xe-0.5*r6i_Xe);
          *vir_Xe/=(3*N/2);
	}
     }
   }
//calcolo energia per Argon
   for (i=0;i<=(N-2)/2;i++) {
     for (j=i+1;j<N/2;j++) {
	dx_Ar  = (rx_Ar[i]-rx_Ar[j]);
	dy_Ar  = (ry_Ar[i]-ry_Ar[j]);
	dz_Ar  = (rz_Ar[i]-rz_Ar[j]);
	//Convenzione dell'immagine minima
	if (dx_Ar>hL)       dx_Ar-=L;
	else if (dx_Ar<-hL) dx_Ar+=L;
	if (dy_Ar>hL)       dy_Ar-=L;
	else if (dy_Ar<-hL) dy_Ar+=L;
	if (dz_Ar>hL)       dz_Ar-=L;
	else if (dz_Ar<-hL) dz_Ar+=L;
	r2_Ar = dx_Ar*dx_Ar + dy_Ar*dy_Ar + dz_Ar*dz_Ar;
	if (r2_Ar<rc2) {
	  r6i_Ar   = (pow(sigma_Ar,3))/(r2_Ar*r2_Ar*r2_Ar);
	  e_Ar  += 4*(epsi_Ar)*(r6i_Ar*r6i_Ar - r6i_Ar);
	  *vir_Ar += 48*epsi_Ar*(r6i_Ar*r6i_Ar-0.5*r6i_Ar);
          *vir_Ar/=(3*N/2); 
	}
     }
   }
//calcolo energia per Argon-Xenon 
  for (i=0;i<N/2;i++) {
     for (j=0;j<N/2;j++) {
	dx_XA  = (rx_Ar[i]-rx_Xe[j]);
	dy_XA  = (ry_Ar[i]-ry_Xe[j]);
	dz_XA  = (rz_Ar[i]-rz_Xe[j]);
	//Convenzione dell'immagine minima
	if (dx_XA>hL)       dx_XA-=L;
	else if (dx_XA<-hL) dx_XA+=L;
	if (dy_XA>hL)       dy_XA-=L;
	else if (dy_XA<-hL) dy_XA+=L;
	if (dz_XA>hL)       dz_XA-=L;
	else if (dz_XA<-hL) dz_XA+=L;
	r2_XA = dx_XA*dx_XA + dy_XA*dy_XA + dz_XA*dz_XA;
	if (r2_XA<rc2) {
	  r6i_XA   = (pow(sigma_XA,3))/(r2_XA*r2_XA*r2_XA);
	  e_XA += 4*(epsi_XA)*(r6i_XA*r6i_XA - r6i_XA);
	  *vir_XA += 48*epsi_XA*(r6i_XA*r6i_XA-0.5*r6i_XA);
          *vir_XA/=(3*N/2); 
	}
     }
   }
   e_tot=e_Xe+e_Ar+e_XA;
   return e_tot-ecut;
}

/*____ subroutine inizializzazione posizione_____ */
void init ( double * rx_Xe, double * ry_Xe, double * rz_Xe,double * rx_Ar, double * ry_Ar, double * rz_Ar, double L){
  int i,ix,iy,iz;
  double n=N/2;
  int n3=2;
  FILE *new;
  new=fopen("conf_init.data","w");
  //Trovo il più piccolo numero n3 il cui cubo è maggiore o uguale a n 
  while ((n3*n3*n3)<n) n3++;

  ix=iy=iz=0;
  //posizioni per Xe
  for (i=0;i<n;i++) {
    rx_Xe[i] = ((double)ix+0.3)*L/n3;            
    ry_Xe[i] = ((double)iy+0.3)*L/n3;
    rz_Xe[i] = ((double)iz+0.3)*L/n3;
    ix++;
    if (ix==n3) {
      ix=0;
      iy++;
      if (iy==n3) {
	iy=0;
	iz++;
      }
    }
  }
 //posizioni per l'Ar
  ix==0;
  iy==0;
  iz==0;
  for (i=0;i<n;i++) {
    rx_Ar[i] = ((double)ix+0.1)*L/n3;            
    ry_Ar[i] = ((double)iy+0.1)*L/n3;
    rz_Ar[i] = ((double)iz+0.1)*L/n3;
    ix++;
    if (ix==n3) {
      ix=0;
      iy++;
      if (iy==n3) {
	iy=0;
	iz++;
      }
    }
  }
 for(i=0;i<n;i++){
 fprintf(new, "%12.4f  %12.4f  %12.4f   %12.4f   %12.4f   %12.4f\n", rx_Xe[i],ry_Xe[i],rz_Xe[i],rx_Ar[i],ry_Ar[i],rz_Ar[i]);
 }
 fclose(new);
}



