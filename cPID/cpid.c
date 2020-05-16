 // cpid.c ----------------------------------------------------------

 #include <stdio.h>

 // Struttura per contenere i dati relativi a un PID: questo
 // consente di usare una sola istanza del codice per implementare
 // più regolatori dato che ognuno agisce sulla sua zona dati
 typedef struct {
         float K,Ti,Td,N,b,c,CSmax,CSmin,Ts; // parametri
         float SP,SPo,PV,PVo,Do,CS,CSo,TR;   // variabili; o=old --> 4 elementi da "ricordare"
																 // (2 stati del PID e 2 per i delta)
         char  TS; // TR=track reference, TS=track switch
         } PID_DATA;
               
 // Struttura analoga alla precedente per contenere i dati di un
 // sistema dinamico SISO con funzione di trasferimento mu/(1+sT)
 // discretizzato a passo Ts con il metodo di Eulero implicito
 typedef struct {
         float mu,T,Ts; // parametri
         float u,y,xo;  // variabili (ingresso u, uscita y)
         } TF1P0Z_DATA; // funzione di trasferimento 1 polo, nessuno zero

 //------------------------------------------------------------------
 // Funzione che implementa il PID, da chiamarsi ad ogni passo
 // Si noti che il _solo_ argomento è il puntatore alla zona dati
 // del PID: la funzione è void per semplicità, si potrebbe usare
 // un valore int di ritorno per eventuale diagnostica
 void isapid(PID_DATA *pd)
      {
      float DSP,DPV,DP,DI,D,DD,DCS;
      if (pd->TS==0) // modalità Automatico
         {
         DSP    = pd->SP-pd->SPo;
         DPV    = pd->PV-pd->PVo;
         DP     = pd->K*(pd->b*DSP-DPV);
         DI     = pd->K*pd->Ts/pd->Ti*(pd->SP-pd->PV);
         D      = (pd->Td*pd->Do+pd->K*pd->N*pd->Td
                  *(pd->c*DSP-DPV))
                  /(pd->Td+pd->N*pd->Ts);
         DD     = D-pd->Do;
         DCS    = DP+DI+DD;
         pd->CS = pd->CSo+DCS;
         }
      else           // modalità Tracking
         {
         pd->CS = pd->TR;
         D      = 0; // Assegnammento a rigore arbitrario
         }
      if (pd->CS>pd->CSmax) pd->CS = pd->CSmax; // antiwindup
      if (pd->CS<pd->CSmin) pd->CS = pd->CSmin;
      pd->SPo = pd->SP;                         // memorizzazione
      pd->PVo = pd->PV;                         // delle var. di
      pd->CSo = pd->CS;                         // stato per il
      pd->Do  = D;                              // prossimo passo
      }

 //------------------------------------------------------------------
 // Funzione che implementa un sistema del prim'ordine senza zeri,
 // da chiamarsi ad ogni passo
 void tf1p0z(TF1P0Z_DATA *td)
      {
      float x;
      x      = (td->T*(td->xo+td->u))/(td->T+td->Ts);
      td->y  = td->mu*td->Ts/(td->T+td->Ts)*(x+td->u);
      td->xo = x;
      }

 //------------------------------------------------------------------
 // Funzioni scalino e rampa
 float sca(float t) { return t>=0?1:0; }
 float ram(float t) { return t*sca(t); }

 //------------------------------------------------------------------
 // Programma principale che simula un sistema in retroazione in
 // cui il regolatore è un PID ISA e il processo un sistema del
 // prim'ordine senza zeri
 int main(void)
     {
     TF1P0Z_DATA  dataTF1;  // struttura dati per il processo
     PID_DATA     dataPID1; // struttura dati per il PID
     float        Ts;       // passo di campionamento
     float        t;        // tempo continuo (per plotting)
     int          nSteps;   // numero di passi da simulare
     int          k;        // tempo discreto (per plotting)
     FILE        *h;        // file per i risultati

     dataTF1.mu     = 1;    // assegnamento dei parametri
     dataTF1.T      = 10;   // del processo
    
     dataPID1.K     = 20;   // assegnamento dei parametri
     dataPID1.Ti    = 5;    // del PID
     dataPID1.Td    = 0;
     dataPID1.N     = 1;
     dataPID1.b     = 1;
     dataPID1.c     = 0;
     dataPID1.CSmax = 5;
     dataPID1.CSmin = -5;
    
     Ts             = 0.25; // assegnamento di Ts
     nSteps         = 800;  // e del numero di passi da simulare
    
     dataTF1.Ts     = Ts;
     dataPID1.Ts    = Ts;

     dataTF1.xo     =  0;   // inizializzazione delle variabili
     dataTF1.y      =  0;   // di stato (e dell'uscita del
     dataPID1.SPo   =  0;   // processo dato che al primo passo
     dataPID1.PVo   =  0;   // essa non è stata ancora calcolata
     dataPID1.CSo   =  0;   // ma il PID ne ha bisogno)
     dataPID1.Do    =  0;

     h = fopen("data.txt","w");  // apertura del file dati
     for(k=0;k<nSteps;k++)       // loop di simulazione
        {
        t           = k*Ts;
       
        dataPID1.TS = sca(t-120)-sca(t-150);  // assegnamento
        dataPID1.TR = 1;                      // ingressi esogeni
        dataPID1.SP = ram(t-1)-ram(t-5)
                      -0.1*ram(t-50)
                      +0.1*ram(t-80)
                      +sca(t-100);

        // Le istruzioni seguenti realizzano le connessioni dello
        // schema a blocchi dicendo "chi scrive e legge cosa nelle
        // strutture dati di chi"
        dataPID1.PV = dataTF1.y;              // lettura di y
        isapid(&dataPID1);                    // calcolo di u
        dataTF1.u   = dataPID1.CS;            // applicazione di u
        tf1p0z(&dataTF1);                     // calcolo di y
       
        fprintf(h,"%f\t%f\t%f\t%f\t%d\t%f\t%f\n", // scrittura dei dati
               t,dataPID1.SP,dataPID1.PV,         // sul file
               dataPID1.CS,dataPID1.TS,
               dataPID1.TR,dataPID1.CSo);
        } // termine del loop di simulazione
    fclose(h); // chiusura del file dati
    }
