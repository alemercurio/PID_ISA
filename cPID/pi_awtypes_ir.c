 // pi_awtypes.c ----------------------------------------------------

 #include <stdio.h>

 // Struttura per contenere i dati relativi a un PI a 1gdl con aw
 // ottenuto tramite ricalcolo dell'azione integrale
 typedef struct {
         float K,Ti,CSmax,CSmin,Ts; // parametri
         float SP,PV,Io,CS,TR;   // variabili
         char  TS;
         } PI1dof_aw_ir_DATA;

 // Struttura per contenere i dati relativi a un PI a 1gdl con aw
 // ottenuto tramite ricalcolo dell'azione integrale
 typedef struct {
         float K,Ti,CSmax,CSmin,Ts; // parametri
         float SP,PV,xo,CS,TR;   // variabili
         char  TS;
         } PI1dof_aw_fb_DATA;
               
 // Struttura analoga alla precedente per contenere i dati di un
 // sistema dinamico SISO con funzione di trasferimento mu/(1+sT)
 // discretizzato a passo Ts con il metodo di Eulero implicito
 typedef struct {
         float mu,T,Ts; // parametri
         float u,y,xo;  // variabili (ingresso u, uscita y)
         } TF1P0Z_DATA;

 //------------------------------------------------------------------
 // Funzione che implementa il PI a 1gdl con aw
 // ottenuto tramite ricalcolo dell'azione integrale
 void pi1dof_aw_ir(PI1dof_aw_ir_DATA *pd)
      {
      float P,I;
      if (pd->TS==0) // modalità Automatico
         {
         P      = pd->K*(pd->SP-pd->PV);
         I      = pd->Io+pd->K*pd->Ts/pd->Ti*(pd->SP-pd->PV);
         pd->CS = P+I;
         }
      else           // modalità Tracking
         {
         pd->CS = pd->TR;
         }
      if (pd->CS>pd->CSmax) pd->CS = pd->CSmax; // antiwindup
      if (pd->CS<pd->CSmin) pd->CS = pd->CSmin;
      pd->Io  = pd->CS-P;
      }

 //------------------------------------------------------------------
 // Funzione che implementa il PI a 1gdl con aw
 // ottenuto tramite feedback
 void pi1dof_aw_fb(PI1dof_aw_fb_DATA *pd)
      {
      float x;
      if (pd->TS==0) // modalità Automatico
         {
         pd->CS = pd->xo+pd->K*(pd->SP-pd->PV);
         }
      else           // modalità Tracking
         {
         pd->CS = pd->TR;
         }
      if (pd->CS>pd->CSmax) pd->CS = pd->CSmax; // antiwindup
      if (pd->CS<pd->CSmin) pd->CS = pd->CSmin;
      x       = (pd->Ti*pd->xo+pd->Ts*pd->CS)/(pd->Ts+pd->Ti);
      pd->xo  = x;
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
     PI1dof_aw_ir_DATA dataPI1;  // struttura dati per il PI
     float        Ts;       // passo di campionamento
     float        t;        // tempo continuo
     int          nSteps;   // numero di passi da simulare
     int          k;        // tempo discreto
     FILE        *h;        // file per i risultati

     dataTF1.mu     = 1;    // assegnamento dei parametri
     dataTF1.T      = 10;   // del processo
    
     dataPI1.K      = 10;   // assegnamento dei parametri
     dataPI1.Ti     = 5;    // del PID
     dataPI1.CSmax  = 5;
     dataPI1.CSmin  = -5;
    
     Ts             = 0.25; // assegnamento di Ts
     nSteps         = 800;  // e del numero di passi da simulare
    
     dataTF1.Ts     = Ts;
     dataPI1.Ts     = Ts;

     dataTF1.xo     =  0;   // inizializzazione
     dataTF1.y      =  0;   
     dataPI1.Io    =  0;

     h = fopen("data.txt","w");  // apertura del file dati
     for(k=0;k<nSteps;k++)       // loop di simulazione
        {
        t          = k*Ts;
       
        dataPI1.TS = sca(t-120)-sca(t-150);  // assegnamento
        dataPI1.TR = 1;                      // ingressi esogeni
        dataPI1.SP = ram(t-1)-ram(t-5)
                     -0.1*ram(t-50)
                     +0.1*ram(t-80)
                     +sca(t-100);

        // Le istruzioni seguenti realizzano le connessioni dello
        // schema a blocchi dicendo "chi scrive e legge cosa nelle
        // strutture dati di chi"
        dataPI1.PV = dataTF1.y;              // lettura di y
        pi1dof_aw_ir(&dataPI1);                    // calcolo di u
        dataTF1.u   = dataPI1.CS;            // applicazione di u
        tf1p0z(&dataTF1);                     // calcolo di y
       
        fprintf(h,"%f\t%f\t%f\t%f\t%d\t%f\n", // scrittura dei dati
               t,dataPI1.SP,dataPI1.PV,         // sul file
               dataPI1.CS,dataPI1.TS,
               dataPI1.TR);
        } // termine del loop di simulazione
    fclose(h); // chiusura del file dati
    }
