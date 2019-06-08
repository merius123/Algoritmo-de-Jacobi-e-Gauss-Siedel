#include <iostream>
#include <cstdlib>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>



using namespace std;
#define MAX 1000
float A[MAX][MAX];
float V[MAX][MAX];
float C[MAX][MAX];
float X[MAX];
float Xw[MAX];
int linhasA, colunasA, linhasV, colunasV, linhasC, colunasC;
//-------------------------------------------------------------------

double resid[MAX];
double prov_max;
//_______________________

void DefTamA(){
    cout << "Quantidade de linhas e colunhas de A" <<endl;
    cin >> linhasA;
    colunasA = linhasA;
    if (linhasA>1000){
        cout<< "\n Valor excede o limite.";
        linhasA=0;
        return ;

    }
}

void DefTamV(){
    linhasV = colunasA;
    colunasV = 1;
    if (linhasV>1000){
        cout<< "\n Valor excede o limite.";
        linhasV=0;
        return ;

    }
}

void PreMatrizA(){
    for (int i=0; i<linhasA;i++){
        for (int j=0; j<colunasA; j++){
            if(i==j){

                A[i][j]=4;
            }
            else if ((i-j==1)||(i-j==-1)){
                A[i][j]=-1;
            }
            else{
                A[i][j]=0;
            }
        }
    }
  cout<<"\n VALORES GERADOS COM SUCESSO\n";
}
void PreVetorV(){
    srand(time(NULL));
    for (int i=0; i<linhasV;i++){
        for (int j=0; j<colunasV; j++){
            V[i][j]=rand()%2;
            if(V[i][j]==0) V[i][j]=-1;
        }
    }

}

void Produto(float M1[MAX][MAX], float M2[MAX][MAX], float M3[MAX][MAX]){
    linhasC=linhasA;
    colunasC=colunasV;
    //Testa as matrizes para ver se elas tem o mesmo numero de linhas e colunnas
    if(colunasA==linhasV){
        for(int i=0; i<linhasV; i++){
            for (int j=0; j<linhasV; j++){
                for (int k=0; k<colunasA;k++){
                    M3[i][0]=M3[i][0] +( M1[i][k]*M2[k][j] );}
            }
        }

    }
    //caso as matrizes tenham um numero diferente de linhasA e colunasB ou o contrario, o programa nao efetua o produto e mostra essa mensagem
    if((colunasA != linhasV)){
       cout<<"NAO EH POSSIVEL FAZER O PRODUTO DESSAS MATRIZES ESSAS MATRIZES\n)";
    }

}

void ImprA(){
    cout<<"\n MATRIZ A:\n";
    for (int i=0; i<linhasA;i++){
        for (int j=0; j<colunasA; j++){
            cout<< A[i][j];
            cout<<' ';
        }
        cout<<"\n";
    }
}

void ImprV(){
    cout<<"\n VETOR V:\n";
    for (int i=0; i<linhasV;i++){
        for (int j=0; j<colunasV; j++){
           cout<< V[i][j];
        }
        printf("\n");
    }
}

void ImprC(){
    cout<<"\n MATRIZ C:\n";
    for (int i=0; i<linhasC;i++){
        for (int j=0; j<colunasC; j++){
            cout << C[i][j];
        }
        printf("\n");
    }
}
//---------------------------------------------------------------------------------------------------

void residuo(float Xw[], float C[MAX][MAX],double resid[], int linhasA)
{
    float maximo = 0.0;

    for(int i=0;i<linhasA;i++){
        resid[i] = C[i][0];
        for(int j=0;j<linhasA;j++){
        resid[i] = resid[i] - (A[i][j]*Xw[j]);
        }
    }
}

double modulo (double resid[], int k) //funcao para achar o modulo do residuo para usar no normalinha
{
    if (resid[k] < 0) resid[k] = -resid[k];
    else resid[k] = resid[k];
    return resid[k];
}

double normaLinha(double resid[] ,int linhasA)
{
  double prov_max = 0.0;
   for (int i = 0; i < linhasA; i++)
   {
       resid[i] = modulo(resid,i);
   }
   for (int j = 0; j < linhasA; j++)
   {
        if (resid[j] > prov_max) prov_max = resid[j];
   }
   return prov_max;
}
//+________________________-



void Jacobi(){
    int iteracao;
    cout<<"Digite o numero de iteracoes:";
    cin>>iteracao;
//--------------------------------------
    float prov;
//______________
    for(int p =0;p<linhasA;p++)X[p]=0.0;
    for(int k =0;k<iteracao;k++)
    {
        for(int i =0;i<linhasA;i++)
        {
            Xw[i]=0.0;

            for(int j=0;j<linhasA;j++)
            {
                if(i!=j)
                {

                    Xw[i]=Xw[i]+A[i][j]*X[j];

                }

            }

            Xw[i]=(C[i][0]-Xw[i])/A[i][i];


        }
        for(int gol=0;gol<linhasA;gol++)
        {
            X[gol]=Xw[gol];

           // cout<<X[gol]<<endl;
        }


        //cout<<"iteracao numero:"<<k+1<<endl;
//----------------------------------------------------------------------------

       for (int i; i < linhasA; i++) //calculo do residuo por cada iteracao
       {
           for (int j; j < linhasA; j++)
           {
               resid[i] = resid[i] + (A[i][j] * X[j]);
           }
       }
    residuo(Xw,C,resid,linhasA);
    prov = normaLinha(resid,linhasA);
    cout << "Residuo numero " << k + 1 << " = " <<prov << endl;

//_________________________


   }

}


void GaussSeidel() {
int iteracao,j;
float sig;
    for(int gol=0;gol<linhasA;gol++)X[gol]=0;
    for(int gol=0;gol<linhasA;gol++)Xw[gol]=0;
    cout<<"Digite o numero de iteracoes:";
    cin>>iteracao;
    float prov;

    for(int k =0;k<iteracao;k++)
    {
        for(int i =0;i<linhasA;i++)
        {
            sig = 0;
            for (int j=0;j<i;j++){
                sig=sig+A[i][j]*Xw[j];

            }
            for (j=i+1;j<linhasA;j++){
                sig = sig + A[i][j]*X[j];
            }
            Xw[i] = (C[i][0]-sig)/A[i][i];

        }
        for(int gol=0;gol<linhasA;gol++)
        {
            X[gol]=Xw[gol];

            //cout<<X[gol]<<endl;
        }
    //cout<<"interacao numero:"<<k+1<<endl;

    //----------------------------------------------------------------------------

           for (int i; i < linhasA; i++) //calculo do residuo por cada iteracao
           {
               for (int j; j < linhasA; j++)
               {
                   resid[i] = resid[i] + (A[i][j] * X[j]);
               }
           }
        residuo(Xw,C,resid,linhasA);
        prov = normaLinha(resid,linhasA);
        cout << "Residuo numero " << k + 1 << " = " <<prov<< endl;

    //_________________________
   }

}



int main()
{
  DefTamA();
  DefTamV();
  PreMatrizA();
  PreVetorV();
  //ImprA();
  //ImprV();
  Produto(A,V,C);
  //ImprC();

  int opcao;
  do
  {
      cout<<"\n"<<"=============== MENU ===============\n";
      cout<<"(1) Fazer por Jacobi\n";
      cout<<"(2) Fazer por Gauss-Seidel\n";
      cout<<"(3) Sair\n";
      cout<<"Digite o numero do metodo que quer fazer:";
      cin>>opcao;

      switch(opcao)
      {

      case 1:
          Jacobi();
          cout<<endl;
          break;
      case 2:
          GaussSeidel();
          cout<<endl;
          break;
      case 3:

          break;
      default:
          cout<<"Digite uma opcao valida!";
          cout<<endl;
          break;
      }


  }while(opcao!=3);
  return 0;


}
