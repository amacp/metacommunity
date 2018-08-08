#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <new>
#include <memory>
#include <time.h>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <stdio.h>
#include <iomanip>
#include <math.h>
#include <complex>
#include <vector>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/variance.hpp>
// #include <boost/random.hpp>
// #include <boost/math/distributions/gamma.hpp>
#include <boost/filesystem.hpp>


using namespace std;

//global variables
char filename0[200];
char filename1[200];
char filename2[200];
char filename3[200];

//output files
ofstream out_Pars;
ofstream out_mcom;
//classes
class ind;
class deme;
class mcom
{
    public:
        //Variables
        int gmax,gen,xmax,ymax,nb,nh,np,hTot0,pTot0,rad;
        double**LMtrx,**alphaMtrx,**betaMtrx,**gammaMtrx,r,K,**Amtrx,theta,sigma,ch,cp,rhoh,rhop,xih,xip,*psiWVec,m,dh,dp,epsilonhMax,epsilonpMax,bTot0,energy0,ETot;
        double***censusMtrx,***FherbMtrx,***FpredMtrx;
        vector<ind*>**herbTemp;
        vector<ind*>**predTemp;
        deme** landscape;
        //functions
        mcom();
        void input();
        void initialize();
        void calculateF();
        void census(ofstream* strem);
        void biomassGrowth();
        void herbivory();
        void predation();
        void reproduction();
        void death();
        void movement();
        void movement2();
};
class deme
{
    public:
        //variables
        int x,y;
        double E;
        vector<ind*> herb;
        vector<ind*> pred;
        //functions
        deme();
        void initialize(mcom metaCom,int xin, int yin);
        void census(mcom metaCom,ofstream* stream);
        void predation(mcom metaCom);
        void reproduction(mcom metaCom);
        //void calculatePsi(mcom metaCom);
        void death();
};
class ind
{
    public:
        //Variables
        int age,x,y,species;
            //Trophic: 1 (Herbivores), 2 (Preditors)
        double energy,epsilon,psiMax;
        //functions
        ind();
        ind(mcom metaCom, int xIn, int yIn, int speciesIn,double epsilonMax);
        //ind(int ageIn, int xIn, int yIn, int speciesIn,double epsilonIn, double energyIn);
        void calculatePsi(mcom metaCom,int trophic);
        void calculatePsi2(mcom metaCom,int trophic);
};
int* shuffleList(int length);
int main()
{
    //random seed
    srand((unsigned int)time(NULL));//Initializing the random seed.
    //Setting up files
    sprintf(filename0,"/Users/ailenemacpherson/Documents/VisualStudio/metacommunity");
    //sprintf(filename0,"/Users/lmguzman/Documents/UBC/TrophicMetacommunities/Code/Trophic_metacom_model/C_runs");
    int run=1;
    sprintf(filename1,"%s/run_%d",filename0,run);
    boost::filesystem::create_directories(filename1);
    sprintf(filename2,"%s/out.csv",filename1);
    out_Pars.open(filename2);
    //Variables
    mcom metaCom;
    metaCom.input();
    metaCom.initialize();
    while(metaCom.gen<=metaCom.gmax)
    {
        cout<<"generation: "<<metaCom.gen<<endl;
        metaCom.census(&out_mcom);
        //cout<<"biomass growth:"<<endl;
        metaCom.biomassGrowth();
        //cout<<"herbivory:"<<endl;
        metaCom.herbivory();
        //cout<<"predation:"<<endl;
        metaCom.predation();
        //cout<<"reproduction:"<<endl;
        metaCom.reproduction();
        //cout<<"movement:"<<endl;
        metaCom.movement2();
        //cout<<"death:"<<endl;
        metaCom.death();
        metaCom.gen++;
    }
    metaCom.census(&out_mcom);
    cout<<"DONE!"<<endl;
    return 0;
}
int* shuffleList(int length)
{
	int* list=new int[length];
	int r,temp;
	for(int e=0;e<length;e++) list[e]=e;
	for(int e=0;e<length;e++)
	{
		r=e+rand()%(length-e);
		temp=list[r];
		list[r]=list[e];list[e]=temp;
	}
	return list;
}
mcom::mcom()
{
    LMtrx=NULL;
    gmax=0;xmax=0;ymax=0;nb=0;nh=0;np=0;
}
void mcom::input()
{
    //This function reads parameters in from the input file Input.txt
    int**hInitMtrx,**pInitMtrx;//Temporary matricies for initializing hMtrx and pMtrxx
    ifstream inFile;
	sprintf(filename2,"%s/Input.txt",filename0);//Makes the parameter file name.
	inFile.open(filename2);
	string myString;
	myString="Nothing";   
	while(myString != "(gmax):" && inFile.good())	{inFile>>myString;}
	inFile>>gmax;
	out_Pars<<"gmax:, "<<gmax<<endl; 
	while(myString != "(xmax):" && inFile.good())	{inFile>>myString;}
	inFile>>xmax;
	out_Pars<<"xmax:, "<<xmax<<endl; 
	while(myString != "(ymax):" && inFile.good())	{inFile>>myString;}
	inFile>>ymax;
	out_Pars<<"ymax:, "<<ymax<<endl; 
    while(myString != "(LMtrx):" && inFile.good())	{inFile>>myString;}
    LMtrx=new double*[ymax];
    out_Pars<<"LMtrx: "<<endl;
    for(int y=0;y<ymax;y++)
    {
        LMtrx[y]=new double[xmax];
        for(int x=0;x<xmax;x++)
        {
            inFile>>LMtrx[y][x];
            out_Pars<<LMtrx[y][x]<<",";
        }
        out_Pars<<endl;
    }
	while(myString != "(nb):" && inFile.good())	{inFile>>myString;}
	inFile>>nb;
	out_Pars<<"nb:, "<<nb<<endl; 
	while(myString != "(nh):" && inFile.good())	{inFile>>myString;}
	inFile>>nh;
	out_Pars<<"nh:, "<<nh<<endl; 
	while(myString != "(np):" && inFile.good())	{inFile>>myString;}
	inFile>>np;
	out_Pars<<"np:, "<<np<<endl; 
	while(myString != "(bTot0):" && inFile.good())	{inFile>>myString;}
	inFile>>bTot0;
	out_Pars<<"bTot0:, "<<bTot0<<endl; 
	while(myString != "(hTot0):" && inFile.good())	{inFile>>myString;}
	inFile>>hTot0;
	out_Pars<<"hTot0:, "<<hTot0<<endl; 
 	while(myString != "(hInitMtrx):" && inFile.good())	{inFile>>myString;}
    hInitMtrx=new int*[3];
    out_Pars<<"hInitMtrx: "<<endl;
    for(int r=0;r<3;r++)
    {
        hInitMtrx[r]=new int[hTot0];
        for(int c=0;c<hTot0;c++)
        {
            inFile>>hInitMtrx[r][c];
            out_Pars<<hInitMtrx[r][c]<<",";
        }
        out_Pars<<endl;
    }  
	while(myString != "(pTot0):" && inFile.good())	{inFile>>myString;}
	inFile>>pTot0;
	out_Pars<<"pTot0:, "<<pTot0<<endl; 
 	while(myString != "(pInitMtrx):" && inFile.good())	{inFile>>myString;}
    pInitMtrx=new int*[3];
    out_Pars<<"pInitMtrx: "<<endl;
    for(int r=0;r<3;r++)
    {
        pInitMtrx[r]=new int[pTot0];
        for(int c=0;c<pTot0;c++)
        {
            inFile>>pInitMtrx[r][c];
            out_Pars<<pInitMtrx[r][c]<<",";
        }
        out_Pars<<endl;
    }  
	while(myString != "(energy0):" && inFile.good())	{inFile>>myString;}
	inFile>>energy0;
	out_Pars<<"energy0:, "<<energy0<<endl; 
    while(myString != "(alphaMtrx):" && inFile.good())	{inFile>>myString;}
    alphaMtrx=new double*[nb];
    out_Pars<<"alphaMtrx: "<<endl;
    for(int r=0;r<nb;r++)
    {
        alphaMtrx[r]=new double[nb];
        for(int c=0;c<nb;c++)
        {
            inFile>>alphaMtrx[r][c];
            out_Pars<<alphaMtrx[r][c]<<",";
        }
        out_Pars<<endl;
    }
    while(myString != "(betaMtrx):" && inFile.good())	{inFile>>myString;}
    betaMtrx=new double*[nb];
    out_Pars<<"betaMtrx: "<<endl;
    for(int r=0;r<nb;r++)
    {
        betaMtrx[r]=new double[nh];
        for(int c=0;c<nh;c++)
        {
            inFile>>betaMtrx[r][c];
            out_Pars<<betaMtrx[r][c]<<",";
        }
        out_Pars<<endl;
    }
    while(myString != "(gammaMtrx):" && inFile.good())	{inFile>>myString;}
    gammaMtrx=new double*[nh];
    out_Pars<<"gammaMtrx: "<<endl;
    for(int r=0;r<nh;r++)
    {
        gammaMtrx[r]=new double[np];
        for(int c=0;c<np;c++)
        {
            inFile>>gammaMtrx[r][c];
            out_Pars<<gammaMtrx[r][c]<<",";
        }
        out_Pars<<endl;
    }
	while(myString != "(r):" && inFile.good())	{inFile>>myString;}
	inFile>>r;
	out_Pars<<"r:, "<<r<<endl; 
	while(myString != "(K):" && inFile.good())	{inFile>>myString;}
	inFile>>K;
	out_Pars<<"K:, "<<K<<endl; 
	while(myString != "(theta):" && inFile.good())	{inFile>>myString;}
	inFile>>theta;
	out_Pars<<"theta:, "<<theta<<endl; 
	while(myString != "(sigma):" && inFile.good())	{inFile>>myString;}
	inFile>>sigma;
	out_Pars<<"sigma:, "<<sigma<<endl; 
	while(myString != "(ch):" && inFile.good())	{inFile>>myString;}
	inFile>>ch;
	out_Pars<<"ch:, "<<ch<<endl; 
	while(myString != "(cp):" && inFile.good())	{inFile>>myString;}
	inFile>>cp;
	out_Pars<<"cp:, "<<cp<<endl; 
	while(myString != "(rhoh):" && inFile.good())	{inFile>>myString;}
	inFile>>rhoh;
	out_Pars<<"rhoh:, "<<rhoh<<endl; 
	while(myString != "(rhop):" && inFile.good())	{inFile>>myString;}
	inFile>>rhop;
	out_Pars<<"rhop:, "<<rhop<<endl; 
	while(myString != "(xih):" && inFile.good())	{inFile>>myString;}
	inFile>>xih;
	out_Pars<<"xih:, "<<xih<<endl; 
	while(myString != "(xip):" && inFile.good())	{inFile>>myString;}
	inFile>>xip;
	out_Pars<<"xip:, "<<xip<<endl;
	while(myString != "(m):" && inFile.good())	{inFile>>myString;}
	inFile>>m;
	out_Pars<<"m:, "<<m<<endl; 
	while(myString != "(rad):" && inFile.good())	{inFile>>myString;}
	inFile>>rad;
	out_Pars<<"rad:, "<<rad<<endl; 
	while(myString != "(epsilonhMax):" && inFile.good())	{inFile>>myString;}
	inFile>>epsilonhMax;
	out_Pars<<"epsilonhMax:, "<<epsilonhMax<<endl; 
	while(myString != "(epsilonpMax):" && inFile.good())	{inFile>>myString;}
	inFile>>epsilonpMax;
	out_Pars<<"epsilonpMax:, "<<epsilonpMax<<endl; 
	while(myString != "(dh):" && inFile.good())	{inFile>>myString;}
	inFile>>dh;
	out_Pars<<"dh:, "<<dh<<endl; 
	while(myString != "(dp):" && inFile.good())	{inFile>>myString;}
	inFile>>dp;
	out_Pars<<"dp:, "<<dp<<endl;
    psiWVec=new double[3];
 	while(myString != "(psiWVec):" && inFile.good())	{inFile>>myString;}
    out_Pars<<"psiWVec:, ";
    for(int c=0;c<3;c++)
    {
       inFile>>psiWVec[c];
       out_Pars<<psiWVec[c]<<",";
    }
    out_Pars<<endl; 
    //making census matrix
    censusMtrx=new double**[ymax];
    for(int y=0;y<ymax;y++)
    {
        censusMtrx[y]=new double*[xmax];
        for(int x=0;x<xmax;x++)
        {
            censusMtrx[y][x]=new double[nb+nh+np];
            for(int col=0;col<nb+nh+np;col++)
            {
                censusMtrx[y][x][col]=0;
            }
        }
        
    }
    for(int c=0;c<hTot0;c++) censusMtrx[hInitMtrx[1][c]][hInitMtrx[0][c]][nb+hInitMtrx[2][c]]++;
    for(int c=0;c<pTot0;c++) censusMtrx[pInitMtrx[1][c]][pInitMtrx[0][c]][nb+nh+pInitMtrx[2][c]]++;
    for(int r=0;r<3;r++)
    {
        delete[] hInitMtrx[r];
        delete[] pInitMtrx[r];
    }
}
void mcom::initialize()
{
    //Setting up landscape (also initializes demes).  Initializes biomass, herbivoures and predators.
    gen=0;
    //Making demes and calculating ETot
    landscape=new deme*[ymax];
    ETot=0;
    for(int r=0;r<ymax;r++)
    {
        landscape[r]=new deme[xmax];
        for(int c=0;c<xmax;c++)
        {
            ETot+=LMtrx[r][c];
            landscape[r][c].E=LMtrx[r][c];
            landscape[r][c].x=c;
            landscape[r][c].y=r;
        }
    }
    herbTemp=new vector<ind*>*[ymax];
    predTemp=new vector<ind*>*[ymax];
    for(int y=0;y<ymax;y++)
    {
        herbTemp[y]=new vector<ind*>[xmax];
        predTemp[y]=new vector<ind*>[xmax];
    }
   //Seeding biomass herbivores and preditors
    for(int y=0;y<ymax;y++)
    {
        for(int x=0;x<xmax;x++)
        {
            //seeding biomass
            landscape[y][x].initialize((*this),x,y);
            for(int col=0;col<nh;col++)
            {
                for(int h=0;h<censusMtrx[y][x][nb+col];h++)
                {
                    landscape[y][x].herb.push_back(new ind((*this),x,y,col,epsilonhMax));
                }
            }
            for(int col=0;col<np;col++)
            {
                for(int p=0;p<censusMtrx[y][x][nb+nh+col];p++)
                {
                    landscape[y][x].pred.push_back(new ind((*this),x,y,col,epsilonpMax));
                }
            }
        }
    } 
    out_Pars<<"AMtrx: "<<endl;
    Amtrx=new double*[ymax];
    for(int y=0;y<ymax;y++)
    {
        Amtrx[y]=new double[xmax];
        for(int x=0;x<xmax;x++)
        {
            Amtrx[y][x]=exp(-1.0*pow((LMtrx[y][x]-theta)/sigma,2));
            out_Pars<<Amtrx[y][x]<<",";
        }
        out_Pars<<endl;
    }
    FherbMtrx=new double**[ymax];
    FpredMtrx=new double**[ymax];
    for(int y=0;y<ymax;y++)
    {
        FherbMtrx[y]=new double*[xmax];
        FpredMtrx[y]=new double*[xmax];
        for(int x=0;x<xmax;x++)
        {
            FherbMtrx[y][x]=new double[nh];
            FpredMtrx[y][x]=new double[np];
        }
    }
}
void mcom::calculateF()
{
    double FsquAvg,FAvg,FSD;
    //Herbivore
    for(int h=0;h<nh;h++)
    {
        FsquAvg=0;FAvg=0;
        for(int y=0;y<ymax;y++)
        {
            for(int x=0;x<xmax;x++)
            {

                FherbMtrx[y][x][h]=psiWVec[0]*Amtrx[y][x];
                for(int b=0;b<nb;b++) FherbMtrx[y][x][h]+=psiWVec[1]*betaMtrx[b][h]*censusMtrx[y][x][b];
                for(int p=0;p<landscape[y][x].pred.size();p++) FherbMtrx[y][x][h]-=psiWVec[2]*gammaMtrx[h][landscape[y][x].pred[p]->species];
                FsquAvg+=pow(FherbMtrx[y][x][h],2);
                FAvg+=FherbMtrx[y][x][h];
            }
        }
        FsquAvg/=((double)(ymax*xmax));
        FAvg/=((double)(ymax*xmax));
        FSD=sqrt(FsquAvg-pow(FAvg,2.0));
        for(int y=0;y<ymax;y++)
        {
            for(int x=0;x<xmax;x++)
            {

                FherbMtrx[y][x][h]/=FSD;
            }
        }        
    }
    //Preditor
    for(int p=0;p<np;p++)
    {
        FsquAvg=0;FAvg=0;
        for(int y=0;y<ymax;y++)
        {
            for(int x=0;x<xmax;x++)
            {

                FpredMtrx[y][x][p]=psiWVec[0]*Amtrx[y][x];
                for(int h=0;h<landscape[y][x].herb.size();h++) FpredMtrx[y][x][p]+=(psiWVec[1]+psiWVec[2])*gammaMtrx[landscape[y][x].herb[h]->species][p];
                FsquAvg+=pow(FpredMtrx[y][x][p],2);
                FAvg+=FpredMtrx[y][x][p];
            }
        }
        FsquAvg/=((double)(ymax*xmax));
        FAvg/=((double)(ymax*xmax));
        FSD=sqrt(FsquAvg-pow(FAvg,2.0));
        for(int y=0;y<ymax;y++)
        {
            for(int x=0;x<xmax;x++)
            {

                FpredMtrx[y][x][p]/=FSD;
            }
        }        
    }
}
void mcom::census(ofstream* stream)
{
    //calculating censusMtrx
    for(int y=0;y<ymax;y++)
    {
        for(int x=0;x<xmax;x++)
        {
            for(int col=0;col<nh+np;col++) censusMtrx[y][x][nb+col]=0;
            for(int h=0;h<landscape[y][x].herb.size();h++) censusMtrx[y][x][nb+landscape[y][x].herb[h]->species]++;
            for(int p=0;p<landscape[y][x].pred.size();p++) censusMtrx[y][x][nb+nh+landscape[y][x].pred[p]->species]++;
        }
    }
    //Printing
    sprintf(filename2,"%s/mcom_%d.csv",filename1,gen);
    (*stream).open(filename2);
    for(int y=0;y<ymax;y++)
    {
        for(int x=0;x<xmax;x++)
        {
            (*stream)<<x<<","<<y<<",";
            for(int b=0;b<nb;b++) (*stream)<<censusMtrx[y][x][b]<<",";
            for(int h=0;h<nh;h++) (*stream)<<censusMtrx[y][x][nb+h]<<",";
            for(int p=0;p<np;p++) (*stream)<<censusMtrx[y][x][nb+nh+p]<<",";
            (*stream)<<endl;
        }
        //(*stream)<<endl;
    }
    (*stream).close();
}
void mcom::biomassGrowth()
{
    double* compVec=new double[nb];
    for(int y=0;y<ymax;y++)
    {
        for(int x=0;x<xmax;x++)
        {
            //calculating competition
            for(int b=0;b<nb;b++)
            {
                compVec[b]=0;
                for(int b2=0;b2<nb;b2++)
                {
                    compVec[b]+=alphaMtrx[b][b2]*censusMtrx[y][x][b2];
                }
                compVec[b]/=K;
            }
            for(int b=0;b<nb;b++)
            {
                censusMtrx[y][x][b]=censusMtrx[y][x][b]*(1.0+r*Amtrx[y][x]*(1.0-compVec[b]));
            }
        }
    }
    delete[] compVec;
}
void mcom::herbivory()
{
    double* temp1;
    temp1=new double[nb]; 
    for(int y=0;y<ymax;y++)
    {
        for(int x=0;x<xmax;x++)
        {
            //Calculating consumption
            for(int b=0;b<nb;b++) temp1[b]=0;
            for(int b=0;b<nb;b++)
            {
                for(int h=0;h<landscape[y][x].herb.size();h++)
                {
                    temp1[b]+=betaMtrx[b][landscape[y][x].herb[h]->species];
                    landscape[y][x].herb[h]->energy+=ch*betaMtrx[b][landscape[y][x].herb[h]->species]*censusMtrx[y][x][b]; 
                }
            }
            //Change in biomass 
            for(int b=0;b<nb;b++)
            {
                censusMtrx[y][x][b]-=temp1[b];
                censusMtrx[y][x][b]=max(censusMtrx[y][x][b],0.0);
            }
        }
    }
}
void mcom::predation()
{
    for(int y=0;y<ymax;y++)
    {
        for(int x=0;x<xmax;x++)
        {
            landscape[y][x].predation((*this));
        }
    }
}
void mcom::reproduction()
{
    for(int y=0;y<ymax;y++)
    {
        for(int x=0;x<xmax;x++)
        {
            landscape[y][x].reproduction((*this));
        }
    }    
}
void mcom::death()
{
    for(int y=0;y<ymax;y++)
    {
        for(int x=0;x<xmax;x++)
        {
            landscape[y][x].death();
        }
    }    
}
void mcom::movement()
{
    //ORIGINAL MOVEMENT FUNCTION.  CONTINUS DIMENSION PERCEPTION.
    //calculating where individudals move
    calculateF();
    for(int y=0;y<ymax;y++)
    {
        for(int x=0;x<xmax;x++)
        {
            //cout<<"y "<<y<<" x "<<x<<endl;
            for(int h=0;h<landscape[y][x].herb.size();h++)
            {
                landscape[y][x].herb[h]->calculatePsi((*this),1);
                herbTemp[landscape[y][x].herb[h]->y][landscape[y][x].herb[h]->x].push_back(landscape[y][x].herb[h]);
            }
            for(int p=0;p<landscape[y][x].pred.size();p++)
            {
                landscape[y][x].pred[p]->calculatePsi((*this),2);
                predTemp[landscape[y][x].pred[p]->y][landscape[y][x].pred[p]->x].push_back(landscape[y][x].pred[p]);
            }
            landscape[y][x].herb.clear();
            landscape[y][x].pred.clear();
        }
    }
    //moving individuals
    for(int y=0;y<ymax;y++)
    {
        for(int x=0;x<xmax;x++)
        {
            for(int h=0;h<herbTemp[y][x].size();h++) landscape[y][x].herb.push_back(herbTemp[y][x][h]);
            for(int p=0;p<predTemp[y][x].size();p++) landscape[y][x].pred.push_back(predTemp[y][x][p]);
            herbTemp[y][x].clear();
            predTemp[y][x].clear();
        }
    }
}
void mcom::movement2()
{
    //Movement with a fixed radius
    //calculating where individudals move
    calculateF();
    for(int y=0;y<ymax;y++)
    {
        for(int x=0;x<xmax;x++)
        {
            //cout<<"y "<<y<<" x "<<x<<endl;
            for(int h=0;h<landscape[y][x].herb.size();h++)
            {
                landscape[y][x].herb[h]->calculatePsi2((*this),1);
                herbTemp[landscape[y][x].herb[h]->y][landscape[y][x].herb[h]->x].push_back(landscape[y][x].herb[h]);
            }
            for(int p=0;p<landscape[y][x].pred.size();p++)
            {
                landscape[y][x].pred[p]->calculatePsi2((*this),2);
                predTemp[landscape[y][x].pred[p]->y][landscape[y][x].pred[p]->x].push_back(landscape[y][x].pred[p]);
            }
            landscape[y][x].herb.clear();
            landscape[y][x].pred.clear();
        }
    }
    //moving individuals
    for(int y=0;y<ymax;y++)
    {
        for(int x=0;x<xmax;x++)
        {
            for(int h=0;h<herbTemp[y][x].size();h++) landscape[y][x].herb.push_back(herbTemp[y][x][h]);
            for(int p=0;p<predTemp[y][x].size();p++) landscape[y][x].pred.push_back(predTemp[y][x][p]);
            herbTemp[y][x].clear();
            predTemp[y][x].clear();
        }
    }
}
deme::deme()
{
    y=0;E=0;
}
void deme::initialize(mcom metaCom,int xin, int yin)
{
    //Envionment
    x=xin;y=yin;
    E=metaCom.LMtrx[yin][xin];
    //Filling demes with biomass
    for(int b=0;b<metaCom.nb;b++)
    {
        metaCom.censusMtrx[yin][xin][b]=(metaCom.bTot0/(double)metaCom.nb)*(E/metaCom.ETot);
    }
    //Making arrays
}
void deme::predation(mcom metaCom)
{
    int* temp;
    for(int h=0;h<herb.size();h++)
    {
        temp=shuffleList(pred.size());
        for(int p=0;p<pred.size();p++)
        {
            if(rand()/(double)RAND_MAX<metaCom.gammaMtrx[herb[h]->species][pred[temp[p]]->species])
            {
                //cout<<"here"<<endl;
                //Update predator energy
                pred[temp[p]]->energy+=metaCom.cp;
                //Herbivore is consumed
                //herb[h]->remove(metaCom.ymax);
                herb.erase(herb.begin()+h);
                h--;
                p=pred.size();
            }
        }
        delete[] temp;
    }
}
void deme::reproduction(mcom metaCom)
{
    int temp,num;
    num=herb.size();//Size changes but the number of parents does not
    for(int h=0;h<num;h++)
    {
        temp=herb[h]->energy/metaCom.rhoh;
        for(int i=0;i<temp;i++)
        {
            if(rand()/(double)RAND_MAX<metaCom.xih)
            {
                herb.push_back(new ind(metaCom,x,y,herb[h]->species,metaCom.epsilonhMax));
                herb[h]->energy-=metaCom.rhoh;
            }
        }
    }
    num=pred.size();
    for(int p=0;p<num;p++)
    {
        temp=pred[p]->energy/metaCom.rhop;
        for(int i=0;i<temp;i++)
        {
            if(rand()/(double)RAND_MAX<metaCom.xip)
            {
                pred.push_back(new ind(metaCom,x,y,pred[p]->species,metaCom.epsilonpMax));
                pred[p]->energy-=metaCom.rhop;
            }
        }
    }
}
void deme::death()
{
    for(int h=0;h<herb.size();h++)
    {
        herb[h]->energy--;
        //cout<<"energy: "<<herb[h]->energy;getchar();
        if(herb[h]->energy<=0)
        {
            herb.erase(herb.begin()+h);
            h--;
        }
    }
    for(int p=0;p<pred.size();p++)
    {
        pred[p]->energy--;
        if(pred[p]->energy<=0)
        {
            pred.erase(pred.begin()+p);
            p--;
        }
    }
}
ind::ind(mcom metaCom, int xIn, int yIn, int speciesIn,double epsilonMax)
{
    age=0;
    energy=metaCom.energy0;
    x=xIn;
    y=yIn;
    species=speciesIn;
    epsilon= rand()/(double)RAND_MAX*epsilonMax;
    //A=exp(-1.0*pow((metaCom.LMtrx[y][x]-metaCom.theta)/metaCom.sigma,2));
}
void ind::calculatePsi(mcom metaCom,int trophic)
{
    int xout, yout;
    //cout<<"yin "<<y<<" xin "<<x<<endl;
    double temp;
    if(rand()/(double)RAND_MAX<metaCom.m)
    {
        psiMax=0;
        if(trophic==1)
        {
            for(int yy=0;yy<metaCom.ymax;yy++)
            {
                for(int xx=0;xx<metaCom.xmax;xx++)
                {
                    temp=(metaCom.FherbMtrx[yy][xx][species]+epsilon)*exp(-1.0*metaCom.dh*sqrt(pow(y-yy,2.0)+pow(x-xx,2.0)));
                    //cout<<"ytest "<<yy<<" xtest "<<xx<<" psi "<<temp<<endl;getchar();
                    if(psiMax<temp)
                    {
                        psiMax=temp;
                        xout=xx;yout=yy;
                    }
                }
            }
        }
        else
        {
            for(int yy=0;yy<metaCom.ymax;yy++)
            {
                for(int xx=0;xx<metaCom.xmax;xx++)
                {
                    temp=(metaCom.FpredMtrx[yy][xx][species]+epsilon)*exp(-1.0*metaCom.dp*sqrt(pow(y-yy,2.0)+pow(x-xx,2.0)));
                    //cout<<"ytest "<<yy<<" xtest "<<xx<<" psi "<<temp<<endl;getchar();
                    if(psiMax<temp)
                    {
                        psiMax=temp;
                        xout=xx;yout=yy;
                    }
                }
            }
        }
        x=xout;y=yout;
    }
    // cout<<"yout "<<y<<" xout "<<x<<endl;
    // getchar();
}
void ind::calculatePsi2(mcom metaCom,int trophic)
{
    int xout, yout;
    //cout<<"yin "<<y<<" xin "<<x<<endl;
    double temp,dist;
    if(rand()/(double)RAND_MAX<metaCom.m)
    {
        psiMax=0;
        if(trophic==1)
        {
            for(int yy=0;yy<metaCom.ymax;yy++)
            {
                for(int xx=0;xx<metaCom.xmax;xx++)
                {
                    dist=sqrt(pow(y-yy,2.0)+pow(x-xx,2.0));
                    if(dist<=metaCom.rad)
                    {
                        temp=(metaCom.FherbMtrx[yy][xx][species]+epsilon)*exp(-1.0*metaCom.dh*dist);
                        //cout<<"ytest "<<yy<<" xtest "<<xx<<" psi "<<temp<<endl;getchar();
                        if(psiMax<temp)
                        {
                            psiMax=temp;
                            xout=xx;yout=yy;
                        }
                    }
                }
            }
        }
        else
        {
            for(int yy=0;yy<metaCom.ymax;yy++)
            {
                for(int xx=0;xx<metaCom.xmax;xx++)
                {
                    dist=sqrt(pow(y-yy,2.0)+pow(x-xx,2.0));
                    if(dist<=metaCom.rad)
                    {
                        temp=(metaCom.FpredMtrx[yy][xx][species]+epsilon)*exp(-1.0*metaCom.dp*dist);
                        //cout<<"ytest "<<yy<<" xtest "<<xx<<" psi "<<temp<<endl;getchar();
                        if(psiMax<temp)
                        {
                            psiMax=temp;
                            xout=xx;yout=yy;
                        }
                    }
                }
            }
        }
        x=xout;y=yout;
    }
    // cout<<"yout "<<y<<" xout "<<x<<endl;
    // getchar();
}