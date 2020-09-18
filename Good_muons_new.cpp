#include <iostream>
#include "Riostream.h"
#include <string>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <math.h>
#include "TTree.h"
#include "TBranch.h"
#include "TFrame.h"
#include "TCanvas.h"
#include "TPaveLabel.h"
#include "TPaveText.h"
#include "TFile.h"
#include "TString.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include "math.h"
#include "time.h"
#include "TRandom.h"
#include "TSpectrum.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include <algorithm>
#include "TLine.h"

using namespace std;

int Good_muons_new(){   
    


    Float_t selMu0_eta;
    Float_t genMu0_eta;
    Float_t selMu0_phi;
    Float_t genMu0_phi;

    Float_t selMu1_eta;
    Float_t genMu1_eta;
    Float_t selMu1_phi;
    Float_t genMu1_phi;

    Float_t selMu2_eta;
    Float_t genMu2_eta;
    Float_t selMu2_phi;
    Float_t genMu2_phi;

    Float_t selMu3_eta;
    Float_t genMu3_eta;
    Float_t selMu3_phi;
    Float_t genMu3_phi;

   
    Float_t genMu0_pT;
    Float_t genMu1_pT;
    Float_t genMu2_pT;
    Float_t genMu3_pT;

    Float_t selMu0_pT;
    Float_t selMu1_pT;
    Float_t selMu2_pT;
    Float_t selMu3_pT;

    Bool_t is1GenMuHighPt;
    Bool_t is2GenMuHighPt;
    Bool_t is3GenMuLowPt;
    Bool_t is4GenMuLowPt;


    Int_t counter1;
    counter1 = 0;
    Int_t counter2;
    counter2 =0;
    Int_t counter3;
    counter3 =0;
    Int_t counter4;
    counter4 = 0;
    Int_t counter5;
    counter5 = 0;   

    Int_t counter6;
    counter6 = 0; 
    Int_t counter7;
    counter7 = 0; 
    Int_t counter8;
    counter8 = 0; 
    Int_t counter9;
    counter9 = 0; 
    Int_t counter10;
    counter10 = 0;
    Int_t counter11;
    counter11 = 0; 
    Int_t counter12;
    counter12 = 0; 
    Int_t counter13;
    counter13 = 0; 
    Int_t counter14;
    counter14 = 0; 
    Int_t counter15;
    counter15 = 0; 
    Int_t counter16;
    counter16 = 0; 
    Int_t counter17;
    counter17 = 0; 
    Int_t counter18;
    counter18 = 0; 
    Int_t counter19;
    counter19 = 0;
    Int_t counter20;
    counter20 = 0;
    Int_t counter21;
    counter21 = 0;

    TFile *myFile = new TFile("ALL_ana_700.root");

    // TCanvas *cnv1 = new TCanvas();
    // TH1F *sD_Lxy1 = new TH1F("#DeltaR_{0}", "#DeltaR_{0}", 100, -0.01, 0.01);

    // TCanvas *cnv2 = new TCanvas();
    // TH1F *sD_Lxy2 = new TH1F("#DeltaR_{1}", "#DeltaR_{1}", 100, -0.01, 0.01);

    // TCanvas *cnv3 = new TCanvas();
    // TH1F *sD_Lxy3 = new TH1F("#DeltaR_{2}", "#DeltaR_{2}", 100, -0.01, 0.01);

    // TCanvas *cnv4 = new TCanvas();
    // TH1F *sD_Lxy4 = new TH1F("#DeltaR_{3}", "#DeltaR_{3}", 100, -0.01, 0.01);

    // sD_Lxy1->SetXTitle("#DeltaR (best match for sel#mu_{0})");
    // sD_Lxy1->SetYTitle("Number of Events");
    // sD_Lxy2->SetXTitle("#DeltaR (best match for sel#mu_{1})");
    // sD_Lxy2->SetYTitle("Number of Events ");
    // sD_Lxy3->SetXTitle("#DeltaR (best match for sel#mu_{2})");
    // sD_Lxy3->SetYTitle("Number of Events");
    // sD_Lxy4->SetXTitle("#DeltaR (best match for sel#mu_{3})");
    // sD_Lxy4->SetYTitle("Number of Events");

    gDirectory->cd("cutFlowAnalyzerPXBL4PXFL3;1");
    gDirectory->pwd();
    TTree *myTree1 = nullptr;
    gDirectory->GetObject("Events;1",myTree1);

    //Number of entries
    int N = myTree1->GetEntries();
    cout << "Number of entries: " << N << endl;

   
    



    myTree1->SetBranchAddress("genMu0_eta", &genMu0_eta);
    myTree1->SetBranchAddress("selMu0_eta", &selMu0_eta);
    myTree1->SetBranchAddress("genMu1_eta", &genMu1_eta);
    myTree1->SetBranchAddress("selMu1_eta", &selMu1_eta);
    myTree1->SetBranchAddress("genMu2_eta", &genMu2_eta);
    myTree1->SetBranchAddress("selMu2_eta", &selMu2_eta);
    myTree1->SetBranchAddress("genMu3_eta", &genMu3_eta);
    myTree1->SetBranchAddress("selMu3_eta", &selMu3_eta);

    myTree1->SetBranchAddress("genMu0_phi", &genMu0_phi);
    myTree1->SetBranchAddress("selMu0_phi", &selMu0_phi);
    myTree1->SetBranchAddress("genMu1_phi", &genMu1_phi);
    myTree1->SetBranchAddress("selMu1_phi", &selMu1_phi);
    myTree1->SetBranchAddress("genMu2_phi", &genMu2_phi);
    myTree1->SetBranchAddress("selMu2_phi", &selMu2_phi);
    myTree1->SetBranchAddress("genMu3_phi", &genMu3_phi);
    myTree1->SetBranchAddress("selMu3_phi", &selMu3_phi);

    myTree1->SetBranchAddress("genMu0_pT", &genMu0_pT);
    myTree1->SetBranchAddress("genMu1_pT", &genMu1_pT);
    myTree1->SetBranchAddress("genMu2_pT", &genMu2_pT);
    myTree1->SetBranchAddress("genMu3_pT", &genMu3_pT);

    myTree1->SetBranchAddress("selMu0_pT", &selMu0_pT);
    myTree1->SetBranchAddress("selMu1_pT", &selMu1_pT);
    myTree1->SetBranchAddress("selMu2_pT", &selMu2_pT);
    myTree1->SetBranchAddress("selMu3_pT", &selMu3_pT);

    myTree1->SetBranchAddress("is1GenMuHighPt", &is1GenMuHighPt);
    myTree1->SetBranchAddress("is2GenMuHighPt", &is2GenMuHighPt);
    myTree1->SetBranchAddress("is3GenMuLowPt", &is3GenMuLowPt);
    myTree1->SetBranchAddress("is4GenMuLowPt", &is4GenMuLowPt);



for(int ii = 0; ii < N; ii++){
    myTree1->GetEntry(ii);

// if ( is1GenMuHighPt == true ) { counter1=counter1+1; }
// if ( is1GenMuHighPt == true ) { counter2=counter2+1; }
// if ( is3GenMuLowPt == true ) { counter3=counter3+1; }
if ( is1GenMuHighPt == true && is2GenMuHighPt == true && is3GenMuLowPt == true && is4GenMuLowPt == true ) { 

            counter21=counter21+1;

            Float_t a, b, c, d, min1;
            a = ( (selMu0_eta - genMu0_eta) * (selMu0_eta - genMu0_eta) ) +  ( (selMu0_phi - genMu0_phi) * (selMu0_phi - genMu0_phi) ) ;
            b = ( (selMu0_eta - genMu1_eta) * (selMu0_eta - genMu1_eta) ) +  ( (selMu0_phi - genMu1_phi) * (selMu0_phi - genMu1_phi) ) ;
            c = ( (selMu0_eta - genMu2_eta) * (selMu0_eta - genMu2_eta) ) +  ( (selMu0_phi - genMu2_phi) * (selMu0_phi - genMu2_phi) ) ;
            d = ( (selMu0_eta - genMu3_eta) * (selMu0_eta - genMu3_eta) ) +  ( (selMu0_phi - genMu3_phi) * (selMu0_phi - genMu3_phi) ) ;
            min1= min(min(min(a,b), c),d); 
            //sD_Lxy1->Fill(sqrt(min1));



            Float_t e, f, g, h, min2;
            e = ( (selMu1_eta - genMu0_eta) * (selMu1_eta - genMu0_eta) ) +  ( (selMu1_phi - genMu0_phi) * (selMu1_phi - genMu0_phi) ) ;
            f = ( (selMu1_eta - genMu1_eta) * (selMu1_eta - genMu1_eta) ) +  ( (selMu1_phi - genMu1_phi) * (selMu1_phi - genMu1_phi) ) ;
            g = ( (selMu1_eta - genMu2_eta) * (selMu1_eta - genMu2_eta) ) +  ( (selMu1_phi - genMu2_phi) * (selMu1_phi - genMu2_phi) ) ;
            h = ( (selMu1_eta - genMu3_eta) * (selMu1_eta - genMu3_eta) ) +  ( (selMu1_phi - genMu3_phi) * (selMu1_phi - genMu3_phi) ) ;
            min2 = std::min(std::min(std::min(e,f), g),h);
            //sD_Lxy2->Fill(sqrt(min2));            


            Float_t i, j, k, l, min3;
            i = ( (selMu2_eta - genMu0_eta) * (selMu2_eta - genMu0_eta) ) +  ( (selMu2_phi - genMu0_phi) * (selMu2_phi - genMu0_phi) ) ;
            j = ( (selMu2_eta - genMu1_eta) * (selMu2_eta - genMu1_eta) ) +  ( (selMu2_phi - genMu1_phi) * (selMu2_phi - genMu1_phi) ) ;
            k = ( (selMu2_eta - genMu2_eta) * (selMu2_eta - genMu2_eta) ) +  ( (selMu2_phi - genMu2_phi) * (selMu2_phi - genMu2_phi) ) ;
            l = ( (selMu2_eta - genMu3_eta) * (selMu2_eta - genMu3_eta) ) +  ( (selMu2_phi - genMu3_phi) * (selMu2_phi - genMu3_phi) ) ;
            min3 = std::min(std::min(std::min(i,j), k),l);
            //sD_Lxy3->Fill(sqrt(min3));            

 

            Float_t m, n, o, p, min4;
            m = ( (selMu3_eta - genMu0_eta) * (selMu3_eta - genMu0_eta) ) +  ( (selMu3_phi - genMu0_phi) * (selMu3_phi - genMu0_phi) ) ;
            n = ( (selMu3_eta - genMu1_eta) * (selMu3_eta - genMu1_eta) ) +  ( (selMu3_phi - genMu1_phi) * (selMu3_phi - genMu1_phi) ) ;
            o = ( (selMu3_eta - genMu2_eta) * (selMu3_eta - genMu2_eta) ) +  ( (selMu3_phi - genMu2_phi) * (selMu3_phi - genMu2_phi) ) ;
            p = ( (selMu3_eta - genMu3_eta) * (selMu3_eta - genMu3_eta) ) +  ( (selMu3_phi - genMu3_phi) * (selMu3_phi - genMu3_phi) ) ;
            min4 = std::min(std::min(std::min(m,n), o),p);
            //sD_Lxy4->Fill(sqrt(min4));           


            if ( min1 == a && min1 < e && min1 < i && min1 < m) {
            counter1 = counter1+1;
            }

            if (min1 == b && min1 < f &&  min1 < j &&  min1 < n) {
            counter2 = counter2+1;
            // cout<<"genMu0_pT - selMu0_pT: "<<genMu0_pT - selMu0_pT <<endl;
            // cout<<"genMu0_pT - selMu1_pT: "<<genMu0_pT - selMu1_pT <<endl;
            // cout<<"genMu0_eta - selMu0_eta: "<<genMu0_eta - selMu0_eta <<endl;
            // cout<<"genMu0_eta - selMu1_eta: "<<genMu0_eta - selMu1_eta <<endl;
            // cout<<"genMu0_phi - selMu0_phi: "<<genMu0_phi - selMu0_phi <<endl;
            // cout<<"genMu0_phi - selMu1_phi: "<<genMu0_phi - selMu1_phi <<endl<<endl;           

            }

            if (min1 == c  && min1 < g  && min1 < k  && min1 < o) {
            counter3 = counter3+1;
            
            }

            if (min1 == d  && min1 < h  && min1 < l  && min1 < p) {
            counter4 = counter4+1;
            
            }

            // else if ( min1 != a && min1 != b && min1 != c && min1 != d) { 
            // counter5 = counter5+1;
            // }



            
            if ( min2 == e && min2 < a && min2 < i && min2 < m) {
            counter6 = counter6+1;
            // cout<<"genMu0_pT - selMu1_pT: "<<genMu0_pT - selMu1_pT <<endl;
            // cout<<"genMu1_pT - selMu1_pT: "<<genMu1_pT - selMu1_pT <<endl;
            // cout<<"genMu0_eta - selMu1_eta: "<<genMu0_eta - selMu1_eta <<endl;
            // cout<<"genMu1_eta - selMu1_eta: "<<genMu1_eta - selMu1_eta <<endl;
            // cout<<"genMu0_phi - selMu1_phi: "<<genMu0_phi - selMu1_phi <<endl;
            // cout<<"genMu1_phi - selMu1_phi: "<<genMu1_phi - selMu1_phi <<endl<<endl;  
            }

            if (min2 == f  && min2 < b &&  min2 < j &&  min2 < n){
            counter7 = counter7+1;
            }

            if (min2 == g && min2 < c  && min2 < k  && min2 < o) {
            counter8 = counter8+1;
            
            }

            if (min2 == h && min2 < d  && min2 < l  && min2 < p) {
            counter9 = counter9+1;
            
            }

            // else if ( min2 != e && min2 != f && min2 != g && min2 != h) { 
            // counter10 = counter10+1;
            // }    


            if ( min3 == i  && min3 < a && min3 < e && min3 < m) {
            counter11 = counter11+1;
            }

            if (min3 == j && min3 < b &&  min3 < f &&  min3 < n) {
            counter12 = counter12+1;
            }

            if (min3 == k && min3 < c  && min3 < g  && min3 < o) {
            counter13 = counter13+1;
            
            }

            if (min3 == l && min3 < d  && min3 < h  && min3 < p) {
            counter14 = counter14+1;
            
            }

            // else if ( min3 != i && min3 != j && min3 != k && min3 != l) { 
            // counter15 = counter15+1;
            // }  


            if ( min4 == m && min4 < a && min4 < e && min4 < i) {
            counter16 = counter16+1;
            }

            if (min4 == n && min4 < b &&  min4 < f &&  min4 < j) {
            counter17 = counter17+1;
            }

            if (min4 == o && min4 < c  && min4 < g  && min4 < k) {
            counter18 = counter18+1;
            
            }

            if (min4 == p && min4 < d  && min4 < h  && min4 < l) {
            counter19 = counter19+1;
            
            }

            // else if ( min4 != m && min4 != n && min4 != o && min4 != p) { 
            // counter20 = counter20+1;
            // } 
            else {counter20=counter20+1;} 
     }
}
    cout << "Number of events with 4 GOOD muons: " << counter21 << endl;

    cout << "leading sel muon matched to leading gen moun:  " << counter1 << endl;
    cout << "leading sel muon matched to sub-leading gen moun:  " << counter2 << endl;
    cout << "leading sel muon matched to 3rd gen moun:  " << counter3 << endl;
    cout << "leading sel muon matched to 4th gen moun:  " << counter4 << endl<< endl;
  //  cout << "leading sel moun not matched: " << counter5 << endl;

    cout << "sub-leading sel muon matched to leading gen moun:  " << counter6 << endl;
    cout << "sub-leading sel muon matched to sub-leading gen moun:  " << counter7 << endl;
    cout << "sub-leading sel muon matched to 3rd gen moun:  " << counter8 << endl;
    cout << "sub-leading sel muon matched to 4th gen moun:  " << counter9 << endl << endl;
 //   cout << "sub-leading sel moun not matched: " << counter10 << endl;


    cout << "3rd sel muon matched to leading gen moun:  " << counter11 << endl;
    cout << "3rd sel muon matched to sub-leading gen moun:  " << counter12 << endl;
    cout << "3rd sel muon matched to 3rd gen moun:  " << counter13 << endl;
    cout << "3rd sel muon matched to 4th gen moun:  " << counter14 << endl<<endl;
  //  cout << "3rd sel moun not matched: " << counter15 << endl;


    cout << "4th sel muon matched to leading gen moun:  " << counter16 << endl;
    cout << "4th sel muon matched to sub-leading gen moun:  " << counter17 << endl;
    cout << "4th sel muon matched to 3rd gen moun:  " << counter18 << endl;
    cout << "4th sel muon matched to 4th gen moun:  " << counter19 << endl<<endl;
    cout << "Number of un-matched mouns: " << counter20 << endl<<endl;

 return 0;

}
