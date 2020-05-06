//2 : mass fitting
#include <sstream>
#include <iostream>
#include <iomanip>
#include <string>

using namespace std;

TH1D* h_fit_ori = 0;
TH1D* h_fit = 0;

void fcn(int &npar, double *gin, double &f, double *par, int iflag);
double model(double x, double* par);
std::string to_string_with_precision(const double& a_value, const int n = 6);

//TH1F* h_fit_signal = NULL;

void mass_fit()
{
    string opt = "h_Mjjb";
    
    TFile* f1 = new TFile("Mass_Hist.root");
    f1->GetObject((char*)opt.c_str(),h_fit_ori);

    double X_M = 500.;
    double X_m = 0.;

    h_fit = new TH1D("h_fit","",int(X_M-X_m),X_m,X_M);

    for(int i=1;i<(int)h_fit->GetNbinsX();++i)
    {
        h_fit->SetBinContent( i, h_fit_ori->GetBinContent( (int)X_m + i ) );
    }
    
    TMinuit* gMinuit = new TMinuit(5);
    gMinuit->SetFCN(fcn);
    
    int NPar = 24;

    gMinuit->DefineParameter(0,"g1_no",500000.,1.,1.,1000000000.);
    gMinuit->DefineParameter(1,"g1_mean",170,1.,0.,600.);
    gMinuit->DefineParameter(2,"g1_sigma",10,0.5,0.,200.);

    gMinuit->DefineParameter(3,"g2_no",50000.,1.,1.,1000000000.);
    gMinuit->DefineParameter(4,"g2_mean",300,1.,0.,600.);
    gMinuit->DefineParameter(5,"g2_sigma",20,1.,0.,200.);

    gMinuit->DefineParameter(6,"g3_no",5000.,1.,1.,1000000000.);
    gMinuit->DefineParameter(7,"g3_mean",200,1.,0.,600.);
    gMinuit->DefineParameter(8,"g3_sigma",20,1.,0.,200.);

    gMinuit->DefineParameter(9,"g4_no",5000.,1.,1.,1000000000.);
    gMinuit->DefineParameter(10,"g4_mean",150,1.,0.,600.);
    gMinuit->DefineParameter(11,"g4_sigma",20,1.,0.,200.);

    gMinuit->DefineParameter(12,"g5_no",5000.,1.,1.,1000000000.);
    gMinuit->DefineParameter(13,"g5_mean",200,1.,0.,600.);
    gMinuit->DefineParameter(14,"g5_sigma",20,1.,0.,200.);

    gMinuit->DefineParameter(15,"g6_no",5000.,1.,1.,1000000000.);
    gMinuit->DefineParameter(16,"g6_mean",150,1.,0.,600.);
    gMinuit->DefineParameter(17,"g6_sigma",20,1.,0.,200.);

    gMinuit->DefineParameter(18,"g7_no",5000.,1.,1.,1000000000.);
    gMinuit->DefineParameter(19,"g7_mean",150,1.,0.,600.);
    gMinuit->DefineParameter(20,"g7_sigma",20,1.,0.,200.);

    gMinuit->DefineParameter(21,"g8_no",5000.,1.,1.,1000000000.);
    gMinuit->DefineParameter(22,"g8_mean",150,1.,0.,600.);
    gMinuit->DefineParameter(23,"g8_sigma",20,1.,0.,200.);

  
    gMinuit->Command("MIGRAD");
    gMinuit->Command("MIGRAD");
    //gMinuit->Command("MINOS");        
    //calling MINOS , the assymmetric error will be obtained
    
    double par[NPar],err[NPar];
    for(int i=0;i<NPar;i++)
    {   gMinuit->GetParameter(i,par[i],err[i]); }        
    //take the parameter and error after the fitting
    
    TH1D* curve = new TH1D("curve","curve",(h_fit->GetNbinsX())*5,\
                           h_fit->GetXaxis()->GetXmin(),h_fit->GetXaxis()->GetXmax());
    
    //cout << par[0] << " " << par[1] << " " << par[2] << endl;

    for(int i=1;i<curve->GetNbinsX();i++)
    {
        double x = curve->GetBinCenter(i);
        double f = model(x,par);        //now the par is the parameter after fitting
        double BinWidth = curve->GetBinWidth(1);
        //curve->SetBinContent(i,f*BinWidth);
        //cout << "Bin " << curve->GetBinCenter(i) << ", " << f << endl;
        curve->SetBinContent(i,f);
    }
    
    //--- bkg ---//

    TH1D* h_bkg = new TH1D("h_bkg","",int(X_M-X_m)*10,X_m,X_M);

    for(int i=1;i<(int)h_bkg->GetNbinsX();++i)
    {
        double x = h_bkg->GetBinCenter(i);

        double norm3 = 1./sqrt(2.*TMath::Pi())/par[8];
        double G3 = norm3*exp(-0.5 * (x-par[7])*(x-par[7])/par[8]/par[8]);

        double norm5 = 1./sqrt(2.*TMath::Pi())/par[14];
        double G5 = norm5*exp(-0.5 * (x-par[13])*(x-par[13])/par[14]/par[14]);

        double norm6 = 1./sqrt(2.*TMath::Pi())/par[17];
        double G6 = norm6*exp(-0.5 * pow((x-par[16])/par[17],2));
        
        double norm8 = 1./sqrt(2.*TMath::Pi())/par[23];
        double G8 = norm8*exp(-0.5 * pow((x-par[22])/par[23],2));
        
        double f = par[6] * G3 + par[12] * G5 +  par[15] * G6 + par[21] * G8;
        h_bkg->SetBinContent(i,f);
    }

    //--- signal ---//

    TH1D* h_sig = new TH1D("h_sig","",int(X_M-X_m)*10,X_m,X_M);

    for(int i=1;i<(int)h_sig->GetNbinsX();++i)
    {
        double x = h_sig->GetBinCenter(i);
        
        double norm3 = 1./sqrt(2.*TMath::Pi())/par[8];
        double G3 = norm3*exp(-0.5 * (x-par[7])*(x-par[7])/par[8]/par[8]);
        double norm5 = 1./sqrt(2.*TMath::Pi())/par[14];
        double G5 = norm5*exp(-0.5 * (x-par[13])*(x-par[13])/par[14]/par[14]);
        double norm6 = 1./sqrt(2.*TMath::Pi())/par[17];
        double G6 = norm6*exp(-0.5 * pow((x-par[16])/par[17],2));
        double norm8 = 1./sqrt(2.*TMath::Pi())/par[23];
        double G8 = norm8*exp(-0.5 * pow((x-par[22])/par[23],2));
        
        double f = model(x,par) - (par[6] * G3 + par[12] * G5 + par[15] * G6 + par[21] * G8);
    
        h_sig->SetBinContent(i,f);
    }

    //--- Do the mean and width calculation ---//

    int bin_M = h_sig->GetMaximumBin();
    double mean_signal = h_sig->GetBinCenter(bin_M);
    double width_signal = 0.;

    cout << "Mean" << h_sig->GetMean() << endl;
    
    double binW = h_sig->GetBinWidth(1);
    //double all_integral = h_sig->Integral();
    double all_integral = 0.;
    for(int i=1;i<(int)h_sig->GetNbinsX();++i)
    {
        all_integral += h_sig->GetBinContent(i)*binW;
    }

    //cout << "all " << all_integral << endl;
    double acc_integral = h_sig->GetBinContent(bin_M)*binW;
    for(int i=1;i<(int)h_sig->GetNbinsX();++i)
    {

        acc_integral += h_sig->GetBinContent(bin_M+i)*binW;
        acc_integral += h_sig->GetBinContent(bin_M-i)*binW;
        //cout << "acc " << acc_integral << endl;
        if( acc_integral/all_integral >= 0.682 ) {
            width_signal = h_sig->GetBinCenter(bin_M+i) - mean_signal;
            break;
        }
    }
    

    //--- Draw ---//

    TGaxis::SetMaxDigits(4);

    h_fit->SetMaximum( h_fit->GetMaximum()*1.25 );
    h_fit->GetXaxis()->SetTitle("M_{jjb} [GeV]");
    h_fit->GetYaxis()->SetTitle("Events No.");
    
    curve->SetLineWidth(2);
    
    //h_bkg->SetLineStyle(9);      // TAttLine.h
    h_bkg->SetLineColor(50);  //red
    h_bkg->SetFillColor(50);  //red
    h_bkg->SetFillStyle(3005);
    h_bkg->SetLineWidth(1);

    h_sig->SetLineColor(9);
    h_sig->SetFillColor(9);
    h_sig->SetFillStyle(3004);
    h_sig->SetLineWidth(1);

    string text_mean = "#mu : " + to_string_with_precision(mean_signal,2);
    string text_width = "#sigma : " + to_string_with_precision(width_signal,2);

    auto leg = new TLegend(0.6,0.6,0.8,0.8);
    leg->SetNColumns(2);
    leg->AddEntry(curve,"Fitting line","l");
    leg->AddEntry((TObject*)0,"","");
    leg->AddEntry(h_sig,"Reco top","F");
    leg->AddEntry((TObject*)0,text_mean.c_str(),"");
    leg->AddEntry((TObject*)0,"","");
    leg->AddEntry((TObject*)0,text_width.c_str(),"");
    leg->AddEntry(h_bkg,"Reco Bkg.","F");
    leg->AddEntry((TObject*)0,"","");

    leg->SetTextFont( 43 );
    leg->SetTextSize( 14 );
    leg->SetBorderSize( 0 );
    

    TCanvas* c = new TCanvas();
    h_fit->SetStats(false);
    h_fit->Draw();
    
    curve->Draw( "csame" );
    h_bkg->Draw("SAME");
    h_sig->Draw("SAME");
    leg->Draw("SAME");

    string pdf_name = "Plot_" + opt + ".pdf"; 
    c->Print( pdf_name.c_str() );

}

double model(double x, double* par)
{
    double g1_n = par[0];
    double g1_mean = par[1];
    double g1_sigma = par[2];
    double norm1 = 1./sqrt(2.*TMath::Pi())/g1_sigma;
    double G1 = norm1*exp(-0.5 * pow((x-g1_mean)/g1_sigma, 2));

//    return g1_n * G1;         //don't times binwidth because it's an unbinned fit

    double g2_n = par[3];
    double g2_mean = par[4];
    double g2_sigma = par[5];
    double norm2 = 1./sqrt(2.*TMath::Pi())/g2_sigma;
    double G2 = norm2*exp(-0.5 * (x-g2_mean)*(x-g2_mean)/g2_sigma/g2_sigma);
//    return g1_n * G1 + g2_n * G2;

    double g3_n = par[6];
    double g3_mean = par[7];
    double g3_sigma = par[8];
    double norm3 = 1./sqrt(2.*TMath::Pi())/g3_sigma;
    double G3 = norm3*exp(-0.5 * (x-g3_mean)*(x-g3_mean)/g3_sigma/g3_sigma);
//    return g1_n * G1 + g2_n * G2 + g3_n * G3;

    double g4_n = par[9];
    double g4_mean = par[10];
    double g4_sigma = par[11];
    double norm4 = 1./sqrt(2.*TMath::Pi())/g4_sigma;
    double G4 = norm4*exp(-0.5 * (x-g4_mean)*(x-g4_mean)/g4_sigma/g4_sigma);
//    return g1_n * G1 + g2_n * G2 + g3_n * G3 + g4_n * G4;

    double g5_n = par[12];
    double g5_mean = par[13];
    double g5_sigma = par[14];
    double norm5 = 1./sqrt(2.*TMath::Pi())/g5_sigma;
    double G5 = norm5*exp(-0.5 * (x-g5_mean)*(x-g5_mean)/g5_sigma/g5_sigma);
//    return g1_n * G1 + g2_n * G2 + g3_n * G3 + g4_n * G4 + g5_n * G5;

    double g6_n = par[15];
    double g6_mean = par[16];
    double g6_sigma = par[17];
    double norm6 = 1./sqrt(2.*TMath::Pi())/g6_sigma;
    double G6 = norm6*exp(-0.5 * (x-g6_mean)*(x-g6_mean)/g6_sigma/g6_sigma);
//    return g1_n * G1 + g2_n * G2 + g3_n * G3 + g4_n * G4 + g5_n * G5 + g6_n * G6;

    double g7_n = par[18];
    double g7_mean = par[19];
    double g7_sigma = par[20];
    double norm7 = 1./sqrt(2.*TMath::Pi())/g7_sigma;
    double G7 = norm7*exp(-0.5 * (x-g7_mean)*(x-g7_mean)/g7_sigma/g7_sigma);
//    return g1_n * G1 + g2_n * G2 + g3_n * G3 + g4_n * G4 + g5_n * G5 + g6_n * G6 + + g7_n * G7;

    double g8_n = par[21];
    double g8_mean = par[22];
    double g8_sigma = par[23];
    double norm8 = 1./sqrt(2.*TMath::Pi())/g8_sigma;
    double G8 = norm8*exp(-0.5 * (x-g8_mean)*(x-g8_mean)/g8_sigma/g8_sigma);

    return g1_n * G1 + g2_n * G2 + g3_n * G3 + g4_n * G4 + g5_n * G5 + g6_n * G6 + + g7_n * G7 + g8_n * G8;

}


void fcn(int &npar, double *gin, double &f, double *par, int iflag)
{
    f = 0.;
    
    double g1_n = par[0];
    //f = 2. * (g1_n);


    double g2_n = par[3];
    //f = 2. * (g1_n + g2_n);
    double g3_n = par[6];
    //f = 2. * (g1_n + g2_n + g3_n);
    double g4_n = par[9];
    //f = 2. * (g1_n + g2_n + g3_n + g4_n);
    double g5_n = par[12];
    //f = 2. * (g1_n + g2_n + g3_n + g4_n + g5_n);
    double g6_n = par[15];
    //f = 2. * (g1_n + g2_n + g3_n + g4_n + g5_n + g6_n);
    double g7_n = par[18];
    //f = 2. * (g1_n + g2_n + g3_n + g4_n + g5_n + g6_n + g7_n);
    double g8_n = par[21];
    f = 2. * (g1_n + g2_n + g3_n + g4_n + g5_n + g6_n + g7_n + g8_n);

    
    for(int i=1;i<=h_fit->GetNbinsX();i++)
    {
        double x = h_fit->GetBinCenter(i);
        double weight = h_fit->GetBinContent(i);    //binned likelihood fit
        double L = model(x,par);
        if(L>0.){ f -= 2.*log(L)*weight;}
        else {f = HUGE; return;}
    }
}

std::string to_string_with_precision(const double& a_value, const int n = 6)
{
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return out.str();
}
