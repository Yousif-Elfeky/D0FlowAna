#ifndef MYFUNCTION_H_INCLUDED
#define MYFUNCTION_H_INCLUDED
void globalSetting()
{
    gStyle->Reset("plain");
    gStyle->SetPalette(1,0);  //1--for colored TH2F
    gStyle->SetFillColor(0); //
    gStyle->SetOptStat(0);   //0--not show histograph para, 1--show histo para
    gStyle->SetOptFit(000);  //000--not show fit para, 111--show fit para.
    gStyle->SetEndErrorSize(0.01);
    TGaxis::SetMaxDigits(4);
    
    gStyle->SetTitleX(0.55); //center title
    gStyle->SetTitleW(0.6);
    gStyle->SetTitleAlign(23);
}
//__________________________________________________________________________________________________________________
// divide each bin by bin width
template <typename HIST>    //TH1, TH1F, TH1D, TH1I
void densityHist(HIST* h) {
    if(!h) return;
    int nbinsX = h->GetNbinsX();
    if(nbinsX<=0) return;
    for(int ibinX=1; ibinX<=nbinsX; ibinX++) {
        float num = h->GetBinContent(ibinX);
        float err = h->GetBinError(ibinX);
        float binWidth = h->GetBinWidth(ibinX);
        //cout << "before: " << "num = " << num << "\t" << "err = " << err << endl;
        num /= binWidth;
        err /= binWidth;
        //cout << "after: " << "num = " << num << "\t" << "err = " << err << endl;
        h->SetBinContent(ibinX,num);
        h->SetBinError(ibinX,err);
    }
}
template <typename HIST2D>    //TH2, TH2F, TH2D, TH2I
void densityHist2D(HIST2D* h) {
    if(!h) return;
    int nbinsX = h->GetNbinsX();
    int nbinsY = h->GetNbinsY();
    if(nbinsX<=0) return;
    for(int ibinX=1; ibinX<=nbinsX; ibinX++) {
        for(int ibinY=1; ibinY<=nbinsY; ibinY++) {
            float num = h->GetBinContent(ibinX,ibinY);
            float err = h->GetBinError(ibinX,ibinY);
            float binWidthX = h->GetXaxis()->GetBinWidth(ibinX);
            float binWidthY = h->GetYaxis()->GetBinWidth(ibinY);
            //cout << "before: " << "num = " << num << "\t" << "err = " << err << endl;
            num = num / binWidthX / binWidthY;
            err = err / binWidthX / binWidthY;
            //cout << "after: " << "num = " << num << "\t" << "err = " << err << endl;
            h->SetBinContent(ibinX,ibinY,num);
            h->SetBinError(ibinX,ibinY,err);
        }
    }
}

//__________________________________________________________________________________________________________________
// Scale TGraph, TGraphErrors, TGraphAsymmErrors
TGraph* ScaleGraph(TGraph* gr, float scale) {
    for(int i=0; i<gr->GetN(); i++) {
        gr->GetY()[i] *= scale;
    }
    return gr;
}
TGraphErrors* ScaleGraph(TGraphErrors* gr, float scale) {
    for(int i=0; i<gr->GetN(); i++) {
        gr->GetY()[i] *= scale;
        gr->GetEY()[i] *= scale;
    }
    return gr;
}
TGraphAsymmErrors* ScaleGraph(TGraphAsymmErrors* gr, float scale) {
    for(int i=0; i<gr->GetN(); i++) {
        gr->GetY()[i] *= scale;
        gr->GetEYlow()[i] *= scale;
        gr->GetEYhigh()[i] *= scale;
    }
    return gr;
}

//__________________________________________________________________________________________________________________
template <typename HIST3>
HIST3* setHisto(HIST3* h, const char* titlename, const char* xtitle, const char* ytitle,
               int optSize=1,
               double titleSize=0.05, double labelSize=0.05,
               int optOffset=1,
               double x_Offset=1.0, double y_Offset=1.73,
               int optRange=0,
               double x_low=0, double x_up=1, double y_low=0, double y_up=1)
{
    if(h==0) h = new HIST3();
    h->SetTitle(titlename);
    h->GetXaxis()->SetTitle(xtitle);
    h->GetYaxis()->SetTitle(ytitle);
    h->GetXaxis()->SetTitleFont(22);
    h->GetYaxis()->SetTitleFont(22);
    h->GetXaxis()->SetLabelFont(42);
    h->GetYaxis()->SetLabelFont(42);
    h->GetZaxis()->SetTitleFont(22);
    h->GetZaxis()->SetLabelFont(42);
    h->GetYaxis()->SetLabelOffset(0.021);

    h->SetLineWidth(2);

    h->GetXaxis()->SetNdivisions(505);
    h->GetYaxis()->SetNdivisions(505);
    h->GetZaxis()->SetNdivisions(505);

    if(optSize) {
        h->GetXaxis()->SetTitleSize(titleSize);
        h->GetXaxis()->SetLabelSize(labelSize);  //axis number size
        h->GetYaxis()->SetTitleSize(titleSize);
        h->GetYaxis()->SetLabelSize(labelSize);
        //h->GetZaxis()->SetTitleSize(titleSize);
        //h->GetZaxis()->SetLabelSize(labelSize);
        }

    if(optOffset) {
        h->GetXaxis()->SetTitleOffset(x_Offset);
        h->GetYaxis()->SetTitleOffset(y_Offset);}

    if(optRange) {
        h->GetXaxis()->SetRangeUser(x_low, x_up);
        h->GetYaxis()->SetRangeUser(y_low, y_up);}
    return h;
}
//__________________________________________________________________________________________________________________
template <typename HIST1>
TH1F* getHisto(char name[], Double_t xlow, Double_t xup, Double_t ylow, Double_t yup, char xTitle[], char yTitle[])
{
    TH1F* dd = new TH1F(name,"",100,xlow,xup);
    dd->SetMinimum(ylow);
    dd->SetMaximum(yup);
    dd->GetXaxis()->SetTitle(xTitle);
    dd->GetYaxis()->SetTitle(yTitle);
    dd->GetXaxis()->SetTitleFont(22);
    dd->GetYaxis()->SetTitleFont(22);
    dd->GetXaxis()->SetLabelFont(62);
    dd->GetYaxis()->SetLabelFont(62);

    dd->GetXaxis()->SetTitleSize(0.05);
    dd->GetXaxis()->SetTitleOffset(0.9);
    dd->GetXaxis()->SetLabelSize(0.05);
    dd->GetYaxis()->SetTitleSize(0.05);
    dd->GetYaxis()->SetTitleOffset(1);
    dd->GetYaxis()->SetLabelSize(0.05);
    //dd->GetXaxis()->CenterTitle(kTRUE);
    //dd->GetYaxis()->CenterTitle(kTRUE);
    dd->GetXaxis()->SetNdivisions(512);
    return dd;
}
//__________________________________________________________________________________________________________________
void setPad(TPad* pad, float left=0.17, float right=0.05, float top=0.05, float bottom=0.12, int nlogy=0, int nlogz=1)
{
    pad->SetLogy(nlogy);
    pad->SetLogz(nlogz);
    pad->SetFillColor(10);
    pad->SetBorderMode(0);
    pad->SetBorderSize(0);
    pad->SetFrameFillColor(10);
    pad->SetFrameBorderMode(0);
    pad->SetFrameBorderSize(0);

    #if 1
    pad->SetGridx(1);
    pad->SetGridy(1);
    pad->SetTickx(1);
    pad->SetTicky(1);
    #endif

    #if 1
    pad->SetLeftMargin(left);
    pad->SetRightMargin(right);
    pad->SetTopMargin(top);
    pad->SetBottomMargin(bottom);
    #endif
}
//__________________________________________________________________________________________________________________

TLatex* drawLatex(Double_t x, Double_t y, char* text, Int_t textFont, Double_t textSize, Int_t colorIndex)
{
    TLatex *latex = new TLatex(x,y,text);
    latex->SetNDC();
    latex->SetTextFont(textFont);
    latex->SetTextSize(textSize);
    latex->SetTextColor(colorIndex);
    latex->Draw("same");
    return latex;
}
//__________________________________________________________________________________________________________________
template <typename HIST2>
TPaveStats* setPaveStats(HIST2* h, float x_low = 0.15, float y_low = 0.74, float x_up = 0.55, float y_up = 0.94, int n = 0, char name[][250] = NULL)
{
    TPaveStats* ptstats = new TPaveStats(x_low,y_low,x_up,y_up,"brNDC");
    ptstats->SetName("stats");
    ptstats->SetBorderSize(2);
    //ptstats->SetFillColor(10);
    ptstats->SetTextAlign(12);
    ptstats->SetTextSize(0.04);
    ptstats->SetTextFont(22);
    ptstats->SetShadowColor(0);
    //ptstats->SetLineColor(0);
    ptstats->SetFitFormat("5.3g");
    ptstats->SetStatFormat("6.3g");

    if(name)
    {
        gStyle->SetOptStat(0);   //0--not show histograph para, 1--show histo para
        gStyle->SetOptFit(000);  //000--not show fit para, 111--show fit para.
    }
    else
    {
        gStyle->SetOptStat(0);   //0--not show histograph para, 1--show histo para
        gStyle->SetOptFit(111);  //000--not show fit para, 111--show fit para.
    }

    TText* text;
    for(int i=0; i<n; i++) text = ptstats->AddText(name[i]);

    ptstats->SetOptStat(0);
    ptstats->SetOptFit(000);
    h->GetListOfFunctions()->Add(ptstats);
    ptstats->SetParent(h);

    /*#if 0
    TList *list = ps->GetListOfLines();
    // Remove the RMS line
    TText *tconst = ps->GetLineWith("RMS");
    list->Remove(tconst);
    // Add a new line in the stat box.
    // Note that "=" is a control character
    TLatex *myt = new TLatex(0,0,"Test = 10");
    myt ->SetTextFont(42);
    myt ->SetTextSize(0.04);
    myt ->SetTextColor(kRed);
    list->Add(myt);
    #endif 0*/

    ptstats->Draw();

    //the following line is needed to avoid that the automatic redrawing of stats
    //h->SetStats(0);
    return ptstats;
}

//__________________________________________________________________________________________________________________
// e purity------------one two five Gauss-------------------//

Double_t oneGaus(Double_t* x, Double_t* par)
{
    return par[0]*TMath::Gaus(x[0], par[1], par[2], kTRUE)*par[3];
}

Double_t twoGaus(Double_t* x, Double_t* par)
{
    Double_t sum = 0;
    for(Int_t i=0; i<2; i++) sum += par[3*i+0]*TMath::Gaus(x[0], par[3*i+1], par[3*i+2], kTRUE)*par[6];
    return sum;
}

Double_t threeGaus(Double_t* x, Double_t* par)
{
    Double_t sum = 0;
    for(Int_t i=0; i<3; i++) sum += par[3*i+0]*TMath::Gaus(x[0], par[3*i+1], par[3*i+2], kTRUE)*par[9];
    return sum;
}

Double_t fourGaus(Double_t* x, Double_t* par)
{
    Double_t sum = 0;
    for(Int_t i=0; i<4; i++) sum += par[3*i+0]*TMath::Gaus(x[0], par[3*i+1], par[3*i+2], kTRUE)*par[12];
    return sum;
}

Double_t fiveGaus(Double_t* x, Double_t* par)
{
    Double_t sum = 0;
    for(Int_t i=0; i<5; i++) sum += par[3*i+0]*TMath::Gaus(x[0], par[3*i+1], par[3*i+2], kTRUE)*par[15];
    return sum;
}

Double_t fInvBeta(Double_t* x, Double_t* par)  //par[0]---m, par[1]---per
{
    Double_t p = x[0];
    Double_t mass = par[0];
    Double_t percent = par[1];
    return sqrt(pow(mass,2)/pow(p,2)+1) * (1+percent);
}
// ------------possion function-------------------//
//Double_t myPossion(Double_t* x, Double_t* par) {
    
//}

/*Double_t fMass2(Double_t* x, Double* par)
 {
 Double_t p = x[0];
 Double_t mass = par[0];
 Double_t percent = par[1];
 Double_t mass_c = mass * (1+percent);
 return pow(p,2) * sqrt(pow)
 }*/

Double_t oneExp(Double_t* x, Double_t* par)
{
    return par[0]*pow(10,par[1]*x[0])+par[2];
}

//__________________________________________________________________________________________________________________
double fk_lw(double* x, double* par)
{
    if(x[0]<1.372653) return par[0] + par[1]*x[0] +par[2]*x[0]*x[0];
    else return par[0] + par[1]*1.372653 +par[2]*1.372653*1.372653;
}
double fk_up(double* x, double* par)
{
    if(x[0]<1.899554) return par[0] + par[1]*x[0] +par[2]*x[0]*x[0];
    else return par[0] + par[1]*1.899554 +par[2]*1.899554*1.899554;
}
double fpi_lw(double* x, double* par)
{
    return -2.0;
}
double fpi_up(double* x, double* par)
{
    if(x[0]<1.600093) return par[0] + par[1]*x[0];
    else return par[0] + par[1]*1.600093;
}

//__________________________________________________________________________________________________________________
//double HtoF(double x, TH1F* h)
//{ 
//  int bin = (int)h->FindBin(x);
//  return h->GetBinContent(bin);
//} 
  

//__________________________________________________________________________________________________________________
double round(double val)
{
    return (val> 0.0) ? floor(val+ 0.5) : ceil(val- 0.5);
}
//__________________________________________________________________________________________________________________
//function used in centrality define
double gRefTail(double* x, double* par)
{
    double A = par[0];
    double h = par[1];
    double sigma = par[2];
    return A*TMath::Erf(-sigma*(x[0]-h))+A;
}

double fVzDepand(double* x, double* par)  //polN
{
    double sum = par[0];
    double pol = 1;
    int n=1;
    for(int i=0; i<n; i++) {
        pol *= x[0];
        sum += par[i+1]*pol;
    }
    return sum;
}
//__________________________________________________________________________________________________________________
//__________________________________________________________________________________________________________________
//__________________________________________________________________________________________________________________
#endif // MYFUNCTION_H_INCLUDED
