/** @file Dispersion.cxx
@brief implement Dispersion and Dispersion::Hist

$Header$
*/

#include "Dispersion.h"

#include "IRF.h"
#include "IrfAnalysis.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPaveStats.h"

#include <cmath>
#include <iomanip>
#include <cassert>

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//                Disperion::Hist
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
// specify fit function -- this devised by Riccardo Rando
const char* Dispersion::Hist::function="[0]*pow(1+x,[1])/(1+exp(x/[2]))";
const char* Dispersion::Hist::pnames[]={"dnorm", "ltail", "rwidth"};
double      Dispersion::Hist::pinit[] ={1,      20,      0.02};
double      Dispersion::Hist::pmin[]  ={0.,      2,      0.005};
double      Dispersion::Hist::pmax[]  ={1e6,    50,      0.5};
double      Dispersion::Hist::fitrange[]={-0.4, 0.5};
bool        Dispersion::Hist::s_logy = true;
int         Dispersion::Hist::s_minEntries( 10);
int Dispersion::Hist::npars(){return sizeof(pnames)/sizeof(const char*);}

double Dispersion::function(double * x, double * p)
{
    return p[0]*pow(1+*x,p[1])/(1+exp(*x/p[2]));

}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Dispersion::Hist::Hist(std::string id, std::string title)
{
    static std::string axistitles(";fit/generated-1");
    m_h = new TH1F(id.c_str(), (title+axistitles).c_str(),  50, -1., 1.);
    m_h->Sumw2(); // just to make error bars

}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void Dispersion::Hist::fill(double diff){
    m_h->Fill(diff);
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void Dispersion::Hist::fit(std::string opts)
{
    TH1F & h = *m_h;; 

    // add overflow to last bin for display
    int nbins ( h.GetNbinsX() );
    h.SetBinContent(nbins, h.GetBinContent(nbins) +h.GetBinContent(nbins+1));

    TF1* f1 = new TF1("f1",function,fitrange[0],fitrange[1]);
    for (unsigned int i = 0; i < sizeof(pmin)/sizeof(double); i++) {
        f1->SetParLimits(i, pmin[i], pmax[i]);
        f1->SetParName(i, pnames[i]);
    }
    f1->SetParameters(pinit);

    if( h.GetEntries()>= s_minEntries ){
        h.Fit(f1, opts.c_str());

        // measure tails
        double binsize( 2./ nbins ),
            predlow ( f1->Integral(-1.,0.)/binsize ),
            predhigh( f1->Integral(0, 1.)/binsize ),
            measlow ( h.Integral(1,nbins/2) ),
            meashigh( h.Integral(nbins/2+1, nbins) );

        m_lowTail = (measlow-predlow)/(measlow+meashigh);
        m_highTail= (meashigh-predhigh)/(measlow+meashigh);
    }
    h.Write(); // save fit info 
}
double Dispersion::Hist::entries()const{return m_h->GetEntries();}
double Dispersion::Hist::parameter(int n)const
{
    TF1* f1 = m_h->GetFunction("f1");
    if( f1 ==0 ) return 99; // no fit.
    if( n>=0 )  return  f1->GetParameter(n);
    return f1->GetChisquare();
}
double Dispersion::Hist::ltail()const{

    return  parameter(1);
}
double Dispersion::Hist::rwidth()const{
    return  parameter(2);
}
double Dispersion::Hist::chisq()const{
    return  parameter(-1);
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void Dispersion::Hist::summaryTitle(std::ostream& out)
{
    using std::setw;
    out << "\n\t\t Dispersion fit summary\n" 
        << setw(30)<< std::right << "title" 
        << setw(10) << "count"
        << std::left<<"  "
        << setw(10) << "chisq"
        << setw(10) << "ltail"
        << setw(10) << "rwidth"
        << setw(10) << "lotail(%)"
        << setw(10) << "hitail(%)"
        << std::endl;

}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void Dispersion::Hist::summarize(std::ostream & out)
{
    using std::setw; using std::setprecision;

    int n= static_cast<int>( entries() );
    out  
        << setprecision(3) << std::right << std::fixed
        << setw(30) << m_h->GetTitle() 
        << setw(10) << n;
    if( n < s_minEntries) {
        out << "    -- " ;
    }else{
        out << std::left<<"  " << setprecision(1)
            << setw(10) << chisq()
            << setw(10) << ltail()
            << setprecision(3)
            << setw(10) << rwidth()
            << std::setprecision(1)
            << setw(10) << m_lowTail*100.
            << setw(10) << m_highTail*100.;
    }
    out   << std::endl;

}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void Dispersion::Hist::draw()const
{//
    // set up appearance, and draw to current pad
    if( s_logy) gPad->SetLogy();

    TH1F & h = *m_h;
    h.SetStats(true);
    h.SetLineColor(kRed);
    h.GetXaxis()->CenterTitle(true);

    TList * list = h.GetListOfFunctions();
    TPaveStats  *s = (TPaveStats*)list->FindObject("stats");
    if( s!= 0 ){
        s->SetY1NDC(0.6);
    }
    h.Draw();
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void Dispersion::Hist::getFitPars(std::vector<double>& params)const
{
    params.resize(npars());
    TF1* f1 = m_h->GetFunction("f1");
    if( f1==0 ){ // no fit, just return default
        copy(pinit, pinit+npars(), params.begin());
    }else {
        f1->GetParameters(&params[0]);
    }
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//               Dispersion
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Dispersion::Dispersion( IrfAnalysis& irf, std::ostream& log)
: m_irf(irf)
, m_log(&log)
{
    m_hists.resize(irf.size());
    for (int ebin = 0; ebin < irf.energy_bins; ++ebin) {
        for (int abin = 0; abin <= irf.angle_bins; ++abin) {
            int id = irf.ident(ebin,abin);
            std::ostringstream title;
            title << (int)(irf.eCenter(ebin)+0.5) << " MeV," ;
            if ( abin < irf.angle_bins ) {
                title << irf.angles[abin] << "-"<< irf.angles[abin+1] << " degrees";
            }else {
                title <<  irf.angles[0] << "-"<< irf.angles[irf.angle_bins-2] << " degrees";
            }
            m_hists[id] = Hist(IRF::hist_name(abin,ebin,"dsp"),title.str());
        }
    }
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Dispersion::~Dispersion()
{}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void Dispersion::fill(double diff, double energy, double costheta, bool/* front*/)
{
    int z_bin = m_irf.angle_bin( costheta );     if( z_bin>= m_irf.angle_bins) return;
    int e_bin = m_irf.energy_bin(energy);        if( e_bin<0 || e_bin>= m_irf.energy_bins )return;

    int id =  m_irf.ident(e_bin, z_bin);
    m_hists[id].fill(diff);

    // set special combined hist, accumulate all but last bins of angles
    if( z_bin< m_irf.angle_bins-2) {
        m_hists[m_irf.ident(e_bin, m_irf.angle_bins)].fill(diff);
    }
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void Dispersion::summarize()
{
    // summarize psf plot contents, fit 
    Hist::summaryTitle( out());
    for( HistList::iterator it = m_hists.begin(); it!=m_hists.end(); ++it){
        (*it).summarize(out());
    }
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void Dispersion::fit(std::string opts)
{
    for( HistList::iterator it = m_hists.begin(); it!=m_hists.end(); ++it){
        (*it).fit(opts);
    }
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void Dispersion::draw(const std::string &ps_filename ) {
    gStyle->SetOptFit(111);

    TCanvas c;

    for( int abin=0; abin<= m_irf.angle_bins; ++abin){
        int rows=2;
        m_irf.divideCanvas(c,m_irf.energy_bins/rows,rows, 
            std::string("Plots from ") +m_irf.summary_filename());
        for(int ebin=0; ebin<m_irf.energy_bins; ++ebin){
            c.cd(ebin+1);
            gPad->SetRightMargin(0.02);
            gPad->SetTopMargin(0.03);
            m_hists[m_irf.ident(ebin,abin)].draw();
        }
        std::cout << "Printing page #" << (abin+1) << std::endl; 
        if( abin==0) c.Print( (ps_filename+"(").c_str());
        else if (abin<m_irf.angle_bins) c.Print(ps_filename.c_str());
        else c.Print( (ps_filename+")").c_str());
    }
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void Dispersion::fillParameterTables()
{
    // make a set of 2-d histograms with values of the fit parameters
    // binning according to energy and costheta bins

    for( int i = 0; i< Hist::npars(); ++i){
        std::string name(Hist::pnames[i]);

        TH2F* h2 = new TH2F(name.c_str(), (name+";log energy; costheta").c_str() 
        ,IRF::energy_bins, IRF::logemin, IRF::logemin+IRF::energy_bins*IRF::logedelta
        ,IRF::angle_bins,  0.3, 1.0
        );
        std::vector<double> pars;

        int index(0);
        for( HistList::iterator it = m_hists.begin(); it!=m_hists.end(); ++it, ++index){
            it->getFitPars(pars);
            double costheta = 0.95 - 0.1*(index/IRF::angle_bins);
            double logenergy = log10(IRF::eCenter(index));

            h2->Fill(logenergy, costheta, pars[i]);
        }
        h2->GetXaxis()->CenterTitle();
        h2->GetYaxis()->CenterTitle();
        h2->Write();
    }
}
