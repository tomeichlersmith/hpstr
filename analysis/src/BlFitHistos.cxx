#include "BlFitHistos.h"
#include <sstream>

//Globally used fit window
double fitmin;
double fitmax;

BlFitHistos::BlFitHistos() {
    //ModuleMapper used to translate between hw and sw names
    mmapper_ = new ModuleMapper();

    //Global svtIDMap built in ModuleMapper. Required to output baselines in database format
    svtIDMap = mmapper_->buildChannelSvtIDMap();
}

BlFitHistos::~BlFitHistos() {
}

void BlFitHistos::getHistosFromFile(TFile* inFile, std::string layer){
    
    for (auto hist: _h_configs.items()) {
        std::string h_name = hist.key();
        std::cout << "h_name: " << h_name << std::endl;

        //Check for substring "L<layer>" in input histograms
        std::size_t found = (h_name).find_first_of("_");
        std::string h_layer = h_name.substr(0,found);
        if (h_layer.find(layer) == std::string::npos)
            continue;
        if(debug_)
            std::cout << "found layer " << layer << " in " << h_name << std::endl;

        found = (h_name).find_last_of("_");
        std::string sample = h_name.substr(found+1);
        if(debug_)
            std::cout << "time sample is " << sample << std::endl;

        TIter next(inFile->GetListOfKeys());
        TKey *key;
        bool got = false;
        std::string histoname;
        while ((key = (TKey*)next())) {
            std::string name(key->GetName());
            if(debug_)
                std::cout << "Checking if histokey " << name << " matches " << h_name << std::endl;
            if (name.find(h_name) != std::string::npos){
                TH2F *hh = (TH2F*) inFile-> Get(key->GetName());
                histos2d[key->GetName()] = hh;
                std::cout << "Adding histo " << key->GetName() << " to list of histos to fit" << std::endl;
            }
        }
    }
}

void BlFitHistos::iterativeGausFit(TH1D* hist, double min, double max, double sigmaRange, double hardminimum, double hardmaximum) {

    double minthresh = hardminimum;
    double threshold = hardmaximum;
    //perform single Gaus fit across full range of histo
    if (min < 0.0 )
        min = hist->GetBinLowEdge(hist->FindFirstBinAbove(0.0,1));

    if (max < 0.0 )
        max = hist->GetBinLowEdge(hist->FindLastBinAbove(0.0,1));

    //Initial fit to establish rough mean and sigma
    if (debug_)
        std::cout << "initial min: " << min << " | initial max: " << max << std::endl;

    TF1 *fitA = new TF1("fitA", "gaus", min, max);
    hist->Fit("fitA","ORQN","");
    double fitAMean = fitA->GetParameter(1);
    double fitASig = fitA->GetParameter(2);

    delete fitA;

    if(fitAMean + fitASig*sigmaRange < threshold)
        max = fitAMean + fitASig*sigmaRange;
    if(fitAMean - fitASig*sigmaRange > minthresh)
        min = fitAMean - fitASig*sigmaRange;

    if (debug_){
        std::cout << "meanA " << fitAMean << " | sigmaA: " << fitASig << std::endl;
        std::cout << "minA: " << min << " | maxA: " << max << std::endl;
    }

    //Fit second time using updated min and max
    TF1 *fitB = new TF1("fitB", "gaus", min, max);
    hist->Fit("fitB","ORQN","");
    double fitMean = fitB->GetParameter(1);
    double fitSig = fitB->GetParameter(2);
    if (debug_)
        std::cout << "meanB " << fitMean << " | sigmaB: " << fitSig << std::endl;

    if(fitMean + fitSig*sigmaRange < threshold)
        max = fitMean + fitSig*sigmaRange;
    if(fitMean - fitSig*sigmaRange > minthresh)
        min = fitMean - fitSig*sigmaRange;

    if (debug_)
        std::cout << "minB: " << min << " | maxB: " << max << std::endl;

    TF1 *fit = new TF1("fit", "gaus", min, max);
    hist->Fit("fit","ORQN","");

    double newFitSig = 99999;
    double newFitMean = 99999;
    int i = 0;
    while (std::abs(fitSig - newFitSig) > 0.0005 || std::abs(fitMean - newFitMean) > 0.0005) {
        double itermax = max;
        double itermin = min;

        if(i > 0){
            fitMean = newFitMean;
            fitSig = newFitSig;
        }

        if(fitMean + fitSig*sigmaRange < threshold)
            max = fitMean + fitSig*sigmaRange;
        if(fitMean - fitSig*sigmaRange > minthresh)
            min = fitMean - fitSig*sigmaRange;
        fit->SetRange(min,max);
        hist->Fit("fit","ORQN","");

        newFitMean = fit->GetParameter(1);
        newFitSig = fit->GetParameter(2);

        if(i > 20)
            break;
        i = i + 1;
    }

    fitmin = min;
    fitmax = max;
}

void BlFitHistos::fit2DHistoChannelBaselines(std::map<std::string,TH2F*> histos2d, int rebin_, int minStats_, int deadRMS_, std::string thresholdsFileIn_, FlatTupleMaker* flat_tuple_) {
     
    //Get half module string names 
    std::vector<std::string> halfmodule_strings;
    mmapper_->getStrings(halfmodule_strings);

    //Build FebHybAPV thresholds map in modulemapper
    mmapper_->buildApvChannelMap();
    //Read apv channel thresholds from file for closest run
    mmapper_->ReadThresholdsFile(thresholdsFileIn_);

    //Loop over rawsvthit 2D histograms, one for each selected halfmodule
    for(std::map<std::string, TH2F*>::iterator it = histos2d.begin(); it != histos2d.end(); ++it)
    {
        TH2F* halfmodule_hh = it->second; 
        halfmodule_hh->RebinY(rebin_);
        halfmodule_hh->Write();
        std::string hh_name = it->first;

        //get the hardware tag for this F<n>H<M>. Required for svtid mapping
        std::string hwTag;

        for(std::vector<std::string>::iterator it = halfmodule_strings.begin(); it != halfmodule_strings.end(); ++it){
            if(hh_name.find(*it) != std::string::npos){
                hwTag = mmapper_->getHwFromString(*it);
                if(debug_)
                    std::cout << "hwTag for " << hh_name << " is " << hwTag << std::endl;
                break;
            }
        }

        //Feb and Hybrid numbers
        std::string feb = (hwTag.substr(1,1));
        std::string hyb = (hwTag.substr(3,1));

        //Perform fitting procedure over all channels on a sensor
        for(int cc=0; cc < 640 ; ++cc) 
        {
            if(debug_){
                std::cout << hh_name << " " << cc << std::endl;
                std::cout << "CHANNEL " << cc << std::endl;
            }
            //std::cout << hh_name << " " << cc << std::endl;
            if (cc%100 == 0)
                std::cout << "CHANNEL " << cc << std::endl;

            //get the global svt_id for channel
            int svt_id = mmapper_->getSvtIDFromHWChannel(cc, hwTag, svtIDMap);
            if(debug_)
                std::cout << "Global SVT ID: " << svt_id << std::endl;

            //Feb 0-1 have max_channel = 512.
            //svt_id = 99999 means max_channel reached. Skip cc > 512 
            if(svt_id == 99999)  
                continue;

            //load apv channel readout threshold value from run_thresholds.dat file
            double threshold = (double) mmapper_->getThresholdValue(feb, hyb, cc); 
            if(debug_)
                std::cout << "THRESHOLD F" <<feb << "H" <<hyb << "channel " << cc << ": " << threshold << std::endl;

            //Set Channel and Hybrid information and paramaters in the flat tuple
            flat_tuple_->setVariableValue("halfmodule_hh", hh_name);
            flat_tuple_->setVariableValue("channel", cc);
            flat_tuple_->setVariableValue("svt_id", svt_id);
            flat_tuple_->setVariableValue("minStats", (double)minStats_);
            flat_tuple_->setVariableValue("rebin", (double)rebin_);

            //Get YProjection (1D Channel Histo) from 2D Histogram 
            TH1D* projy_h = halfmodule_hh->ProjectionY(Form("%s_proY_ch%i",hh_name.c_str(),cc),
                    cc+1,cc+1,"e");
            projy_h->Smooth(1);
            projy_h->SetTitle(Form("%s_proY_ch%i",hh_name.c_str(),cc));

            //Check number of entries and RMS of channel
            flat_tuple_->setVariableValue("n_entries", projy_h->GetEntries());
            double chRMS = projy_h->GetRMS();
            flat_tuple_->setVariableValue("rms", chRMS);

            //If the channel RMS is under some threshold, flag it as "dead"
            if(chRMS < deadRMS_ || projy_h->GetEntries() == 0)
                flat_tuple_->setVariableValue("dead",1.0);

            //Fit window max set by threshold value loaded from file
            //Minum value of x set to first bin with fraction of maximum value
            double maxbin = projy_h->GetBinContent(projy_h->GetMaximumBin());
            //double frac = 0.15;
            double frac = 0.25;
            int minbin = projy_h->FindFirstBinAbove((double)frac*maxbin,1);
            double minx = projy_h->GetBinLowEdge(minbin);
            double binwidth = projy_h->GetBinWidth(minbin);
            double minxVal = projy_h->GetBinContent(minbin);

            double maxx = threshold - binwidth*1;
            flat_tuple_->setVariableValue("threshold",maxx);

            //If channel does not have the minimum statistics required, set all variables to -9999.9
            //and skip the fit procedure on this channel
            if (minbin == -1 || projy_h->GetEntries() < minStats_ ) 
            {
                flat_tuple_->setVariableValue("BlFitMean", -9999.9);
                flat_tuple_->setVariableValue("BlFitSigma", -9999.9);
                flat_tuple_->setVariableValue("BlFitNorm", -9999.9);
                flat_tuple_->setVariableValue("BlFitRangeLower", -9999.9);
                flat_tuple_->setVariableValue("BlFitRangeUpper", -9999.9);
                flat_tuple_->setVariableValue("BlFitChi2", -9999.9);
                flat_tuple_->setVariableValue("BlFitNdf", -9999.9);
                flat_tuple_->setVariableValue("lowdaq", 0.0);
                flat_tuple_->setVariableValue("suplowDaq", 0.0);
                flat_tuple_->setVariableValue("lowStats",1.0);
                flat_tuple_->setVariableValue("badfit",0.0);
                flat_tuple_->fill();
                continue;
            }
            flat_tuple_->setVariableValue("lowStats",0.0);

            /*
            //If baseline fitting an online baseline, must set simpleGauseFit_ to true!
            if(simple_ == true){
                TF1* simpleFit = singleGausIterative(projy_h, 1.5,minx);
                const double* parameters;
                parameters = simpleFit->GetParameters();
                flat_tuple_->setVariableValue("BlFitMean", parameters[1]);
                flat_tuple_->setVariableValue("BlFitSigma", parameters[2]);
                flat_tuple_->setVariableValue("BlFitNorm", parameters[0]);
                flat_tuple_->setVariableValue("BlFitChi2", simpleFit->GetChisquare());
                flat_tuple_->setVariableValue("BlFitNdf", simpleFit->GetNDF());
                flat_tuple_->setVariableValue("BlFitRangeLower", fitmin);
                flat_tuple_->setVariableValue("BlFitRangeUpper", fitmax);

                flat_tuple_->fill();
                delete simpleFit;
                delete projy_h;
                continue;
            }
            */

            

            /*
            int iter = 0;
            int itermaxbin = projy_h->FindBin(maxx);
            itermaxbin = itermaxbin - 1;
            double itermaxX = projy_h->GetBinLowEdge(itermaxbin);
            double itermaxVal = projy_h->GetBinContent(itermaxbin);

            int iterminbin = itermaxbin - 20;
            double iterminX = projy_h->GetBinLowEdge(iterminbin);
            double iterminVal = projy_h->GetBinContent(iterminbin);

            if(debug_){
                std::cout << "xmax: " << itermaxX << " value: " << itermaxVal << std::endl;
                std::cout << "xmin: " << minx << " value: " << minxVal << std::endl;
            }

            double tallestVal = 0;
            double tallestX = 0;
            bool peakfound = false;

            double maxPeak = 0.0;
            double peakminbin = iterminbin;
            double peakmaxbin = itermaxbin;
            int j = 0;
            while(iterminbin > minbin) {

                peakfound = false;
                double peakval = 0;
                iterminbin--;
                iterminVal = projy_h->GetBinContent(iterminbin);
                iterminX = projy_h->GetBinLowEdge(iterminbin);
                itermaxVal = projy_h->GetBinContent(itermaxbin);
                itermaxX = projy_h->GetBinLowEdge(itermaxbin);

                if(debug_){
                    std::cout << "iterminX: " << iterminX << " has value: " << iterminVal << std::endl; 
                    std::cout << "itermaX: " << itermaxX << " has value: " << itermaxVal << std::endl; 
                }

                if(iterminVal > itermaxVal){
                    j = 0;
                }
                else if(iterminVal <= itermaxVal){
                    iterminbin--; 
                    j++;
                }

                if(debug_)
                    std::cout << "j = " << j << std::endl;

                if(j > 10){
                    iterminbin = iterminbin + j;
                    peakminbin = iterminbin;
                    peakmaxbin = itermaxbin;
                    int peakbin = projy_h->GetMaximumBin();

                    projy_h->GetXaxis()->SetRange(peakminbin, peakmaxbin);
                    peakval = projy_h->GetBinContent(peakbin);
                    projy_h->GetXaxis()->SetRange(0,projy_h->FindLastBinAbove(0));
                    peakfound = true;
                    j = 0;
                    if(debug_){
                        std::cout << "local peak found at " << projy_h->GetBinLowEdge(peakbin) << std::endl; 
                        std::cout << "peak value = " << peakval << std::endl;
                    }
                }

                if(peakfound){
                    peakfound = false;
                    int k = 0;
                    //check if peak value is found again to the left
                    while(iterminbin > minbin){
                        iterminbin--;
                        iterminVal = projy_h->GetBinContent(iterminbin);
                        if(debug_){
                            std::cout << "iterminx = " << projy_h->GetBinLowEdge(iterminbin) << std::endl;
                            std::cout << "iterminVal = " << iterminVal << std::endl;
                        }
                        if(iterminVal >= peakval){
                            k++;    
                        }
                        else{
                            k=0;
                        }
                        if(k > 10){
                            itermaxbin = iterminbin + k;
                            iterminbin = itermaxbin - 3;
                            k = 0;
                            peakfound = true;
                            if(debug_){
                                std::cout << "New large peak found beyond " << projy_h->GetBinLowEdge(iterminbin) << std::endl;
                                std::cout << "Set iterxmax to : " << projy_h->GetBinLowEdge(itermaxbin) << std::endl;
                            }
                            break;
                        }
                    }
                    if(!peakfound)
                    {
                        iterminbin = peakminbin; 
                        break;
                    }
                }
            }

            fitmax = projy_h->GetBinLowEdge(itermaxbin); 
            fitmin = projy_h->GetBinLowEdge(iterminbin); 
            */

            fitmin = minx;
            fitmax = maxx;
            iterativeGausFit(projy_h, fitmin, fitmax, 1, minx, threshold);

            TF1 *fit = new TF1("fit", "gaus", fitmin, fitmax);
            projy_h->Fit("fit","ORQN","");
            double fitmean = fit->GetParameter(1);
            double fitsigma = fit->GetParameter(2);
            double fitnorm = fit->GetParameter(0);
            double fitchi2 = fit->GetChisquare();
            double fitndf = fit->GetNDF();

            if(debug_){
                std::cout << "Fit Mean: " << fitmean << std::endl;
                std::cout << "Fit sigma: " << fitsigma << std::endl;
                std::cout << "Fit norm: " << fitnorm << std::endl;
                std::cout << "Fit min: " << fitmin << std::endl;
                std::cout << "Fit max: " << fitmax << std::endl;
            }
            delete fit;

            //If fit mean is less than fitmax, channel does not have full gaussian shape.
            //Mark channel as super_low_Daq
            bool badfit = false;
            bool suplowDaq = false;
            if(fitmean <= fitmin || fitmin > fitmax){
                badfit = true;
                if(debug_)
                    std::cout << "bad fit!" << std::endl;
            }

            if(!badfit && fitmean >= fitmax){
                suplowDaq = true;
                if(debug_)
                    std::cout << "Super low Daq threshold" << std::endl;
            }

            /*
            if(suplowDaq){
                //refit with xmax at threshold value and xmin NSigma lower
                fitmax = maxx;
                fitmin = fitmin - 1.5*fitsigma;
                if(debug_)
                    std::cout << "suplowdaq max and min: " << fitmax << " " << fitmin << std::endl;
                if(fitmin < minx)
                    fitmin = minx;
                TF1 *fitL = new TF1("fitL", "gaus", fitmin, fitmax);
                projy_h->Fit("fitL","ORQN","");
                fitmean = fitL->GetParameter(1);
                fitsigma = fitL->GetParameter(2);
                fitnorm = fitL->GetParameter(0);
                fitchi2 = fitL->GetChisquare();
                fitndf = fitL->GetNDF();

                if(fitmean <= fitmin || fitmin > fitmax){
                    badfit = true;
                    if(debug_)
                        std::cout << "bad fit!" << std::endl;
                }
            }
            */

            bool lowdaq = false;
            if(!badfit && !suplowDaq){
                //Check channel to see if it has a low DAQ threshold, where landau shape interferes with baseline
                //If maxbin occurs outside of fit mean by NSigma...flag
                double maxbinx = projy_h->GetBinLowEdge(projy_h->GetMaximumBin());
                if ((std::abs(maxbinx - fitmean) > fitsigma)){
                    lowdaq = true;
                }
                //If fitmean > fitmax or fitmean < fitmin...flag
                if(fitmean > fitmax || fitmean < fitmin)
                    lowdaq = true;

                //If bins after fitmax averaged to the right are greater than the fitmean...flag
                double maxavg = 0;
                int fitmaxbin = projy_h->FindBin(fitmax);
                for (int i = 1; i < 6; i++){
                    maxavg = maxavg + projy_h->GetBinContent(fitmaxbin + i); 
                }
                maxavg = maxavg/5;
                if(maxavg > fitnorm)
                    lowdaq = true;
            }

            if(debug_)
                if(lowdaq)
                    std::cout << "Low daq threshold" << std::endl;

            flat_tuple_->setVariableValue("lowdaq", lowdaq);
            flat_tuple_->setVariableValue("suplowDaq", suplowDaq);
            flat_tuple_->setVariableValue("badfit", badfit);
            flat_tuple_->setVariableValue("BlFitMean", fitmean);
            flat_tuple_->setVariableValue("BlFitSigma", fitsigma);
            flat_tuple_->setVariableValue("BlFitNorm", fitnorm);
            flat_tuple_->setVariableValue("BlFitChi2", fitchi2);
            flat_tuple_->setVariableValue("BlFitNdf", fitndf);
            flat_tuple_->setVariableValue("BlFitRangeLower", fitmin);
            flat_tuple_->setVariableValue("BlFitRangeUpper", fitmax);

            flat_tuple_->fill();

            delete projy_h;
            continue;
                        
        }
    }
}

