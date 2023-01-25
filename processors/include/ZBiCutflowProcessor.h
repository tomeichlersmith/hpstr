#ifndef __ZBICUTFLOW_ANAPROCESSOR_H__
#define __ZBICUTFLOW_ANAPROCESSOR_H__

// HPSTR
#include "Processor.h"
#include "ZBiHistos.h"
#include "IterativeCutSelector.h"
#include "SimpEquations.h"

// ROOT 
#include "TFile.h"
#include "TTree.h"
#include "TRefArray.h"
#include "TBranch.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TH1F.h"
#include "TEfficiency.h"

// C++ 
#include <memory>

class ZBiCutflowProcessor : public Processor {

    public:

        ZBiCutflowProcessor(const std::string& name, Process& process);

        ~ZBiCutflowProcessor();

        virtual void configure(const ParameterSet& parameters);

        virtual void initialize(TTree* tree) {};

        virtual bool process(IEvent* event) {};

        virtual void finalize();

        virtual bool process();

        virtual void initialize(std::string inFilename, std::string outFilename);

        void initializeFlatTuple(TTree* tree, std::map<std::string, double*> &tuple_map);

        double calculateZBi(double n_on, double n_off, double tau);

        bool doesCutVariableExist(std::string cutvariable);

        void printZBiMatrix();


    private:

        int debug_{0}; 
        std::string outFileName_{"zbi_out.root"};
        std::string cuts_cfgFile_{""};
        std::string radSlicFilename_{""};
        std::string vdSimFilename_{""};
        double vdMassMeV_;
        double ApMassMeV_;

        std::map<std::string,double> mcScale_;

        std::vector<std::string> cutlist_strings_{};
        std::vector<std::string> cut_vars_{};

        TFile* outFile_{nullptr};

        //cuts
        typedef std::map<std::string, std::pair<double,int>>::iterator cut_iter_;
        std::map<std::string, std::pair<double,int>>* cuts_;
        std::vector<std::string> cutVariables_;
        std::map<std::string,double> initialIntegrals_;
        std::map<std::string,std::vector<std::pair<double,double>>> global_ZBi_map_;

        //cut selector
        IterativeCutSelector *cutSelector_{nullptr};
        ZBiHistos* cutHistos_{nullptr};

        //signal
        std::string signalHistCfgFilename_{""};
        std::string signalFilename_{""};
        ZBiHistos* signalHistos_{nullptr};
        std::map<std::string,double*> signal_tuple_;
        TTree* signalTree_{nullptr};
        //signal pretrigger vtx distribution
        TH1F* vdSimZ_h_{nullptr};

        //background
        std::string tritrigFilename_{""};
        ZBiHistos* tritrigHistos_{nullptr};
        std::map<std::string,double*> tritrig_tuple_;
        TTree* tritrigTree_{nullptr};

        //simp equations
        SimpEquations* simpEqs_{nullptr};
};



#endif
