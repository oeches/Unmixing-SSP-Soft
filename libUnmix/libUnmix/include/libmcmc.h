#ifndef LIBMCMC_H
#define LIBMCMC_H
#include "libUnmix_global.h"
#include <iostream>
#include <armadillo>


class LIBUNMIX_EXPORT MCMCChains
{
public:
    MCMCChains(const unsigned int length, const unsigned int dimension, const std::string name);
    MCMCChains(const MCMCChains &ChainToCopy);
    ~MCMCChains();
    void fillInitialSamples(const arma::mat samplesInit);
    void fillSamples(const arma::mat samples, const unsigned int index);
    void setNameMCMC(std::string nameMCMC);
    arma::mat getNSamples(const unsigned int Nsamples, const unsigned int Nbegin) const;

private:
    unsigned int m_lengthMCMC;
    arma::mat m_SampleValues;
    std::string m_nameMCMC;
};

class LIBUNMIX_EXPORT Samples
{
public:
    Samples(const unsigned int dimension, const std::string name, const arma::mat sampleInit);
    ~Samples();
    arma::mat getPreviousSamples() const;
    void setPreviousSamples(arma::mat sampleValue);
    arma::mat getCurrentSamples() const;
    void setCurrentSamples(arma::mat sampleValue);
    unsigned int getDimension() const;

private:
    unsigned int m_dimension;
    arma::mat m_previousSample;
    arma::mat m_currentSample;
    std::string m_name;
};

class LIBUNMIX_EXPORT SamplerMCMC
{
public:
    SamplerMCMC(const arma::mat spectraPix, const arma::mat endmSpectra);
    virtual ~SamplerMCMC();
    virtual void generateSamples(const Samples samples, Samples &samplesExit) = 0;

protected:
    arma::mat m_spectraPix;
    arma::mat m_EndmSpectra;
};

class LIBUNMIX_EXPORT MetropolisWithinGibbs : public SamplerMCMC
{
public:
    MetropolisWithinGibbs(const arma::mat spectraPix, const arma::mat endmSpectra, const double acceptanceProb);
    virtual ~MetropolisWithinGibbs();
    virtual void generateSamples(const Samples samples, Samples &samplesExit) = 0;

protected:
    double m_acceptanceProb;
};

class LIBUNMIX_EXPORT MwG_abundances : public MetropolisWithinGibbs
{
public:
    MwG_abundances(const arma::mat spectraPix, const arma::mat endmSpectra, const arma::vec Tx, const arma::vec varRandomWalk, const double acceptanceProb,
                   const Samples samplesRho, const unsigned int m_compt);
    virtual ~MwG_abundances();
    arma::vec getTx() const ;
    void setTx(const arma::vec TxVal);
    arma::vec getRandomWalk() const ;
    void setRandomWalk(const arma::vec randomWalkVar);
    void setmCompt(const unsigned int m_compt);
    void setSamplesRho(const Samples samplesRho);
    Samples getSamplesRho();
    virtual void generateSamples(const Samples samplesSig, Samples &samplesAb);

private:
    arma::vec m_Tx;
    arma::vec m_varRandomWalk;
    Samples m_samplesRho;
    unsigned int m_mCompt;

};

class LIBUNMIX_EXPORT GibbsForVariance : public SamplerMCMC
{
public:
    GibbsForVariance(const arma::mat spectraPix, const arma::mat endmSpectra, const double delta);
    virtual ~GibbsForVariance();
    void setDelta(const double delta);
    virtual void generateSamples(const Samples samplesAb, Samples &samplesSig);

private:
    double m_delta;
};

class LIBUNMIX_EXPORT GibbsForHyperparameter : public SamplerMCMC
{
public:
    GibbsForHyperparameter(const arma::mat spectraPix, const arma::mat endmSpectra);
    virtual ~GibbsForHyperparameter();
    virtual void generateSamples(const Samples samplesSig, Samples &samplesDelta);

};

#endif // LIBMCMC_H
