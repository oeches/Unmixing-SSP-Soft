#ifndef UNMIXING_H
#define UNMIXING_H

#include <armadillo>
#include <iostream>
#include "libmcmc.h"

template<class T>
void unmixing(T imagePixel, int chainLength, MCMCChains &TAlphaPlus, MCMCChains &TSigma2r);

template<class T>
void unmixing(T imagePixel, int chainLength, MCMCChains &TAlphaPlus, MCMCChains &TSigma2r)
{
    //initialization quantities
    unsigned int numEnd = imagePixel.getEndmNumb();
    arma::Mat<double> endmMat = imagePixel.getEmdmMat();
    arma::Mat<double> dataVal = imagePixel.getSpectralValues();

    // Hyperparameter
    double delta = 0.001;
    Samples deltaSamples(1, "Hyperparameter", delta*arma::ones(1,1));

    // endmember variance
    arma::mat sigma2r(1, numEnd);
    sigma2r = 0.001*arma::ones(1, numEnd);
    Samples sigma2rSamples(numEnd, "Endmember variance", sigma2r);

    // abundance vector
    arma::mat abundance(1, numEnd, arma::fill::ones);
    abundance /= numEnd;
    Samples abundanceSamples(numEnd, "Abundances", abundance);

    // boolean
    arma::mat rho(1, numEnd, arma::fill::zeros);
    Samples rhoSamples(numEnd, "Boolean", rho);

    // useful quantities
    arma::vec sig(numEnd);
    sig = 0.05 * arma::ones<arma::vec>(numEnd);
    arma::vec Tx(numEnd);
    Tx = 0.3*arma::ones<arma::vec>(numEnd);
    arma::mat vecRhoSamples(numEnd, 100, arma::fill::zeros);

    // Initializing chains
    MCMCChains TRho(chainLength, numEnd, "Acceptation rate chain");
    TRho.fillInitialSamples(rho);
    TAlphaPlus.fillInitialSamples(abundance);
    TSigma2r.fillInitialSamples(sigma2r);

    //Initializing samplers
    MwG_abundances sampleAbund(dataVal, endmMat,
                               Tx, sig, 0.5, rhoSamples, 0);
    GibbsForVariance sampleVar(dataVal, endmMat,
                               delta);
    GibbsForHyperparameter sampleHyper(dataVal, endmMat);

    // for loop
    for (int m_compt = 0; m_compt < chainLength; ++m_compt)
    {
        // updating the abundance coefficients
        sampleAbund.setmCompt(m_compt);
        sampleAbund.setSamplesRho(rhoSamples);
        sampleAbund.generateSamples(sigma2rSamples, abundanceSamples);

        TAlphaPlus.fillSamples(abundanceSamples.getCurrentSamples(), m_compt);
        rhoSamples = sampleAbund.getSamplesRho();
        TRho.fillSamples(rhoSamples.getCurrentSamples(), m_compt);

        if ((m_compt > 0) && (m_compt % 99 == 0))
        {
            // abundance vector acceptance rate for every 100 iterations
            vecRhoSamples = arma::trans(TRho.getNSamples(99, m_compt - 99));
            Tx = arma::vectorise(arma::mean(vecRhoSamples, 1));
            sampleAbund.setTx(Tx);
        }

        // sampling variance
        sampleVar.setDelta(deltaSamples.getCurrentSamples()(0,0));
        sampleVar.generateSamples(abundanceSamples, sigma2rSamples);
        TSigma2r.fillSamples(sigma2rSamples.getCurrentSamples(), m_compt);

        // updating hyperparameter
        sampleHyper.generateSamples(sigma2rSamples, deltaSamples);
    }

}



#endif // UNMIXING_H
