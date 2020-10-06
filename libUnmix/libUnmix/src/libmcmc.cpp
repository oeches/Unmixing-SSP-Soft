#include "include/libmcmc.h"
#include <armadillo>

MCMCChains::MCMCChains(const unsigned int length, const unsigned int dimension, const std::string name) : m_lengthMCMC(length), m_nameMCMC(name)
{
    m_SampleValues.set_size(dimension, length);
    m_SampleValues.fill(0);
}

MCMCChains::MCMCChains(const MCMCChains &ChainToCopy)
{
    m_lengthMCMC = ChainToCopy.m_lengthMCMC;
    m_SampleValues = ChainToCopy.m_SampleValues;
    m_nameMCMC = ChainToCopy.m_nameMCMC;
}

MCMCChains::~MCMCChains()
{

}

void MCMCChains::fillInitialSamples(const arma::mat samplesInit)
{
    m_SampleValues.head_cols(1) = arma::vectorise(samplesInit);
}

void MCMCChains::fillSamples(const arma::mat samples, const unsigned int index)
{
    m_SampleValues.col(index) = arma::vectorise(samples);
}

void MCMCChains::setNameMCMC(std::string nameMCMC)
{
    m_nameMCMC = nameMCMC;
}

arma::mat MCMCChains::getNSamples(const unsigned int Nsamples, const unsigned int Nbegin) const
{
    return arma::trans(m_SampleValues.cols(Nbegin, Nsamples+Nbegin));
}

Samples::Samples(const unsigned int dimension, const std::string name, const arma::mat sampleInit)
    : m_dimension(dimension), m_name(name)
{
    m_previousSample.set_size(dimension, 1);
    m_currentSample.set_size(dimension, 1);
    m_currentSample = sampleInit;
}

Samples::~Samples()
{

}

arma::mat Samples::getPreviousSamples() const
{
    return m_previousSample;
}

void Samples::setPreviousSamples(arma::mat sampleValue)
{
    m_previousSample = sampleValue;
}

arma::mat Samples::getCurrentSamples() const
{
    return m_currentSample;
}

void Samples::setCurrentSamples(arma::mat sampleValue)
{
    m_currentSample = sampleValue;
}

unsigned int Samples::getDimension() const
{
    return m_dimension;
}

SamplerMCMC::SamplerMCMC(const arma::mat spectraPix, const arma::mat endmSpectra) : m_spectraPix(spectraPix), m_EndmSpectra(endmSpectra)
{

}

SamplerMCMC::~SamplerMCMC()
{

}

MetropolisWithinGibbs::MetropolisWithinGibbs(const arma::mat spectraPix, const arma::mat endmSpectra, const double acceptanceProb) :
    SamplerMCMC(spectraPix, endmSpectra), m_acceptanceProb(acceptanceProb)
{

}

MetropolisWithinGibbs::~MetropolisWithinGibbs()
{

}

MwG_abundances::MwG_abundances(const arma::mat spectraPix, const arma::mat endmSpectra, const arma::vec Tx, const arma::vec varRandomWalk, const double acceptanceProb,
                               const Samples samplesRho, const unsigned int m_compt)
    : MetropolisWithinGibbs(spectraPix, endmSpectra, acceptanceProb), m_Tx(Tx), m_varRandomWalk(varRandomWalk),
    m_samplesRho(samplesRho), m_mCompt(m_compt)
{

}

MwG_abundances::~MwG_abundances()
{

}

arma::vec MwG_abundances::getTx() const
{
    return m_Tx;
}

void MwG_abundances::setTx(const arma::vec TxVal)
{
    m_Tx = TxVal;
}
arma::vec MwG_abundances::getRandomWalk() const
{
    return m_varRandomWalk;
}

void MwG_abundances::setRandomWalk(const arma::vec randomWalkVar)
{
    m_varRandomWalk = randomWalkVar;
}

void MwG_abundances::setmCompt(const unsigned int m_compt)
{
    m_mCompt = m_compt;
}

void MwG_abundances::setSamplesRho(const Samples samplesRho)
{
    m_samplesRho = samplesRho;
}

Samples MwG_abundances::getSamplesRho()
{
    return m_samplesRho;
}

void MwG_abundances::generateSamples(const Samples samplesSig, Samples &samplesAb)
{
    unsigned int bandsNumber = m_EndmSpectra.n_rows;
    unsigned int numEnd = samplesAb.getDimension();
    arma::mat rho = m_samplesRho.getCurrentSamples();
    arma::mat alphaCand = samplesAb.getCurrentSamples();

    //selecting the first abundance coefficient to be sampled
    std::vector<unsigned int> listAbArr;
    for (unsigned int i=0; i<numEnd; ++i)
    {
        listAbArr.push_back(i);
    }
    std::random_shuffle(listAbArr.begin(), listAbArr.end());

    unsigned int firstAb = listAbArr.at(0);
    listAbArr.erase(listAbArr.begin());
    unsigned int ind_k = 0;
    arma::mat spectra_k(bandsNumber, numEnd-1);
    arma::mat Sk(1, numEnd-1);
    arma::mat alpha_k(1, numEnd-1);

    for (const auto &itList : listAbArr)
    {

        if (numEnd < 3)
        {
            spectra_k = arma::zeros(bandsNumber,1);
            Sk = arma::zeros(1,1);
            alpha_k = arma::zeros(1,1);
        }
        else
        {
            unsigned int indL = 0;
            spectra_k.set_size(bandsNumber, numEnd-1);
            Sk.set_size(1, numEnd-1);
            alpha_k.set_size(numEnd-1,1);
            for (const auto &itListBuild : listAbArr)
            {
                spectra_k.col(indL) = m_EndmSpectra.col(itListBuild);
                Sk.col(indL) = samplesSig.getCurrentSamples().col(itListBuild);
                alpha_k.row(indL) = arma::trans(alphaCand.col(itListBuild));
                ++indL;
            }
            spectra_k.shed_col(ind_k);
            Sk.shed_col(ind_k);
            alpha_k.shed_row(ind_k);

        }

        // random walk

        //determining the Gaussian variance in function of the acceptance rate every 100 iterations
        if ((m_Tx(itList) > 0.4) && (m_mCompt % 99 == 0))
        {
            m_varRandomWalk(itList) *= 5;
        }
        if ((m_Tx(itList) < 0.3) && (m_mCompt % 99 == 0))
        {
            m_varRandomWalk(itList) /= 5;
        }

        // proposed samples
        arma::vec valRandn = arma::randn<arma::vec>(1);
        double alpha = alphaCand(itList) + std::sqrt(m_varRandomWalk(itList))*valRandn(0);
        double sum_ak = arma::accu(alpha_k);
        double varDiff;
        double alpha_star = 0.;

        if ((alpha > 0) && ( alpha < 1 - sum_ak))
        {
            alpha_star= alpha;
            double sum_a = arma::accu(arma::sum(alphaCand.cols(0,numEnd-2)));
            arma::vec a_sq = arma::trans(arma::square(alphaCand.cols(0,numEnd-2)));

            arma::mat mu_alpha = m_EndmSpectra.cols(0,numEnd-2)*arma::trans(alphaCand.cols(0,numEnd-2)) + m_EndmSpectra.col(numEnd-1)*(1 - sum_a);
            arma::mat C_alpha = samplesSig.getCurrentSamples().cols(0,numEnd-2)*a_sq +
                    samplesSig.getCurrentSamples().col(numEnd-1)*((1-sum_a)*(1-sum_a));
            arma::mat mu_alpha_star = m_EndmSpectra.col(itList)*alpha_star + spectra_k*alpha_k
                    + m_EndmSpectra.col(firstAb)*(1-(alpha_star+sum_ak));
            arma::mat C_alpha_star = samplesSig.getCurrentSamples().col(itList)*(alpha_star*alpha_star)
                    + Sk*arma::square(alpha_k) + samplesSig.getCurrentSamples().col(firstAb) * ((1-alpha_star-sum_ak)*(1-alpha_star-sum_ak));

            //difference between the logarithms of the distributions
            varDiff = 0.5*(std::pow(arma::norm(m_spectraPix - mu_alpha_star),2.)/C_alpha_star(0)
                           - std::pow(arma::norm(m_spectraPix - mu_alpha),2.)/C_alpha(0))
                           + (bandsNumber/2)*std::log(C_alpha_star(0)/C_alpha(0));
        }
        else
        {
            varDiff = arma::datum::inf;
        }

        if ((varDiff < 0.) || (std::exp(-varDiff) > std::rand()))
        {
            alphaCand.col(itList) = alpha_star;
            rho.col(itList) = 1;
        }
        else
        {
            rho.col(itList) = 0;
        }

        ++ind_k;
    }

    samplesAb.setPreviousSamples(samplesAb.getCurrentSamples());
    arma::mat alphaCandDis = alphaCand;
    alphaCandDis.shed_col(firstAb);
    alphaCand.col(firstAb) = 1 - arma::accu(alphaCandDis);
    m_samplesRho.setPreviousSamples(m_samplesRho.getCurrentSamples());
    samplesAb.setCurrentSamples(alphaCand);
    rho.col(firstAb) = 0;
    m_samplesRho.setCurrentSamples(rho);

}

GibbsForVariance::GibbsForVariance(const arma::mat data, const arma::mat spectra, const double delta)
    : SamplerMCMC(data, spectra), m_delta(delta)
{

}

GibbsForVariance::~GibbsForVariance()
{

}

void GibbsForVariance::setDelta(const double delta)
{
    m_delta = delta;
}

void GibbsForVariance::generateSamples(const Samples samplesAb, Samples &samplesSig)
{
    unsigned int L = m_EndmSpectra.n_rows;
    unsigned int R = m_EndmSpectra.n_cols;
    arma::mat matOnes(1,R, arma::fill::ones);

    double normOfModel = arma::norm(m_spectraPix - m_EndmSpectra * arma::trans(samplesAb.getCurrentSamples()));
    double normOfModelSq = normOfModel*normOfModel;

    double B = normOfModelSq / (2*arma::accu(arma::sum(arma::square(samplesAb.getCurrentSamples()),1))) + m_delta;

    double sig2inv_scal = arma::randg(arma::distr_param(L/2 + 1, 1/B));
    arma::mat sig2inv = sig2inv_scal * matOnes;

    samplesSig.setPreviousSamples(samplesSig.getCurrentSamples());
    samplesSig.setCurrentSamples(1/sig2inv);

}

GibbsForHyperparameter::GibbsForHyperparameter(const arma::mat data, const arma::mat spectra)
    : SamplerMCMC(data, spectra)
{

}

GibbsForHyperparameter::~GibbsForHyperparameter()
{

}

void GibbsForHyperparameter::generateSamples(const Samples samplesSig, Samples &samplesDelta)
{
    double sig2 = samplesSig.getCurrentSamples()(0,0);

    samplesDelta.setCurrentSamples(arma::randg(arma::distr_param(1., sig2))*arma::ones(1,1));

}
