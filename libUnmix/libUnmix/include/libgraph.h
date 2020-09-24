#ifndef LIBGRAPH_H
#define LIBGRAPH_H

#include "libUnmix_global.h"
#include <memory>
#include <iostream>
#include <armadillo>
#include "libmcmc.h"
#include "qcustomplot.h"
#include <QVector>

class LIBUNMIX_EXPORT GraphBase
{
public:
    GraphBase();
    virtual ~GraphBase();
    virtual QCustomPlot *buildPlot(const arma::vec xvalues, const arma::vec yvalues);
    static void convertToVectorDouble(const arma::vec Avec, std::vector<double> &Aout);

protected:
    std::string m_title;
    std::string m_xlabel;
    std::string m_ylabel;
    QCustomPlot *m_customPlot;
};

class LIBUNMIX_EXPORT HistogramSamples : public GraphBase
{
public:
    HistogramSamples(const unsigned int numDisplay, const unsigned int Nclass);
    virtual ~HistogramSamples();
    QCustomPlot *getCustomPlot();
    virtual QCustomPlot *buildPlot(const unsigned int Nmc, const unsigned int Nbi, const MCMCChains SamplesToAnalyze, const std::string name, arma::vec &meanOfSamples);

private:
    unsigned int m_Nclass;
    unsigned int m_numDisplay;
    QVector<QCPGraph*> m_graphVec;
};

class LIBUNMIX_EXPORT AbundanceMaps
{
public:
    AbundanceMaps(const unsigned int NumMaps, const unsigned int sizeX, const unsigned int sizeY);
    ~AbundanceMaps();
    QVector<QCustomPlot*> showMap(const arma::mat meanSamples);

private:
    unsigned int m_NumMaps;
    unsigned int m_SizeX;
    unsigned int m_SizeY;
    QVector<QCustomPlot*> m_custplotVec;
    QVector<QCPColorMap*> m_mapVec;
};

#endif // LIBGRAPH_H
