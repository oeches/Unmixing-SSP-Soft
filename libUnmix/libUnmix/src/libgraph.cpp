#include "include/libgraph.h"
#include "include/libUnmix_global.h"
#include <iostream>
#include <armadillo>
#include <QVector>
#include "include/libmcmc.h"
#include "include/qcustomplot.h"

GraphBase::GraphBase()
{
    m_title = "Your title";
    m_xlabel = "Your x_label here";
    m_ylabel = "Your y_label here";
    m_customPlot = new QCustomPlot();
    m_customPlot->plotLayout()->clear();
}

GraphBase::~GraphBase()
{

}

QCustomPlot *GraphBase::buildPlot(const arma::vec xvalues, const arma::vec yvalues)
{
    std::vector<double> xvalues_vec, yvalues_vec;
    GraphBase::convertToVectorDouble(xvalues, xvalues_vec);
    GraphBase::convertToVectorDouble(yvalues, yvalues_vec);

    QVector<double> x_vec(xvalues_vec.size()), y_vec(yvalues_vec.size());
    x_vec = QVector<double>::fromStdVector(xvalues_vec);
    y_vec = QVector<double>::fromStdVector(yvalues_vec);

    QCPAxisRect *oneAxis = new QCPAxisRect(m_customPlot);
    QCPGraph *oneGraph = m_customPlot->addGraph(oneAxis->axis(QCPAxis::atBottom),oneAxis->axis(QCPAxis::atLeft));

    oneGraph->setData(x_vec, y_vec);
    // first we create and prepare a text layout element:
    QCPTextElement *title = new QCPTextElement(m_customPlot);
    title->setText(QString::fromStdString(m_title));
    title->setFont(QFont("sans", 12, QFont::Bold));
    // then we add it to the main plot layout:
    m_customPlot->plotLayout()->insertRow(0); // insert an empty row above the axis rect
    m_customPlot->plotLayout()->addElement(0, 0, title); // place the title in the empty cell we've just created

    oneGraph->keyAxis()->setLabel(QString::fromStdString(m_xlabel));
    oneGraph->valueAxis()->setLabel(QString::fromStdString(m_ylabel));
    oneGraph->rescaleAxes();

    return m_customPlot;

}

void GraphBase::convertToVectorDouble(const arma::vec Avec, std::vector<double> &Aout)
{
    unsigned int nElem = Avec.n_elem;
    for (unsigned int i = 0; i < nElem; ++i)
    {
        Aout.at(i) = Avec(i);
    }
}

HistogramSamples::HistogramSamples(const unsigned int numDisplay, const unsigned int Nclass) : GraphBase()
{
    m_title = "Posterior distribution of MCMC samples";
    m_xlabel = "Sample values";
    m_ylabel = "Posterior distribution";
    m_Nclass = Nclass;
    m_numDisplay = numDisplay;

    QCPLayoutGrid *subLayout = new QCPLayoutGrid;
    m_customPlot->plotLayout()->addElement(0,0,subLayout);
    unsigned int j = 0, k = 0;
    for (unsigned int i = 0; i < m_numDisplay; ++i)
    {
        QCPAxisRect *testAxis = new QCPAxisRect(m_customPlot);
        subLayout->addElement(j, k, testAxis);

        m_graphVec.push_back(m_customPlot->addGraph(testAxis->axis(QCPAxis::atBottom),
                                        testAxis->axis(QCPAxis::atLeft)));
        ++k;
        if (k == 3)
        {
            j++;
            k = 0;
        }
    }

    QCPTextElement *title = new QCPTextElement(m_customPlot);
    title->setText(QString::fromStdString(m_title));
    title->setFont(QFont("sans", 12, QFont::Bold));

    m_customPlot->plotLayout()->insertRow(0); // insert an empty row above the axis rect
    m_customPlot->plotLayout()->addElement(0, 0, title); // place the title in the empty cell we've just created
}

HistogramSamples::~HistogramSamples()
{
}

QCustomPlot *HistogramSamples::getCustomPlot()
{
    return m_customPlot;
}

QCustomPlot *HistogramSamples::buildPlot(const unsigned int Nmc, const unsigned int Nbi, const MCMCChains SamplesToAnalyze, const std::string name, arma::vec &meanOfSamples)
{
    arma::mat TSamplesR = SamplesToAnalyze.getNSamples( (Nmc-1)-Nbi,Nbi);
    m_Nclass = (unsigned int) std::sqrt(Nmc-Nbi+1);
    arma::vec Tsamples(TSamplesR.n_rows);
    arma::uvec histVal(m_Nclass);
    arma::vec histValD(m_Nclass);
    std::vector<double> histValVec(m_Nclass);
    std::vector<double> abcissesVec(m_Nclass);

    for (unsigned int i = 0; i < m_numDisplay; ++i)
    {
        Tsamples = TSamplesR.col(i);

        double xmin = Tsamples.min();
        double xmax = Tsamples.max();
        double step = (xmax - xmin)/m_Nclass;
        histVal = arma::hist(Tsamples, m_Nclass);
        histValD = arma::conv_to<arma::vec>::from(histVal);
        histValD /= step*(Nmc - Nbi);
        meanOfSamples(i) = arma::mean(Tsamples);
        GraphBase::convertToVectorDouble(histValD, histValVec);
        int l = 0;
        double step_2 = (xmax - xmin)/(m_Nclass-1);
        std::generate(abcissesVec.begin(), abcissesVec.end(), [&] () -> double {return xmin+(l++)*step_2;} );
        abcissesVec.back() = xmax;
        QVector<double> x_vec(abcissesVec.size()), y_vec(histValVec.size());
        x_vec = QVector<double>::fromStdVector(abcissesVec);
        y_vec = QVector<double>::fromStdVector(histValVec);

        if (m_numDisplay > 1)
        {
            m_ylabel = "Post. distri. of " + name + "_" + std::to_string(i+1);
        }
        else
        {
            m_ylabel = "Post. distri. of " + name;
        }
        m_graphVec.at(i)->setData(x_vec, y_vec);
        m_graphVec.at(i)->keyAxis()->setLabel(QString::fromStdString(m_xlabel));
        m_graphVec.at(i)->valueAxis()->setLabel(QString::fromStdString(m_ylabel));
        m_graphVec.at(i)->rescaleAxes();

    }
    return m_customPlot;
}

AbundanceMaps::AbundanceMaps(const unsigned int NumMaps, const unsigned int sizeX, const unsigned int sizeY)
{
    m_NumMaps = NumMaps;
    m_SizeX = sizeX;
    m_SizeY = sizeY;

    unsigned int j = 0, k = 0;
    for (unsigned int i = 0; i < m_NumMaps; ++i)
    {
        m_custplotVec.push_back(new QCustomPlot);
        m_mapVec.push_back(new QCPColorMap(m_custplotVec.at(i)->xAxis,
                                          m_custplotVec.at(i)->yAxis));
        ++k;
        if (k == 3)
        {
            j++;
            k = 0;
        }
        QString title_text = "Abundance map of elem. " + QString::number(i+1);
        QCPTextElement *title = new QCPTextElement(m_custplotVec.at(i));
        title->setText(title_text);
        title->setFont(QFont("sans", 12, QFont::Bold));
        // then we add it to the main plot layout:
        m_custplotVec.at(i)->plotLayout()->insertRow(0); // insert an empty row above the axis rect
        m_custplotVec.at(i)->plotLayout()->addElement(0, 0, title); // place the title in the empty cell we've just created
    }
}

AbundanceMaps::~AbundanceMaps()
{

}

QVector<QCustomPlot*> AbundanceMaps::showMap(const arma::mat meanSamples)
{
    std::vector<double> im2show(m_SizeX*m_SizeY);
    arma::vec meanA(m_SizeX*m_SizeY);

    for (unsigned int indR = 0; indR < m_NumMaps; ++indR)
    {
        meanA = meanSamples.col(indR);
        // conversion to std
        GraphBase::convertToVectorDouble(meanA, im2show);
        std::transform(im2show.begin(), im2show.end(), im2show.begin(),
                       [](double d) -> double {return (255*d);});
        const std::vector<unsigned char> im2showInt(im2show.begin(), im2show.end());

        m_mapVec.at(indR)->data()->setSize(m_SizeX, m_SizeY);
        m_mapVec.at(indR)->data()->setRange(QCPRange(0, m_SizeX), QCPRange(0, m_SizeY));
        for (unsigned int i =0; i < m_SizeX; ++i)
        {
            for (unsigned int j = 0; j < m_SizeY; ++j)
            {
                m_mapVec.at(indR)->data()->setCell(j,m_SizeX-1-i,im2show.at(j+i*m_SizeY));
            }
        }
        m_mapVec.at(indR)->setGradient(QCPColorGradient::gpGrayscale);
        m_mapVec.at(indR)->setInterpolate(false);
        m_mapVec.at(indR)->setTightBoundary(false);
        m_mapVec.at(indR)->rescaleDataRange(true);
        m_custplotVec.at(indR)->rescaleAxes();
        m_custplotVec.at(indR)->replot();
    }

    return m_custplotVec;
}
