#include "fenresult.h"

FenSynthResult::FenSynthResult(arma::vec meanVal, int numbEndm, QWidget *parent, QCustomPlot *customPlot)
    : QDialog(parent)
{
    m_buttonQuit = new QPushButton("Quit", this);
    m_mainLayout = new QVBoxLayout();

    setWindowTitle("Abundance results (means)");
    m_resultsDispText.setAlignment(Qt::AlignCenter);

    m_mainLayout->addWidget(&m_resultsDispText);
    m_mainLayout->addWidget(m_buttonQuit);
    m_mainLayout->addWidget(customPlot);
    setLayout(m_mainLayout);
    displayResults(meanVal, numbEndm);

    connect(m_buttonQuit, SIGNAL(clicked()), this, SLOT(closeIt()));

}

FenSynthResult::~FenSynthResult()
{
    std::cout << "Debug : close" << std::endl;
}

void FenSynthResult::closeIt()
{
    this->close();
}

void FenSynthResult::displayResults(arma::vec &meanVal, int &numbEndm)
{
    QString text2Display = "";
    for (int i = 0; i < numbEndm; ++i)
    {
        text2Display += "a_" + QString::number(i+1) + " = " + QString::number(meanVal(i)) + "\n";
    }
    m_resultsDispText.setText(text2Display);

}

FenImageResult::FenImageResult(unsigned int numbEndm, QWidget *parent, QVector<QCustomPlot*> customPlotVec)
    : QDialog(parent)
{
    m_buttonQuit = new QPushButton("Quit");
    m_mainLayout = new QVBoxLayout();
    m_subGridLayout = new QGridLayout();
    m_mainLayout->addWidget(m_buttonQuit);
    m_mainLayout->addLayout(m_subGridLayout);

    unsigned int j = 0, k = 0;
    for (unsigned int i = 0; i < numbEndm; ++i)
    {
        m_subGridLayout->addWidget(customPlotVec.at(i), j, k);
        ++k;
        if (k == 3)
        {
            j++;
            k = 0;
        }
    }

    setLayout(m_mainLayout);

    connect(m_buttonQuit, SIGNAL(clicked()), this, SLOT(closeIt()));
}

FenImageResult::~FenImageResult()
{

}

void FenImageResult::closeIt()
{
    this->close();
}
