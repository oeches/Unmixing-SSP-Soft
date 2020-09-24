#ifndef FENSYNTHRESULT_H
#define FENSYNTHRESULT_H

#include <QWidget>
#include <QPushButton>
#include <QLabel>
#include <QDialog>
#include <QVBoxLayout>
#include <armadillo>
#include <memory>
#include "qcustomplot.h"

class FenSynthResult : public QDialog
{
    Q_OBJECT

public:
    explicit FenSynthResult(arma::vec meanVal, int numbEndm, QWidget *parent, QCustomPlot *customPlot);
    ~FenSynthResult();
    void displayResults(arma::vec &meanVal, int &numbEndm);

public slots:
    void closeIt();

private:
    QLabel m_resultsDispText;
    QPushButton *m_buttonQuit;
    QVBoxLayout *m_mainLayout;

};

class FenImageResult : public QDialog
{
    Q_OBJECT

public:
    explicit FenImageResult(unsigned int numbEndm, QWidget *parent, QVector<QCustomPlot*> customPlotVec);
    ~FenImageResult();

public slots:
    void closeIt();

private:
    QPushButton *m_buttonQuit;
    QVBoxLayout *m_mainLayout;
    QGridLayout *m_subGridLayout;

};

#endif // FENSYNTHRESULT_H
