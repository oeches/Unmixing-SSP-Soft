#include "fensynthpixel.h"
#include "ui_fensynthpixel.h"
#include "fenresult.h"
#include "libpixelsynth.h"
#include "unmixing.h"
#include "libmcmc.h"
#include "libgraph.h"
#include <memory>
#include <qglobal.h>
#include <QMessageBox>
#include <QTime>
#include <armadillo>

FenSynthPixel::FenSynthPixel(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::FenSynthPixel)
{
    ui->setupUi(this);

    m_numbEndm = 2;
    m_maxNumbEndm = 6;
    m_spectraLoaded = false;
    m_pxSynth = std::make_unique<libPixelSynth>();

    ui->RadioButtonGroup->setId(ui->radioButtR2, 0);
    ui->RadioButtonGroup->setId(ui->radioButtR3, 1);
    ui->RadioButtonGroup->setId(ui->radioButtR4, 2);
    ui->RadioButtonGroup->setId(ui->radioButtR5, 3);
    ui->RadioButtonGroup->setId(ui->radioButtR6, 4);

    connect(ui->lineEditLib, SIGNAL(editingFinished()), this, SLOT(textHasChanged()));
}

FenSynthPixel::~FenSynthPixel()
{
    if (ui != nullptr)
    {
        delete ui;
        ui = nullptr;
    }
}

void FenSynthPixel::on_ButtonBrowseLib_clicked()
{
    QString libFile = QFileDialog::getOpenFileName(this, "Open spectral library", QString(), "Data files (*.h5)");
    ui->lineEditLib->setText(libFile);
    checkLoadSpectraLib();
}

void FenSynthPixel::textHasChanged()
{
    checkLoadSpectraLib();
}

void FenSynthPixel::on_RadioButtonGroup_buttonClicked(int id)
{
    switch(id)
    {
        case 0:
            ui->a_3DoubleSpinBox->setEnabled(false);
            ui->a_4DoubleSpinBox->setEnabled(false);
            ui->a_5DoubleSpinBox->setEnabled(false);
            ui->a_6DoubleSpinBox->setEnabled(false);
            m_numbEndm = 2;
        break;

        case 1:
            ui->a_3DoubleSpinBox->setEnabled(true);
            ui->a_4DoubleSpinBox->setEnabled(false);
            ui->a_5DoubleSpinBox->setEnabled(false);
            ui->a_6DoubleSpinBox->setEnabled(false);
            m_numbEndm = 3;
        break;

        case 2:
            ui->a_3DoubleSpinBox->setEnabled(true);
            ui->a_4DoubleSpinBox->setEnabled(true);
            ui->a_5DoubleSpinBox->setEnabled(false);
            ui->a_6DoubleSpinBox->setEnabled(false);
            m_numbEndm = 4;
        break;

        case 3:
            ui->a_3DoubleSpinBox->setEnabled(true);
            ui->a_4DoubleSpinBox->setEnabled(true);
            ui->a_5DoubleSpinBox->setEnabled(true);
            ui->a_6DoubleSpinBox->setEnabled(false);
            m_numbEndm = 5;
        break;

        case 4:
            ui->a_3DoubleSpinBox->setEnabled(true);
            ui->a_4DoubleSpinBox->setEnabled(true);
            ui->a_5DoubleSpinBox->setEnabled(true);
            ui->a_6DoubleSpinBox->setEnabled(true);
            m_numbEndm = 6;
        break;

        default:
        break;
    }
}

void FenSynthPixel::on_ButtonGenerate_clicked()
{
    arma::mat rMat(m_maxNumbEndm,1,arma::fill::zeros);
    rMat.rows(0,m_numbEndm-1) = arma::randu(m_numbEndm,1);
    double sumAb = arma::accu(rMat);
    rMat /= sumAb;

    // precision : 2 digits
    arma::Mat<int> rMat_int = arma::conv_to<arma::Mat<int>>::from(100*rMat);
    rMat_int(m_numbEndm-1,0) = 100 - arma::accu(rMat_int.rows(0,m_numbEndm-2));
    rMat = arma::conv_to<arma::mat>::from(rMat_int)/100;
    ui->a_1DoubleSpinBox->setValue(rMat(0,0));
    ui->a_2DoubleSpinBox->setValue(rMat(1,0));
    ui->a_3DoubleSpinBox->setValue(rMat(2,0));
    ui->a_4DoubleSpinBox->setValue(rMat(3,0));
    ui->a_5DoubleSpinBox->setValue(rMat(4,0));
    ui->a_6DoubleSpinBox->setValue(rMat(5,0));
}

void FenSynthPixel::on_ButtonUnmix_clicked()
{
    try
    {
        if (!m_spectraLoaded)
        {
            throw std::logic_error("Please load the spectra first");
        }
    }
    catch (std::logic_error &str)
    {
        QString msgError = QString::fromStdString(str.what());
        QMessageBox::warning(this, "No spectra", msgError, QMessageBox::Close);
        return;
    }
    // abundances
    arma::mat AbundancesRead(m_maxNumbEndm,1);
    arma::mat AbundancesMat(m_numbEndm,1);
    AbundancesRead(0,0) = ui->a_1DoubleSpinBox->value();
    AbundancesRead(1,0) = ui->a_2DoubleSpinBox->value();
    AbundancesRead(2,0) = ui->a_3DoubleSpinBox->value();
    AbundancesRead(3,0) = ui->a_4DoubleSpinBox->value();
    AbundancesRead(4,0) = ui->a_5DoubleSpinBox->value();
    AbundancesRead(5,0) = ui->a_6DoubleSpinBox->value();
    AbundancesMat = AbundancesRead.rows(0,m_numbEndm-1);

    double sumRes = arma::accu(AbundancesMat);
    try
    {
        if (sumRes != 1.)
        {
            std::ostringstream strs;
            strs << "Problem with generated abundances : sum equal to " << sumRes;
            throw (std::out_of_range(strs.str()));
        }
    } catch (std::out_of_range &str_e)
    {
        QString msgError = QString::fromStdString(str_e.what());
        QMessageBox::warning(this, "Abundances coefficients !", msgError, QMessageBox::Close);
        return;
    }

    // sampling step
    unsigned int samplingStep = 2;

    // noise variance
    double noiseVar = ui->noiseVarianceDoubleSpinBox->value();

    try
    {
        m_pxSynth->generateSynthPixel(m_numbEndm, samplingStep, AbundancesMat, noiseVar);
    }
    catch(std::string &str_e)
    {
        QString msgError = QString::fromStdString(str_e);
        QMessageBox::warning(this, "Wrong number of endmember", msgError, QMessageBox::Close);
        m_pxSynth->clearPixel();
        return;
    }

    unsigned int Nmc = ui->iterationsSpinBox->value();
    unsigned int Nbi = ui->burnInSpinBox->value();

    std::shared_ptr<MCMCChains> MCMCAbund = std::make_shared<MCMCChains>(Nmc, m_numbEndm, "MCMC Abundances");
    std::shared_ptr<MCMCChains> MCMCVar = std::make_shared<MCMCChains>(Nmc, m_numbEndm, "MCMC Variance");

    unmixing<libPixelSynth>(*m_pxSynth, Nmc, *MCMCAbund, *MCMCVar);
    m_pxSynth->clearPixel();

    // Histograms results
    std::unique_ptr<HistogramSamples> histoAbund = std::make_unique<HistogramSamples>
            (m_numbEndm, 50);
    arma::vec meanValAb(m_numbEndm);
    //std::shared_ptr<QCustomPlot> customPlotHistoAb = std::make_shared<QCustomPlot>();
    QCustomPlot *customPlotHistoAb = new QCustomPlot();
    customPlotHistoAb = histoAbund->buildPlot(Nmc, Nbi, *MCMCAbund, "Ab", meanValAb);
    std::unique_ptr<FenSynthResult> windowResult = std::make_unique<FenSynthResult>(meanValAb, m_numbEndm, this, customPlotHistoAb);
    windowResult->setWindowTitle("Synthetic pixel results");
    windowResult->setGeometry(100, 100, 500, 400);
    windowResult->exec();
    
}

void FenSynthPixel::checkLoadSpectraLib()
{
    try
    {
        int retRead = 0;
        std::string libFileRead = ui->lineEditLib->text().toStdString();
        retRead = m_pxSynth->loadSpectraLib(libFileRead);
        if (retRead != 0)
        {
            throw std::string ("Please select a correct hdf5 file that contains spectra data");
        }
        if (m_pxSynth->getSpectralLibrarySize() <= 1)
        {
            throw std::out_of_range ("Please select a correct hdf5 file that contains at least 2 spectra");
        }
    } catch (std::string &str)
    {
        QString msgError = QString::fromStdString(str);
        QMessageBox::warning(this, "Wrong hdf5 file", msgError, QMessageBox::Close);
        m_spectraLoaded = false;
        return;
    }
    catch (std::out_of_range &str_e)
    {
        QString msgError = QString::fromStdString(str_e.what());
        QMessageBox::warning(this, "Too few spectra", msgError, QMessageBox::Close);
        m_spectraLoaded = false;
        return;
    }

    if (m_pxSynth->getSpectralLibrarySize() > m_maxNumbEndm)
    {
        m_pxSynth->setSpectralLibrarySize(m_maxNumbEndm);
        // No more than 6 endmember spectra for the moment
    }
    enableRadioButtonGroup(m_pxSynth->getSpectralLibrarySize());
    m_spectraLoaded = true;
}

void FenSynthPixel::enableRadioButtonGroup(unsigned int maxEndm)
{
    ui->radioButtR2->setEnabled(true);
    switch(maxEndm)
    {
        case 2:
            ui->radioButtR3->setEnabled(false);
            ui->radioButtR4->setEnabled(false);
            ui->radioButtR5->setEnabled(false);
            ui->radioButtR6->setEnabled(false);
            break;
        case 3:
            ui->radioButtR3->setEnabled(true);
            ui->radioButtR4->setEnabled(false);
            ui->radioButtR5->setEnabled(false);
            ui->radioButtR6->setEnabled(false);
            break;
        case 4:
            ui->radioButtR3->setEnabled(true);
            ui->radioButtR4->setEnabled(true);
            ui->radioButtR5->setEnabled(false);
            ui->radioButtR6->setEnabled(false);
            break;
        case 5:
            ui->radioButtR3->setEnabled(true);
            ui->radioButtR4->setEnabled(true);
            ui->radioButtR5->setEnabled(true);
            ui->radioButtR6->setEnabled(false);
            break;
        case 6:
            ui->radioButtR3->setEnabled(true);
            ui->radioButtR4->setEnabled(true);
            ui->radioButtR5->setEnabled(true);
            ui->radioButtR6->setEnabled(true);
            break;
        default:
            ui->radioButtR2->setEnabled(false);
            ui->radioButtR3->setEnabled(false);
            ui->radioButtR4->setEnabled(false);
            ui->radioButtR5->setEnabled(false);
            ui->radioButtR6->setEnabled(false);
            break;
    }
}
