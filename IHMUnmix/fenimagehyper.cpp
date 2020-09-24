#include "fenimagehyper.h"
#include "ui_FenImageHyper.h"
#include "libmcmc.h"
#include "libgraph.h"
#include "unmixing.h"
#include "fenresult.h"
#include <QFileDialog>
#include <QProgressDialog>
#include <QMessageBox>
#include <memory>

FenImageHyper::FenImageHyper(QWidget *parent) : QWidget(parent), ui(new Ui::FenImageHyper)
{
    ui->setupUi(this);
    m_spectraLoaded = false;
    m_imageLoaded = false;
    m_imageHyper = std::make_unique<LibImage>();
}

FenImageHyper::~FenImageHyper()
{
    if (ui != nullptr)
    {
        delete ui;
        ui = nullptr;
    }
}

void FenImageHyper::on_ButtonBrowseLib_clicked()
{
    QString libFile = QFileDialog::getOpenFileName(this, "Open spectral library", QString(), "Data files (*.h5)");
    ui->lineEditSpectralLib->setText(libFile);
}

void FenImageHyper::on_ButtonBrowseIm_clicked()
{
    QString imageFile = QFileDialog::getOpenFileName(this, "Open image file", QString(), "Data files (*.tif)");
    ui->lineEditImage->setText(imageFile);
}

void FenImageHyper::on_ButtonUnmix_clicked()
{
    try
    {
        checkLoad();
        if (!m_spectraLoaded || !m_imageLoaded)
        {
            throw std::logic_error("Please load the spectra and image first");
        }
    }
    catch (std::logic_error &str)
    {
        QString msgError = QString::fromStdString(str.what());
        QMessageBox warnSpectra(QMessageBox::Warning, "No Data", msgError, QMessageBox::Close, this);
        return;
    }

    unsigned int Nmc = ui->iterationsSpinBox->value();
    unsigned int Nbi = ui->burnInSpinBox->value();

    unsigned int numPixels = m_imageHyper->getNPixels();
    unsigned int numEnd = m_imageHyper->getEndmNumb();

    std::shared_ptr<MCMCChains> MCMCAbund = std::make_shared<MCMCChains>(Nmc, numEnd, "MCMC Abundances");
    std::shared_ptr<MCMCChains> MCMCVar = std::make_shared<MCMCChains>(Nmc, numEnd, "MCMC Variance");

    arma::mat valMeanAbund(numPixels,numEnd);
    bool wasCanceled = false;

    // create QDialog for progress bar
    QProgressDialog windowProgress("Unmixing image", "Cancel", 0, numPixels);
    windowProgress.setModal(true);
    windowProgress.show();
    for (unsigned int indPx = 0; indPx < numPixels; ++indPx)
    {
        m_imageHyper->getPixel(indPx);

        unmixing<LibImage>(*m_imageHyper, Nmc, *MCMCAbund, *MCMCVar);

        valMeanAbund.row(indPx) = arma::mean(MCMCAbund->getNSamples(Nmc-Nbi-1, Nbi));
        windowProgress.setValue(indPx+1);
        if (windowProgress.wasCanceled())
        {
            wasCanceled = true;
            break;
        }

    }
    windowProgress.close();

    if (!wasCanceled)
    {
        std::unique_ptr<AbundanceMaps> map_result = std::make_unique<AbundanceMaps>(numEnd, m_imageHyper->getDimX(), m_imageHyper->getDimY());
        QVector<QCustomPlot*> customPlotVec;
        customPlotVec = map_result->showMap(valMeanAbund);
        std::unique_ptr<FenImageResult> fen_map_result = std::make_unique<FenImageResult>(numEnd, this, customPlotVec);
        fen_map_result->setWindowTitle("Abundance maps");
        fen_map_result->setGeometry(100, 100, 500, 400);
        fen_map_result->exec();
    }
}

void FenImageHyper::checkLoad()
{
    try
    {
        int retRead = 0;
        std::string libFileReadSpectra = ui->lineEditSpectralLib->text().toStdString();
        std::string libFileReadImage = ui->lineEditImage->text().toStdString();
        try
        {
            retRead = m_imageHyper->loadImage(libFileReadImage, libFileReadSpectra);
        }
        catch(std::string &str_e)
        {
            QString msgError = QString::fromStdString(str_e);
            QMessageBox::warning(this, "Image file error", msgError, QMessageBox::Close);
            m_imageLoaded = false;
            return;
        }

        if (retRead != 0)
        {
            throw std::string ("Please select a correct hdf5 file that contains spectra data");
        }
        if (m_imageHyper->getEndmNumb() <= 1)
        {
            throw std::out_of_range ("Please select a correct hdf5 file that contains at least 2 spectra");
        }
        if (m_imageHyper->getDimX() <= 1 || m_imageHyper->getDimY() <= 1)
        {
            throw std::out_of_range ("Please select a correct image file");
        }
    }
    catch (std::string &str)
    {
        QString msgError = QString::fromStdString(str);
        QMessageBox::warning(this, "Wrong hdf5 file", msgError, QMessageBox::Close);
        m_spectraLoaded = false;
        return;
    }
    catch (std::out_of_range &str_e)
    {
        QString msgError = QString::fromStdString(str_e.what());
        QMessageBox::warning(this, "Wrong spectra or image", msgError, QMessageBox::Close);
        m_spectraLoaded = false;
        m_imageLoaded = false;
        return;
    }

    m_spectraLoaded = true;
    m_imageLoaded = true;
}
