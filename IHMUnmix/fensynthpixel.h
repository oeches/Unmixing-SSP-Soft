#ifndef FENSYNTHPIXEL_H
#define FENSYNTHPIXEL_H

#include <QDialog>
#include <QFileDialog>
#include <memory>
#include "libpixelsynth.h"

namespace Ui {
class FenSynthPixel;
}

class FenSynthPixel : public QWidget
{
    Q_OBJECT

public:
    explicit FenSynthPixel(QWidget *parent = nullptr);
    ~FenSynthPixel();
    void checkLoadSpectraLib();
    void enableRadioButtonGroup(unsigned int maxEndm);

public slots:
    void on_ButtonBrowseLib_clicked();
    void on_RadioButtonGroup_buttonClicked(int id);
    void on_ButtonGenerate_clicked();
    void on_ButtonUnmix_clicked();
    void textHasChanged();

private:
    Ui::FenSynthPixel *ui;
    std::unique_ptr<libPixelSynth> m_pxSynth;
    unsigned int m_numbEndm;
    unsigned int m_maxNumbEndm;
    bool m_spectraLoaded;

};

#endif // FENSYNTHPIXEL_H
