#ifndef FENPRIMAIRE_H
#define FENPRIMAIRE_H

#include <QWidget>
#include <QPushButton>
#include <QHBoxLayout>
#include "fensynthpixel.h"
#include "fenimagehyper.h"

class FenPrimaire : public QWidget
{
    Q_OBJECT
public:
    explicit FenPrimaire(QWidget *parent = nullptr);
    ~FenPrimaire();

public slots:
    void generateSynthPixel();
    void unmixRealImage();

private:
    std::unique_ptr<QPushButton> m_buttonImage;
    std::unique_ptr<QPushButton> m_buttonSynthPixel;
    std::unique_ptr<FenSynthPixel> m_winSynth;
    std::unique_ptr<FenImageHyper>  m_winImage;

};

#endif // FENPRIMAIRE_H
