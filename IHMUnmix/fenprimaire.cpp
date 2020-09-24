#include "fenprimaire.h"

FenPrimaire::FenPrimaire(QWidget *parent) : QWidget(parent)
{
    setWindowTitle("What do you want to Unmix ?");
    QHBoxLayout *mainLayout = new QHBoxLayout;

    m_buttonImage = std::make_unique<QPushButton>("Unmix image", this);
    m_buttonSynthPixel = std::make_unique<QPushButton>("Unmix synthetic pixel", this);

    mainLayout->addWidget(m_buttonImage.get());
    mainLayout->addWidget(m_buttonSynthPixel.get());

    this->setLayout(mainLayout);

    QObject::connect(m_buttonSynthPixel.get(), SIGNAL(clicked()), this, SLOT(generateSynthPixel()));
    QObject::connect(m_buttonImage.get(), SIGNAL(clicked()), this, SLOT(unmixRealImage()));
}

FenPrimaire::~FenPrimaire()
{

}

void FenPrimaire::generateSynthPixel()
{
    m_winSynth = std::make_unique<FenSynthPixel>();
    m_winSynth->show();
}

void FenPrimaire::unmixRealImage()
{
    m_winImage = std::make_unique<FenImageHyper>();
    m_winImage->show();
}
