#ifndef FENIMAGEHYPER_H
#define FENIMAGEHYPER_H

#include <QWidget>
#include <memory>
#include "libImage.h"

namespace Ui {
class FenImageHyper;
}

class FenImageHyper : public QWidget
{
    Q_OBJECT
public:
    explicit FenImageHyper(QWidget *parent = nullptr);
    ~FenImageHyper();
    void checkLoad();

public slots:
    void on_ButtonBrowseLib_clicked();
    void on_ButtonBrowseIm_clicked();
    void on_ButtonUnmix_clicked();

private:
    Ui::FenImageHyper *ui;
    bool m_spectraLoaded;
    bool m_imageLoaded;
    std::unique_ptr<LibImage> m_imageHyper;
};

#endif // FENIMAGEHYPER_H
