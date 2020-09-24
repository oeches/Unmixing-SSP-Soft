#include <QApplication>
#include "fenprimaire.h"

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);

    FenPrimaire fenetre;
    fenetre.show();

    return app.exec();
}
