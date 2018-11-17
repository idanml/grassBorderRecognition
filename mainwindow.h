#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();
    QString data() const;
    QString image() const;
    QString persents() const;
    QString feedback() const;
    QString fileName;
private slots:
    void on_pushButton_history_clicked();

    void on_pushButton_main_clicked();

    void on_browse_btn_clicked();

    void on_analize_btn_clicked();

    void on_correct_btn_clicked();

    void on_incorrect_btn_clicked();

private:
    Ui::MainWindow *ui;

    enum colum {
        DATE, IMAGE, PERSENT, FEEDBACK
    };
};

#endif // MAINWINDOW_H
