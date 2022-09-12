import sys
from collections import defaultdict
import sqlite3
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import QDir
from PyQt5.QtGui import QStandardItemModel, QIcon
from PyQt5.QtWidgets import *
from pathlib import Path
import pyqtgraph as pg
from pyqtgraph import PlotWidget

import bloom_filter_algorithm
import smith_waterman_algorithm
from interface import Ui_MainWindow
import fm_index_algorithm
import kmer_algorithm
import aho_corasick_algorithm


# TO DO
# Veritabanı
# Arayüzdeki küçük hatalar başlık hataları vs.

class MainWindow(QMainWindow, Ui_MainWindow):
    def __init__(self):
        QMainWindow.__init__(self)
        self.setupUi(self)
        self.main_window = QMainWindow()  # QMainWindow instance
        self.ui = Ui_MainWindow()  # Arayüzü kullanabilmek için instance
        self.ui.setupUi(self.main_window)

        # self.setWindowIcon(QIcon(":/icon/icons/Logo.jpg"))

        self.ui.stackedWidget.setCurrentWidget(
            self.ui.main_page)  # Stacked widget'ın current widget'ı mainwindow yapıldı.

        # butonların fonksiyonları atandı
        self.ui.btn_main.clicked.connect(lambda: self.goToPage(self.ui.main_page))
        self.ui.btn_bloom_aho.clicked.connect(lambda: self.goToPage(self.ui.bloom_aho_page))
        self.ui.btn_fm_kmer.clicked.connect(lambda: self.goToPage(self.ui.fm_page))
        self.ui.btn_bloom_aho.clicked.connect(lambda: self.goToPage(self.ui.bloom_aho_page))
        self.ui.btn_smith.clicked.connect(lambda: self.goToPage(self.ui.smith_page))

        # FM INDEX BUTTONS
        self.ui.fm_btn_upload_txt.clicked.connect(self.fm_load_txt)
        self.ui.fm_btn_upload_ptrn.clicked.connect(self.fm_load_ptrn)
        self.ui.fm_btn_search.clicked.connect(self.fm_search)
        self.ui.fm_btn_clear.clicked.connect(self.fm_clear)
        self.ui.fm_btn_calc_kmer.clicked.connect(self.calc_kmer)
        self.ui.fm_btn_calc_kmer.clicked.connect(self.btn_kmer)
        self.ui.fm_btn_calc_kmer.clicked.connect(self.kmer_freq_top)
        self.ui.fm_btn_calc_kmer.clicked.connect(lambda: self.kmer_draw_graph_1())
        self.ui.fm_btn_calc_kmer.clicked.connect(lambda: self.kmer_draw_graph_2())
        self.ui.fm_btn_calc_kmer.clicked.connect(lambda: self.kmer_draw_graph_3())

        # AHO CORASICK BUTTONS
        self.ui.aho_btn_upload_txt.clicked.connect(self.aho_load_txt)
        self.ui.aho_btn_upload_ptrn.clicked.connect(self.aho_load_ptrn)
        self.ui.aho_btn_search.clicked.connect(self.aho_search)
        self.ui.aho_btn_clear.clicked.connect(self.aho_clear)

        # BLOOM FILTER BUTTONS
        self.ui.bloom_btn_upload_txt.clicked.connect(self.bloom_load_txt)
        self.ui.bloom_btn_upload_ptrn.clicked.connect(self.bloom_load_ptrn)
        self.ui.bloom_btn_search.clicked.connect(self.bloom_search)
        self.ui.bloom_btn_clear.clicked.connect(self.bloom_clear)

        # SMITH WATERMAN BUTTONS
        self.ui.smith_btn_upload_row.clicked.connect(self.smith_load_row)
        self.ui.smith_btn_upload_col.clicked.connect(self.smith_load_col)
        self.ui.smith_btn_search.clicked.connect(self.smith_waterman_search)
        self.ui.smith_btn_search.clicked.connect(self.smith_draw_matrix)
        self.ui.smith_btn_clear.clicked.connect(self.smith_clear)

    def show(self):
        self.main_window.show()

    def goToPage(self, widget):
        self.ui.stackedWidget.setCurrentWidget(widget)

    # FM INDEX İÇİN GEREKLİ KODLAR BÖLÜMÜ
    txt = ""
    ptrn = ""
    sonuc = []
    k = 0
    f_name_txt = ""
    f_name_ptrn = ""

    def fm_clear(self):
        self.ui.fm_output.clear()
        self.ui.fm_kmer_output.clear()
        self.ui.fm_kmer_output_top5.clear()
        self.ui.fm_lbl_txt.clear()
        self.ui.fm_lbl_ptrn.clear()
        self.ui.kmer_graph1.clear()
        self.ui.kmer_graph2.clear()
        self.ui.kmer_graph3.clear()

    # ALIGNMENT
    def fm_search(self):
        t = self.txt
        p = self.ptrn
        count = 0
        output = []
        run_time = 0.0

        if self.f_name_txt == "":
            QMessageBox.about(main_window, "ERROR", "No FASTA file for text sequence selected. Please select one.")
            QMessageBox.show()

        elif self.f_name_ptrn == "":
            QMessageBox.about(main_window, "ERROR", "No FASTA file for pattern sequence selected. Please select one.")
            QMessageBox.show()
        else:
            result = fm_index_algorithm.align(t, p)
            output = result[0]
            run_time = result[1]
            # print(result)

            self.ui.fm_output.setText("The start position of pattern sequence is: ")

            for i in output:
                # self.textBrowser_3.append("i: {}".format(i))
                self.ui.fm_output.append("- {}".format(i))
                count += 1
            fm_index_algorithm.insertVaribleIntoTable(p, count)
            self.ui.fm_output.append(
                "\n-----------------------------------------------------------------------------------------------")
            self.ui.fm_output.append("\n {} start position found".format(count))
            self.ui.fm_output.append(
                "\n-----------------------------------------------------------------------------------------------")
            self.ui.fm_output.append("\n Run time is: {}".format(run_time))

    def fm_load_txt(self):
        dialog = QFileDialog(directory="D:\Masaüstü\FASTA files")
        dialog.setWindowTitle("Choose FASTA file for text sequence.")
        dialog.setNameFilter("FNA files (*.fna)")
        dialog.setFilter(QDir.Files)

        str_fname = ""

        if dialog.exec_():
            fm_txt_file_name = dialog.selectedFiles()
            self.f_name_txt = fm_txt_file_name

            if fm_txt_file_name[0].endswith('.fna'):
                with open(fm_txt_file_name[0], 'r') as f:
                    text = f.read()
                    self.txt = text
                    f.close()
            else:
                pass

        str_fname = str_fname.join(fm_txt_file_name)

        if fm_txt_file_name:
            fname = Path(str_fname)
            self.ui.fm_lbl_txt.setText("Selected file: {}".format(fname.name))
        else:
            pass
            self.ui.fm_lbl_txt.setText("No file is chosen. Please choose a file.")

            return text

    def fm_load_ptrn(self):
        dialog = QFileDialog(directory="D:\Masaüstü\FASTA files")
        dialog.setWindowTitle("Choose FASTA file for pattern sequence.")
        dialog.setNameFilter("FNA files (*.fna)")
        dialog.setFilter(QDir.Files)

        str_fname = ""

        if dialog.exec_():
            fm_ptrn_file_name = dialog.selectedFiles()
            self.f_name_ptrn = fm_ptrn_file_name

            if fm_ptrn_file_name[0].endswith('.fna'):
                with open(fm_ptrn_file_name[0], 'r') as f:
                    pattern = f.read()
                    self.ptrn = pattern
                    f.close()
            else:
                pass

        str_fname = str_fname.join(fm_ptrn_file_name)

        if fm_ptrn_file_name:
            fname = Path(str_fname)
            self.ui.fm_lbl_ptrn.setText("Selected file: {}".format(fname.name))
        else:
            self.ui.fm_lbl_ptrn.setText("No file was choosed. Please choose a file.")

            return pattern

    # KMER
    def btn_kmer(self):
        sequence = self.ptrn
        get_k = self.ui.fm_input_kmer_value.text()

        if self.f_name_ptrn == "":
            QMessageBox.about(main_window, "ERROR", "No FASTA file for pattern sequence selected. Please select one.")
            QMessageBox.show()

        elif get_k == "":
            QMessageBox.about(main_window, "ERROR",
                              "K value for K-Mer algorithm can't be 0. Please enter a valid value.")
            QMessageBox.show()
        else:
            self.calc_kmer()
            self.kmer_freq_top()
            self.kmer_draw_graph_1()
            self.kmer_draw_graph_2()
            self.kmer_draw_graph_3()

    def calc_kmer(self):
        sequence = self.ptrn
        get_k = self.ui.fm_input_kmer_value.text()

        ksize = int(get_k)
        self.k = ksize
        kmers = kmer_algorithm.build_kmers(sequence, ksize)

        for kmer, freq in kmers.items():
            txt = "kmer: {} - freq: {}".format(kmer, freq)
            self.ui.fm_kmer_output.append(txt)

        # verScrollBar = self.txt_kmer_freq.setVerticalScrollBar()
        # verScrollBar.setValue(verScrollBar.minimum())

    def kmer_freq_top(self):

        sequence = self.ptrn
        get_k = self.ui.fm_input_kmer_value.text()
        ksize = int(get_k)
        self.k = ksize

        mc = kmer_algorithm.kmer_freq_top(sequence, ksize)

        for kmer, freq in mc:
            txt = "kmer: {} - freq: {}".format(kmer, freq)
            self.ui.fm_kmer_output_top5.append("{}".format(txt))
            sonuc = kmer + " freq " + str(freq)
            kmer_algorithm.insertVaribleIntoTable(ksize, sonuc)

    def kmer_draw_graph_1(self):
        sequence = self.ptrn
        get_k = self.ui.fm_input_kmer_value.text()
        ksize = int(get_k)
        self.k = ksize

        x = []
        y = []

        sonuc = kmer_algorithm.kmer_graph_1(sequence, ksize)
        x = sonuc[0]
        y = sonuc[1]

        pen = pg.mkPen(color=(255, 0, 0), width=1, style=QtCore.Qt.DotLine)
        styles = {'color': 'r', 'font-size': '15px'}
        self.ui.kmer_graph1.setLabel('left', '1/freq', **styles)
        self.ui.kmer_graph1.setLabel('bottom', 'freqs', **styles)
        self.ui.kmer_graph1.setTitle("K-Mer Period Graph", color="red", size="10pt")
        self.ui.kmer_graph1.plot(x, y, pen=pen)
        self.ui.kmer_graph1.updateMatrix()

    def kmer_draw_graph_2(self):
        sequence = self.ptrn
        get_k = self.ui.fm_input_kmer_value.text()
        ksize = int(get_k)
        self.k = ksize

        x = []
        y = []

        sonuc2 = kmer_algorithm.kmer_graph_2(sequence, ksize)
        x = sonuc2[0]
        y = sonuc2[1]

        # pen = pg.mkPen(color="blue", width=1, style=QtCore.Qt.SolidLine)
        styles = {'color': 'r', 'font-size': '15px'}
        self.ui.kmer_graph2.setLabel('left', 'freq', **styles)
        self.ui.kmer_graph2.setLabel('bottom', 'k-mer no', **styles)
        self.ui.kmer_graph2.setTitle("K-Mer Frequence Graph", color="b", size="10pt")
        # self.graph_2.plot(x, y,pen=pen, symbol='+', symbolSize=5, symbolBrush=('b'))
        self.ui.kmer_graph2.plot(x, y, pen=None, symbol='o')
        self.ui.kmer_graph2.updateMatrix()

    def kmer_draw_graph_3(self):
        sequence = self.ptrn
        get_k = self.ui.fm_input_kmer_value.text()
        ksize = int(get_k)
        self.k = ksize

        x = []
        y = []

        sonuc3 = kmer_algorithm.kmer_graph_3(sequence, ksize)
        x = sonuc3[0]
        y = sonuc3[1]

        # pen = pg.mkPen(color="white", width=1, style=QtCore.Qt.SolidLine)
        bargraph = pg.BarGraphItem(x=x, height=y, width=0.6, brush='g')
        styles = {'color': 'r', 'font-size': '15px'}
        self.ui.kmer_graph3.setLabel('left', 'freq', **styles)
        self.ui.kmer_graph3.setLabel('bottom', 'kmer no', **styles)
        self.ui.kmer_graph3.setTitle("Top 5 Frequence Graph", color="g", size="10pt")
        self.ui.kmer_graph3.plot(x, y)
        self.ui.kmer_graph3.addItem(bargraph)
        self.ui.kmer_graph3.updateMatrix()

    # AHO CORASICK KODLARI
    file_path = ""
    a_ptrn_fnames = []
    aho_txt = ""
    aho_ptrns = []
    aho_fname_txt = ""
    aho_fname_ptrn = ""

    def aho_load_txt(self):

        dialog = QFileDialog(directory="D:\Masaüstü\FASTA files")
        dialog.setWindowTitle("Choose FASTA file for text sequence.")
        dialog.setNameFilter("FNA files (*.fna)")
        dialog.setFilter(QDir.Files)
        dialog.setFileMode(QFileDialog.ExistingFiles)

        str_fname = ""

        if dialog.exec_():
            aho_txt_file_name = dialog.selectedFiles()
            self.aho_fname_txt = aho_txt_file_name

            if aho_txt_file_name[0].endswith('.fna'):
                with open(aho_txt_file_name[0], 'r') as f:
                    text = f.read()
                    self.aho_txt = text
                    f.close()
            else:
                pass

        str_fname = str_fname.join(aho_txt_file_name)

        if aho_txt_file_name:
            fname = Path(str_fname)
            self.ui.aho_lbl_txt.setText("Selected file: {}".format(fname.name))
        else:
            pass
            self.ui.aho_lbl_txt.setText("No file is chosen. Please choose a file.")

            return text

    def aho_load_ptrn(self):

        caption = "Choose FASTA file for pattern sequence."
        dir = 'D:\Masaüstü\FASTA files'
        filter = "FNA files (*.fna)"
        fnames = []

        aho_ptrn_file_names = QFileDialog.getOpenFileNames(None, caption, dir, filter)[0]
        for file in aho_ptrn_file_names:
            fname = Path(file)  # file yolunu aldı
            self.file_path = fname
            fnames.append(fname.name)  # file adını fnames listesine attı
            self.a_ptrn_fnames = fnames
            aho_ptrn_fnames = ' '.join(map(str, fnames))  # listeyi stringe dönüştürdü
            self.aho_fname_ptrn = aho_ptrn_fnames

            if file.endswith(".fna"):
                with open(file, 'r') as f:
                    pattern = f.read()
                    self.aho_ptrns.append(pattern)
                    f.close()
            else:
                pass

        if aho_ptrn_file_names:
            size_aho_patterns = len(fnames)
            self.ui.aho_lbl_ptrn.setText("{} patterns selected".format(size_aho_patterns))
        else:
            self.ui.aho_lbl_ptrn.setText("No file was choosed. Please choose a file.")
            return pattern

    def aho_search(self):
        t = self.aho_txt
        p = []
        sonuclar = {}

        for ptrn in self.aho_ptrns:
            p.append(ptrn)

        if self.aho_fname_txt == "":
            QMessageBox.about(main_window, "ERROR", "No FASTA file for text sequence selected. Please select one.")
            QMessageBox.show()

        elif self.aho_fname_ptrn == "":
            QMessageBox.about(main_window, "ERROR", "No FASTA file for pattern sequences selected. Please select one.")
            QMessageBox.show()
        else:
            aho = aho_corasick_algorithm.AhoCorasick(p)
            result = aho.run_aho_corasick(t, p)

            j = 0
            for pattern in result:
                # pattern = self.a_ptrn_fnames[j]
                pat = self.a_ptrn_fnames[j]
                j = j + 1
                for i in result[pattern]:
                    self.ui.aho_output.append(
                        "Pattern {} appears from {}.index to {}.index".format(pat, i, i + len(pattern) - 1))
                    conn = sqlite3.connect("database.db")
                    cursor = conn.cursor()
                    sqlcommand = """
                               CREATE TABLE IF NOT EXISTS CAROSICK
                               (
                               WORD TEXT PRIMARY KEY,
                               START INTEGER,
                               FINISH INTEGER 
                               )      
                               """
                    cursor.execute(sqlcommand)
                    aho_corasick_algorithm.AhoCorasick.insertVaribleIntoTable(pattern, i, i + len(pattern) - 1)
                    conn.commit()

    def aho_clear(self):
        self.ui.aho_output.clear()
        self.ui.aho_lbl_txt.clear()
        self.ui.aho_lbl_ptrn.clear()

    # BLOOM KODLAR
    probability = 0
    bloom_ptrn_fnames = []
    bloom_txt = ""
    bloom_ptrns = []
    bloom_fname_txt = ""
    bloom_fname_ptrn = ""

    def bloom_load_txt(self):
        dialog = QFileDialog(directory="D:\Masaüstü\FASTA files")
        dialog.setWindowTitle("Choose FASTA file for text sequence.")
        dialog.setNameFilter("FNA files (*.fna)")
        dialog.setFilter(QDir.Files)
        dialog.setFileMode(QFileDialog.ExistingFiles)

        str_fname = ""

        if dialog.exec_():
            bloom_txt_file_name = dialog.selectedFiles()
            self.bloom_fname_txt = bloom_txt_file_name

            if bloom_txt_file_name[0].endswith('.fna'):
                with open(bloom_txt_file_name[0], 'r') as f:
                    text = f.read()
                    self.bloom_txt = text
                    f.close()
            else:
                pass

        str_fname = str_fname.join(bloom_txt_file_name)

        if bloom_txt_file_name:
            fname = Path(str_fname)
            self.ui.bloom_lbl_txt.setText("Selected file: {}".format(fname.name))
        else:
            pass
            self.ui.bloom_lbl_txt.setText("No file is chosen. Please choose a file.")

            return text

    def bloom_load_ptrn(self):
        caption = "Choose FASTA file for pattern sequence."
        dir = 'D:\Masaüstü\FASTA files'
        filter = "FNA files (*.fna)"
        fnames = []

        bloom_ptrn_file_names = QFileDialog.getOpenFileNames(None, caption, dir, filter)[0]
        for file in bloom_ptrn_file_names:
            fname = Path(file)  # file yolunu aldı
            fnames.append(fname.name)  # file adını fnames listesine attı
            bloom_ptrn_fnames = ' '.join(map(str, fnames))  # listeyi stringe dönüştürdü
            self.bloom_fname_ptrn = bloom_ptrn_fnames

            if file.endswith(".fna"):
                with open(file, 'r') as f:
                    pattern = f.read()
                    self.bloom_ptrns.append(pattern)
                    f.close()
            else:
                pass

        if bloom_ptrn_file_names:
            size_bloom_patterns = len(fnames)
            # self.ui.bloom_lbl_ptrn.setText(bloom_ptrn_fnames)
            self.ui.bloom_lbl_ptrn.setText("{} patterns selected".format(size_bloom_patterns))
        else:
            self.ui.bloom_lbl_ptrn.setText("No file was choosed. Please choose a file.")
            return pattern

    def bloom_search(self):
        t = self.txt
        p = self.ptrn
        n = 20
        t = self.bloom_txt
        p = []

        get_prob = self.ui.bloom_prob_value.text()
        prob = float(get_prob)
        self.probability = prob

        for ptrn in self.bloom_ptrns:
            p.append(ptrn)

        if self.bloom_fname_txt == "":
            QMessageBox.about(main_window, "ERROR", "No FASTA file for text sequence selected. Please select one.")
            QMessageBox.show()

        elif self.bloom_fname_ptrn == "":
            QMessageBox.about(main_window, "ERROR", "No FASTA file for pattern sequences selected. Please select one.")
            QMessageBox.show()
        else:
            bloomf = bloom_filter_algorithm.BloomFilter(n, prob)
            test_words = bloomf.run_bloom_filter(t, p)

            self.ui.bloom_output.append("Size of bit array:{}".format(bloomf.size))
            self.ui.bloom_output.append("False positive Probability:{}".format(bloomf.fp_prob))
            self.ui.bloom_output.append("Number of hash functions:{}".format(bloomf.hash_count))

            for word in test_words:
                if bloomf.check(word):
                    if word in p:
                        self.ui.bloom_output.append("'{}' is a false positive!".format(word))
                        output = 'false positive'
                    else:
                        self.ui.bloom_output.append("'{}' is probably present!".format(word))
                        output = 'probably present'
                else:
                    output = 'definitely not present'
                    self.ui.bloom_output.append("'{}' is definitely not present!".format(word))
                bloom_filter_algorithm.BloomFilter.insertVaribleIntoTable(word, output)
            # print(result)

    def bloom_clear(self):
        self.ui.bloom_lbl_txt.clear()
        self.ui.bloom_lbl_ptrn.clear()
        self.ui.bloom_output.clear()

    # SMITH-WATERMAN KODLARI
    smith_row = ""
    smith_col = ""
    match = 0
    mismatch = 0
    gap = 0
    f_smithname_row = ""
    f_smithname_col = ""

    # bu stringin her bir karakteri bir satırın başına gelir
    def smith_load_row(self):
        dialog = QFileDialog(directory="D:\Masaüstü\FASTA files")
        dialog.setWindowTitle("Choose FASTA file for text sequence.")
        dialog.setNameFilter("FNA files (*.fna)")
        dialog.setFilter(QDir.Files)

        str_fname = ""

        if dialog.exec_():
            smith_row_file_name = dialog.selectedFiles()
            self.f_smithname_row = smith_row_file_name

            if smith_row_file_name[0].endswith('.fna'):
                with open(smith_row_file_name[0], 'r') as f:
                    text = f.read()
                    self.smith_row = text
                    f.close()
            else:
                pass

        str_fname = str_fname.join(smith_row_file_name)

        if smith_row_file_name:
            fname = Path(str_fname)
            self.ui.smith_lbl_output_row.setText("Selected file: {}".format(fname.name))
        else:
            pass
            self.ui.smith_lbl_output_row.setText("No file is chosen. Please choose a file.")

            return text

    # bu stringin her bir karakteri bir kolonun üstüne gelir
    def smith_load_col(self):
        dialog = QFileDialog(directory="D:\Masaüstü\FASTA files")
        dialog.setWindowTitle("Choose FASTA file for pattern sequence.")
        dialog.setNameFilter("FNA files (*.fna)")
        dialog.setFilter(QDir.Files)

        str_fname = ""

        if dialog.exec_():
            smith_col_file_name = dialog.selectedFiles()
            self.f_smithname_col = smith_col_file_name

            if smith_col_file_name[0].endswith('.fna'):
                with open(smith_col_file_name[0], 'r') as f:
                    pattern = f.read()
                    self.smith_col = pattern
                    f.close()
            else:
                pass

        str_fname = str_fname.join(smith_col_file_name)

        if smith_col_file_name:
            fname = Path(str_fname)
            self.ui.smith_lbl_output_col.setText("Selected file: {}".format(fname.name))
        else:
            self.ui.smith_lbl_output_col.setText("No file was choosed. Please choose a file.")

            return pattern

    def smith_waterman_search(self):

        if self.f_smithname_col == "":
            QMessageBox.about(main_window, "ERROR", "No FASTA file for text sequence selected. Please select one.")
            QMessageBox.show()

        elif self.f_smithname_col == "":
            QMessageBox.about(main_window, "ERROR", "No FASTA file for pattern sequence selected. Please select one.")
            QMessageBox.show()
        else:

            get_match = self.ui.smith_match_value.text()
            match = int(get_match)
            self.match = match

            get_mismatch = self.ui.smith_mismatch_value.text()
            mismatch = int(get_mismatch)
            self.mismatch = mismatch

            get_gap = self.ui.smith_gap_value.text()
            gap = int(get_gap)
            self.gap = gap

            score_matrix, start_pos = smith_waterman_algorithm.create_score_matrix(self.smith_row, self.smith_col,
                                                                                   self.match, self.mismatch, self.gap)
            seq1_aligned, seq2_aligned = smith_waterman_algorithm.traceback(score_matrix, start_pos, self.smith_row,
                                                                            self.smith_col)
            alignment_str, idents, gaps, mismatches = smith_waterman_algorithm.alignment_string(seq1_aligned,
                                                                                                seq2_aligned)
            alength = len(seq1_aligned)



            self.ui.smith_output.append("Identities = {0}/{1} ({2:.1%}), Gaps = {3}/{4} ({5:.1%})".format(idents,
                                                                                                          alength,
                                                                                                          idents / alength,
                                                                                                          gaps,
                                                                                                          alength,
                                                                                                          gaps / alength))
            smith_waterman_algorithm.insertVaribleIntoTable(self.smith_row,self.smith_col, match, mismatch, gap, idents / alength, gaps / alength)

            for i in range(0, alength, 60):
                seq1_slice = seq1_aligned[i:i + 60]
                seq2_slice = seq2_aligned[i:i + 60]

                self.ui.smith_output.append("Query {0:<4}  {1}  {2:<4}".format(i + 1, seq1_slice, i + len(seq1_slice)))
                self.ui.smith_output.append("                 {0}".format(alignment_str[i:i + 60]))
                self.ui.smith_output.append("Sbjct  {0:<4}  {1}  {2:<4}".format(i + 1, seq2_slice, i + len(seq2_slice)))

    def smith_draw_matrix(self):
        if self.f_smithname_col == "":
            QMessageBox.about(main_window, "ERROR", "No FASTA file for Row Sequence selected. Please select one.")
            QMessageBox.show()

        elif self.f_smithname_col == "":
            QMessageBox.about(main_window, "ERROR", "No FASTA file for Column Sequence selected. Please select one.")
            QMessageBox.show()
        else:
            seq_row = self.smith_row
            seq_col = self.smith_col
            len_row = len(seq_row)
            len_col = len(seq_col)
            match = self.match
            mismatch = self.mismatch
            gap = self.gap

            # matriste en üstte her kolonun üstüne col'un bir karakteri yazılır
            # matriste dikey tarafa (sola) her satır başına row'un bir karakteri yazılır
            # seq_col = ACTG
            # seq_row = ACTGA

            ust = []
            ust.append("0")
            for i in seq_col:
                ust.append(i)
            size_ust = len(ust)

            sol = []
            sol.append("0")
            for i in seq_row:
                sol.append(i)
            size_sol = len(sol)

            # row ve col stringlerine 0 da eklendiği için matrisin satır ve sütun sayıları +1 yapıp ayarlar
            self.ui.smith_table.setRowCount(len_row + 1)
            self.ui.smith_table.setColumnCount(len_col + 1)

            # kolonların genişlik ve yüksekliğini içindeki dataya göre ayarlar
            header = self.ui.smith_table.horizontalHeader()
            header.setSectionResizeMode(QtWidgets.QHeaderView.ResizeToContents)

            # kolon ve satırların başna gerekleri karakterleri atar
            self.ui.smith_table.setHorizontalHeaderLabels(ust)
            self.ui.smith_table.setVerticalHeaderLabels(sol)

            score_matrix, max_pos = smith_waterman_algorithm.create_score_matrix(seq_row, seq_col, match, mismatch, gap)

            row = 0
            col = 0
            for list in score_matrix:
                for element in list:
                    self.ui.smith_table.setItem(row, col, QtWidgets.QTableWidgetItem(str(element)))
                    col = col + 1
            row = row + 1

    def smith_clear(self):
        self.ui.smith_output.clear()
        self.ui.smith_lbl_output_row.clear()
        self.ui.smith_lbl_output_col.clear()
        self.ui.smith_table.clear()


if __name__ == "__main__":
    app = QApplication(sys.argv)
    main_window = MainWindow()
    main_window.show()
    sys.exit(app.exec_())
