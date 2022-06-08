from PySide6 import QtGui
from PySide6.QtGui import QFont
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QApplication,
    QMainWindow,
    QMenuBar,
    QHBoxLayout,
    QWidget,
    QVBoxLayout,
    QComboBox,
    QTableWidgetItem,
    QLineEdit,
    QSplitter,
    QFileDialog,
    QLabel,
    QGridLayout
)

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import sys

from protein import Protein
from builtin_windows import (
    AptamerListTableWidget, 
    HistCanvas, 
    Visualisation
)

class Label(QLabel):
    def __init__(self, parent=None, text=None):
        super().__init__(parent=parent, text=text)
        self.setFont(QFont("Consolas", 12))

class MenuBar(QMenuBar):
    def __init__(self, parent):
        super().__init__()
        self.parent = parent
        self._protein = None
        self._add_menus()

    def _add_menus(self):
        file_menu = self.addMenu("&File")
        open_action = QtGui.QAction(QtGui.QIcon("open.png"), "&Open", self.parent)
        open_action.setShortcut("Ctrl+O")
        open_action.setStatusTip("Open PDB file")
        open_action.triggered.connect(self.read_pdb_file_action)
        exit_action = QtGui.QAction(QtGui.QIcon("exit.png"), "&Exit", self.parent)
        exit_action.setShortcut("Ctrl+Q")
        exit_action.setStatusTip("Exit application")
        exit_action.triggered.connect(self.parent.close)
        file_menu.addAction(open_action)
        file_menu.addAction(exit_action)

    def read_pdb_file_action(self):
        pdb_path = QFileDialog.getOpenFileName(self, "Open Image", "", "PDB Files (*.pdb)")[0]
        self.parent.file_name.setText(f"{Path(pdb_path).stem} protein structure")
        self.parent.protein = Protein(pdb_path)
        self.parent.chain_1_items, self.parent.chain_2_items = [], []
        for chain in self.parent.protein.chains:
            self.parent.chain_1.addItem(chain)
            self.parent.chain_1_items.append(chain)
            self.parent.chain_2.addItem(chain)
            self.parent.chain_2_items.append(chain)
        for method in self.parent.protein.methods:
            self.parent.method.addItem(method)
        self.parent.chain_1_value = self.parent.chain_1.currentText()
        self.parent.chain_2_value = self.parent.chain_2.currentText()
        self.parent.chain_1.model().item(0).setEnabled(True)
        self.parent.chain_2.model().item(0).setEnabled(False)
        self.parent.chain_2.setCurrentIndex(1)
        self.parent.chain_2_value = self.parent.chain_2.currentText()

        self.parent.aptamer_size_value = "5"
        self.parent.aptamer_size.setText(self.parent.aptamer_size_value)

        self.parent.show_protein_structure()
        self.parent.visualisation()

class OptionLayout(QGridLayout):
    def __init__(self, parent=None):
        super().__init__(parent)

    def addWidget(self, widget, row, column, alignment=Qt.Alignment()):
        widget.setMaximumWidth(130)
        widget.setMaximumHeight(20)
        super().addWidget(widget, row, column, alignment)

class MainWindow(QMainWindow):
    def __init__(self):
        super(MainWindow, self).__init__()
        self.init_ui()

    def init_ui(self):
        self.widget = QWidget()
        self.setCentralWidget(self.widget)
        self.setWindowTitle("Protein Viewer")

        self.table = AptamerListTableWidget(self)
        self.canvas = HistCanvas(self, dpi=100)
        self.visual = Visualisation()

        self.splitter = QSplitter(Qt.Vertical)
        self.splitter.addWidget(self.visual)
        self.splitter.addWidget(self.canvas)
        
        self.chain_1 = QComboBox()
        self.chain_2 = QComboBox()
        self.aptamer_size = QLineEdit()
        self.method = QComboBox()

        self.options_layout = OptionLayout()
        self.options_layout.addWidget(Label(self, "Chain 1"), 0, 0)
        self.options_layout.addWidget(Label(self, "Chain 2"), 0, 1)
        self.options_layout.addWidget(Label(self, "Size"), 0, 2)
        self.options_layout.addWidget(Label(self, "Method"), 0, 3)
        self.options_layout.addWidget(self.chain_1, 1, 0)
        self.options_layout.addWidget(self.chain_2, 1, 1)
        self.options_layout.addWidget(self.aptamer_size, 1, 2)
        self.options_layout.addWidget(self.method, 1, 3)

        self.right_layout = QVBoxLayout()
        self.right_layout.addLayout(self.options_layout, 0)
        self.right_layout.addWidget(self.table, 1)

        self.left_layout = QVBoxLayout()
        self.left_layout.addWidget(self.splitter, 0)

        self.vert_layout = QHBoxLayout()
        self.vert_layout.addLayout(self.left_layout, 0)
        self.vert_layout.addLayout(self.right_layout, 1)

        self.file_name = Label(self)
        self.file_name.setAlignment(Qt.AlignCenter)

        self.layout = QVBoxLayout()
        self.layout.addWidget(self.file_name, 0)
        self.layout.addLayout(self.vert_layout, 1)

        self.widget.setLayout(self.layout)

        self.setMenuBar(MenuBar(self))

        # Init slots
        self.chain_1.currentTextChanged.connect(self.chain_1_changed)
        self.chain_2.currentTextChanged.connect(self.chain_2_changed)
        self.aptamer_size.textChanged.connect(self.aptamer_size_changed)

    def chain_1_changed(self, value):
        self.chain_2.model().item(self.chain_1_items.index(self.chain_1_value)).setEnabled(True)
        for i, chain in enumerate(self.protein.chains):
            if chain == value:
                self.chain_2.model().item(i).setEnabled(False)
        self.chain_1_value = value
        self.show_protein_structure()

    def chain_2_changed(self, value):
        self.chain_1.model().item(self.chain_1_items.index(self.chain_2_value)).setEnabled(True)
        for i, chain in enumerate(self.protein.chains):
            if chain == value:
                self.chain_1.model().item(i).setEnabled(False)
        self.chain_2_value = value
        self.show_protein_structure()

    def aptamer_size_changed(self, value):
        self.aptamer_size_value = value
        self.show_protein_structure()

    def show_protein_structure(self):
        kmer_strings_with_distances = sorted(
            self.protein.get_kmer_strings_with_distances(
                self.chain_1_value,
                self.chain_2_value,
                int(self.aptamer_size_value),
            ),
            key=lambda x: x[-1],
        )

        self.table_data, self.hist_data, self.hist_labels = (
            [],
            [],
            [],
        )

        for kmer_pair in kmer_strings_with_distances:
            self.table_data.append([f"{kmer_pair[0]}\n{kmer_pair[1]}", str(kmer_pair[3])])
            self.hist_data.append(kmer_pair[2])
            self.hist_labels.append([f"{char}\n{kmer_pair[1][i]}" for i, char in enumerate(kmer_pair[0])])

        self.table.setRowCount(len(self.table_data))
        self.table.setColumnCount(len(self.table_data[0]))
        for i, table_row in enumerate(self.table_data):
            for j, value in enumerate(self.table_data[i]):
                item = QTableWidgetItem(value)
                item.setTextAlignment(Qt.AlignCenter)
                item.setFlags(Qt.ItemIsSelectable | Qt.ItemIsEnabled)
                self.table.setItem(i, j, item)

        self.table.setCurrentItem(self.table.item(0, 0))
        self.table.draw_hist(0, 0)

    def visualisation(self):
        v = self.protein.get_alpha()
        for idx, chain_id in enumerate(self.protein.chains):
            vect = v[chain_id]
            res = []
            prev = vect[0, 3]
            cur = [vect[0, :]]
            for el in vect[1:, :]:
                if el[3] - prev > 1:
                    cur = np.array(cur)
                    res.append(cur)
                    cur = [el]
                else:
                    cur.append(el)
                prev = el[3]
            cur = np.array(cur)
            res.append(cur)

            for el in res:
                self.visual.ax.plot3D(el[:, 0], el[:, 1], el[:, 2], color=self.visual.colors[idx], linewidth=0.5)
                self.visual.ax.scatter3D(vect[:, 0], vect[:, 1], vect[:, 2], color=self.visual.colors[idx], s=3)
                
        plt.show()

if __name__ == "__main__":
    app = QApplication([])

    widget = MainWindow()
    widget.resize(1280, 720)
    widget.show()

    sys.exit(app.exec())