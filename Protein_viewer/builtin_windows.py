from PySide6 import QtGui
from PySide6.QtCore import Qt
from PySide6.QtWidgets import QAbstractItemView, QTableWidget

import matplotlib
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure

matplotlib.use("Qt5Agg")

class HistCanvas(FigureCanvasQTAgg):
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(dpi=dpi)
        self.axes = fig.add_subplot(111)
        self.set_text()
        super(HistCanvas, self).__init__(fig)

    def set_text(self):
        self.axes.set_ylabel("Distance")
        self.axes.set_xlabel("Residues")
        self.axes.set_title("CA - CA distance distribution")

class AptamerListTableWidget(QTableWidget):
    def __init__(self, parent):
        super().__init__()
        self.parent = parent
        self._configure_table()

    def _configure_table(self):
        self.setMaximumWidth(545)
        self.table_font = QtGui.QFont("Consolas", 12)
        self.setFont(self.table_font)
        self.verticalHeader().setDefaultSectionSize(80)
        self.horizontalHeader().setDefaultSectionSize(240)
        self.verticalHeader().setDefaultAlignment(Qt.AlignCenter)
        self.setShowGrid(False)
        self.horizontalHeader().hide()
        self.setSelectionMode(QAbstractItemView.SingleSelection)
        self.setFocusPolicy(Qt.NoFocus)
        self.setStyleSheet(
            "QTableWidget::item:selected{background-color: white; color: black; border: 1px solid black;}"
        )

        # Init signals
        self.cellClicked.connect(self.draw_hist)

    def draw_hist(self, row, column):
        self.parent.canvas.axes.cla()
        self.parent.canvas.set_text()
        self.parent.canvas.axes.bar(self.parent.hist_labels[row], self.parent.hist_data[row])
        self.parent.canvas.draw()

class Visualisation(FigureCanvasQTAgg):
    def __init__(self):
        self.colors = [
    (0.964, 0.211, 0.192),
    (0.192, 0.376, 0.964),
    (0.192, 0.964, 0.713),
    (0.964, 0.917, 0.192),
    (0.964, 0.592, 0.192),
    (0.807, 0.274, 0.764),
]   

        fig = Figure()
        self.ax = fig.add_subplot(projection = '3d')
        self.ax.axis('off')
        super(Visualisation, self).__init__(fig)
