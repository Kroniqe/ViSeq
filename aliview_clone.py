import sys
import os
import shutil
import tempfile
import subprocess
import re
import math
import csv
import numpy as np
from Bio import AlignIO, SeqIO, Phylo
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

# Optional Import for Tree Plotting
try:
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

# Biopython compatibility check for MeltingTemp
try:
    from Bio.SeqUtils import MeltingTemp as mt
except ImportError:
    try:
        import Bio.SeqUtils.MeltingTemp as mt
    except ImportError:
        mt = None

from PyQt6.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, 
                             QHBoxLayout, QFileDialog, QMessageBox, QSplitter,
                             QScrollBar, QAbstractScrollArea, QFrame, QMenu,
                             QTableWidget, QTableWidgetItem, QHeaderView, QDialog,
                             QLabel, QPushButton, QTabWidget, QListWidget, 
                             QAbstractItemView, QSpinBox, QWidgetAction, QLineEdit,
                             QFormLayout, QDialogButtonBox, QGroupBox, QToolBar, 
                             QSizePolicy, QTextEdit, QListWidgetItem, QGridLayout)
from PyQt6.QtCore import Qt, QRect, QSize, QPoint, pyqtSignal, QEvent, QTimer
from PyQt6.QtGui import (QPainter, QColor, QFont, QPen, QBrush, QAction, 
                         QFontMetrics, QKeySequence, QIcon)

# --- Configuration ---
APP_NAME = "ViSeq"
APP_TAGLINE = "Sequence visualization, reimagined"

# --- Colors ---
NUCLEOTIDE_COLORS = {
    b'A': QColor(230, 100, 100),
    b'C': QColor(100, 230, 100),
    b'G': QColor(230, 230, 100),
    b'T': QColor(100, 100, 230),
    b'U': QColor(100, 100, 230),
    b'-': QColor(240, 240, 240),
    b'N': QColor(180, 180, 180),
    b'?': QColor(180, 180, 180)
}
DEFAULT_BG = QColor(255, 255, 255)
TEXT_COLOR = QColor(0, 0, 0)
SELECTION_COLOR = QColor(0, 0, 255, 50)
PRIMER_F_COLOR = QColor(0, 255, 0, 80) # Semi-transparent Green
PRIMER_R_COLOR = QColor(255, 0, 0, 80) # Semi-transparent Red

# Pre-calculate byte-to-string lookup for rendering performance
BYTE_TO_CHAR = {bytes([b]): bytes([b]).decode('ascii', errors='replace') for b in range(256)}

# --- Primer Logic ---
class PrimerDesigner:
    def __init__(self, sequence_str):
        self.seq = sequence_str.upper()
        self.length = len(self.seq)

    def calculate_gc(self, seq):
        if len(seq) == 0: return 0
        return (seq.count('G') + seq.count('C')) / len(seq) * 100

    def check_runs(self, seq):
        # Optimum: no runs >3
        if re.search(r'(A{4,}|C{4,}|G{4,}|T{4,})', seq): return False 
        return True

    def check_3_prime(self, seq):
        if len(seq) < 3: return "Poor"
        end = seq[-3:]
        gc_count = end.count('G') + end.count('C')
        has_clamp = seq[-1] in 'GC'
        # Avoid >2 consecutive G/C at 3'
        if end in ['GGG', 'CCC', 'GCG', 'CGC']: return "Poor"
        
        # Optimum: 1 to 2 G or C at 3' end
        if 1 <= gc_count <= 2 and has_clamp: return "Optimal"
        # Okayish: no more than 3 G/Cs at 3' end
        elif gc_count <= 3: return "Okay"
        return "Poor"

    def analyze_candidate(self, seq_fragment, start_index, direction="F"):
        length = len(seq_fragment)
        gc = self.calculate_gc(seq_fragment)
        tm = 0
        if mt:
            try: tm = mt.Tm_NN(Seq(seq_fragment))
            except: tm = 0 
        
        score = 0
        notes = []
        
        # 1. Length Check (18-22 optimal, 17-25 okay)
        if 18 <= length <= 22: score += 2 
        elif 17 <= length <= 25: score += 1 
        else: return None # Strict Fail
        
        # 2. Tm Check (58-62 optimal, 55-65 okay)
        if 58 <= tm <= 62: score += 2 
        elif 55 <= tm <= 65: score += 1 
        else: return None # Strict Fail
        
        # 3. GC Content (45-55 optimal, 40-60 okay)
        if 45 <= gc <= 55: score += 2 
        elif 40 <= gc <= 60: score += 1 
        else: return None # Strict Fail
        
        # 4. Runs
        if not self.check_runs(seq_fragment):
            score -= 1 
            notes.append("Runs")
            
        # 5. 3' Stability
        end_qual = self.check_3_prime(seq_fragment)
        if end_qual == "Optimal": score += 2
        elif end_qual == "Okay": score += 1
        else: 
            score -= 2
            notes.append("3'End")
            
        return {
            "seq": seq_fragment, "start": start_index, "length": length,
            "tm": round(tm, 2), "gc": round(gc, 1), "score": score,
            "notes": ", ".join(notes), "direction": direction
        }

    def find_primers(self, min_amp=100, max_amp=1000):
        candidates_f = []
        candidates_r = []
        search_seq = self.seq
        
        # Scan for Forward and Reverse candidates
        for l in range(17, 26): 
            for i in range(len(search_seq) - l + 1):
                sub = search_seq[i : i+l]
                res = self.analyze_candidate(sub, i, "F")
                if res: candidates_f.append(res)
                
                sub_sense = search_seq[i : i+l]
                primer_seq = str(Seq(sub_sense).reverse_complement())
                res_r = self.analyze_candidate(primer_seq, i, "R")
                if res_r: 
                    res_r['sense_start'] = i
                    res_r['sense_end'] = i + l
                    candidates_r.append(res_r)
        
        pairs = []
        for f in candidates_f:
            for r in candidates_r:
                dist = r['sense_end'] - f['start']
                if min_amp <= dist <= max_amp:
                    tm_diff = abs(f['tm'] - r['tm'])
                    pair_score = f['score'] + r['score']
                    
                    # Penalty for Tm difference
                    if tm_diff > 2: pair_score -= 1
                    if tm_diff > 5: pair_score -= 2
                    
                    # Annealing Temp Estimate (min Tm - 5)
                    ta = min(f['tm'], r['tm']) - 5
                    
                    pairs.append({
                        "f": f, "r": r, "amp_len": dist,
                        "score": pair_score, "tm_diff": round(tm_diff, 1),
                        "ta": round(ta, 1)
                    })
        
        pairs.sort(key=lambda x: x['score'], reverse=True)
        return candidates_f, candidates_r, pairs

# --- Dialogs ---
class TreeResultDialog(QDialog):
    def __init__(self, newick_str, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Phylogenetic Tree Result")
        self.resize(600, 400)
        self.newick_str = newick_str
        layout = QVBoxLayout()
        layout.addWidget(QLabel("Newick Tree Format:"))
        self.text_area = QTextEdit()
        self.text_area.setPlainText(newick_str)
        layout.addWidget(self.text_area)
        btn_layout = QHBoxLayout()
        save_btn = QPushButton("Save .nwk File")
        save_btn.clicked.connect(self.save_tree)
        btn_layout.addWidget(save_btn)
        if HAS_MATPLOTLIB:
            plot_btn = QPushButton("Show/Save Graphical Tree")
            plot_btn.clicked.connect(self.plot_tree)
            btn_layout.addWidget(plot_btn)
        close_btn = QPushButton("Close")
        close_btn.clicked.connect(self.accept)
        btn_layout.addWidget(close_btn)
        layout.addLayout(btn_layout)
        self.setLayout(layout)
        
    def save_tree(self):
        fname, _ = QFileDialog.getSaveFileName(self, "Save Tree", "", "Newick Tree (*.nwk);;All Files (*)")
        if fname:
            with open(fname, "w") as f:
                f.write(self.newick_str)
                
    def plot_tree(self):
        if not HAS_MATPLOTLIB: return
        try:
            from io import StringIO
            tree = Phylo.read(StringIO(self.newick_str), "newick")
            fig = plt.figure(figsize=(10, 8))
            ax = fig.add_subplot(1, 1, 1)
            Phylo.draw(tree, axes=ax, do_show=False)
            fname, _ = QFileDialog.getSaveFileName(self, "Save Tree Image", "", "PNG Image (*.png);;PDF Document (*.pdf)")
            if fname:
                plt.savefig(fname)
                QMessageBox.information(self, "Saved", f"Image saved to {fname}")
            plt.show()
        except Exception as e:
            QMessageBox.critical(self, "Plot Error", str(e))

class PrimerParamsDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Design Parameters")
        layout = QFormLayout()
        self.min_amp = QSpinBox()
        self.min_amp.setRange(50, 10000)
        self.min_amp.setValue(100)
        self.max_amp = QSpinBox()
        self.max_amp.setRange(50, 10000)
        self.max_amp.setValue(500)
        layout.addRow("Min Amplicon Length:", self.min_amp)
        layout.addRow("Max Amplicon Length:", self.max_amp)
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)
        self.setLayout(layout)
    def get_values(self):
        return self.min_amp.value(), self.max_amp.value()

class PrimerResultsDialog(QDialog):
    def __init__(self, f_primers, r_primers, pairs, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Primer Design Results")
        self.resize(1100, 600)
        self.parent_ref = parent
        layout = QVBoxLayout()
        tabs = QTabWidget()
        
        self.pair_table = self.create_pair_table(pairs)
        self.f_table = self.create_table(f_primers, "Forward")
        self.r_table = self.create_table(r_primers, "Reverse")
        
        tabs.addTab(self.pair_table, f"Valid Pairs ({len(pairs)})")
        tabs.addTab(self.f_table, f"Individual F ({len(f_primers)})")
        tabs.addTab(self.r_table, f"Individual R ({len(r_primers)})")
        
        layout.addWidget(tabs)
        
        btn_layout = QHBoxLayout()
        export_btn = QPushButton("Export Results to Excel/CSV")
        export_btn.clicked.connect(self.export_csv)
        btn_layout.addWidget(export_btn)
        close_btn = QPushButton("Close")
        close_btn.clicked.connect(self.accept)
        btn_layout.addWidget(close_btn)
        
        layout.addLayout(btn_layout)
        layout.addWidget(QLabel("Click a row to navigate and highlight."))
        self.setLayout(layout)

    def export_csv(self):
        fname, _ = QFileDialog.getSaveFileName(self, "Export Primers", "", "CSV Files (*.csv)")
        if fname:
            try:
                with open(fname, 'w', newline='') as csvfile:
                    writer = csv.writer(csvfile)
                    writer.writerow(["--- PAIRED PRIMERS ---"])
                    headers = [self.pair_table.horizontalHeaderItem(c).text() for c in range(self.pair_table.columnCount())]
                    writer.writerow(headers)
                    for i in range(self.pair_table.rowCount()):
                        row = [self.pair_table.item(i, c).text() for c in range(self.pair_table.columnCount())]
                        writer.writerow(row)
                    writer.writerow([])
                    writer.writerow(["--- INDIVIDUAL FORWARD ---"])
                    headers_f = [self.f_table.horizontalHeaderItem(c).text() for c in range(self.f_table.columnCount())]
                    writer.writerow(headers_f)
                    for i in range(self.f_table.rowCount()):
                        row = [self.f_table.item(i, c).text() for c in range(self.f_table.columnCount())]
                        writer.writerow(row)
                    writer.writerow([])
                    writer.writerow(["--- INDIVIDUAL REVERSE ---"])
                    headers_r = [self.r_table.horizontalHeaderItem(c).text() for c in range(self.r_table.columnCount())]
                    writer.writerow(headers_r)
                    for i in range(self.r_table.rowCount()):
                        row = [self.r_table.item(i, c).text() for c in range(self.r_table.columnCount())]
                        writer.writerow(row)
                QMessageBox.information(self, "Export Successful", f"Saved to {fname}")
            except Exception as e:
                QMessageBox.critical(self, "Export Error", str(e))

    def create_pair_table(self, pairs):
        table = QTableWidget()
        table.setColumnCount(9)
        table.setHorizontalHeaderLabels([
            "Score", "Amp Len", "Anneal T", 
            "F Pos", "F Tm", "F GC%", 
            "R Pos", "R Tm", "R GC%"
        ])
        table.setRowCount(len(pairs))
        table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        table.cellClicked.connect(lambda r, c: self.highlight_pair(r, table))
        
        for i, p in enumerate(pairs):
            table.setItem(i, 0, QTableWidgetItem(str(p['score'])))
            table.setItem(i, 1, QTableWidgetItem(str(p['amp_len'])))
            table.setItem(i, 2, QTableWidgetItem(str(p['ta'])))
            f_pos = f"{p['f']['start']}-{p['f']['start']+p['f']['length']}"
            table.setItem(i, 3, QTableWidgetItem(f_pos))
            table.setItem(i, 4, QTableWidgetItem(str(p['f']['tm'])))
            table.setItem(i, 5, QTableWidgetItem(str(p['f']['gc'])))
            r_pos = f"{p['r']['sense_start']}-{p['r']['sense_end']}"
            table.setItem(i, 6, QTableWidgetItem(r_pos))
            table.setItem(i, 7, QTableWidgetItem(str(p['r']['tm'])))
            table.setItem(i, 8, QTableWidgetItem(str(p['r']['gc'])))
            table.item(i, 0).setData(Qt.ItemDataRole.UserRole, p)
        return table

    def create_table(self, data, p_type):
        table = QTableWidget()
        table.setColumnCount(6)
        table.setHorizontalHeaderLabels(["Score", "Sequence", "Len", "Tm", "GC", "Note"])
        table.setRowCount(len(data))
        table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        table.cellClicked.connect(lambda r, c, t=table, pt=p_type: self.highlight_single(r, t, pt))
        for i, p in enumerate(data):
            table.setItem(i, 0, QTableWidgetItem(str(p['score'])))
            seq = QTableWidgetItem(p['seq'])
            seq.setFont(QFont("Consolas", 9))
            table.setItem(i, 1, seq)
            table.setItem(i, 2, QTableWidgetItem(str(p['length'])))
            table.setItem(i, 3, QTableWidgetItem(str(p['tm'])))
            table.setItem(i, 4, QTableWidgetItem(str(p['gc'])))
            table.setItem(i, 5, QTableWidgetItem(p['notes']))
            table.item(i, 0).setData(Qt.ItemDataRole.UserRole, p)
        return table

    def highlight_single(self, row, table, p_type):
        item = table.item(row, 0)
        data = item.data(Qt.ItemDataRole.UserRole)
        if self.parent_ref:
            start = data['start'] if p_type == "Forward" else data['sense_start']
            self.parent_ref.highlight_primer_region(start, data['length'], p_type)

    def highlight_pair(self, row, table):
        item = table.item(row, 0)
        p = item.data(Qt.ItemDataRole.UserRole)
        if self.parent_ref:
            self.parent_ref.highlight_primer_pair(p['f'], p['r'])

# --- Core Data Logic ---
class AlignmentData:
    def __init__(self):
        self.data_matrix = None 
        self.ids = []
        self.history = [] 

    def save_state(self):
        if self.data_matrix is not None:
            snapshot = (self.data_matrix.copy(), list(self.ids))
            self.history.append(snapshot)
            if len(self.history) > 50: self.history.pop(0)

    def undo(self):
        if not self.history: return False
        matrix, ids = self.history.pop()
        self.data_matrix = matrix
        self.ids = ids
        return True

    def clear_history(self):
        self.history = []

    def load_fasta(self, filepath):
        self.clear_history()
        try:
            raw_records = list(SeqIO.parse(filepath, "fasta"))
        except Exception as e:
            raise ValueError(f"Failed to parse FASTA: {str(e)}")

        if not raw_records:
            raise ValueError("File is empty or not a valid FASTA.")

        self.ids = [rec.id for rec in raw_records]
        max_len = max(len(rec.seq) for rec in raw_records)
        
        data_list = []
        for rec in raw_records:
            s_str = str(rec.seq).upper()
            if len(s_str) < max_len:
                s_str = s_str.ljust(max_len, '-')
            # Safe Decode with 'replace'
            safe_str = s_str.encode('ascii', errors='replace').decode('ascii')
            data_list.append(list(safe_str)) 

        self.data_matrix = np.array(data_list, dtype='|S1')
            
    def get_slice(self, row_start, row_end, col_start, col_end):
        if self.data_matrix is None: return None
        return self.data_matrix[row_start:row_end, col_start:col_end]

    def update_base(self, row, col, new_char):
        if self.data_matrix is not None:
            self.save_state() 
            self.data_matrix[row, col] = new_char.upper()

    def delete_columns(self, start_col, end_col):
        if self.data_matrix is None: return
        self.save_state()
        indices = list(range(start_col, end_col + 1))
        self.data_matrix = np.delete(self.data_matrix, indices, axis=1)

    def delete_rows(self, row_indices):
        if self.data_matrix is None: return
        self.save_state()
        self.data_matrix = np.delete(self.data_matrix, row_indices, axis=0)
        for index in sorted(row_indices, reverse=True):
            if 0 <= index < len(self.ids):
                del self.ids[index]

    def to_biopython_obj(self):
        if self.data_matrix is None: return None
        new_records = []
        rows, cols = self.data_matrix.shape
        for i in range(rows):
            seq_chars = self.data_matrix[i].tobytes().decode('utf-8', errors='replace')
            new_records.append(SeqRecord(Seq(seq_chars), id=self.ids[i], description=""))
        return MultipleSeqAlignment(new_records)
    
    def save_fasta(self, filepath):
        if self.data_matrix is None: return
        records = []
        rows, cols = self.data_matrix.shape
        for r in range(rows):
            seq_chars = self.data_matrix[r].tobytes().decode('ascii', errors='replace')
            records.append(SeqRecord(Seq(seq_chars), id=self.ids[r], description=""))
        with open(filepath, "w") as f:
            SeqIO.write(records, f, "fasta")
        
    def reorder_rows(self, new_order_indices):
        if self.data_matrix is None: return
        self.save_state()
        self.data_matrix = self.data_matrix[new_order_indices]
        self.ids = [self.ids[i] for i in new_order_indices]

# --- UI Components ---

class RulerWidget(QWidget):
    def __init__(self, viewer_ref):
        super().__init__()
        self.viewer = viewer_ref
        self.setFixedHeight(30)
        self.setStyleSheet("background-color: #E0E0E0; border-bottom: 1px solid #CCC;")
    
    def paintEvent(self, event):
        if self.viewer.data.data_matrix is None: return
        painter = QPainter(self)
        font = QFont("Arial", 9, QFont.Weight.Bold)
        painter.setFont(font)
        painter.setPen(TEXT_COLOR)
        scroll_x = self.viewer.horizontalScrollBar().value()
        char_w = self.viewer.char_width
        start_col = max(0, scroll_x // char_w)
        visible_cols = (self.width() // char_w) + 2
        
        for c in range(start_col, start_col + visible_cols):
            x_pos = (c * char_w) - scroll_x
            idx = c + 1
            if idx % 10 == 0:
                painter.drawLine(int(x_pos), 10, int(x_pos), 30)
                painter.drawText(int(x_pos) + 2, 12, str(idx))
            elif idx % 5 == 0 and char_w > 5:
                painter.drawLine(int(x_pos), 20, int(x_pos), 30)
            elif char_w > 15:
                painter.drawLine(int(x_pos), 25, int(x_pos), 30)

class DraggableListWidget(QListWidget):
    reorder_signal = pyqtSignal() 

    def __init__(self, data_ref, viewer_ref):
        super().__init__()
        self.data = data_ref
        self.viewer = viewer_ref
        self.setFixedWidth(250)
        self.setStyleSheet("""
            QListWidget { 
                background-color: #F0F0F0; 
                color: black; 
                font-family: Consolas; 
                font-weight: bold;
                border-right: 1px solid #999; 
            }
            QListWidget::item:selected { background-color: #3399FF; color: white; }
        """)
        self.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        self.setDragEnabled(True)
        self.setAcceptDrops(True)
        self.setDragDropMode(QAbstractItemView.DragDropMode.InternalMove)
        self.setDefaultDropAction(Qt.DropAction.MoveAction)
        self.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)

    def populate(self):
        self.clear()
        if self.data.ids:
            for i, seq_id in enumerate(self.data.ids):
                self.addItem(f"{i+1}. {seq_id}")
            self.update_font_size(self.viewer.font_size, self.viewer.char_height)

    def update_font_size(self, size, row_height):
        font = QFont("Consolas", size, QFont.Weight.Bold)
        self.setFont(font)
        for i in range(self.count()):
            item = self.item(i)
            item.setSizeHint(QSize(self.width(), row_height))

    def dropEvent(self, event):
        super().dropEvent(event)
        new_order_ids_clean = []
        for i in range(self.count()):
            text = self.item(i).text()
            if ". " in text:
                clean = text.split(". ", 1)[1]
                new_order_ids_clean.append(clean)
            else:
                new_order_ids_clean.append(text)
        old_ids = self.data.ids
        from collections import defaultdict
        id_map = defaultdict(list)
        for idx, id_ in enumerate(old_ids): id_map[id_].append(idx)
        new_indices = []
        for id_ in new_order_ids_clean:
            if id_map[id_]: new_indices.append(id_map[id_].pop(0))
        self.data.reorder_rows(new_indices)
        self.populate()
        self.reorder_signal.emit()

    def selectionChanged(self, selected, deselected):
        super().selectionChanged(selected, deselected)
        rows = [self.row(item) for item in self.selectedItems()]
        if rows: self.viewer.select_rows(rows)
            
    def set_selected_row(self, row_index):
        self.clearSelection()
        if row_index < self.count():
            item = self.item(row_index)
            item.setSelected(True)
            self.scrollToItem(item)

class SequenceViewer(QAbstractScrollArea):
    row_height_changed = pyqtSignal(int, int) 
    grid_clicked = pyqtSignal(int) 

    def __init__(self, data_ref):
        super().__init__()
        self.data = data_ref
        self.font_size = 10
        self.char_width = 12
        self.char_height = 20
        self.font = QFont("Consolas", self.font_size)
        self.draw_text = True
        self.selected_cells = set()
        self.highlight_f_cells = set()
        self.highlight_r_cells = set()
        self.selection_start = None
        self.setMouseTracking(True)
        self.viewport().setMouseTracking(True)
        self.color_mode = 0 

    def set_global_font_size(self, size):
        self.font_size = size
        self.char_width = int(size * 1.0) 
        self.char_height = int(size * 1.8)
        self.font = QFont("Consolas", self.font_size)
        self.update_scrollbars()
        self.viewport().update()
        self.row_height_changed.emit(self.font_size, self.char_height)

    def set_zoom(self, level):
        new_size = int(10 * level)
        self.set_global_font_size(max(6, new_size))
    
    def set_color_mode(self, mode):
        self.color_mode = mode
        self.viewport().update()

    def update_scrollbars(self):
        if self.data.data_matrix is None: return
        rows, cols = self.data.data_matrix.shape
        total_width = cols * self.char_width
        total_height = rows * self.char_height
        self.horizontalScrollBar().setRange(0, total_width - self.viewport().width())
        self.verticalScrollBar().setRange(0, total_height - self.viewport().height())
        self.horizontalScrollBar().setPageStep(self.viewport().width())
        self.verticalScrollBar().setPageStep(self.viewport().height())

    def resizeEvent(self, event):
        self.update_scrollbars()
        super().resizeEvent(event)

    def paintEvent(self, event):
        if self.data.data_matrix is None: return
        painter = QPainter(self.viewport())
        painter.setFont(self.font)
        scroll_x = self.horizontalScrollBar().value()
        scroll_y = self.verticalScrollBar().value()
        start_col = max(0, scroll_x // self.char_width)
        start_row = max(0, scroll_y // self.char_height)
        rows, cols = self.data.data_matrix.shape
        grid_slice = self.data.get_slice(start_row, min(rows, start_row + (self.viewport().height()//self.char_height)+2), start_col, min(cols, start_col + (self.viewport().width()//self.char_width)+2))
        if grid_slice is None: return
        metrics = QFontMetrics(self.font)
        text_y_offset = (self.char_height + metrics.ascent() - metrics.descent()) // 2
        
        selection_pen = QPen(Qt.GlobalColor.black)
        selection_pen.setStyle(Qt.PenStyle.DashLine)
        selection_pen.setWidth(1)

        for r in range(grid_slice.shape[0]):
            abs_row = start_row + r
            y_pos = (abs_row * self.char_height) - scroll_y
            for c in range(grid_slice.shape[1]):
                abs_col = start_col + c
                x_pos = (abs_col * self.char_width) - scroll_x
                char_byte = grid_slice[r, c]
                
                # Base Color
                if self.color_mode == 0: 
                    bg_color = NUCLEOTIDE_COLORS.get(char_byte, DEFAULT_BG)
                    txt_color = TEXT_COLOR
                elif self.color_mode == 1: 
                    bg_color = DEFAULT_BG
                    txt_color = NUCLEOTIDE_COLORS.get(char_byte, TEXT_COLOR)
                else: 
                    bg_color = DEFAULT_BG
                    txt_color = TEXT_COLOR
                
                if (abs_row, abs_col) in self.highlight_f_cells:
                    bg_color = PRIMER_F_COLOR
                elif (abs_row, abs_col) in self.highlight_r_cells:
                    bg_color = PRIMER_R_COLOR
                elif (abs_row, abs_col) in self.selected_cells:
                    # Keep underlying color but prepare for dashed border
                    pass

                painter.fillRect(x_pos, y_pos, self.char_width, self.char_height, bg_color)
                
                if (abs_row, abs_col) in self.selected_cells:
                    painter.save()
                    painter.setPen(selection_pen)
                    painter.setBrush(Qt.BrushStyle.NoBrush)
                    painter.drawRect(int(x_pos), int(y_pos), self.char_width, self.char_height)
                    painter.restore()
                
                if self.draw_text and self.char_width > 6:
                    # Use lookup table
                    char_str = BYTE_TO_CHAR.get(char_byte, '?')
                    painter.setPen(txt_color)
                    painter.drawText(x_pos, y_pos + text_y_offset, char_str)

    def mousePressEvent(self, event):
        x = event.position().x() + self.horizontalScrollBar().value()
        y = event.position().y() + self.verticalScrollBar().value()
        col = int(x // self.char_width)
        row = int(y // self.char_height)
        self.selection_start = (row, col)
        
        self.highlight_f_cells = set()
        self.highlight_r_cells = set()
        
        self.selected_cells = {(row, col)}
        self.grid_clicked.emit(row)
        self.viewport().update()

    def mouseMoveEvent(self, event):
        if not event.buttons() & Qt.MouseButton.LeftButton: return
        if self.selection_start is None: return
        x = event.position().x() + self.horizontalScrollBar().value()
        y = event.position().y() + self.verticalScrollBar().value()
        col = int(x // self.char_width)
        row = int(y // self.char_height)
        r1, c1 = self.selection_start
        r_min, r_max = min(r1, row), max(r1, row)
        c_min, c_max = min(c1, col), max(c1, col)
        new_selection = set()
        for r in range(r_min, r_max + 1):
            for c in range(c_min, c_max + 1):
                new_selection.add((r, c))
        self.selected_cells = new_selection
        self.viewport().update()

    def wheelEvent(self, event):
        if event.modifiers() == Qt.KeyboardModifier.ControlModifier:
            delta = event.angleDelta().y()
            new_size = self.font_size + (1 if delta > 0 else -1)
            self.set_global_font_size(max(6, new_size))
        else:
            super().wheelEvent(event)

    def keyPressEvent(self, event):
        if event.matches(QKeySequence.StandardKey.Copy):
            self.copy_selection_to_clipboard()
            return
        if not self.selected_cells: return
        key = event.key()
        text = event.text()
        if key == Qt.Key.Key_Delete:
            for r, c in self.selected_cells: self.data.update_base(r, c, "-")
            self.viewport().update()
        elif text and text.upper() in "ACGTU-N":
            for r, c in self.selected_cells: self.data.update_base(r, c, text)
            self.viewport().update()
            
    def copy_selection_to_clipboard(self):
        if not self.selected_cells or self.data.data_matrix is None: return
        rows = sorted(list(set(r for r, c in self.selected_cells)))
        cols = sorted(list(set(c for r, c in self.selected_cells)))
        if not rows or not cols: return
        start_col = min(cols)
        end_col = max(cols)
        fasta_output = []
        for r in rows:
            seq_id = self.data.ids[r]
            row_bytes = self.data.data_matrix[r, start_col:end_col+1]
            seq_str = row_bytes.tobytes().decode('ascii', errors='replace')
            fasta_output.append(f">{seq_id}")
            fasta_output.append(seq_str)
        text_data = "\n".join(fasta_output)
        QApplication.clipboard().setText(text_data)

    def get_selected_column_range(self):
        if not self.selected_cells: return None
        cols = [c for r, c in self.selected_cells]
        return min(cols), max(cols)

    def get_selected_rows(self):
        if not self.selected_cells: return []
        return sorted(list(set(r for r, c in self.selected_cells)))
    
    def select_rows(self, row_indices):
        if self.data.data_matrix is None: return
        rows, cols = self.data.data_matrix.shape
        self.selected_cells = set()
        for r in row_indices:
            if 0 <= r < rows:
                for c in range(cols): self.selected_cells.add((r, c))
        self.viewport().update()
        
    def scroll_to_cell(self, row, col):
        x = col * self.char_width
        y = row * self.char_height
        self.horizontalScrollBar().setValue(x - self.viewport().width() // 2)
        self.verticalScrollBar().setValue(y - self.viewport().height() // 2)

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle(f"{APP_NAME} - {APP_TAGLINE}")
        self.resize(1200, 800)
        self.data = AlignmentData()
        
        # STATUS BAR for selection info
        self.status_label = QLabel("Ready")
        self.statusBar().addWidget(self.status_label)
        
        self.toolbar = QToolBar("Main Toolbar")
        self.addToolBar(self.toolbar)
        self.toolbar.setMovable(False)
        spacer = QWidget()
        spacer.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        self.toolbar.addWidget(spacer)
        
        search_container = QWidget()
        search_layout = QHBoxLayout(search_container)
        search_layout.setContentsMargins(0, 0, 0, 0)
        self.search_input = QLineEdit()
        self.search_input.setPlaceholderText("Search sequence/name...")
        self.search_input.setFixedWidth(200)
        self.search_input.returnPressed.connect(self.find_next)
        find_btn = QPushButton("Find")
        find_btn.clicked.connect(self.find_next)
        search_layout.addWidget(QLabel("Find:"))
        search_layout.addWidget(self.search_input)
        search_layout.addWidget(find_btn)
        self.toolbar.addWidget(search_container)

        self.viewer = SequenceViewer(self.data)
        self.name_list = DraggableListWidget(self.data, self.viewer)
        self.ruler = RulerWidget(self.viewer)
        
        # --- SCROLLBAR UNIFICATION ---
        self.viewer.verticalScrollBar().valueChanged.connect(self.name_list.verticalScrollBar().setValue)
        self.name_list.verticalScrollBar().valueChanged.connect(self.viewer.verticalScrollBar().setValue)
        self.viewer.horizontalScrollBar().valueChanged.connect(self.ruler.update)
        self.viewer.row_height_changed.connect(lambda s, h: self.ruler.update())
        self.name_list.reorder_signal.connect(self.viewer.viewport().update)
        self.viewer.row_height_changed.connect(self.name_list.update_font_size)
        self.viewer.grid_clicked.connect(self.name_list.set_selected_row)
        
        # Update Status Bar on clicks
        self.viewer.grid_clicked.connect(self.update_status_from_grid)
        self.name_list.itemClicked.connect(self.update_status_from_list)

        splitter = QSplitter(Qt.Orientation.Horizontal)
        left_widget = QWidget()
        left_layout = QVBoxLayout(left_widget)
        left_layout.setContentsMargins(0,0,0,0)
        left_layout.setSpacing(0)
        list_header = QLabel(" Sequences")
        list_header.setFixedHeight(30)
        list_header.setStyleSheet("background-color: #DDD; border-bottom: 1px solid #CCC; border-right: 1px solid #999; font-weight: bold;")
        left_layout.addWidget(list_header)
        left_layout.addWidget(self.name_list)
        right_widget = QWidget()
        right_layout = QVBoxLayout(right_widget)
        right_layout.setContentsMargins(0,0,0,0)
        right_layout.setSpacing(0)
        right_layout.addWidget(self.ruler)
        right_layout.addWidget(self.viewer)
        splitter.addWidget(left_widget)
        splitter.addWidget(right_widget)
        splitter.setStretchFactor(1, 4)
        
        central_widget = QWidget()
        layout = QVBoxLayout(central_widget)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        layout.addWidget(splitter)
        
        self.setCentralWidget(central_widget)
        self.create_menu()
        
        self.primer_design_row = None
        self.primer_design_col_start = 0
        self.primer_design_raw_gapped = ""
        self.last_search_idx = -1
    
    def update_status_from_grid(self, row):
        if 0 <= row < len(self.data.ids):
            name = self.data.ids[row]
            self.status_label.setText(f"Selected: {name} (Row {row+1})")
            
    def update_status_from_list(self, item):
        self.status_label.setText(f"Selected: {item.text()}")

    def create_menu(self):
        menu = self.menuBar()
        file_menu = menu.addMenu("&File")
        open_action = QAction("&Open Sequence File...", self)
        open_action.setShortcut(QKeySequence.StandardKey.Open)
        open_action.triggered.connect(self.open_file)
        file_menu.addAction(open_action)
        save_action = QAction("&Save Alignment As...", self)
        save_action.setShortcut("Ctrl+S")
        save_action.triggered.connect(self.save_alignment)
        file_menu.addAction(save_action)
        exit_action = QAction("E&xit", self)
        exit_action.triggered.connect(self.close)
        file_menu.addAction(exit_action)
        
        edit_menu = menu.addMenu("&Edit")
        undo_action = QAction("&Undo", self)
        undo_action.setShortcut("Ctrl+Z")
        undo_action.triggered.connect(self.undo_last_action)
        edit_menu.addAction(undo_action)
        copy_action = QAction("&Copy Selection (FASTA)", self)
        copy_action.setShortcut("Ctrl+C")
        copy_action.triggered.connect(lambda: self.viewer.copy_selection_to_clipboard())
        edit_menu.addAction(copy_action)
        edit_menu.addSeparator()
        del_cols_action = QAction("Delete Selected &Columns (Shorten)", self)
        del_cols_action.setShortcut("Ctrl+Delete")
        del_cols_action.triggered.connect(self.delete_selected_columns)
        edit_menu.addAction(del_cols_action)
        del_rows_action = QAction("Delete Selected &Sequences", self)
        del_rows_action.triggered.connect(self.delete_selected_rows)
        edit_menu.addAction(del_rows_action)
        fill_gap_action = QAction("&Fill Selection with Gaps (-)", self)
        fill_gap_action.setShortcut("Delete")
        fill_gap_action.triggered.connect(self.fill_selection_gaps)
        edit_menu.addAction(fill_gap_action)
        
        view_menu = menu.addMenu("&View")
        zoom_in = QAction("Zoom &In", self)
        zoom_in.setShortcut(QKeySequence.StandardKey.ZoomIn)
        zoom_in.triggered.connect(lambda: self.viewer.set_zoom(self.viewer.char_width/10.0 + 0.2))
        view_menu.addAction(zoom_in)
        zoom_out = QAction("Zoom &Out", self)
        zoom_out.setShortcut(QKeySequence.StandardKey.ZoomOut)
        zoom_out.triggered.connect(lambda: self.viewer.set_zoom(self.viewer.char_width/10.0 - 0.2))
        view_menu.addAction(zoom_out)
        toggle_text = QAction("Toggle &Text Visibility", self)
        toggle_text.triggered.connect(self.toggle_text_mode)
        view_menu.addAction(toggle_text)

        settings_menu = menu.addMenu("&Settings")
        font_action = QAction("Increase Font Size", self)
        font_action.triggered.connect(lambda: self.viewer.set_global_font_size(self.viewer.font_size + 2))
        settings_menu.addAction(font_action)
        font_down = QAction("Decrease Font Size", self)
        font_down.triggered.connect(lambda: self.viewer.set_global_font_size(max(6, self.viewer.font_size - 2)))
        settings_menu.addAction(font_down)

        colors_menu = menu.addMenu("&Colors")
        c_blocks = QAction("&Color Blocks (Default)", self)
        c_blocks.triggered.connect(lambda: self.viewer.set_color_mode(0))
        colors_menu.addAction(c_blocks)
        c_text = QAction("Color &Text (White BG)", self)
        c_text.triggered.connect(lambda: self.viewer.set_color_mode(1))
        colors_menu.addAction(c_text)
        c_bw = QAction("&Black & White", self)
        c_bw.triggered.connect(lambda: self.viewer.set_color_mode(2))
        colors_menu.addAction(c_bw)

        align_menu = menu.addMenu("&Align")
        muscle_action = QAction("Realign with &MUSCLE", self)
        muscle_action.triggered.connect(self.run_muscle)
        align_menu.addAction(muscle_action)
        
        phylo_menu = menu.addMenu("&Phylogeny")
        fasttree_action = QAction("Build Tree with &FastTree (Fast)", self)
        fasttree_action.triggered.connect(self.run_fasttree)
        phylo_menu.addAction(fasttree_action)
        nj_action = QAction("Build Tree with &Neighbor-Joining (Slow/Built-in)", self)
        nj_action.triggered.connect(self.run_nj_tree)
        phylo_menu.addAction(nj_action)
        
        tools_menu = menu.addMenu("&Tools")
        primer_action = QAction("&Design Primers (from selection)", self)
        primer_action.triggered.connect(self.design_primers)
        tools_menu.addAction(primer_action)

    def open_file(self):
        fname, _ = QFileDialog.getOpenFileName(self, "Open Sequence File", "", "Sequence Files (*.fasta *.fa *.fas *.txt);;All Files (*)")
        if fname:
            try:
                self.data.load_fasta(fname)
                self.refresh_views()
                self.setWindowTitle(f"{APP_NAME} - {os.path.basename(fname)}")
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to load file:\n{str(e)}")
    
    def save_alignment(self):
        fname, _ = QFileDialog.getSaveFileName(self, "Save Alignment As", "", "FASTA Files (*.fasta);;All Files (*)")
        if fname:
            try:
                self.data.save_fasta(fname)
                QMessageBox.information(self, "Success", f"Saved alignment to:\n{fname}")
            except Exception as e:
                QMessageBox.critical(self, "Save Error", str(e))

    def refresh_views(self):
        self.viewer.selected_cells = set()
        self.viewer.highlight_f_cells = set()
        self.viewer.highlight_r_cells = set()
        self.viewer.update_scrollbars()
        self.viewer.viewport().update()
        self.ruler.update() 
        self.name_list.populate() 

    def toggle_text_mode(self):
        self.viewer.draw_text = not self.viewer.draw_text
        self.viewer.viewport().update()
    
    def undo_last_action(self):
        if self.data.undo(): self.refresh_views()
        else: QMessageBox.information(self, "Undo", "Nothing left to undo!")

    def delete_selected_columns(self):
        col_range = self.viewer.get_selected_column_range()
        if not col_range: return
        start, end = col_range
        reply = QMessageBox.question(self, "Confirm Delete", 
                                     f"Are you sure you want to delete columns {start} to {end}?",
                                     QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No)
        if reply == QMessageBox.StandardButton.Yes:
            self.data.delete_columns(start, end)
            self.refresh_views()

    def delete_selected_rows(self):
        rows = self.viewer.get_selected_rows()
        if not rows: return
        reply = QMessageBox.question(self, "Confirm Delete", 
                                     f"Are you sure you want to delete {len(rows)} sequence(s)?",
                                     QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No)
        if reply == QMessageBox.StandardButton.Yes:
            self.data.delete_rows(rows)
            self.refresh_views()

    def fill_selection_gaps(self):
        if not self.viewer.selected_cells: return
        for r, c in self.viewer.selected_cells:
            self.data.update_base(r, c, "-")
        self.viewer.viewport().update()

    def run_muscle(self):
        if self.data.data_matrix is None: return
        self.data.save_state()
        if getattr(sys, 'frozen', False): base_dir = os.path.dirname(sys.executable)
        else: base_dir = os.path.dirname(os.path.abspath(__file__))
        muscle_exe = None
        possible_local = os.path.join(base_dir, "muscle.exe")
        if os.path.exists(possible_local): muscle_exe = possible_local
        else: muscle_exe = shutil.which("muscle")
        if not muscle_exe:
            QMessageBox.warning(self, "Missing Tool", f"MUSCLE executable not found.\nPlace 'muscle.exe' in: {base_dir}")
            return
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fa') as input_tmp:
            input_path = input_tmp.name
            bio_obj = self.data.to_biopython_obj()
            AlignIO.write(bio_obj, input_tmp, "fasta")
        output_path = input_path + ".aligned"
        try:
            cmd = [muscle_exe, "-in", input_path, "-out", output_path]
            startupinfo = None
            if os.name == 'nt':
                startupinfo = subprocess.STARTUPINFO()
                startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW
            result = subprocess.run(cmd, capture_output=True, text=True, startupinfo=startupinfo)
            if result.returncode != 0: raise Exception(f"Muscle Error: {result.stderr}")
            self.data.load_fasta(output_path)
            self.refresh_views()
            QMessageBox.information(self, "Success", "Alignment completed!")
        except Exception as e: QMessageBox.critical(self, "Alignment Failed", str(e))
        finally:
            if os.path.exists(input_path): os.remove(input_path)
            if os.path.exists(output_path): os.remove(output_path)
            
    def run_fasttree(self):
        if self.data.data_matrix is None: return
        if getattr(sys, 'frozen', False): base_dir = os.path.dirname(sys.executable)
        else: base_dir = os.path.dirname(os.path.abspath(__file__))
        ft_exe = None
        possible = os.path.join(base_dir, "FastTree.exe")
        if os.path.exists(possible): ft_exe = possible
        else: ft_exe = shutil.which("FastTree")
        if not ft_exe:
            QMessageBox.warning(self, "Missing Tool", f"FastTree executable not found.\nPlace 'FastTree.exe' in: {base_dir}")
            return
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fa') as input_tmp:
            input_path = input_tmp.name
            bio_obj = self.data.to_biopython_obj()
            AlignIO.write(bio_obj, input_tmp, "fasta")
        try:
            cmd = [ft_exe, "-nt", input_path]
            startupinfo = None
            if os.name == 'nt':
                startupinfo = subprocess.STARTUPINFO()
                startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW
            result = subprocess.run(cmd, capture_output=True, text=True, startupinfo=startupinfo)
            if result.returncode != 0: raise Exception(result.stderr)
            tree_str = result.stdout
            dlg = TreeResultDialog(tree_str, self)
            dlg.exec()
        except Exception as e: QMessageBox.critical(self, "FastTree Error", str(e))
        finally:
            if os.path.exists(input_path): os.remove(input_path)

    def run_nj_tree(self):
        if self.data.data_matrix is None: return
        rows, cols = self.data.data_matrix.shape
        if rows > 100:
            reply = QMessageBox.question(self, "Large Alignment", 
                "Calculating Neighbor-Joining for >100 sequences in Python is very slow.\nContinue?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No)
            if reply == QMessageBox.StandardButton.No: return
        try:
            aln = self.data.to_biopython_obj()
            calc = DistanceCalculator('identity')
            dm = calc.get_distance(aln)
            constructor = DistanceTreeConstructor()
            tree = constructor.nj(dm)
            from io import StringIO
            handle = StringIO()
            Phylo.write(tree, handle, "newick")
            tree_str = handle.getvalue()
            dlg = TreeResultDialog(tree_str, self)
            dlg.exec()
        except Exception as e: QMessageBox.critical(self, "Tree Error", str(e))
    
    def find_next(self):
        if self.data.data_matrix is None: return
        query = self.search_input.text().upper()
        if not query: return
        rows, cols = self.data.data_matrix.shape
        start_search_idx = self.last_search_idx + 1
        found = False
        for i, seq_id in enumerate(self.data.ids):
            if i <= self.last_search_idx and self.last_search_idx < rows: continue 
            if query in seq_id.upper():
                self.name_list.set_selected_row(i)
                self.last_search_idx = i
                return
        query_bytes = query.encode('ascii', errors='replace')
        for r in range(rows):
            row_bytes = self.data.data_matrix[r].tobytes()
            match_col = row_bytes.find(query_bytes)
            if match_col != -1:
                current_linear = r * cols + match_col
                if current_linear > self.last_search_idx:
                    self.viewer.selected_cells = set()
                    for k in range(len(query)): self.viewer.selected_cells.add((r, match_col + k))
                    self.viewer.scroll_to_cell(r, match_col)
                    self.viewer.viewport().update()
                    self.name_list.set_selected_row(r)
                    self.last_search_idx = current_linear
                    found = True
                    return
        if not found:
            self.last_search_idx = -1 
            QMessageBox.information(self, "Find", f"'{query}' not found (or wrapped to start).")

    def design_primers(self):
        if self.data.data_matrix is None: return
        rows = self.viewer.get_selected_rows()
        col_range = self.viewer.get_selected_column_range()
        if not rows or not col_range:
            QMessageBox.warning(self, "Selection Required", "Please select a region to analyze.")
            return
        param_dialog = PrimerParamsDialog(self)
        if param_dialog.exec() != QDialog.DialogCode.Accepted: return
        min_amp, max_amp = param_dialog.get_values()
        start, end = col_range
        row_idx = rows[0]
        if len(rows) > 1:
            consensus_seq = []
            for c in range(start, end + 1):
                col_bytes = self.data.data_matrix[rows, c] 
                unique, counts = np.unique(col_bytes, return_counts=True)
                best_char = unique[np.argmax(counts)]
                consensus_seq.append(best_char.decode('ascii'))
            self.primer_design_raw_gapped = "".join(consensus_seq)
            self.primer_design_row = rows[0] 
        else:
            self.primer_design_row = row_idx
            raw_bytes = self.data.data_matrix[row_idx, start:end+1]
            self.primer_design_raw_gapped = raw_bytes.tobytes().decode('ascii', errors='replace')
        self.primer_design_col_start = start
        seq_fragment = self.primer_design_raw_gapped.replace('-', '')
        if len(seq_fragment) < 20:
             QMessageBox.warning(self, "Too Short", "Selected region is too short.")
             return
        designer = PrimerDesigner(seq_fragment)
        f_cands, r_cands, pairs = designer.find_primers(min_amp, max_amp)
        if not f_cands and not r_cands:
            QMessageBox.information(self, "No Primers", "No valid primers found.")
            return
        dialog = PrimerResultsDialog(f_cands, r_cands, pairs, self)
        dialog.exec()
    
    def map_coords(self, gapless_start, gapless_length):
        if self.primer_design_row is None: return None
        visual_start = -1
        visual_end = -1
        current_gapless_idx = 0
        for i, char in enumerate(self.primer_design_raw_gapped):
            if char != '-':
                if current_gapless_idx == gapless_start:
                    visual_start = i
                    break
                current_gapless_idx += 1
        count = 0
        for i in range(visual_start, len(self.primer_design_raw_gapped)):
            if self.primer_design_raw_gapped[i] != '-': count += 1
            if count == gapless_length:
                visual_end = i + 1 
                break
        return visual_start + self.primer_design_col_start, visual_end + self.primer_design_col_start

    def highlight_primer_region(self, start, length, p_type):
        res = self.map_coords(start, length)
        if not res: return
        abs_start, abs_end = res
        self.viewer.selected_cells = set() 
        self.viewer.highlight_f_cells = set()
        self.viewer.highlight_r_cells = set()
        target_set = self.viewer.highlight_f_cells if p_type == "Forward" else self.viewer.highlight_r_cells
        for c in range(abs_start, abs_end): target_set.add((self.primer_design_row, c))
        self.viewer.scroll_to_cell(self.primer_design_row, abs_start)
        self.viewer.viewport().update()

    def highlight_primer_pair(self, f, r):
        f_res = self.map_coords(f['start'], f['length'])
        r_res = self.map_coords(r['sense_start'], r['length'])
        if not f_res or not r_res: return
        self.viewer.selected_cells = set()
        self.viewer.highlight_f_cells = set()
        self.viewer.highlight_r_cells = set()
        for c in range(f_res[0], f_res[1]): self.viewer.highlight_f_cells.add((self.primer_design_row, c))
        for c in range(r_res[0], r_res[1]): self.viewer.highlight_r_cells.add((self.primer_design_row, c))
        self.viewer.scroll_to_cell(self.primer_design_row, f_res[0])
        self.viewer.viewport().update()

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())