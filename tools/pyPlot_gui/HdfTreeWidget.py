# hdftreewidget.py ---
#
# Filename: hdftreewidget.py
# Description:
# Author: subha
# Maintainer:
# Created: Fri Jul 24 20:54:11 2015 (-0400)
# Version:
# Last-Updated: Fri Sep 18 21:08:59 2015 (-0400)
#           By: subha
#     Update #: 539
# URL:
# Keywords:
# Compatibility:
#
#

# Commentary:
#
#
#
#

# Change log:
#
#
#
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation; either version 3, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; see the file COPYING.  If not, write to
# the Free Software Foundation, Inc., 51 Franklin Street, Fifth
# Floor, Boston, MA 02110-1301, USA.
#
#

# Code:
"""Defines class to display HDF5 file tree.

"""

import numpy as np
import h5py as h5

from PyQt5 import QtCore
from PyQt5.QtWidgets import QApplication, QMainWindow, QWidget, QTreeView

from HdfTreeModel import *

class HDFTreeWidget(QTreeView):
    """Convenience class to display HDF file trees.

    HDFTreeWidget wraps an HDFTreeModel and displays it.
    In addition it allows opening multiple files from a list.

    Signals
    -------

    sigDatasetWidgetCreated(QWidget): emitted by createDatasetWidget slot
    after a new dataset widget is created with the created
    widget. This provides a way to send it back to the DataViz widget
    (top level widget and caller of the slot) for incorporation into
    the main window.

    sigDatasetWidgetClosed(QWidget): this is emitted for each of the
    widgets showing datasets under a given file tree when the file is
    closed. This allows the toplevel widget to close the corresponding
    mdi child window.

    sigAttributeWidgetCreated(QWidget): same as sigDatasetWidgetCreated but
    for the widget displaying HDF5 attributes.

    sigAttributeWidgetClosed(QWidget): same as sigDatasetWidgetClosed but
    for the widget displaying HDF5 attributes.

    sigPlotWidgetCreated(QWidget): same as sigDatasetWidgetCreated but
    for the widget displaying HDF5 dataset plots.

    sigPlotWidgetClosed(QWidget): same as sigDatasetWidgetClosed but
    for the widget displaying HDF5 dataset plots.

    """

    def __init__(self, parent=None):
        super(HDFTreeWidget, self).__init__(parent)
        model = HDFTreeModel([])
        self.setModel(model)


    def openFiles(self, files, mode='r+'):
        """Open the files listed in argument.

        files: list of file paths. For example, output of
               QFileDialog::getOpenFileNames

        mode: string specifying open mode for h5py

        """
        for fname in files:
            self.model().openFile(fname, mode)

    def closeFiles(self):
        """Close the files selected in the model.

        If there are datasets in the file that are being displayed via
        a `HDFDatasetWidget`, then it emits a sigDatasetWidgetClosed
        signal with the HDFDatasetWidget as a parameter for each of
        them. Same for `HDFAttributeWidget`s.

        """
        indices = self.selectedIndexes()
        for index in indices:
            item = self.model().getItem(index)
            filename = item.h5node.file.filename
            if self.model().closeFile(index):
                 print( "one file is closed!" )


if __name__ == '__main__':
    import sys
    app = QApplication(sys.argv)
    window = QMainWindow()
    widget = HDFTreeWidget()
    window.setCentralWidget(widget)
    widget.openFiles(['test.h5'], 'r+')
    window.show()
    sys.exit(app.exec_())

#
# hdftreewidget.py ends here
