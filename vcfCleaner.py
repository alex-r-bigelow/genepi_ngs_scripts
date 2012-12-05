#!/usr/bin/env python
import sys
from PySide.QtCore import QFile
from PySide.QtUiTools import QUiLoader
from PySide.QtGui import QApplication

class gui:
    def __init__(self):
        self.loader = QUiLoader()
        infile = QFile("ui/Main.ui")
        infile.open(QFile.ReadOnly)
        self.window = self.loader.load(infile, None)
        infile.close()
        
        self.window.quitButton.clicked.connect(self.quit)
        self.window.browseVcfButton.clicked.connect(self.browseVcf)
        self.window.browseOutputButton.clicked.connect(self.browseOutput)
        self.window.calculateFieldsButton.clicked.connect(self.calculateFields)
        self.window.removeFieldButton.clicked.connect(self.removeField)
        self.window.addRemovedFieldButton.clicked.connect(self.addRemovedField)
        self.window.browseCsvButton.clicked.connect(self.browseCsv)
        self.window.addCsvFieldButton.clicked.connect(self.addCsvField)
        
        self.window.show()
    
    def quit(self):
        self.window.reject()
    
    def browseVcf(self):
        pass
    
    def browseOutput(self):
        pass
    
    def calculateFields(self):
        infile = QFile("ui/Calculate.ui")
        infile.open(QFile.ReadOnly)
        calcWindow = self.loader.load(infile, None)
        infile.close()
        
        calcWindow.show()
    
    def removeField(self):
        pass
    
    def addRemovedField(self):
        pass
    
    def browseCsv(self):
        pass
    
    def addCsvField(self):
        pass


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = gui()
    sys.exit(app.exec_())