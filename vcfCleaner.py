#!/usr/bin/env python
import sys, os, webbrowser
from PySide.QtCore import QFile
from PySide.QtUiTools import QUiLoader
from PySide.QtGui import QApplication

REPR_CHAR_LIMIT = 30

class attribute:
    def __init__(self, name, sourcePath, additionalText=""):
        self.name = name
        self.sourcePath = sourcePath
        self.sourceFile = os.path.split(sourcePath)[1]
        self.additionalText = additionalText
        self.globalName = name
    
    def __repr__(self):
        outstr = "%s (%s" % (self.globalName,self.sourceFile)
        if self.additionalText != "":
            outstr += ")"
        else:
            outstr += ", %s)" % self.additionalText
        return outstr

class statistic:
    def __init__(self, function, population):
        self.function = function
        self.population = population
        self.globalName = self.population + "_" + self.function
    
    def __repr__(self):
        return self.globalName

class filterExpression:
    def __init__(self, expression):
        self.expression = expression
    
    def __repr__(self):
        temp = self.expression
        if len(temp) > REPR_CHAR_LIMIT-15:
            temp = temp[:REPR_CHAR_LIMIT-15] + "..."
        return "EXPRESSION: %s" % temp

class filterBed:
    def __init__(self, path):
        self.path = path
    
    def __repr__(self):
        temp = self.path
        if len(temp) > REPR_CHAR_LIMIT-8:
            temp = "..." + temp[-REPR_CHAR_LIMIT-8:]
        return "BED: %s" % temp

class gui:
    def __init__(self):
        self.loader = QUiLoader()
        infile = QFile("ui/Main.ui")
        infile.open(QFile.ReadOnly)
        self.window = self.loader.load(infile, None)
        self.childWindow = None
        infile.close()
        
        # data structures
        self.populationOrder = []
        self.populations = {}
        self.loadPopulations('KGP_Populations.txt')
        
        self.variantFilters = []
        
        self.includeAttributes = [] # contains actual attribute objects
        self.removeAttributes = [] # contains actual attribute objects
        
        self.activeCsvFile = None
        self.csvAttributes = [] # really just contains strings
        self.activeBedFile = None
        self.bedAttributes = [] # really just contains strings
        
        self.calculatedAttributes = []
        
        # set up GUI
        self.window.browseVcfButton.clicked.connect(self.browseVcf)
        self.window.browseOutputButton.clicked.connect(self.browseOutput)
        
        self.window.populationComboBox.currentIndexChanged.connect(self.changePopulation)
        
        self.window.filterList.itemSelectionChanged.connect(self.changeFilter)
        self.window.removeFilterButton.clicked.connect(self.removeFilter)
        self.window.browseBedFilterButton.clicked.connect(self.addBedFilter)
        self.window.addExpressionFilterButton.clicked.connect(self.addExpressionFilter)
        
        self.window.includeAttributeList.itemSelectionChanged.connect(self.changeIncludeAttribute)
        #self.window.removeAttributeButton.clicked.connect(self.removeAttributes)
        
        #self.window.removeInfoList.itemSelectionChanged.connect(self.changeRemoveInfo)
        self.window.addInfoButton.clicked.connect(self.addInfo)
        
        #self.window.csvColumnList.itemSelectionChanged.connect(self.changeCsvColumn)
        self.window.addCsvColumnButton.clicked.connect(self.addCsvAttribute)
        
        #self.window.bedScoreList.itemSelectionChanged.connect(self.changeBedScore)
        self.window.addBedScoreButton.clicked.connect(self.addBedAttribute)
        
        self.window.calculateAdditionalButton.clicked.connect(self.addAdditionalAttribute)
        
        self.window.helpButton.clicked.connect(self.openHelpPage)
        self.window.quitButton.clicked.connect(self.quit)
        self.window.runButton.clicked.connect(self.run)
        
        # Initially disable/hide lots of crap
        self.window.vcfAlleleSettings.hide()
        self.window.csvAlleleSettings.disable()
        
        self.window.populationComboBox.disable()
        self.window.removeNonBiallelicCheckBox.disable()
        
        self.window.removeFilterButton.disable()
        self.window.addExpressionFilterButton.disable()
        self.window.filterSettingsGroupBox.disable()
        
        self.window.removeAttributeButton.disable()
        self.window.addCsvAttributeButton.disable()
        self.window.addBedAttributeButton.disable()
        self.window.calculateAdditionalButton.disable()
        self.window.attributeSettingsGroupBox.disable()
        
        self.window.runButton.disable()
        
        self.window.show()
    
    def loadPopulations(self, path):
        pass
    
    # GUI Event handlers
    def browseVcf(self):
        pass
    
    def browseOutput(self):
        pass
    
    def changePopulation(self):
        pass
    
    def changeFilter(self):
        pass
    
    def removeFilter(self):
        pass
    
    def addBedFilter(self):
        pass
    
    def addExpressionFilter(self):
        pass
    
    def changeIncludeAttribute(self):
        pass
    
    def removeAttributes(self):
        pass
    
    def changeRemoveInfo(self):
        pass
    
    def addInfo(self):
        pass
    
    def changeCsvColumn(self):
        pass
    
    def addCsvAttribute(self):
        pass
    
    def changeBedScore(self):
        pass
    
    def addBedAttribute(self):
        pass
    
    def addAdditionalAttribute(self):
        infile = QFile("ui/Calculate.ui")
        infile.open(QFile.ReadOnly)
        self.childWindow = self.loader.load(infile, None)
        infile.close()
        
        # TODO: populate with existing calculated attributes
        
        def addFields():
            self.childWindow.accept()
            self.childWindow = None
        
        def browsePopulation():
            pass
        
        def addStat():
            pass
        
        def removeStat():
            pass
        
        def calculateHelp():
            pass
        
        self.childWindow.calculateOkButton.clicked.connect(addFields)
        self.childWindow.browsePopulationButton.clicked.connect(browsePopulation)
        self.childWindow.addStatButton.clicked.connect(addStat)
        self.childWindow.removeStatButton.clicked.connect(removeStat)
        self.childWindow.calculateHelpButton.clicked.connect(calculateHelp)
        
        self.childWindow.show()
    
    def openHelpPage(self):
        webbrowser.open('http://sci.utah.edu/~abigelow/vcfCleanerHelp.php',new=2)
    
    def quit(self):
        self.window.reject()
    
    def run(self):
        self.window.accept()

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = gui()
    sys.exit(app.exec_())