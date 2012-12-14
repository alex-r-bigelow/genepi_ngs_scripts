#!/usr/bin/env python
import sys, os, webbrowser
from PySide.QtCore import Qt, QFile
from PySide.QtUiTools import QUiLoader
from PySide.QtGui import QApplication, QFileDialog, QProgressDialog, QMessageBox, QComboBox, QTableWidgetItem
from genome_utils import genomeException, parsePopulations, MAX_INFO_STRINGS
from cleanVCF import extractInfoFields
from addCSVtoVCF import sniffCsv
from addBEDtoVCF import sniffBed
from calcStats import allStats

RESERVED_TAGS = ["CHROM","POS","ID"]

class attribute:
    def __init__(self, name, sourcePath, additionalText=""):
        self.name = name
        self.sourcePath = sourcePath
        self.sourceFile = os.path.split(sourcePath)[1]
        self.additionalText = additionalText
        self.globalName = name
    
    def __repr__(self):
        outstr = "%s (%s" % (self.globalName,self.sourceFile)
        if self.additionalText == "":
            outstr += ")"
        else:
            outstr += ", %s)" % self.additionalText
        return outstr
    
    def __eq__(self, other):
        return isinstance(other,attribute) and self.name == other.name and self.sourcePath == other.sourcePath and self.additionalText == other.additionalText
    
    def __ne__(self, other):
        return not self.__eq__(other)

class statistic:
    def __init__(self, function, population,note=""):
        self.function = function
        self.population = population
        self.globalName = self.population + "_" + self.function
        self.note=note
    
    def __repr__(self):
        return self.globalName
    
    def __eq__(self, other):
        return isinstance(other,statistic) and self.function == other.function and self.population == other.population and self.note == other.note
    
    def __ne__(self, other):
        return not self.__eq__(other)

class filterExpression:
    def __init__(self, expression, columns):
        self.expression = expression
        self.columns = columns
    
    def __repr__(self):
        temp = self.expression % tuple(self.columns)
        return "EXPRESSION: %s" % temp
    
    def __eq__(self, other):
        return isinstance(other,filterExpression) and self.expression == other.expression and self.columns == other.columns
    
    def __ne__(self, other):
        return not self.__eq__(other)

class filterBed:
    def __init__(self, path):
        self.path = path
    
    def __repr__(self):
        temp = os.path.split(self.path)[1]
        return "BED: %s" % temp
    
    def __eq__(self, other):
        return isinstance(other,filterBed) and self.path == other.path
    
    def __ne__(self, other):
        return not self.__eq__(other)

class fieldBit:
    def __init__(self, field):
        self.field = field

class gui:
    def __init__(self):
        self.loader = QUiLoader()
        infile = QFile("ui/Main.ui")
        infile.open(QFile.ReadOnly)
        self.window = self.loader.load(infile, None)
        infile.close()
        
        self.childWindow = None
        self.childWidgets = []
        
        # data structures
        self.populationFiles = set()
        self.populationOrder = []
        self.populations = {}
        self.lastBackgroundPopIndex = 0
        self.loadPopulations(sys.path[0] + "/KGP_populations.txt")
        
        self.variantFilters = []
        
        self.infoFields = {}
        self.hasValidInput = False
        self.hasValidOutput = False
        
        self.includedAttributes = [] # contains actual attribute objects
        self.removedAttributes = [] # contains actual attribute objects
        
        self.activeCsvFile = None
        self.csvAttributes = [] # really just contains strings
        self.activeBedFile = None
        self.bedAttributes = [] # really just contains strings
        
        self.calculatedAttributes = []
        
        self.window.allAttributeComboBox.addItem("")
        self.window.allAttributeComboBox.addItems(RESERVED_TAGS)
        
        # set up GUI
        self.window.browseVcfButton.clicked.connect(self.browseVcf)
        self.window.browseOutputButton.clicked.connect(self.browseOutput)
        
        self.window.populationComboBox.currentIndexChanged.connect(self.changePopulation)
        
        self.window.includeAttributeList.itemSelectionChanged.connect(self.changeIncludeAttribute)
        self.window.removeAttributeButton.clicked.connect(self.removeAttributes)
        
        self.window.removeInfoList.itemSelectionChanged.connect(self.changeRemoveInfo)
        self.window.addInfoButton.clicked.connect(self.addInfo)
        
        self.window.browseCsvButton.clicked.connect(self.browseCsv)
        self.window.csvColumnList.itemSelectionChanged.connect(self.changeCsvColumn)
        self.window.addCsvColumnButton.clicked.connect(self.addCsvAttribute)
        
        self.window.browseBedButton.clicked.connect(self.browseBed)
        self.window.bedScoreList.itemSelectionChanged.connect(self.changeBedScore)
        self.window.addBedScoreButton.clicked.connect(self.addBedAttribute)
        
        self.window.calculateAdditionalButton.clicked.connect(self.addAdditionalAttribute)
        
        self.window.filterList.itemSelectionChanged.connect(self.changeFilter)
        self.window.removeFilterButton.clicked.connect(self.removeFilter)
        self.window.browseBedFilterButton.clicked.connect(self.addBedFilter)
        self.window.allAttributeComboBox.currentIndexChanged.connect(self.insertAttribute)
        self.window.addExpressionFilterButton.clicked.connect(self.addExpressionFilter)
        
        self.window.helpButton.clicked.connect(self.openHelpPage)
        self.window.quitButton.clicked.connect(self.quit)
        self.window.runButton.clicked.connect(self.run)
        
        self.globalEnable()
        
        self.window.show()
    
    def loadPopulations(self, path):
        if path in self.populationFiles:
            return
        self.populationFiles.add(path)
        
        self.populationOrder.append(None)#(os.path.split(path)[1]))
        newPops,popOrder = parsePopulations(path)
        for k in popOrder:
            p = newPops[k]
            temp = k
            appendDigit = 2
            while self.populations.has_key(temp):
                temp = "%s_#%i" % (k, appendDigit)
                appendDigit += 1
            self.populations[temp] = (p,path)
            self.populationOrder.append(temp)
        
        self.lastBackgroundPopIndex = max(0,self.window.populationComboBox.currentIndex())
        
        self.window.populationComboBox.clear()
        for i,p in enumerate(self.populationOrder):
            if p == None:
                if not i == 0:
                    #self.window.populationComboBox.addItem(p[0])
                    self.window.populationComboBox.insertSeparator(i)   # TODO: put the file name in instead
            else:
                self.window.populationComboBox.addItem(p)
        
        self.window.populationComboBox.insertSeparator(self.window.populationComboBox.count()+1)
        self.window.populationComboBox.addItem('Load population .txt file...')
        self.window.populationComboBox.setCurrentIndex(self.lastBackgroundPopIndex)
        
        for w in self.childWidgets:
            w.clear()
            for i,p in enumerate(self.populationOrder):
                if p == None:
                    if not i == 0:
                        #self.window.populationComboBox.addItem(p[0])
                        self.window.populationComboBox.insertSeparator(i)   # TODO: put the file name in instead
                else:
                    self.window.populationComboBox.addItem(p)
    
    def globalEnable(self):
        if self.hasValidInput:
            self.window.tabWidget.widget(1).setEnabled(True)
        else:
            self.window.tabWidget.widget(1).setEnabled(False)
        
        if self.hasValidOutput:
            if self.window.outputPathField.text().endswith("vcf"):
                self.window.vcfAlleleSettings.show()
                self.window.csvAlleleSettings.hide()
            else:
                self.window.vcfAlleleSettings.hide()
                self.window.csvAlleleSettings.show()
            
            self.window.generalSettingsPanel.setEnabled(True)
            self.window.tabWidget.widget(2).setEnabled(True)
        else:
            self.window.vcfAlleleSettings.hide()
            self.window.csvAlleleSettings.show()
            
            self.window.generalSettingsPanel.setEnabled(False)
            self.window.tabWidget.widget(2).setEnabled(False)
        
        if self.hasValidInput and self.hasValidOutput:
            self.window.runButton.setEnabled(True)
        else:
            self.window.runButton.setEnabled(False)
    
    def updateLists(self):
        self.window.filterList.clear()
        self.window.includeAttributeList.clear()
        self.changeIncludeAttribute()
        self.window.removeInfoList.clear()
        self.changeRemoveInfo()
        self.window.csvColumnList.clear()
        self.changeCsvColumn()
        self.window.bedScoreList.clear()
        self.changeBedScore()
        self.window.allAttributeComboBox.clear()
        
        # TODO: color the > n VALUES! elements red
        self.window.filterList.addItems([str(i) for i in self.variantFilters])
        self.window.includeAttributeList.addItems([str(i) for i in self.includedAttributes])
        self.window.removeInfoList.addItems([str(i) for i in self.removedAttributes])
        self.window.csvColumnList.addItems(self.csvAttributes)
        self.window.bedScoreList.addItems(self.bedAttributes)
        
        self.window.allAttributeComboBox.addItem("")
        self.window.allAttributeComboBox.addItems(RESERVED_TAGS)
        self.window.allAttributeComboBox.insertSeparator(self.window.allAttributeComboBox.count()+1)
        self.window.allAttributeComboBox.addItems([a.globalName for a in self.includedAttributes])
        self.window.allAttributeComboBox.addItems([a.globalName for a in self.removedAttributes])
    
    def includeAttribute(self, attr):
        if attr in self.removedAttributes:
            self.removedAttributes.remove(attr)
        if not attr in self.includedAttributes:
            self.includedAttributes.append(attr)
        self.uniquifyAttributes()
    
    def excludeAttribute(self, attr):
        if attr in self.includedAttributes:
            self.includedAttributes.remove(attr)
        if not attr in self.removedAttributes and attr.name in self.infoFields.iterkeys():  # only list excluded attributes that are already in the .vcf file
            self.removedAttributes.append(attr)
        self.uniquifyAttributes()
    
    def uniquifyAttributes(self):
        takenTags = set(RESERVED_TAGS)
        for attr in self.includedAttributes:
            if isinstance(attr,attribute):
                name = attr.name
            else:
                name = attr.population + "_" + attr.function
            temp = name
            appendDigit = 2
            while temp in takenTags:
                temp = "%s_#%i" % (name,appendDigit)
                appendDigit += 1
            attr.globalName = temp
            takenTags.add(temp)
        
        takenTags = set(RESERVED_TAGS)
        for attr in self.removedAttributes:
            if isinstance(attr,attribute):
                name = attr.name
            else:
                name = attr.population + "_" + attr.function
            temp = name
            appendDigit = 2
            while temp in takenTags:
                temp = "%s_#%i" % (name,appendDigit)
                appendDigit += 1
            attr.globalName = temp
            takenTags.add(temp)
        
    # GUI Event handlers
    def browseVcf(self):
        fileName = QFileDialog.getOpenFileName(caption=u"Open .vcf file", filter=u"Variant Call Files (*.vcf)")[0]
        if fileName == '':
            return
        
        # count number of categorical values - this probably needs a progress bar
        progress = QProgressDialog(u"Loading %s" % os.path.split(fileName)[1], u"Cancel", 0, 1000, parent=None)
        progress.setWindowModality(Qt.WindowModal)
        progress.show()
        
        def tick():
            if progress.wasCanceled():
                return False
            newValue = min(progress.maximum(),progress.value()+1)
            progress.setValue(newValue)
            return True
        
        self.infoFields = extractInfoFields(fileName,tickFunction=tick)
        
        progress.close()
        
        if self.infoFields == None:
            self.window.vcfPathField.setText("")
            self.hasValidInput = False
            
            self.infoFields = {}
            self.includedAttributes = [] # contains actual attribute objects
            self.removedAttributes = [] # contains actual attribute objects
            
            self.activeCsvFile = None
            self.csvAttributes = [] # really just contains strings
            self.activeBedFile = None
            self.bedAttributes = [] # really just contains strings
            
            self.calculatedAttributes = []
            self.updateLists()
        else:
            self.window.vcfPathField.setText(fileName)
            self.hasValidInput = True
            
            self.includedAttributes = [] # contains actual attribute objects
            self.removedAttributes = [] # contains actual attribute objects
            
            for f in self.infoFields.itervalues():
                if f.maxedOut:
                    self.excludeAttribute(attribute(f.id, fileName, ">%i VALUES!" % MAX_INFO_STRINGS))
                else:
                    self.includeAttribute(attribute(f.id, fileName))
            
            self.activeCsvFile = None
            self.csvAttributes = [] # really just contains strings
            self.activeBedFile = None
            self.bedAttributes = [] # really just contains strings
            
            self.calculatedAttributes = []
            self.updateLists()
            
        self.globalEnable()
    
    def browseOutput(self):
        fileName = QFileDialog.getSaveFileName(caption=u"Save file", filter=u"Variant Call File (*.vcf);;Tabular File (*.csv);;CompreheNGSive Variant File (*.cvf)")[0]
        if fileName == '':
            return
        
        # just in case they type an extension that doesn't match the menu item they chose (which can happen), let's prefer what they typed
        f1,e1 = os.path.splitext(fileName)
        f2,e2 = os.path.splitext(f1)
        if e2.lower() in ['.vcf','.csv','.cvf']:
            fileName = f1
        
        self.window.outputPathField.setText(fileName)
        self.hasValidOutput = True
        self.globalEnable()
    
    def askForPopFile(self):
        temp = self.lastBackgroundPopIndex
        fileName = QFileDialog.getOpenFileName(caption=u"Open .txt file", filter=u"Population text file (*.txt)")[0]
        if fileName != '':
            self.loadPopulations(fileName)
        self.lastBackgroundPopIndex = temp
        self.window.populationComboBox.setCurrentIndex(self.lastBackgroundPopIndex)
    
    def changePopulation(self):
        if self.window.populationComboBox.currentText() == "Load population .txt file...":
            self.askForPopFile()
    
    def changeIncludeAttribute(self):
        if len(self.window.includeAttributeList.selectedItems()) > 0:
            self.window.removeAttributeButton.setEnabled(True)
        else:
            self.window.removeAttributeButton.setEnabled(False)
    
    def removeAttributes(self):
        attrs = []
        for a in self.window.includeAttributeList.selectedItems():
            attrs.append(self.includedAttributes[self.window.includeAttributeList.indexFromItem(a).row()])
        
        for a in attrs:
            self.includedAttributes.remove(a)
            if isinstance(a,attribute) and self.window.vcfPathField.text() == a.sourcePath:
                self.removedAttributes.append(a)
        self.uniquifyAttributes()
        self.updateLists()
    
    def changeRemoveInfo(self):
        if len(self.window.removeInfoList.selectedItems()) > 0:
            self.window.addInfoButton.setEnabled(True)
        else:
            self.window.addInfoButton.setEnabled(False)
    
    def addInfo(self):
        attrs = []
        for a in self.window.removeInfoList.selectedItems():
            attrs.append(self.removedAttributes[self.window.removeInfoList.indexFromItem(a).row()])
        
        for a in attrs:
            self.removedAttributes.remove(a)
            self.includedAttributes.append(a)
        self.uniquifyAttributes()
        self.updateLists()
    
    def browseCsv(self):
        fileName = QFileDialog.getOpenFileName(caption=u"Open .csv file", filter=u"Tabular file (*.csv)")[0]
        if fileName != '':
            self.activeCsvFile = fileName
            self.csvAttributes = []
            try:
                headers,chromColumn,posColumn,idColumn = sniffCsv(fileName)[1:]
                for i,a in enumerate(headers):
                    if i != chromColumn and i != posColumn and i != idColumn:
                        self.csvAttributes.append(a)
            except genomeException:
                self.activeCsvFile = None
                self.csvAttributes = []
                m = QMessageBox()
                m.setText("The .csv file was formatted incorrectly; make sure it has a \"CHROM\", \"POS\", and optionally an \"ID\" header.")
                m.setIcon(QMessageBox.Warning)
                m.exec_()
            self.updateLists()
    
    def changeCsvColumn(self):
        if len(self.window.csvColumnList.selectedItems()) > 0:
            self.window.addCsvColumnButton.setEnabled(True)
        else:
            self.window.addCsvColumnButton.setEnabled(False)
    
    def addCsvAttribute(self):
        if self.window.exactRadioButton.isChecked():
            mode = "Exact"
        elif self.window.copyRadioButton.isChecked():
            mode = "Copy"
        else:
            mode = "Interpolate"
        for a in self.window.csvColumnList.selectedItems():
            attr = attribute(a.text(),self.activeCsvFile,mode)
            if attr not in self.includedAttributes:
                self.includedAttributes.append(attr)
        self.uniquifyAttributes()
        self.updateLists()
    
    def browseBed(self):
        fileName = QFileDialog.getOpenFileName(caption=u"Open .bed file", filter=u"BED file (*.bed)")[0]
        if fileName != '':
            self.activeBedFile = fileName
            self.bedAttributes = []
            try:
                scoreNames,bedRegions = sniffBed(fileName)
                for a in scoreNames:
                    self.bedAttributes.append(a)
            except genomeException:
                self.activeBedFile = None
                self.bedAttributes = []
                m = QMessageBox()
                m.setText("The .bed file must have names and scores for every feature in order to attach scores to variants.")
                m.setIcon(QMessageBox.Warning)
                m.exec_()
            self.updateLists()
    
    def changeBedScore(self):
        if len(self.window.bedScoreList.selectedItems()) > 0:
            self.window.addBedScoreButton.setEnabled(True)
        else:
            self.window.addBedScoreButton.setEnabled(False)
    
    def addBedAttribute(self):
        for a in self.window.bedScoreList.selectedItems():
            attr = attribute(a.text(),self.activeBedFile)
            if not attr in self.includedAttributes:
                self.includedAttributes.append(attr)
        self.uniquifyAttributes()
        self.updateLists()
    
    def addAdditionalAttribute(self):
        infile = QFile("ui/Calculate.ui")
        infile.open(QFile.ReadOnly)
        self.childWindow = self.loader.load(infile, None)
        infile.close()
        
        # some helper functions
        def appendRow(function=None,population=None,note=""):
            newRow = self.childWindow.tableWidget.rowCount()
            self.childWindow.tableWidget.setRowCount(newRow+1)
            
            funcWidget = QComboBox()
            funcWidget.addItems(allStats.STAT_NAMES)
            self.childWindow.tableWidget.setCellWidget(newRow,0,funcWidget)
            if function != None:
                funcWidget.setCurrentIndex(allStats.STAT_NAMES.index(function))
            
            popWidget = QComboBox()
            selectedIndex = 0
            for i,p in enumerate(self.populationOrder):
                if p == None:
                    if not i == 0:
                        #self.window.populationComboBox.addItem(p[0])
                        popWidget.insertSeparator(i)   # TODO: put the file name in instead
                else:
                    if p == population:
                        selectedIndex = i-1
                    popWidget.addItem(p)
            self.childWindow.tableWidget.setCellWidget(newRow,1,popWidget)
            popWidget.setCurrentIndex(selectedIndex)
            
            self.childWindow.tableWidget.setItem(newRow,2,QTableWidgetItem(note))
            return popWidget
        
        def addFields():
            statsToAdd = []
            for r in xrange(self.childWindow.tableWidget.rowCount()):
                func = self.childWindow.tableWidget.cellWidget(r,0).currentText()
                pop = self.childWindow.tableWidget.cellWidget(r,1).currentText()
                note = self.childWindow.tableWidget.takeItem(r,2).text()
                
                stat = statistic(func,pop,note)
                if stat not in self.includedAttributes:
                    self.includedAttributes.append(stat)
                statsToAdd.append(stat)
            
            for a in self.removedAttributes:
                if isinstance(a,statistic) and not a in statsToAdd:
                    self.removedAttributes.remove(a)
            for a in self.includedAttributes:
                if isinstance(a,statistic) and not a in statsToAdd:
                    self.includedAttributes.remove(a)
            
            self.childWindow.accept()
            self.childWindow = None
            self.childWidgets = []
            
            self.uniquifyAttributes()
            self.updateLists()
                
        def addStat():
            self.childWidgets.append(appendRow())
        
        def removeStat():
            currentRow = self.childWindow.tableWidget.currentRow()
            self.childWidgets.pop(currentRow)
            self.childWindow.tableWidget.removeRow(currentRow)
        
        def changeField():
            if self.childWindow.tableWidget.currentRow() != -1:
                self.childWindow.removeStatButton.setEnabled(True)
            else:
                self.childWindow.removeStatButton.setEnabled(False)
        
        # okay, populate the table with what already exists
        self.childWidgets = []
        
        for attr in self.includedAttributes:
            if isinstance(attr,statistic):
                self.childWidgets.append(appendRow(attr.function,attr.population,attr.note))
        
        # and link up events in our dialog
        self.childWindow.calculateOkButton.clicked.connect(addFields)
        self.childWindow.browsePopulationButton.clicked.connect(self.askForPopFile)
        self.childWindow.addStatButton.clicked.connect(addStat)
        self.childWindow.removeStatButton.clicked.connect(removeStat)
        self.childWindow.calculateHelpButton.clicked.connect(self.openHelpPage)
        self.childWindow.tableWidget.itemSelectionChanged.connect(changeField)
        changeField()
        
        self.childWindow.show()
    
    def changeFilter(self):
        f = self.variantFilters[self.window.filterList.currentIndex().row()]
        if isinstance(f,filterExpression):
            self.window.expressionField.setPlainText(f.expression % tuple(f.columns))
        else:
            self.window.expressionField.setPlainText("")
    
    def removeFilter(self):
        self.variantFilters.pop(self.window.filterList.currentIndex().row())
        self.updateLists()
    
    def addBedFilter(self):
        fileName = QFileDialog.getOpenFileName(caption=u"Open .bed file", filter=u"BED File (*.bed)")[0]
        if fileName == '':
            return
        self.variantFilters.append(filterBed(fileName))
        self.updateLists()
    
    def insertAttribute(self):
        f = self.window.allAttributeComboBox.currentText()
        if f != "":
            self.window.expressionField.insertPlainText(f)
        self.window.allAttributeComboBox.setCurrentIndex(0)
    
    def addExpressionFilter(self):
        expr = [self.window.expressionField.toPlainText()]
        
        allTags = [""]
        allTags.extend(RESERVED_TAGS)
        allTags.extend([a.globalName for a in self.includedAttributes])
        allTags.extend([a.globalName for a in self.removedAttributes])
        
        # extract out the variables by splitting up the expression
        for a in allTags[1:]:
            newExpr = []
            f = fieldBit(a)
            for e in expr:
                if isinstance(e,fieldBit):
                    newExpr.append(e)
                else:
                    codeBits = e.split(a)
                    firstBit = True
                    for c in codeBits:
                        if firstBit:
                            firstBit = False
                        else:
                            newExpr.append(f)
                        newExpr.append(c)
            expr = newExpr
        
        # now assemble the expression string and the columns that will be fed in
        exprStr = ""
        columns = []
        for e in expr:
            if isinstance(e,fieldBit):
                exprStr += "%s"
                columns.append(e.field)
            else:
                exprStr += e
        
        result = filterExpression(exprStr,columns)
        if result not in self.variantFilters:
            self.variantFilters.append(result)
            self.updateLists()
    
    def openHelpPage(self):
        webbrowser.open('http://sci.utah.edu/~abigelow/vcfCleanerHelp.php',new=2)
    
    def quit(self):
        self.window.reject()
    
    def run(self):
        # TODO: sort all the files, reorder alleles
        self.window.accept()

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = gui()
    sys.exit(app.exec_())