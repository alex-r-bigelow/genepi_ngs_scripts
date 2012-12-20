#!/usr/bin/env python
import sys, os, webbrowser
from PySide.QtCore import Qt, QFile
from PySide.QtUiTools import QUiLoader
from PySide.QtGui import QApplication, QFileDialog, QProgressDialog, QMessageBox, QComboBox, QTableWidgetItem, QTreeWidgetItem, QPushButton
from genome_utils import genomeException, parsePopulations, MAX_INFO_STRINGS
from cleanVCF import extractInfoFields
from addCSVtoVCF import sniffCsv
from addBEDtoVCF import sniffBed
from calcStats import allStats

RESERVED_TAGS = ["CHROM","POS","ID","QUAL","FILTER"]

class sample(object):
    def __init__(self, name, sourcePop):
        self.name = name
        self.sourcePop = sourcePop
        self.globalName = '%s (%s)' % (name, sourcePop.globalName)
    
    def __hash__(self):
        return hash(self.globalName)
    
    def __eq__(self, other):
        return isinstance(other, sample) and self.name == other.name and self.sourcePop == other.sourcePop
    
    def __ne__(self, other):
        return not self.__eq__(other)

class population(object):
    def __init__(self, name, allSamples, source=False):
        self.name = name
        self.samples = []
        self.source = source
        self.globalName = name
        
        self.allSamples = allSamples
    
    def addSample(self, s):
        if self.source:
            assert not isinstance(s,sample)
            s = sample(s,self)
        else:
            assert isinstance(s,sample)
        self.samples.append(s)
        return s
    
    def __repr__(self):
        return self.name
    
    def __eq__(self, other):
        return isinstance(other, population) and self.name == other.name and self.filePath == other.filePath
    
    def __ne__(self, other):
        return not self.__eq__(other)
    
    def addSamplesToPop(self):
        print self.name

class attribute(object):
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

class statistic(object):
    def __init__(self, targetPop, function, backPop, ascending, description=""):
        self.targetPop = targetPop
        self.function = function
        self.backPop = backPop
        self.ascending = ascending
        self.description = description
        self.globalName = self.targetPop + "_" + self.function + "_" + ('ASC' if self.ascending else 'DEC') + '_' + self.backPop
    
    def __repr__(self):
        return self.globalName
    
    def __eq__(self, other):
        return isinstance(other,statistic) and self.targetPop == other.targetPop and self.function == other.function and self.backPop == other.backPop and self.ascending == other.ascending
    
    def __ne__(self, other):
        return not self.__eq__(other)

class filterExpression(object):
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

class filterBed(object):
    def __init__(self, path):
        self.path = path
    
    def __repr__(self):
        temp = os.path.split(self.path)[1]
        return "BED: %s" % temp
    
    def __eq__(self, other):
        return isinstance(other,filterBed) and self.path == other.path
    
    def __ne__(self, other):
        return not self.__eq__(other)

class fieldBit(object):
    def __init__(self, field):
        self.field = field

class gui:
    def __init__(self):
        self.loader = QUiLoader()
        infile = QFile("ui/Main.ui")
        infile.open(QFile.ReadOnly)
        self.window = self.loader.load(infile, None)
        infile.close()
                
        # data structures
        self.popButtonBag = set()   # I get some strange seg fault errors due to some pyside garbage collection bug... I need to keep a pointer to
                                    # buttons I create dynamically to avoid this
        
        self.samples = {}
        self.populations = {}
        self.loadPopulations(sys.path[0] + "/KGP_populations.txt", source=True)
        
        self.variantFilters = []
        
        self.infoFields = {}
        
        self.includedAttributes = [] # contains actual attribute objects
        self.removedAttributes = [] # contains actual attribute objects
        
        self.activeCsvFile = None
        self.csvAttributes = [] # really just contains strings
        self.activeBedFile = None
        self.bedAttributes = [] # really just contains strings
        
        self.calculatedAttributes = []
        
        # set up GUI - Data Sources Tab
        self.window.browseVcfButton.clicked.connect(self.browseVcf)
        self.window.browseKgpButton.clicked.connect(self.browseKGP)
        self.window.browseOutputButton.clicked.connect(self.browseOutput)
        self.window.browseLogButton.clicked.connect(self.browseLog)
        self.window.browseErrorButton.clicked.connect(self.browseError)
        self.window.browseNonBiallelicButton.clicked.connect(self.browseNonBiallelic)
        self.window.compressOutputCheckbox.stateChanged.connect(self.compressOutput)
        
        self.window.browsePopulationButton.clicked.connect(self.browsePopulation)
        self.window.savePopulationButton.clicked.connect(self.savePopulation)
        self.window.createPopulationButton.clicked.connect(self.createPopulation)
        self.window.removePopOrSampleButton.clicked.connect(self.removePopOrSample)
        
        # Because the "Add Sample(s) to Population" and "Remove" buttons are created dynamically, we bind
        # the methods later, but I need to remember to bind to self.addSampleToPop and call self.removePopOrSample
        # appropriately
        
        # set up GUI - Attribute Settings Tab
        self.window.includeAttributeList.itemSelectionChanged.connect(self.changeIncludeAttribute)
        self.window.removeAttributeButton.clicked.connect(self.removeAttributes)
        
        self.window.removeInfoList.itemSelectionChanged.connect(self.changeRemoveInfo)
        self.window.addInfoButton.clicked.connect(self.addInfo)
        
        self.window.functionComboBox.currentIndexChanged.connect(self.calculateEnable)
        self.window.targetPopComboBox.currentIndexChanged.connect(self.calculateEnable)
        self.window.alleleReorderRadioButton.clicked.connect(self.calculateEnable)
        self.window.alleleVcfOrderRadioButton.clicked.connect(self.calculateEnable)
        self.window.backPopComboBox.currentIndexChanged.connect(self.calculateEnable)
        self.window.createAttributeButton.clicked.connect(self.createAttribute)
        
        self.window.browseCsvButton.clicked.connect(self.browseCsv)
        self.window.csvColumnList.itemSelectionChanged.connect(self.changeCsvColumn)
        self.window.addCsvColumnButton.clicked.connect(self.addCsvAttribute)
        
        self.window.browseBedButton.clicked.connect(self.browseBed)
        self.window.bedScoreList.itemSelectionChanged.connect(self.changeBedScore)
        self.window.addBedScoreButton.clicked.connect(self.addBedAttribute)
        
        # set up GUI - Variant Filters Tab
        self.window.filterList.itemSelectionChanged.connect(self.changeFilter)
        self.window.removeFilterButton.clicked.connect(self.removeFilter)
        self.window.browseBedFilterButton.clicked.connect(self.addBedFilter)
        self.window.allAttributeComboBox.currentIndexChanged.connect(self.insertAttribute)
        self.window.presetFiltersComboBox.currentIndexChanged.connect(self.switchToPresetFilter)
        self.window.addExpressionFilterButton.clicked.connect(self.addExpressionFilter)
        
        self.window.helpButton.clicked.connect(self.openHelpPage)
        self.window.quitButton.clicked.connect(self.quit)
        self.window.runButton.clicked.connect(self.run)
        
        self.globalEnable()
        
        self.window.show()
    
    # ***** Helper methods *****
    
    def loadPopulations(self, path, source=False):
        newPops,popOrder = parsePopulations(path)
        
        popsToAdd = []
        
        for k in popOrder:
            newPop = population(k,self.samples,source)
            appendDigit = 2
            while self.populations.has_key(newPop.globalName):
                newPop.globalName = "%s_#%i" % (k, appendDigit)
                appendDigit += 1
            for s in newPops[k]:
                if source:
                    s = newPop.addSample(s)
                    self.samples[s.name] = s    # This could potentially overwrite sample names we already loaded; this is okay (we'll already have
                                                # set up KGP samples, and the natural assumption is that a user-defined group will use the last-loaded
                                                # user data; the user can tweak this if it doesn't import how they expect)
                else:
                    if not self.samples.has_key(s):
                        m = QMessageBox()
                        m.setText("Couldn't load population file; unknown sample: %s" % s)
                        m.setIcon(QMessageBox.Critical)
                        m.exec_()
                        return
                    else:
                        newPop.addSample(self.samples[s])
            popsToAdd.append(newPop)
        # Now that we loaded without an error, we can add the populations
        for p in popsToAdd:
            self.populations[p.globalName] = p
        
        self.updatePopLists()
    
    def updatePopLists(self):
        popOrder = sorted(self.populations.iterkeys())
        for cBox in [self.window.targetPopComboBox,self.window.backPopComboBox]:
            cBox.clear()
            cBox.addItem("Select...")
            for p in popOrder:
                cBox.addItem(p)
            cBox.setCurrentIndex(0)
        
        self.window.treeWidget.clear()
        self.popButtonBag = set()
        
        for p in reversed(popOrder):
            pop = self.populations[p]
            parent = QTreeWidgetItem([p,''])
            self.window.treeWidget.insertTopLevelItem(0,parent)
            if pop.source:
                parent.setFlags(parent.flags() & Qt.ItemIsEnabled)
                addSamplesToPopButton = QPushButton('Add Samples...')
                self.popButtonBag.add(addSamplesToPopButton)    # a hack to keep a pointer to the button around... otherwise we seg fault
                addSamplesToPopButton.setDisabled(True) # these are here just for decoration
                self.window.treeWidget.setItemWidget(parent,1,addSamplesToPopButton)
                for s in pop.samples:
                    child = QTreeWidgetItem(parent)
                    child.setText(0,s.name)
                    child.setFlags(child.flags() & Qt.ItemIsEnabled)
            else:
                parent.setFlags(parent.flags() | Qt.ItemIsEditable)
                addSamplesToPopButton = QPushButton('Add Samples...')
                self.popButtonBag.add(addSamplesToPopButton)    # a hack to keep a pointer to the button around... otherwise we seg fault
                self.window.treeWidget.setItemWidget(parent,1,addSamplesToPopButton)
                addSamplesToPopButton.clicked.connect(pop.addSamplesToPop)
                for s in pop.samples:
                    child = QTreeWidgetItem(parent)
                    child.setText(0,s.globalName)
    
    def globalEnable(self):
        if os.path.exists(self.window.vcfPathField.getText()) and self.window.outputPathField != "":
            self.window.tabWidget.widget(1).setEnabled(True)
            self.window.tabWidget.widget(2).setEnabled(True)
        else:
            self.window.tabWidget.widget(1).setEnabled(False)
            self.window.tabWidget.widget(2).setEnabled(False)
        
        '''if self.hasValidOutput:
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
            self.window.runButton.setEnabled(False)'''
    
    # ***** GUI Bindings *****
    def browseVcf(self):
        pass
    def browseKGP(self):
        pass
    def browseOutput(self):
        pass
    def browseLog(self):
        pass
    def browseError(self):
        pass
    def browseNonBiallelic(self):
        pass
    def compressOutput(self):
        pass
    
    def browsePopulation(self):
        pass
    def savePopulation(self):
        pass
    def createPopulation(self):
        pass
    def removePopOrSample(self):
        pass
    
    
    
    def changeIncludeAttribute(self):
        pass
    def removeAttributes(self):
        pass
    
    def changeRemoveInfo(self):
        pass
    def addInfo(self):
        pass
    
    def calculateEnable(self):
        pass
    def createAttribute(self):
        pass
    
    def browseCsv(self):
        pass
    def changeCsvColumn(self):
        pass
    def addCsvAttribute(self):
        pass
    
    def browseBed(self):
        pass
    def changeBedScore(self):
        pass
    def addBedAttribute(self):
        pass
    
    # set up GUI - Variant Filters Tab
    def changeFilter(self):
        pass
    def removeFilter(self):
        pass
    def addBedFilter(self):
        pass
    def insertAttribute(self):
        pass
    def switchToPresetFilter(self):
        pass
    def addExpressionFilter(self):
        pass
    
    def openHelpPage(self):
        pass
    def quit(self):
        pass
    def run(self):
        pass



    '''
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
        
    # ***** GUI Event handlers ******
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
    
    def compressOutput(self):
        pass
    
    def browseKGP(self):
        pass
    
    def browsePopulation(self):
        pass
    
    def savePopulation(self):
        pass
    
    def createPopulation(self):
        pass
    
    def removePopOrSample(self):
        pass
    
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
        self.window.accept()'''

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = gui()
    sys.exit(app.exec_())