#!/usr/bin/env python
import sys, os, webbrowser
from PySide.QtCore import Qt, QFile
from PySide.QtUiTools import QUiLoader
from PySide.QtGui import QApplication, QFileDialog, QProgressDialog, QMessageBox, QComboBox, QTableWidgetItem, QTreeWidgetItem, QPushButton
from genome_utils import genomeException, parsePopulations, MAX_INFO_STRINGS
from cleanVCF import extractInfoFields
from addCSVtoVCF import sniffCsv
from addBEDtoVCF import sniffBed
from calcStats import allStats, parseVcfHeader

PYTHON_WORDS = ["and","assert","break","class","continue",
                "def","del","elif","else","except",
                "exec","finally","for","from","global",
                "if","import","in","is","lambda",
                "not","or","pass","print","raise",
                "return","try","while",
                "Data","Float","Int","Numeric","Oxphys",
                "array","close","float","int","input",
                "open","range","type","write","zeros",
                "acos","asin","atan","cos","e",
                "exp","fabs","floor","log","log10",
                "pi","sin","sqrt","tan"]
SPECIAL_WORDS = ["keep_variant"]
RESERVED_TAGS = ["CHROM","POS","ID","REF","ALT","QUAL","FILTER"]
PRESET_FILTERS = {'FILTER == PASS': 'keep_variant = FILTER == "PASS"',
                  'QUAL >= 30': 'keep_variant = QUAL >= 30'}

MAX_FILTER_STRING = 50

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
    def __init__(self, name, gui, source=False):
        self.name = name
        self.samples = set()
        self.source = source
        self.globalName = name
        
        self.gui = gui
    
    def addSample(self, s):
        if self.source:
            assert isinstance(s,str)
            s = sample(s,self)
        else:
            assert isinstance(s,sample)
        self.samples.add(s)
        return s
    
    def __repr__(self):
        return self.name
    
    def __eq__(self, other):
        return isinstance(other, population) and self.name == other.name and self.filePath == other.filePath
    
    def __ne__(self, other):
        return not self.__eq__(other)
    
    def addSamplesToPop(self):
        self.gui.addSamplesToPop(self)

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
    def __init__(self, targetPop, function, backPop, ascending, revertHack):
        self.targetPop = targetPop
        self.function = function
        self.backPop = backPop
        self.ascending = ascending
        self.revertHack = revertHack
        self.globalName = self.targetPop + "_" + self.function + "_" + ('ASC' if self.ascending else 'DEC') + '_' + self.backPop + ('_aH' if self.revertHack else '')
        if backPop == "ALT":
            self.description = "%s calculated by vcfCleaner, with allele-specific values listed in the same order as the ALT alleles"
        else:
            self.description = "%s calculated by vcfCleaner, with allele-specific values listed in order of %s allele frequency in the %s population. When %s is missing data, %s" % (self.function,
                                                                                                                                                                                      ('ascending' if self.ascending else 'descending'),
                                                                                                                                                                                      self.backPop,
                                                                                                                                                                                      self.backPop,
                                                                                                                                                                                      ('the original REF/ALT allele order is used (this technically makes "' + self.globalName + '" a misnomer' if self.revertHack else 'this will be "." even if ' + self.function + 'could have been calculated, because the allele order is undefined'))
    
    def __repr__(self):
        return self.globalName
    
    def __eq__(self, other):
        return isinstance(other,statistic) and self.targetPop == other.targetPop and self.function == other.function and self.backPop == other.backPop and self.ascending == other.ascending and self.revertHack == other.revertHack
    
    def __ne__(self, other):
        return not self.__eq__(other)

class filterExpression(object):
    def __init__(self, expression, columns):
        self.expression = expression
        self.columns = columns
    
    def __repr__(self):
        temp = "EXPRESSION: " + self.expression % tuple(self.columns)
        if len(temp) > MAX_FILTER_STRING:
            temp = temp[:MAX_FILTER_STRING-3] + "..."
        return temp
    
    def __eq__(self, other):
        return isinstance(other,filterExpression) and self.expression == other.expression and self.columns == other.columns
    
    def __ne__(self, other):
        return not self.__eq__(other)

class filterBed(object):
    def __init__(self, path):
        self.path = path
    
    def __repr__(self):
        temp = "BED: %s" % os.path.split(self.path)[1]
        if len(temp) > MAX_FILTER_STRING:
            temp = "..." + temp[len(temp)-MAX_FILTER_STRING+3:]
        return temp
    
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
        self.kgpPops = set()
        self.kgpFiles = {}
        for i in xrange(22):
            self.kgpFiles['chr%i' % (i+1)] = None
        self.kgpFiles['chrX'] = None
        self.vcfPop = None
        self.loadPopulations(sys.path[0] + "/KGP_populations.txt", source=True, kgp=True)
                
        self.variantFilters = []
        
        self.infoFields = {}
        
        self.includedAttributes = [] # contains actual attribute objects
        self.removedAttributes = [] # contains actual attribute objects
        
        self.csvAttributes = [] # really just contains strings
        self.bedAttributes = [] # really just contains strings
        self.calculatedAttributes = [] # really just contains strings
        
        # set up GUI - Data Sources Tab
        self.window.browseVcfButton.clicked.connect(self.browseVcf)
        self.window.browseKgpButton.clicked.connect(self.browseKGP)
        self.window.browseOutputButton.clicked.connect(self.browseOutput)
        self.window.browseLogButton.clicked.connect(self.browseLog)
        self.window.browseErrorButton.clicked.connect(self.browseError)
        self.window.browseNonBiallelicButton.clicked.connect(self.browseNonBiallelic)
        self.window.compressOutputCheckbox.stateChanged.connect(self.compressOutput)
        
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
        
        self.window.functionComboBox.addItems(allStats.STAT_NAMES)
        self.window.functionComboBox.currentIndexChanged.connect(self.calculateEnable)
        self.window.targetPopComboBox.currentIndexChanged.connect(self.calculateEnable)
        self.window.alleleReorderRadioButton.clicked.connect(self.calculateEnable)
        self.window.alleleVcfOrderRadioButton.clicked.connect(self.calculateEnable)
        self.window.backPopComboBox.currentIndexChanged.connect(self.calculateEnable)
        self.window.revertToRefAltCheckBox.clicked.connect(self.calculateEnable)
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
        self.window.presetFiltersComboBox.addItems(sorted(PRESET_FILTERS.iterkeys()))
        self.window.presetFiltersComboBox.currentIndexChanged.connect(self.switchToPresetFilter)
        self.window.addExpressionFilterButton.clicked.connect(self.addExpressionFilter)
        
        self.window.helpButton.clicked.connect(self.openHelpPage)
        self.window.quitButton.clicked.connect(self.quit)
        self.window.runButton.clicked.connect(self.run)
        
        self.globalEnable()
        
        self.window.show()
    
    # ***** Helper methods *****
    
    def loadPopulations(self, path, source=False, kgp=False, vcf=False):
        if vcf:
            newPops = parseVcfHeader(path)[2]
            popOrder = [newPops.iterkeys().next()]
        else:
            newPops,popOrder = parsePopulations(path)
        
        popsToAdd = []
        
        for k in popOrder:
            newPop = population(k,self,source)
            appendDigit = 2
            while self.populations.has_key(newPop.globalName):
                newPop.globalName = "%s_%i" % (k, appendDigit)
                appendDigit += 1
            if kgp:
                self.kgpPops.add(newPop.globalName)
            elif vcf:
                self.vcfPop = newPop.globalName
            for s in newPops[k]:
                if source:
                    s = newPop.addSample(s)
                    assert not self.samples.has_key(s.globalName)
                    self.samples[s.globalName] = s
                else:
                    if not self.samples.has_key("%s (%s)" % (s,k)):
                        m = QMessageBox()
                        m.setText("Couldn't load population file; unknown sample: %s (%s)" % (s,k))
                        m.setIcon(QMessageBox.Critical)
                        m.exec_()
                        return
                    else:
                        newPop.addSample(self.samples["%s (%s)" % (s,k)])
            popsToAdd.append(newPop)
        # Now that we loaded without an error, we can add the populations
        for p in popsToAdd:
            self.populations[p.globalName] = p
        
        self.updatePopLists()
    
    def updatePopLists(self):
        popOrder = sorted(self.populations.iterkeys())
        includeKGP = os.path.isdir(self.window.kgpPathField.text())
        for cBox in [self.window.targetPopComboBox,self.window.backPopComboBox]:
            cBox.clear()
            cBox.addItem("Select...")
            for p in popOrder:
                if not includeKGP and p in self.kgpPops:
                    continue
                cBox.addItem(p)
            cBox.setCurrentIndex(0)
        
        self.window.treeWidget.clear()
        self.popButtonBag = set()
        
        for p in reversed(popOrder):
            if not includeKGP and p in self.kgpPops:
                continue
            pop = self.populations[p]
            parent = QTreeWidgetItem([p,''])
            self.window.treeWidget.insertTopLevelItem(0,parent)
            if pop.source:
                parent.setFlags(parent.flags() & Qt.ItemIsEnabled)
                addSamplesToPopButton = QPushButton('Add Samples...')
                self.popButtonBag.add(addSamplesToPopButton)    # a hack to keep a pointer to the button around... otherwise we seg fault
                addSamplesToPopButton.setDisabled(True) # these are here just for decoration
                self.window.treeWidget.setItemWidget(parent,1,addSamplesToPopButton)
                for s in sorted(pop.samples):
                    child = QTreeWidgetItem(parent)
                    child.setText(0,s.name)
                    child.setFlags(child.flags() & Qt.ItemIsEnabled)
            else:
                parent.setFlags(parent.flags() | Qt.ItemIsEditable)
                addSamplesToPopButton = QPushButton('Add Samples...')
                self.popButtonBag.add(addSamplesToPopButton)    # a hack to keep a pointer to the button around... otherwise we seg fault
                self.window.treeWidget.setItemWidget(parent,1,addSamplesToPopButton)
                addSamplesToPopButton.clicked.connect(pop.addSamplesToPop)
                for s in sorted(pop.samples):
                    child = QTreeWidgetItem(parent)
                    child.setText(0,s.globalName)
    
    def globalEnable(self):
        if os.path.exists(self.window.vcfPathField.text()) and self.window.outputPathField.text() != "":
            self.window.tabWidget.widget(1).setEnabled(True)
            self.window.tabWidget.widget(2).setEnabled(True)
            self.window.runButton.setEnabled(True)
        else:
            self.window.tabWidget.widget(1).setEnabled(False)
            self.window.tabWidget.widget(2).setEnabled(False)
            self.window.runButton.setEnabled(False)
    
    def updateAttrLists(self):
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
        
        self.window.allAttributeComboBox.addItem("Paste...")
        self.window.allAttributeComboBox.addItems(SPECIAL_WORDS)
        self.window.allAttributeComboBox.insertSeparator(self.window.allAttributeComboBox.count()+1)
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
        takenTags.update(PYTHON_WORDS)
        takenTags.update(SPECIAL_WORDS)
        for attr in self.includedAttributes:
            if isinstance(attr,attribute):
                name = attr.name
            else:
                name = attr.targetPop + "_" + attr.function + "_" + ('ASC' if attr.ascending else 'DEC') + '_' + attr.backPop + ('_aH' if attr.revertHack else '')
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
    
    def addSamplesToPop(self, pop):
        loader = QUiLoader()
        infile = QFile("ui/AddSamples.ui")
        infile.open(QFile.ReadOnly)
        window = loader.load(infile, None)
        infile.close()
        
        includeKGP = os.path.isdir(self.window.kgpPathField.text())
        
        for s in sorted(self.samples.itervalues()):
            if not includeKGP:
                keep = True
                for p in self.kgpPops:
                    if s in self.populations[p].samples:
                        keep = False
                        break
                if not keep:
                    continue
            window.listWidget.addItem(s.globalName)
        
        def localAddAll():
            for i in window.listWidget.selectedItems():
                pop.addSample(self.samples[i.text()])
            self.updatePopLists()
        
        def updateCount():
            numSamples = len(window.listWidget.selectedItems())
            window.label.setText("%i Samples Selected" % numSamples)
        
        window.listWidget.itemSelectionChanged.connect(updateCount)
        window.accepted.connect(localAddAll)
        
        window.show()
    
    def updateFilterList(self):
        self.window.filterList.clear()
        for f in self.variantFilters:
            self.window.filterList.addItem(str(f))
    
    # ***** GUI Bindings *****
    def browseVcf(self):
        fileName = QFileDialog.getOpenFileName(caption=u"Open .vcf file", filter=u"Variant Call Files (*.vcf);;Gzip-compressed VCF (*.vcf.gz)")[0]
        if fileName == '':
            return
        
        # count number of categorical values - this probably needs a progress bar
        progress = QProgressDialog(u"Scanning %s" % os.path.split(fileName)[1], u"Cancel", 0, 1000, parent=None)
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
        
        if self.vcfPop != None:
            for s in self.populations[self.vcfPop].samples:
                for p,pop in self.populations.iteritems():
                    if p == self.vcfPop:
                        continue
                    pop.samples.discard(s)
            del self.populations[self.vcfPop]
            self.vcfPop = None
        
        if self.infoFields == None:
            self.window.vcfPathField.setText("")
            self.updatePopLists()
            
            self.infoFields = {}
            self.includedAttributes = [] # contains actual attribute objects
            self.removedAttributes = [] # contains actual attribute objects
            
            self.calculatedAttributes = []
            self.updateAttrLists()
        else:
            self.window.vcfPathField.setText(fileName)
            self.loadPopulations(fileName, source=True, kgp=False, vcf=True)
                        
            self.includedAttributes = [] # contains actual attribute objects
            self.removedAttributes = [] # contains actual attribute objects
            
            for f in self.infoFields.itervalues():
                if f.maxedOut:
                    self.excludeAttribute(attribute(f.id, fileName, ">%i VALUES!" % MAX_INFO_STRINGS))
                else:
                    self.includeAttribute(attribute(f.id, fileName))
            
            self.calculatedAttributes = []
            self.updateAttrLists()
            
        self.globalEnable()
    
    def browseKGP(self):
        newDir = QFileDialog.getExistingDirectory(caption=u"Select the folder containing KGP .vcf files",options=QFileDialog.ShowDirsOnly)
        if newDir == '':
            return
        
        for i in xrange(22):
            self.kgpFiles['chr%i' % (i+1)] = None
        self.kgpFiles['chrX'] = None
        
        for f in os.listdir(newDir):
            if 'genotypes.vcf' in f:
                c = f[f.find('chr'):]
                c = f[:f.find('.')]
                self.kgpFiles[c] = os.path.join(newDir,f)
        for c,f in self.kgpFiles.iteritems():
            if f == None:
                newDir = ''
                for i in xrange(22):
                    self.kgpFiles['chr%i' % (i+1)] = None
                self.kgpFiles['chrX'] = None
                m = QMessageBox()
                m.setText("Couldn't find the %s .genotypes.vcf file" % c)
                m.setIcon(QMessageBox.Critical)
                m.exec_()
                break
        self.window.kgpPathField.setText(newDir)
        self.updatePopLists()
    
    def browseOutput(self):
        fileName = QFileDialog.getSaveFileName(caption=u"Save file", filter=u"Variant Call File (*.vcf);;Tabular File (*.csv);;CompreheNGSive Variant File (*.cvf)")[0]
        if fileName == '':
            return
        e = os.path.splitext(fileName)[1].lower()
        
        if e == '.vcf':
            self.window.removeAttrLabel.setText("Remove these INFO attributes:")
        elif e == '.csv':
            self.window.removeAttrLabel.setText("Remove these columns:")
        elif e == '.cvf':
            self.window.removeAttrLabel.setText("Flag these columns as IGNORE:")
        
        if self.window.compressOutputCheckbox.isChecked():
            fileName += ".gz"
        
        if fileName == self.window.errorPathField.text() or fileName == self.window.nonBiallelicField.text() or fileName == self.window.vcfPathField.text():
            m = QMessageBox()
            m.setText("You are already using that file name; please choose another.")
            m.setIcon(QMessageBox.Critical)
            m.exec_()
            return
        
        self.window.outputPathField.setText(fileName)
        self.window.saveErrorWidget.setEnabled(True)
        self.window.saveNonBiallelicWidget.setEnabled(True)
        self.globalEnable()
    def browseLog(self):
        fileName = QFileDialog.getSaveFileName(caption=u"Save Log file", filter=u"Log File (*.log);;Text File (*.txt)")[0]
        if fileName == '':
            return
        self.window.logPathField.setText(fileName)
    def browseError(self):
        e = self.window.outputPathField.text().lower()
        if not self.window.compressOutputCheckbox.isChecked():
            e += ".gz"
        if e.endswith('.vcf.gz'):
            f = u"Variant Call File (*.vcf)"
        elif e.endswith('.csv.gz'):
            f = u"Tabular File (*.csv)"
        elif e.endswith('.cvf.gz'):
            f = u"CompreheNGSive Variant File (*.cvf)"
        
        fileName = QFileDialog.getSaveFileName(caption=u"Save file", filter=f)[0]
        if fileName == '':
            return
        
        if self.window.compressOutputCheckbox.isChecked():
            fileName += ".gz"
        
        if fileName == self.window.outputPathField.text() or fileName == self.window.nonBiallelicField.text() or fileName == self.window.vcfPathField.text():
            m = QMessageBox()
            m.setText("You are already using that file name; please choose another.")
            m.setIcon(QMessageBox.Critical)
            m.exec_()
            return
        
        self.window.errorPathField.setText(fileName)
        self.globalEnable()
    def browseNonBiallelic(self):
        e = self.window.outputPathField.text().lower()
        if not self.window.compressOutputCheckbox.isChecked():
            e += ".gz"
        if e.endswith('.vcf.gz'):
            f = u"Variant Call File (*.vcf)"
        elif e.endswith('.csv.gz'):
            f = u"Tabular File (*.csv)"
        elif e.endswith('.cvf.gz'):
            f = u"CompreheNGSive Variant File (*.cvf)"
        
        fileName = QFileDialog.getSaveFileName(caption=u"Save file", filter=f)[0]
        if fileName == '':
            return
        
        if self.window.compressOutputCheckbox.isChecked():
            fileName += ".gz"
        
        if fileName == self.window.outputPathField.text() or fileName == self.window.errorPathField.text() or fileName == self.window.vcfPathField.text():
            m = QMessageBox()
            m.setText("You are already using that file name; please choose another.")
            m.setIcon(QMessageBox.Critical)
            m.exec_()
            return
        
        self.window.nonBiallelicField.setText(fileName)
        self.globalEnable()
    def compressOutput(self):
        for f in [self.window.outputPathField,self.window.errorPathField,self.window.nonBiallelicField]:
            if f.text() != "":
                if self.window.compressOutputCheckbox.isChecked():
                    f.setText(f.text() + ".gz")
                else:
                    f.setText(f.text()[:-3])
    def createPopulation(self):
        newPop = population("Population",self,source=False)
        appendDigit = 2
        while self.populations.has_key(newPop.globalName):
            newPop.globalName = "%s_%i" % ("Population", appendDigit)
            appendDigit += 1
        self.populations[newPop.globalName] = newPop
        self.updatePopLists()
    def removePopOrSample(self):
        for i in self.window.treeWidget.selectedItems():
            if not isinstance(i.parent(),QTreeWidgetItem):  # really i.parent() == None would make more sense but a PySide bug doesn't allow that comparison
                # delete a whole population
                if not self.populations[i.text(0)].source:
                    del self.populations[i.text(0)]
            else:
                # delete a sample from a population
                if not self.populations[i.parent().text(0)].source:
                    self.populations[i.parent().text(0)].samples.discard(self.samples[i.text(0)])
        self.updatePopLists()
    
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
            elif isinstance(a,statistic):
                self.calculatedAttributes.remove(a.globalName)
        self.uniquifyAttributes()
        self.updateAttrLists()
        self.calculateEnable()
    
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
        self.updateAttrLists()
    
    def calculateEnable(self):
        statSelected = self.window.functionComboBox.currentIndex() != 0 and self.window.targetPopComboBox.currentIndex() != 0
        if self.window.alleleReorderRadioButton.isChecked():
            self.window.ascendingRadioButton.setEnabled(True)
            self.window.descendingRadioButton.setEnabled(True)
            self.window.backPopComboBox.setEnabled(True)
            self.window.revertToRefAltCheckBox.setEnabled(True)
            if self.window.backPopComboBox.currentIndex() != 0:
                tempStat = statistic(self.window.targetPopComboBox.currentText(),
                                     self.window.functionComboBox.currentText(),
                                     self.window.backPopComboBox.currentText(),
                                     self.window.ascendingRadioButton.isChecked(),
                                     self.window.revertToRefAltCheckBox.isChecked())
            else:
                tempStat = None
        else:
            self.window.ascendingRadioButton.setEnabled(False)
            self.window.descendingRadioButton.setEnabled(False)
            self.window.backPopComboBox.setEnabled(False)
            self.window.revertToRefAltCheckBox.setEnabled(False)
            tempStat = statistic(self.window.targetPopComboBox.currentText(),
                                 self.window.functionComboBox.currentText(),
                                 "ALT",
                                 False,
                                 False)
        if statSelected and tempStat != None and tempStat.globalName not in self.calculatedAttributes:
            self.window.createAttributeButton.setEnabled(True)
        else:
            self.window.createAttributeButton.setEnabled(False)
        
    def createAttribute(self):
        if self.window.alleleReorderRadioButton.isChecked():
            newStat = statistic(self.window.targetPopComboBox.currentText(),
                                 self.window.functionComboBox.currentText(),
                                 self.window.backPopComboBox.currentText(),
                                 self.window.ascendingRadioButton.isChecked(),
                                 self.window.revertToRefAltCheckBox.isChecked())
        else:
            newStat = statistic(self.window.targetPopComboBox.currentText(),
                                 self.window.functionComboBox.currentText(),
                                 "ALT",
                                 False,
                                 False)
        assert newStat.globalName not in self.calculatedAttributes
        self.calculatedAttributes.append(newStat.globalName)
        self.includeAttribute(newStat)
        self.updateAttrLists()
        self.calculateEnable()
    
    def browseCsv(self):
        fileName = QFileDialog.getOpenFileName(caption=u"Open .csv file", filter=u"Tabular file (*.csv)")[0]
        if fileName != '':
            self.window.csvPathField.setText(fileName)
            self.csvAttributes = []
            try:
                headers,chromColumn,posColumn,idColumn = sniffCsv(fileName)[1:]
                for i,a in enumerate(headers):
                    if i != chromColumn and i != posColumn and i != idColumn:
                        self.csvAttributes.append(a)
            except genomeException:
                self.window.csvPathField.setText("")
                self.csvAttributes = []
                m = QMessageBox()
                m.setText("The .csv file was formatted incorrectly; make sure it has a \"CHROM\", \"POS\", and optionally an \"ID\" header.")
                m.setIcon(QMessageBox.Warning)
                m.exec_()
            self.updateAttrLists()
    
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
            attr = attribute(a.text(),self.window.csvPathField.text(),mode)
            if attr not in self.includedAttributes:
                self.includedAttributes.append(attr)
        self.uniquifyAttributes()
        self.updateAttrLists()
    
    def browseBed(self):
        fileName = QFileDialog.getOpenFileName(caption=u"Open .bed file", filter=u"BED file (*.bed)")[0]
        if fileName != '':
            self.window.bedPathField.setText(fileName)
            self.bedAttributes = []
            try:
                scoreNames,bedRegions = sniffBed(fileName)
                for a in scoreNames:
                    self.bedAttributes.append(a)
            except genomeException:
                self.window.bedPathField.setText("")
                self.bedAttributes = []
                m = QMessageBox()
                m.setText("The .bed file must have names and scores for every feature in order to attach scores to variants.")
                m.setIcon(QMessageBox.Warning)
                m.exec_()
            self.updateAttrLists()
    
    def changeBedScore(self):
        if len(self.window.bedScoreList.selectedItems()) > 0:
            self.window.addBedScoreButton.setEnabled(True)
        else:
            self.window.addBedScoreButton.setEnabled(False)
    
    def addBedAttribute(self):
        for a in self.window.bedScoreList.selectedItems():
            attr = attribute(a.text(),self.window.bedPathField.text())
            if not attr in self.includedAttributes:
                self.includedAttributes.append(attr)
        self.uniquifyAttributes()
        self.updateAttrLists()
    
    # set up GUI - Variant Filters Tab
    def changeFilter(self):
        f = self.variantFilters[self.window.filterList.currentRow()]
        self.window.expressionField.setPlainText(f.expression % tuple(f.columns))
    def removeFilter(self):
        if self.window.filterList.count() == 0:
            return
        del self.variantFilters[self.window.filterList.currentRow()]
        self.updateFilterList()
    def addBedFilter(self):
        fileName = QFileDialog.getOpenFileName(caption=u"Open .bed file", filter=u"BED file (*.bed)")[0]
        if fileName != '':
            newFilter = filterBed(fileName)
            if newFilter not in self.variantFilters:
                self.variantFilters.append(newFilter)
            self.updateFilterList()
    def insertAttribute(self):
        if self.window.allAttributeComboBox.currentIndex() != 0:
            self.window.expressionField.insertPlainText(self.window.allAttributeComboBox.currentText())
            self.window.allAttributeComboBox.setCurrentIndex(0)
    def switchToPresetFilter(self):
        if self.window.presetFiltersComboBox.currentIndex() != 0:
            self.window.expressionField.setPlainText(PRESET_FILTERS[self.window.presetFiltersComboBox.currentText()])
            self.window.presetFiltersComboBox.setCurrentIndex(0)
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
            self.updateFilterList()
    
    def openHelpPage(self):
        webbrowser.open('http://sci.utah.edu/~abigelow/vcfCleanerHelp.php',new=2)
    
    def quit(self):
        self.window.reject()
    
    def run(self):
        # TODO: sort all the files, reorder alleles
        self.window.accept()



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