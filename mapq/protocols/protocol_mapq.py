# **************************************************************************
# *
# * Authors:     David Herreros Calero (dherreros@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

from os.path import abspath
from glob import glob

from pwem.convert import toCIF, Ccp4Header
from pwem.convert.atom_struct import toPdb, AtomicStructHandler, addScipionAttribute
from pwem.objects import SetOfAtomStructs, AtomStruct
from pwem.protocols import ProtAnalysis3D

from pyworkflow.protocol import PointerParam, FloatParam, MultiPointerParam, IntParam
from pyworkflow import BETA
import pyworkflow.utils as pwutils

import mapq


class ProtMapQ(ProtAnalysis3D):
    """
    Compute Q-Scores using MapQ software.
    """
    _label = 'compute q-scores'
    _devStatus = BETA
    _ATTRNAME = "MapQ_Score"
    _OUTNAME = "scoredStructures"

    # --------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputVol', PointerParam, pointerClass="Volume", label="Input map",
                      help='Volume to compare to PDBs', important=True)
        form.addParam('pdbs', MultiPointerParam, pointerClass="AtomStruct", important=True,
                      label="Input structures",
                      help='PDBs to compare to input map')
        form.addParam('mapRes', FloatParam, allowsNull=True,
                      label = "Map resolution",
                      help = "Optional - Default is 3.0 - Specifies resolution of map; it is used to output perresidue "
                             "statistics along with expected Q-score at this resolution")
        form.addParam('bFactor', IntParam, allowsNull=True,
                      label="B-factor",
                      help="Optional - If specified, a separate pdb file will be "
                           "written where bfactor=N*(1-Qscore) for each atom.")
        form.addParam('sigma', FloatParam, allowsNull=True,
                      label="Sigma",
                      help="Optional – default is 0.6 – specifies width of reference Gaussian in Å ")
        form.addParallelSection(threads=4, mpi=0)

    # --------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.convertInputStep)
        self._insertFunctionStep(self.computeQScoresStep)
        self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions -------------------------------
    def convertInputStep(self):
        volFile = self.inputVol.get().getFileName()
        sampling = self.inputVol.get().getSamplingRate()
        self.volOutFile = abspath(self._getExtraPath('map.mrc'))
        Ccp4Header.fixFile(volFile, self.volOutFile, (0, 0, 0), sampling,
                           Ccp4Header.START)

        self.pdbOutFile = []
        for pdb in self.pdbs:
            pdbFile = pdb.get().getFileName()
            baseName = pwutils.removeBaseExt(pdbFile)
            self.pdbOutFile.append(abspath(self._getExtraPath('%s.pdb' % baseName)))
            if pwutils.getExt(pdbFile) == ".pdb":
                h = AtomicStructHandler()
                h.read(pdbFile)
                h.writeAsPdb(self.pdbOutFile[-1])
            else:
                toPdb(pdbFile, self.pdbOutFile[-1])

    def computeQScoresStep(self):
        args = '%s %s ' % (mapq.Plugin.getChimeraPath(), self.volOutFile)
        args += " ".join(self.pdbOutFile)

        if self.mapRes.get():
            args += " res=%f" % self.mapRes.get()

        if self.bFactor.get():
            args += ' bfactor=%d' % self.bFactor.get()

        if self.sigma.get():
            args += ' sigma=%f' % self.sigma.get()

        args += " np=%d" % self.numberOfThreads.get()

        python_file, mapq_file = mapq.Plugin.getMapQProgram()
        self.runJob(python_file, mapq_file + " " + args)

    def createOutputStep(self):
        outStructFileBase = self._getExtraPath('{}.cif')
        ASH = AtomicStructHandler()
        outSet = SetOfAtomStructs.create(self._getPath())
        for pdbFile in self.pdbOutFile:
            baseName = pwutils.removeBaseExt(pdbFile)
            outStructFileName = outStructFileBase.format(baseName)
            ASH.read(glob(self._getExtraPath(baseName + "*.pdb"))[1])
            mapQ_dict = {'{}:{}'.format(atom.full_id[2], atom.serial_number): str(round(atom.bfactor, 4))
                         for atom in ASH.getStructure().get_atoms()}
            inpAS = toCIF(pdbFile, outStructFileName)
            cifDic = ASH.readLowLevel(inpAS)
            cifDic = addScipionAttribute(cifDic, mapQ_dict, self._ATTRNAME)
            ASH._writeLowLevel(outStructFileName, cifDic)

            outAS = AtomStruct()
            outAS.setFileName(outStructFileName)
            outSet.append(outAS.clone())

        self._defineOutputs(scoredStructures=outSet)
        for pdb in self.pdbs:
            self._defineSourceRelation(pdb, outSet)

    # --------------------------- INFO functions ------------------------------
    def _methods(self):
        methods = []
        methods.append('QScore computation using MapQ')
        return methods

    def _summary(self):
        summary = []
        if not self.isFinished():
            summary.append("QScores not ready yet.")

        if self.getOutputsSize() >= 1:
            summary.append("*Mean Q-Scores:*")
            ASH = AtomicStructHandler()
            for struct in self.scoredStructures:
                fileName = struct.getFileName()
                fields = ASH.readLowLevel(fileName)
                attributes = fields["_scipion_attributes.name"]
                values = fields["_scipion_attributes.value"]
                mapq_scores = [float(value) for attribute, value in zip(attributes, values)
                               if attribute == self._ATTRNAME]
                mean_score = sum(mapq_scores) / len(mapq_scores)
                summary.append("      - %s --> %.4f" % (pwutils.removeBaseExt(fileName), mean_score))
        return summary
