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

from pwem.emlib.image import ImageHandler
from pwem.objects import FSC, Volume
from pwem.protocols import ProtAnalysis3D

from pyworkflow.protocol import PointerParam, FloatParam, MultiPointerParam
from pyworkflow import BETA
import pyworkflow.utils as pwutils

import mapq


class ProtMapQ(ProtAnalysis3D):
    """
    Compute Q-Scores using MapQ software.
    """
    _label = 'compute q-scores'
    _devStatus = BETA

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
        form.addParam('bFactor', FloatParam, allowsNull=True,
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
        self.volOutFile = abspath(self._getExtraPath('map' + pwutils.getExt(volFile)))
        pwutils.copyFile(volFile, self.volOutFile)

        self.pdbOutFile = []
        for idp, pdb in enumerate(self.pdbs):
            pdbFile = pdb.get().getFileName()
            self.pdbOutFile.append(abspath(self._getExtraPath('pdb_%d' % (idp + 1) + pwutils.getExt(pdbFile))))
            pwutils.copyFile(pdbFile, self.pdbOutFile[-1])

    def computeQScoresStep(self):
        args = '%s %s ' % (mapq.Plugin.getChimeraPath(), self.volOutFile)
        args += " ".join(self.pdbOutFile)

        if self.mapRes.get():
            args += " res=%f" % self.mapRes.get()

        if self.bFactor.get():
            args += ' bfactor=%f' % self.bFactor.get()

        if self.sigma.get():
            args += ' sigma=%f' % self.sigma.get()

        args += " np=%d" % self.numberOfThreads.get()

        python_file, mapq_file = mapq.Plugin.getMapQProgram()
        self.runJob(python_file, mapq_file + " " + args)

    def createOutputStep(self):
        pass

    # --------------------------- INFO functions ------------------------------
    def _methods(self):
        methods = []
        methods.append('QScore computation using MapQ')
        return methods

    def _summary(self):
        summary = []
        if not self.isFinished():
            summary.append("QScores not ready yet.")

        # if self.getOutputsSize() >= 1:
        #     stdout = self._getLogsPath('run.stdout')
        #     with open(stdout) as file:
        #         for num, line in enumerate(file, 1):
        #             if 'Resolution at 1 % FDR-FSC' in line:
        #                 res = [float(s) for s in line.split() if s.replace(".", "", 1).isdigit()][1]
        #     summary.append('Resolution at 1 %% FDR-FSC: %.2f Angstrom' % res)
        return summary
