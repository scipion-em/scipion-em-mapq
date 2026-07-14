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
import numpy as np

from pwem.convert import toCIF, Ccp4Header
from pwem.convert.atom_struct import toPdb, AtomicStructHandler, addScipionAttribute
from pwem.objects import SetOfAtomStructs, AtomStruct
from pwem.protocols import ProtAnalysis3D

from pyworkflow.protocol import PointerParam, FloatParam, MultiPointerParam, IntParam, BooleanParam
from pyworkflow import BETA
import pyworkflow.utils as pwutils

import mapq


class ProtMapQ(ProtAnalysis3D):
    """
    Compute Q-Scores using MapQ software.
    """
    _label = 'compute q-scores'
    _devStatus = BETA
    _ATTRNAME = "bfactor"
    _OUTNAME = "scoredStructures"

    # --------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputVol', PointerParam, pointerClass="Volume", label="Input map",
                      help='Volume to compare to PDBs', important=True)
        form.addParam('pdbs', MultiPointerParam, pointerClass="AtomStruct", important=True,
                      label="Input structures",
                      help='PDBs to compare to input map')
        form.addParam('mapRes', FloatParam, allowsNull=False,
                      default=3.0,
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
        form.addParam('autoFit', BooleanParam, default=True, label="Auto fit map and structures?",
                      help="If true, the map and structures will be automatically aligned with Chimera. "
                           "Otherwise, map and structures will be assumed to be aligned")
        form.addParallelSection(threads=4, mpi=0)

    # --------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.convertInputStep)
        # self._insertFunctionStep(self.computeQScoresStep)
        self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions -------------------------------
    def convertInputStep(self):
        volFile = self.inputVol.get().getFileName()
        sampling = self.inputVol.get().getSamplingRate()
        origin = self.inputVol.get().getShiftsFromOrigin()
        self.volOutFile = abspath(self._getExtraPath('map.mrc'))
        Ccp4Header.fixFile(volFile, self.volOutFile, origin, sampling,
                           Ccp4Header.START)

        self.cifOutFile = []
        self.pdbOutFile = []
        for pdb in self.pdbs:
            cifFile = pdb.get().getFileName()
            baseName = pwutils.removeBaseExt(cifFile)
            self.cifOutFile.append(abspath(self._getExtraPath('%s.cif' % baseName)))  #############################
            self.pdbOutFile.append(abspath(self._getExtraPath('%s_qscore.cif' % baseName)))

            h = AtomicStructHandler()
            h.read(cifFile)
            h.writeAsCif(self.cifOutFile[-1])

            self.runChimeraX(baseName)
            

    def runChimeraX(self, baseName):
        cxc_scriptFile = self._getTmpPath(f"{baseName}_fitting.cxc")
        py_scriptFile = self._getTmpPath(f"{baseName}_fitting.py")
        qscore_file = abspath(self._getExtraPath(f"{baseName}.csv"))

        with open(py_scriptFile, 'w') as fhCmd:
            # Open model and map
            fhCmd.write(
"""
from chimerax.atomic import AtomicStructure

structures = [ model for model in session.models.list() if isinstance(model, AtomicStructure) ]

if len(structures) != 1:
    raise RuntimeError("Expected one atomic structure, found more!")

structure = structures[0]
                
for atom in structure.atoms:
    qscore = getattr(atom, "qscore", None)
    if qscore is not None:
        atom.bfactor = float(qscore)
"""
            )
            
        with open(cxc_scriptFile, 'w') as fh:
            # Open inputs
            fh.write(f"open {self.cifOutFile[-1]}\n")
            fh.write(f"open {self.volOutFile}\n")
            # Optional alignment
            if self.autoFit.get():
                fh.write("volume #2 origin 0,0,0\n")
                fh.write(f"fitmap #1 inMap #2 resolution 3.0 metric cam shift true rotate true\n")
            # QScore assignment
            fh.write(f"qscore #1 toVolume #2 useGui false assignAttr true logDetails false outputFile {qscore_file}\n")
            # Copy qscore to bfactor
            fh.write(f"runscript '{abspath(py_scriptFile)}'\n")
            
            fh.write(f"save {self.pdbOutFile[-1]} models #1\n")
            fh.write("exit\n")

        # Tell ChimeraX to run the script
        args = f"--nogui --nocolor --script {cxc_scriptFile}"
        self.runJob(mapq.Plugin.getChimeraXProgram(), args)


    def createOutputStep(self):
        outStructFileBase = self._getExtraPath('{}.cif')
        ASH = AtomicStructHandler()
        outSet = SetOfAtomStructs.create(self._getPath())
        for pdb in self.pdbs:
            pdbFile = pdb.get().getFileName()
            baseName = pwutils.removeBaseExt(pdbFile)
            outStructFileName = outStructFileBase.format(baseName)
            mapq_pdb = self._getExtraPath(baseName + "_qscore.cif")
            ASH.read(mapq_pdb)
            mapQ_dict = self.createMapQDict(mapq_pdb)
            inpAS = toCIF(pdbFile, outStructFileName)
            cifDic = ASH.readLowLevel(inpAS)
            cifDic = addScipionAttribute(cifDic, mapQ_dict, self._ATTRNAME, recipient = 'atoms')
            ASH._writeLowLevel(outStructFileName, cifDic)

            outAS = AtomStruct()
            outAS.setFileName(outStructFileName)
            outSet.append(outAS.clone())

        self._defineOutputs(scoredStructures=outSet)
        for pdb in self.pdbs:
            self._defineSourceRelation(pdb, outSet)

    # --------------------------- UTILS functions -------------------------------
    def moveOriginTo(self, newOrigin, handler):
        centerMass = handler.centerOfMass(geometric=True)
        for atom in handler.getStructure().get_atoms():
            coords = atom.get_coord()
            atom.coord = coords + np.asarray(newOrigin) - np.asarray(centerMass)

    def createMapQDict(self, mapq_pdb):
            mapQ_dict = {}
            
            with open(mapq_pdb, 'r') as f:
                for line in f:
                    # El tipo de registro (ATOM/HETATM) está en los primeros 6 caracteres
                    record_type = line[0:6].strip()
                    
                    if record_type in ['ATOM', 'HETATM']:
                        try:
                            # Extraemos los campos usando los índices definidos en tu colspecs original
                            # chainid: posición 21 (índice 21:22)
                            # serial: posición 6-11 (índice 6:11)
                            # Q_score: posición 60-66 (índice 60:66)
                            
                            chain_id = line[21:22].strip()
                            serial = line[6:11].strip()
                            q_score = line[60:66].strip()
                            
                            # Creamos la llave con el formato "chainid:serial"
                            key = f"{chain_id}:{int(serial)}"
                            mapQ_dict[key] = q_score
                            
                        except (ValueError, IndexError):
                            # En caso de que una línea esté mal formateada o vacía
                            continue
                            
            return mapQ_dict

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
