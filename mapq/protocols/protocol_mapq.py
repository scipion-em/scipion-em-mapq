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
import csv

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
        form.addParam('mapRes', FloatParam, allowsNull=False,
                      default=3.0,
                      label = "Map resolution",
                      help = "Optional - Default is 3.0 - Specifies resolution of map; it is used to output perresidue "
                             "statistics along with expected Q-score at this resolution")
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
        
        # Write the generic scripts that can be reused
        self.generateQtoBScript()
        self.generateTSVScript()

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

            self.generateQScores(baseName)
            

    def generateQScores(self, baseName):
        cxcQscoreFile = self._getChimeraMainScriptFile(baseName)

        # Generate the script that will align (if needed), calculate
        # the Q-Scores and save everything
        self.generateQCoreScript(cxcQscoreFile, self.cifOutFile[-1], baseName)

        # Tell ChimeraX to run the script
        args = f"--nogui --nocolor --script {cxcQscoreFile}"
        self.runJob(mapq.Plugin.getChimeraXProgram(), args)

    def createOutputStep(self):
        outStructFileBase = self._getExtraPath('{}.cif')
        ASH = AtomicStructHandler()
        outSet = SetOfAtomStructs.create(self._getPath())
        for pdb in self.pdbs:
            pdbFile = pdb.get().getFileName()
            baseName = pwutils.removeBaseExt(pdbFile)
            outStructFileName = outStructFileBase.format(baseName)
            mapq_pdb = self._getExtraPath(f"{baseName}_qscore.cif")
            ASH.read(mapq_pdb)
            mapQ_dict = self.createMapQDict(self._getQScoreTSV(baseName))
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
    def generateQtoBScript(self):
        py_scriptFile = self._getChimeraQtoBPythonFile()
        with open(py_scriptFile, 'w') as fh:
            # Open model and map
            fh.write(
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

    def generateTSVScript(self):
        py_scriptFile = self._getChimeraExportPythonFile()
        with open(py_scriptFile, 'w') as fh:
            # TSV dict generation script
            fh.write(
"""
import csv
import sys

from chimerax.atomic import AtomicStructure

def get_structure(session) -> AtomicStructure:
    structures = [
        model
        for model in session.models.list()
        if isinstance(model, AtomicStructure)
    ]

    if len(structures) != 1:
        raise RuntimeError(
            f"Expected exactly one atomic structure, found {len(structures)}"
        )

    return structures[0]

def export_atom_qscores(session, output_path: str) -> None:
    structure = get_structure(session)

    scored = 0
    missing = 0

    with open(output_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")

        writer.writerow([
            "chain_id",
            "serial",
            "residue_number",
            "insertion_code",
            "residue_name",
            "atom_name",
            "alt_loc",
            "qscore",
        ])

        for atom in structure.atoms:
            qscore = getattr(atom, "qscore", None)

            if qscore is None:
                missing += 1
                continue

            residue = atom.residue

            writer.writerow([
                residue.chain_id or "",
                int(atom.serial_number),
                int(residue.number),
                residue.insertion_code or "",
                residue.name,
                atom.name,
                atom.alt_loc or "",
                float(qscore),
            ])

            scored += 1

    session.logger.info(
        f"Exported {scored} atomic Q-scores to {output_path}; "
        f"{missing} atoms had no qscore"
    )

if len(sys.argv) != 2:
    raise RuntimeError(
        "Usage: export_atom_qscores.py OUTPUT_TSV"
    )

export_atom_qscores(session, sys.argv[1])
"""
            )

    def generateQCoreScript(self, fn, inCif, baseName):
        with open(fn, 'w') as fh:
            # Open inputs
            fh.write(f"open {inCif}\n")
            fh.write(f"open {self.volOutFile}\n")
            # Optional alignment
            if self.autoFit.get():
                fh.write("volume #2 origin 0,0,0\n")
                fh.write(f"fitmap #1 inMap #2 resolution {self.mapRes.get()} metric cam shift true rotate true\n")
            # QScore assignment
            fh.write(f"qscore #1 toVolume #2 useGui false assignAttr true logDetails false outputFile {self._getQScoreCSV(baseName)}\n")
            # Generate TSV for Scipion Dictionary generation
            fh.write(f"runscript '{self._getChimeraExportPythonFile()}' '{self._getQScoreTSV(baseName)}'\n")
            # Write the CIF
            fh.write(f"save {self.pdbOutFile[-1]} models #1\n")
            fh.write(f"save {self._getQScoreATTR(baseName)} attrName a:qscore models #1 modelIds false\n")
            fh.write(f"save {self._getChimeraSessionFile(baseName)}\n")
            fh.write("exit\n")

    def moveOriginTo(self, newOrigin, handler):
        centerMass = handler.centerOfMass(geometric=True)
        for atom in handler.getStructure().get_atoms():
            coords = atom.get_coord()
            atom.coord = coords + np.asarray(newOrigin) - np.asarray(centerMass)

    def createMapQDict(self, atom_tsv):
        mapq_dict = {}
        with open(atom_tsv, newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                try:
                    chain_id = row["chain_id"].strip()
                    serial = int(row["serial"])
                    qscore = float(row["qscore"])
                except (KeyError, TypeError, ValueError):
                    continue
                mapq_dict[f"{chain_id}:{serial}"] = qscore
        return mapq_dict


    def createMapQDictOld(self, mapq_pdb):
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
                            # Ojing McOjing: arriba hacemos el truco de guardar
                            # el Q-Score en el campo B-Factor
                            # TODO: in a future, replace this dirty technique with actually parsing
                            # the defattr file saved from Chimera
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

    def _getChimeraMainScriptFile(self, baseName):
        return self._getExtraPath(f"{baseName}_main.cxc")
    
    def _getChimeraQtoBPythonFile(self):
        return self._getExtraPath(f"qscoretobfactor.py")

    def _getChimeraExportPythonFile(self):
        return self._getExtraPath(f"exportTSV.py")
    
    def _getChimeraSessionFile(self, baseName):
        return self._getExtraPath(f"{baseName}.cxs")
    
    def _getQScoreATTR(self, baseName):
        return self._getExtraPath(f"{baseName}.defattr")
    
    def _getQScoreCSV(self, baseName):
        return self._getExtraPath(f"{baseName}_qscores.csv")
    
    def _getQScoreTSV(self, baseName):
        return self._getExtraPath(f"{baseName}_dict.tsv")
