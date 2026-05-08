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
    Computes atomic Q-scores by comparing atomic structures against a
    cryo-EM density map using MapQ.

    AI Generated:

    Compute Q-Scores (ProtMapQ) — User Manual
        Overview

        The Compute Q-Scores protocol evaluates how well one or more
        atomic structures agree locally with an experimental cryo-EM
        density map.

        Its main purpose is to calculate per-atom Q-scores using the
        MapQ software.

        In structural biology workflows, Q-scores provide a local
        measure of resolvability and model-to-map agreement. Higher
        Q-scores generally indicate that atoms are better supported by
        the experimental density.

        From a biological perspective, this protocol helps identify
        well-resolved regions, poorly supported residues, flexible
        domains, and local inconsistencies between the atomic model and
        the experimental map.

        Inputs and General Workflow

        The protocol requires:

            - One input cryo-EM density map.
            - One or more atomic structures.

        Several optional parameters can also be provided:

            - Map resolution
            - B-factor scaling
            - Sigma of the reference Gaussian
            - Automatic fitting of structures into the map

        During execution, the protocol performs three main stages:

            1. Input conversion and preparation
            2. Q-score computation with MapQ
            3. Output generation with embedded per-atom scores

        Input Preparation

        The input map is first converted into a standardized MRC file.

        During this step, the protocol preserves both:

            - Sampling rate
            - Origin coordinates

        This is important because Q-score calculations are highly
        sensitive to spatial consistency between map and atomic model.

        Each input structure is converted into CIF format.

        If automatic fitting is enabled, each structure is also fitted
        into the map using ChimeraX before Q-score computation.

        Automatic Fitting

        When autoFit is enabled, the protocol performs an initial rigid
        fit of each structure into the density map.

        ChimeraX is used to:

            - Open the structure
            - Open the map
            - Run fitmap
            - Save the fitted structure

        This option is useful when the input models are approximately
        correct but not perfectly aligned with the density.

        If the map and structures are already aligned, automatic fitting
        can be disabled.

        From a practical perspective, disabling fitting avoids
        unnecessary geometric perturbation when structures have already
        been carefully positioned.

        Q-Score Computation

        After input preparation, the protocol launches MapQ.

        MapQ evaluates the local density surrounding every atom and
        compares it with an idealized reference Gaussian profile.

        The protocol can optionally use:

            - Map resolution, which helps estimate expected Q-scores at
              that resolution.
            - Sigma, which controls the width of the reference Gaussian.
            - B-factor scaling, which generates an auxiliary PDB where
              atomic B-factors are proportional to 1 - Qscore.

        Computation is parallelized using the number of selected CPU
        threads.

        Interpretation of Q-Scores

        Q-scores are computed per atom.

        In general:

            - High Q-scores indicate strong local density support.
            - Low Q-scores indicate weak density support, disorder, or
              possible local modelling inaccuracies.

        Biologically, low-scoring regions often correspond to flexible
        loops, mobile domains, poorly resolved side chains, or
        heterogeneous conformational states.

        Conversely, highly ordered structural cores usually produce
        higher Q-scores.

        Output Generation

        After MapQ finishes, the protocol reads the resulting scored PDB
        files and extracts the per-atom Q-score values.

        These values are then inserted into the output CIF files as
        Scipion atom-level attributes.

        The final output is a SetOfAtomStructs.

        Each output structure preserves the original atomic model while
        adding per-atom MapQ scores.

        Outputs and Their Interpretation

        The protocol generates one scored atomic structure for each
        input model.

        Each atom contains an associated attribute called:

            MapQ_Score

        This allows direct downstream analysis of local map agreement.

        In addition, the protocol reports the mean Q-score for each
        structure in the protocol summary.

        Biologically, the mean Q-score provides a global estimate of
        model-to-map consistency, while the per-atom values reveal local
        structural reliability.

        Practical Recommendations

        In routine cryo-EM workflows, the protocol performs best when
        structures are already approximately aligned to the density map.

        Automatic fitting is useful for preliminary analyses, but when
        accurate manual fitting has already been performed, disabling
        autoFit may preserve intended structural placement.

        Interpretation of Q-scores should always consider the nominal
        map resolution. Lower local Q-scores do not necessarily imply
        incorrect modelling, especially in flexible or heterogeneous
        regions.

        Visual inspection of low-scoring regions is strongly
        recommended.

        Final Perspective

        For cryo-EM users, this protocol provides a quantitative link
        between atomic modelling and experimental density quality.

        Although computationally simple, Q-score analysis is
        biologically valuable because it highlights which parts of a
        model are strongly supported by the experimental map and which
        parts require more cautious interpretation.
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
        form.addParam('autoFit', BooleanParam, default=True, label="Auto fit map and structures?",
                      help="If true, the map and structures will be automatically aligned with Chimera. "
                           "Otherwise, map and structures will be assumed to be aligned")
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
            self.pdbOutFile.append(abspath(self._getExtraPath('%s.pdb' % baseName)))

            h = AtomicStructHandler()
            h.read(cifFile)
            h.writeAsCif(self.cifOutFile[-1])

            #### CHIMERAX
            if self.autoFit.get():
                print("Fitting %s into map..." % baseName)
                scriptFile = self._getTmpPath("fitting.py")
                with open(scriptFile, 'w') as fhCmd:  # Using 'with' ensures the file is properly closed
                    fhCmd.write("from chimerax.core.commands import run\n")
                    fhCmd.write("run(session, 'open %s')\n" % self.cifOutFile[-1])
                    fhCmd.write("run(session, 'open %s')\n" % self.volOutFile)
                    fhCmd.write("run(session, 'fitmap #1 inMap #2')\n")
                    fhCmd.write("run(session, 'save %s models #1 relModel #2')\n" % self.pdbOutFile[-1])
                    fhCmd.write("run(session, 'exit')\n")  # Ensure ChimeraX exits after running the script
                args = "--nogui --script %s" % scriptFile
                self.runJob(mapq.Plugin.getChimeraXProgram(), args)

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
        for pdb in self.pdbs:
            pdbFile = pdb.get().getFileName()
            baseName = pwutils.removeBaseExt(pdbFile)
            outStructFileName = outStructFileBase.format(baseName)
            mapq_pdb = self._getExtraPath(baseName + ".pdb__Q__map.mrc.pdb")
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
