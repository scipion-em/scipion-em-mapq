# **************************************************************************
# *
# * Authors:    David Herreros Calero (dherreros@cnb.csic.es) [1]
# *
# * [1] Centro Nacional de Biotecnologia, CSIC, Spain
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


from os.path import join
import numpy as np

from pwem.protocols import ProtImportPdb, ProtImportVolumes
from pwem.convert.atom_struct import AtomicStructHandler

from pyworkflow.tests import BaseTest, setupTestProject

from mapq.protocols import ProtMapQ
import mapq


class TestMapQ(BaseTest):
    """This class check if MapQ protocol works properly"""

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.pdb = join(mapq.Plugin.getHome('QScore_Apoferritin_Tutorial'), '3ajo_fitted_refined.pdb')
        cls.map = join(mapq.Plugin.getHome('QScore_Apoferritin_Tutorial'), 'emd20026_prot.mrc')

    def runImportPDBs(cls, label):
        """ Run an Import particles protocol. """
        protImport = cls.newProtocol(ProtImportPdb,
                                     inputPdbData=1,
                                     pdbFile=cls.pdb,
                                     objLabel=label)
        cls.launchProtocol(protImport)
        return protImport.outputPdb

    def runImportVolumes(cls, samplingRate, label):
        """ Run an Import particles protocol. """
        protImport = cls.newProtocol(ProtImportVolumes,
                                     filesPath=cls.map, samplingRate=samplingRate, objLabel=label)
        cls.launchProtocol(protImport)
        return protImport.outputVolume

    def runMapQ(self, pdb, volume):
        prot = self.newProtocol(ProtMapQ, inputVol=volume, pdbs=[pdb])
        self.launchProtocol(prot)
        self.assertIsNotNone(prot.scoredStructures,
                             "There was a problem with MapQ protocol output")
        return prot

    def test_mapq(self):
        pdb = self.runImportPDBs('PDB Struct')
        volume = self.runImportVolumes(0.65, 'Map')
        prot = self.runMapQ(pdb, volume)

        ASH = AtomicStructHandler()
        for struct in prot.scoredStructures:
            fileName = struct.getFileName()
            fields = ASH.readLowLevel(fileName)
            attributes = fields["_scipion_attributes.name"]
            values = fields["_scipion_attributes.value"]
            mapq_scores = [float(value) for attribute, value in zip(attributes, values)
                           if attribute == prot._ATTRNAME]
            mean_score = sum(mapq_scores) / len(mapq_scores)
        self.assertEqual(round(mean_score, 4), -0.0052, "Unexpected score value: mean")

        return prot