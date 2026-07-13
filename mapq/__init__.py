# **************************************************************************
# *
# * Authors:     David Herreros Calero (dherreros@cnb.csic.es) [1]
# *
# * [1] National Center for Biotechnology (CSIC)
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
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

import os
import sys

import pwem
import pyworkflow.utils as pwutils

import mapq.constants as mapqConst
import glob

__version__ = "1.0.0"
_logo = "mapq_logo.png"
_references = ['Pintilie2020']
_url = "https://github.com/scipion-em/scipion-em-mapq"

SCRATCHDIR = pwutils.getEnvVariable('SPOCSCRATCHDIR', default='/tmp/')

class Plugin(pwem.Plugin):
    _supportedVersions = [mapqConst.MAPQ_DEFAULT_VERSION]
    _currentVersion = mapqConst.MAPQ_DEFAULT_VERSION
    _fullversion = f"mapq-{_currentVersion}"

    @classmethod
    def _defineVariables(cls):
        cls._defineEmVar(mapqConst.MAPQ_HOME, cls._fullversion)
        cls._defineEmVar(mapqConst.MAPQ_CHIMERA_HOME, '/usr/bin', "Chimera (OLD) with MAPQ installed")
        cls._defineEmVar(mapqConst.MAPQ_CHIMERAX_HOME, '/usr/bin', "ChimeraX installation")

    @classmethod
    def getEnviron(cls):
        environ = pwutils.Environ(os.environ)
        d = {}
        d['PATH'] = cls.getHome('chimera', 'bin')
        if "REMOTE_MESA_LIB" in os.environ:
            d["LD_LIBRARY_PATH"] = os.environ['REMOTE_MESA_LIB']
        environ.update(d, position=pwutils.Environ.BEGIN)
        return environ

    @classmethod
    def runChimeraProgram(cls, program, args="", cwd=None):
        """ Internal shortcut function to launch chimera program. """
        env = cls.getEnviron()
        pwutils.runJob(None, program, args, env=env, cwd=cwd)

    @classmethod
    def getChimeraPath(cls):
        chimera_home = cls.getVar(mapqConst.MAPQ_CHIMERA_HOME)
        return chimera_home

    @classmethod
    def getChimeraProgram(cls):
        return cls.getChimeraPath() + "/bin/chimera"

    @classmethod
    def getChimeraXPath(cls):
        chimerax_home = cls.getVar(mapqConst.MAPQ_CHIMERAX_HOME)
        return chimerax_home

    @classmethod
    def getChimeraXProgram(cls):
        return cls.getChimeraXPath() + "/bin/ChimeraX"

    @classmethod
    def getMapQProgram(cls):
        """ Return the program binary that will be used. """
        cmd = cls.getHome('mapq', 'mapq_cmd.py')
        return sys.executable, str(cmd)

    @classmethod
    def isVersionActive(cls):
        return cls.getActiveVersion().startswith(mapqConst.MAPQ_DEFAULT_VERSION)

    @classmethod
    def defineBinaries(cls, env):
        # MapQ binaries
        for ver in cls._supportedVersions:
            cls.addMapQPackage(env, ver, default = (ver == mapqConst.MAPQ_DEFAULT_VERSION))

        MAPQ_CHIMERA_INSTALLED = f"mapq_chimera_{cls._currentVersion}_installed"
        chimeraIntegrateCmd  = []
        chimeraIntegrateCmd.append(f" cd {cls.getVar(mapqConst.MAPQ_HOME)} && ")
        chimeraIntegrateCmd.append(f" python install.py {cls.getVar(mapqConst.MAPQ_CHIMERA_HOME)} && ")
        chimeraIntegrateCmd.append(f" && touch ../{MAPQ_CHIMERA_INSTALLED}")
        chimeraIntegrateCmd.append(('wget -c https://github.com/gregdp/mapq/raw/master/data/QScore_Apoferritin_Tutorial.zip',
                             'QScore_Apoferritin_Tutorial.zip'))
        chimeraIntegrateCmd.append(('unzip QScore_Apoferritin_Tutorial.zip', "QScore_Apoferritin_Tutorial"))

        
        chimeraCmds = [(chimeraIntegrateCmd ,MAPQ_CHIMERA_INSTALLED)]

        env.addPackage(mapqConst.MAPQ_CHIMERA, version=cls._currentVersion,
                       commands=chimeraCmds,
                       default=False)

        # Note: ChimeraX and Chimera installation are NOT managed through this package
        # Advise the user to download and install them, and point to them through the 
        # EM Vars instead!

    @classmethod
    def addMapQPackage(cls, env, version, default = False):
        MAPQ_INSTALLED = f"mapq_{version}_installed"
        installCmd = f"git clone https://github.com/gregdp/mapq {cls._fullversion} "
        installCmd += f" touch {MAPQ_INSTALLED}"
        installationCmds = [(installCmd, MAPQ_INSTALLED)]

        env.addPackage(mapqConst.MAPQ,
                       version=version,
                       commands=installationCmds,
                       default=default)
