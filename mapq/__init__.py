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

__version__ = "1.0.0"
_logo = "mapq_logo.png"
_references = ['Pintilie2020']
_url = "https://github.com/scipion-em/scipion-em-mapq"


SCRATCHDIR = pwutils.getEnvVariable('SPOCSCRATCHDIR', default='/tmp/')


class Plugin(pwem.Plugin):
    _homeVar = mapqConst.MAPQ_HOME
    _pathVars = [mapqConst.MAPQ_HOME]
    _supportedVersions = mapqConst.V1_16_1
    _currentVersion = mapqConst.V1_16_1

    @classmethod
    def _defineVariables(cls):
        cls._defineEmVar(mapqConst.MAPQ_HOME, 'chimera-%s' % mapqConst.V1_16_1)

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
        return str(cls.getHome('chimera'))

    @classmethod
    def getMapQProgram(cls):
        """ Return the program binary that will be used. """
        cmd = cls.getHome('mapq', 'mapq_cmd.py')
        return sys.executable, str(cmd)

    @classmethod
    def isVersionActive(cls):
        return cls.getActiveVersion().startswith(mapqConst.V1_16_1)

    @classmethod
    def defineBinaries(cls, env):
        # from scipion.install.funcs import VOID_TGZ

        # cls.defineChimeraXInstallation(env, V1_1, default=True)
        cls.defineChimeraInstallation(env, cls._currentVersion, default=True)

        # # Scipion plugin for chimera. It will depend on the version currently active
        # pathToPlugin = os.path.join(os.path.dirname(__file__),
        #                             "Bundles", "scipion")
        # pathToBinary = cls.getProgram()
        #
        # activeVersion = mapqConst.V1_16_1
        # installationFlagFile = "chimera-%s/installed-%s" % (activeVersion, activeVersion)
        #
        # installPluginsCommand = [("%s --nogui --exit " \
        #                           "--cmd 'devel install %s' && touch %s" % (
        #                           pathToBinary, pathToPlugin, installationFlagFile),
        #                           [installationFlagFile])]
        #
        # env.addPackage('scipionchimera', version='1.3',
        #                tar=VOID_TGZ,
        #                default=True,
        #                commands=installPluginsCommand)

    @classmethod
    def defineChimeraInstallation(cls, env, version, default=False):
        from scipion.install.funcs import VOID_TGZ  # Local import to avoid having scipion-app installed when building the package.

        getchimera_script = os.path.join(os.path.dirname(__file__),
                                         "getchimera.py")

        chimera_cmds = []
        chimera_cmds.append(("cd .. && python %s %s" % (getchimera_script, version),
                             "../chimera-%s-linux_x86_64.bin" % version))
        chimera_cmds.append(("chmod +x ../chimera-%s-linux_x86_64.bin && "
                             "printf './chimera\\nno\\n\\n' | ../chimera-%s-linux_x86_64.bin"
                             % (version, version),
                             "chimera"))

        chimera_cmds.append(('wget -c https://github.com/gregdp/mapq/raw/master/download/mapq_1_8_2.zip',
                             'mapq_1_8_2.zip'))
        chimera_cmds.append(('unzip mapq_1_8_2.zip', "mapq"))
        chimera_cmds.append(("cd mapq && python install.py ../chimera &&"
                             "touch ../mapq_installed", "mapq_installed"))

        env.addPackage('chimera', version=version,
                       tar=VOID_TGZ,
                       default=default,
                       commands=chimera_cmds,
                       )
