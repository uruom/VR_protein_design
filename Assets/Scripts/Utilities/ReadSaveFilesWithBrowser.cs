/*
    ================================================================================
    Copyright Centre National de la Recherche Scientifique (CNRS)
        Contributors and copyright holders :

        Xavier Martinez, 2017-2021
        Marc Baaden, 2010-2021
        baaden@smplinux.de
        http://www.baaden.ibpc.fr

        This software is a computer program based on the Unity3D game engine.
        It is part of UnityMol, a general framework whose purpose is to provide
        a prototype for developing molecular graphics and scientific
        visualisation applications. More details about UnityMol are provided at
        the following URL: "http://unitymol.sourceforge.net". Parts of this
        source code are heavily inspired from the advice provided on the Unity3D
        forums and the Internet.

        This program is free software: you can redistribute it and/or modify
        it under the terms of the GNU General Public License as published by
        the Free Software Foundation, either version 3 of the License, or
        (at your option) any later version.

        This program is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
        GNU General Public License for more details.

        You should have received a copy of the GNU General Public License
        along with this program. If not, see <https://www.gnu.org/licenses/>.

        References : 
        If you use this code, please cite the following reference :         
        Z. Lv, A. Tek, F. Da Silva, C. Empereur-mot, M. Chavent and M. Baaden:
        "Game on, Science - how video game technology may help biologists tackle
        visualization challenges" (2013), PLoS ONE 8(3):e57990.
        doi:10.1371/journal.pone.0057990
       
        If you use the HyperBalls visualization metaphor, please also cite the
        following reference : M. Chavent, A. Vanel, A. Tek, B. Levy, S. Robert,
        B. Raffin and M. Baaden: "GPU-accelerated atom and dynamic bond visualization
        using HyperBalls, a unified algorithm for balls, sticks and hyperboloids",
        J. Comput. Chem., 2011, 32, 2924

    Please contact unitymol@gmail.com
    ================================================================================
*/


using SimpleFileBrowser;
using System.Collections.Generic;
using UnityEngine;
using System.Collections;
using System.IO;
using UnityEngine.XR;
using SFB;

namespace UMol
{

    public class ReadSaveFilesWithBrowser : MonoBehaviour
    {
        public string lastOpenedFolder = "";
        public string initPath = "";
        public string extension = "";

        void loadFileFromPath(string path, bool readHetm)
        {

#if UNITY_STANDALONE_OSX || UNITY_EDITOR_OSX
        path = path.Replace("file:/", "");
#endif
            if (!string.IsNullOrEmpty(path))
            {
                if (path.EndsWith(".xtc"))
                {
                    string lastStructureName = API.APIPython.last().uniqueName;
                    if (lastStructureName != null)
                    {
                        API.APIPython.loadTraj(lastStructureName, path);
                    }
                }
                else if (path.EndsWith(".dx"))
                {
                    string lastStructureName = API.APIPython.last().uniqueName;
                    if (lastStructureName != null)
                    {
                        API.APIPython.loadDXmap(lastStructureName, path);
                    }
                }
                else if (path.EndsWith(".py"))
                {
                    API.APIPython.loadHistoryScript(path);
                }
                else if (path.EndsWith(".itp"))
                {
                    //WARNING Not published yet
                }
                else if (path.EndsWith(".conf"))
                {
                    API.APIPython.StartNAMD(path);
                }
                else
                {
                    API.APIPython.load(path, readHetm,false,true,true);
                }
            }
        }


        public void readFiles(bool readHetm = true, bool forceDesktop = false)
        {
            // Debug.Log("!!!!!path get");
            // return ;
            string[] paths = filesToRead(initPath, extension, readHetm, forceDesktop);
            if (paths != null && paths.Length != 0)
            {
                if (paths[0] != "")
                    initPath = Path.GetDirectoryName(paths[0]);
                for (int i = 0; i < paths.Length; i++)
                {
                    Debug.Log("!!!!!path"+paths[i]);
                    loadFileFromPath(paths[i], readHetm);
                }
            }
        }

        public void saveState()
        {
            string path = stateToRead(initPath);
            if (path != null)
            {
                API.APIPython.saveHistoryScript(path);
            }
        }

        public void readState()
        {
            string[] paths = filesToRead(initPath, "py");
            if (paths != null && paths.Length > 0)
            {
#if UNITY_STANDALONE_OSX || UNITY_EDITOR_OSX
            paths[0] = paths[0].Replace("file:/", "");
#endif
                API.APIPython.loadHistoryScript(paths[0]);
            }
        }

        public string stateToRead(string initPath = "")
        {
#if (UNITY_STANDALONE_OSX || UNITY_EDITOR_OSX || UNITY_STANDALONE_WIN || UNITY_EDITOR_WIN || UNITY_WEBGL)
            var extensions = new[]
            {
            new ExtensionFilter("UnityMol State Files", "py" )
        };

            return StandaloneFileBrowser.SaveFilePanel("Save UnityMol State", initPath, "UMolState.py", extensions);

#else
        StartDialogSave();
#endif
            return null;
        }
        public string[] filesToRead(string initPath = "", string extension = "", bool readHetm = true, bool forceDesktop = false)
        {
            var extensions = new[]
            {
            new ExtensionFilter("Molecule Files", "pdb", "mmcif", "cif", "gro", "mol2", "sdf", "mol", "xyz"),
            new ExtensionFilter("Trajectory Files", "xtc"),
            new ExtensionFilter("Density map Files", "dx"),
            new ExtensionFilter("State/Script Files", "py"),
            new ExtensionFilter("Martini itp Files", "itp"),
            new ExtensionFilter("All Files", "*"),
        };

            string[] paths = null;
            //Use native file browser for Windows and Mac and WebGL (https://github.com/gkngkc/UnityStandaloneFileBrowser)
#if (UNITY_STANDALONE_OSX || UNITY_EDITOR_OSX || UNITY_STANDALONE_WIN || UNITY_EDITOR_WIN || UNITY_WEBGL)
            if (!UnityMolMain.inVR() || forceDesktop)
            {
                if (extension == "")
                {
                    paths = StandaloneFileBrowser.OpenFilePanel("Open File", initPath, extensions, true);
                }
                else if (extension == "*")
                {
                    paths = StandaloneFileBrowser.OpenFilePanel("Open File", initPath, "", true);
                }
                else
                {
                    paths = StandaloneFileBrowser.OpenFilePanel("Open File", initPath, extension, true);
                }
            }
            else
            {
                FileBrowser.SetFilters(true, new FileBrowser.Filter("Supported", ".pdb", ".cif", ".mmcif", ".gro", ".mol2", ".xyz", ".sdf", ".mol", ".py", ".dx", ".xtc", ".itp"));
                FileBrowser.SetDefaultFilter(".pdb");
                StartCoroutine(ShowLoadDialogCoroutine("to be decided by loadFileFromPath()",readHetm));
            }
#else //Use asset based file browser (https://github.com/yasirkula/UnitySimpleFileBrowser)
        //Uses a coroutine
        StartDialog(readHetm);
#endif
            return paths;
        }
        public string pathForConfRead()
        {
            string path = null;
            ExtensionFilter[] _extensions =
            {
                new ExtensionFilter("Config files","conf")
            };
            if (!UnityMolMain.inVR())
            {
                string[] paths = StandaloneFileBrowser.OpenFilePanel("Open File", initPath, _extensions, false);
                if(paths.Length>0)
                    path = paths[0];
            }
            else
            {
                FileBrowser.SetFilters(true, new FileBrowser.Filter("Config files", ".conf"));
                FileBrowser.SetDefaultFilter(".conf");
                StartCoroutine(ShowLoadDialogCoroutine("read namd conf",false));
            }
            return path;
        }

        public string pathForPDBSave()
        {
            string path = null;
            ExtensionFilter[] _extensions =
            {
                new ExtensionFilter("Molecule files","pdb")
            };
            if (!UnityMolMain.inVR())
            {
                path = StandaloneFileBrowser.SaveFilePanel("Save File", initPath, "untitled_molecule", _extensions);
            }
            else
            {
                FileBrowser.SetFilters(true, new FileBrowser.Filter("Supported", ".pdb"));
                FileBrowser.SetDefaultFilter(".pdb");
                StartCoroutine(ShowSaveDialogCoroutine("save pdb"));
            }
            return path;
        }
        void CoroutineReturnHandler(string action)
        {
            switch (action)
            {
                case "save pdb":
                    API.APIPython.saveSelectionToPDB(FileBrowser.Result);
                    break;
                case "read namd conf":
                    API.APIPython.StartNAMD(FileBrowser.Result);
                    break;
                default:
                    foreach (string p in FileBrowser.Results)
                        loadFileFromPath(p, true);
                    break;
            }
        }
        IEnumerator ShowLoadDialogCoroutine(string action,bool readHetm)
        {
            if (initPath == "")
            {
                string savedFolder = PlayerPrefs.GetString("lastOpenedFolderVR");
                if (!string.IsNullOrEmpty(savedFolder))
                {
                    initPath = savedFolder;
                }
            }
            yield return FileBrowser.WaitForLoadDialog(false, initPath, "Load File", "Load");
            if (FileBrowser.Success)
            {
                initPath = Path.GetDirectoryName(FileBrowser.Result);
                CoroutineReturnHandler(action);
                PlayerPrefs.SetString("lastOpenedFolderVR", initPath);
            }
            else
            {
                Debug.LogError("Could not load selected file");
            }
        }
        IEnumerator ShowSaveDialogCoroutine(string action)
        {
            if (initPath == "")
            {
                string savedFolder = PlayerPrefs.GetString("lastOpenedFolderVR");
                if (!string.IsNullOrEmpty(savedFolder))
                {
                    initPath = savedFolder;
                }
            }
            yield return FileBrowser.WaitForSaveDialog(false, initPath, "Save File", "Save");
            if (FileBrowser.Success)
            {
                CoroutineReturnHandler(action);
                initPath = Path.GetDirectoryName(FileBrowser.Result);
                PlayerPrefs.SetString("lastOpenedFolderVR", initPath);
            }
            else
            {
                Debug.LogError("Could not save to selected file");
            }
        }
    }
}