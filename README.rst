======================
PROTAC plugin
======================

**Documentation under development, sorry for the inconvenience**

This is a **Scipion** plugin to model PROTAC-mediated ternary complexes
(protein of interest - PROTAC - E3 ligase). It wraps two independent modeling
pipelines and adds a helper protocol to prepare their inputs.

Current programs implemented:

    - **PROTAC-Model** (`PROTAC-Model <https://github.com/gaoqiweng/PROTAC-Model>`_,
      Zhejiang University): FRODOCK global docking guided by a site on the
      receptor, filtering of the poses by compatibility with the PROTAC, and
      optional RosettaDock refinement.
    - **PRosettaC** (`PRosettaC <https://github.com/LondonLab/PRosettaC>`_,
      Zaidman et al., 2020): PatchDock global docking constrained by the sampled
      distance between the two anchor atoms, Rosetta local docking, constrained
      PROTAC conformations and clustering of the final models.
    - **Warhead transplant**: copies a warhead pose from a structure where it is
      already bound to another structure of the same or a close homologous
      protein, by superposing the binding pocket.

Both modeling protocols need each protein with its warhead already bound, and a
PROTAC SMILES in which each warhead is an exact substructure.

==========================
Install this plugin
==========================

You will need to use `Scipion3 <https://scipion-em.github.io/docs/docs/scipion
-modes/how-to-install.html>`_ to run these protocols.

1. **Binary files**

Installed automatically by the plugin (``scipion3 installb <name>`` or the
plugin manager): FRODOCK, ADFRsuite, Vina, VoroMQA, FCC, the PROTAC-Model and
PRosettaC repositories, and the conda environments both pipelines run in.

Two tools need a manual installation, because their licenses require a
personal registration:

- **Rosetta**, for PRosettaC and for the optional RosettaDock refinement of
  PROTAC-Model. It is shared with the ``scipion-chem-rosetta`` plugin: point
  ``ROSETTA_HOME`` at your installation.
- **PatchDock** (https://bioinfo3d.cs.tau.ac.il/PatchDock/), for PRosettaC:
  point ``PATCHDOCK_HOME`` at your installation.

Both variables can be set in ``scipion.conf`` or in the shell environment.

2. **Install the plugin in Scipion**

- **Developer's version**

    .. code-block::

        git clone https://github.com/danigual/scipion-chem-protac.git
        cd scipion-chem-protac
        scipion3 installp -p . --devel

3. **Tests**

    .. code-block::

        scipion3 test protac.tests.test_protac_model
        scipion3 test protac.tests.test_prosettac

The end-to-end tests skip themselves when the external tools are not installed.
