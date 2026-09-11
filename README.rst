======================
PROTAC plugin
======================

**Documentation under development, sorry for the inconvenience**

This is a **Scipion** plugin that models PROTAC-mediated protein-protein
ternary complexes (POI - linker - E3 ligase) using
`PROTAC-Model <https://github.com/gaoqiweng/PROTAC-Model>`_ (Gao et al.,
Zhejiang University): FRODOCK global docking guided by a receptor-interface
site point, filtering of the resulting poses by compatibility with the PROTAC
linker geometry, and optional RosettaDock refinement of the surviving poses.

Current programs implemented:

    - PROTAC-Model ternary complex modeling

==========================
Install this plugin
==========================

You will need to use `Scipion3 <https://scipion-em.github.io/docs/docs/scipion
-modes/how-to-install.html>`_ to run these protocols.

1. **Binary files**

FRODOCK and Rosetta (used optionally for RosettaDock refinement) are **NOT**
downloaded automatically with the plugin: both require accepting a license
(a click-through agreement for FRODOCK, a personal academic login for
Rosetta) that can't be scripted. Point ``FRODOCK_HOME`` (and ``ROSETTA_HOME``,
via the ``scipion-chem-rosetta`` plugin) at your own downloads in
``scipion.conf`` or your shell environment.

ADFRsuite, Vina, Voromqa and FCC have no such license gate and are installed
automatically by the plugin's own ``defineBinaries()``.

2. **Install the plugin in Scipion**

- **Developer's version**

    .. code-block::

        git clone https://github.com/danigual/scipion-chem-protac.git
        cd scipion-chem-protac
        scipion3 installp -p . --devel
