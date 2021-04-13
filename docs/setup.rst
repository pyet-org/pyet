Setup
=====

To use this template for other projects:

1. Get a copy of this project by cloning_, forking_, copying_, or downloading_.

2. Adjust things as appropriate. In particular ``setup.py`` and all the ``*.rst`` files should be adjusted as needed.

3. Build locally by::

     $ cd doc
     $ make html  (or make live for livehtml refreshing)

4. If all is good, commit back to your git repo

5. For serving through readthedocs_:

   a) Add your repo as a new RTD project

   b) Visit the RTD admin page for the project

   c) Under "Advanced Settings":

      * Install Project: unchecked
      * Requirements file = ``doc/requirements.txt``
      * Python configuration file = ``doc/source/conf/py``

   It should then build. Check the build log output for hints if it doesn't work.



.. _cloning: https://github.com/datadavev/sphinx-doc-template.git
.. _forking: https://github.com/datadavev/sphinx-doc-template#fork-destination-box
.. _copying: https://github.com/datadavev/sphinx-doc-template
.. _downloading: https://github.com/datadavev/sphinx-doc-template/archive/master.zip
.. _readthedocs: https://readthedocs.org/
