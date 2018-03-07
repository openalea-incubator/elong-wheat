===================
README for Elong-Wheat
===================

This is Elong-Wheat model, a mechanistic model of leaf elongation for wheat that accounts for the CN status.


Prerequisites
=============

* To run the model:
    * Python >= 2.7, http://www.python.org/
    * NumPy >= 1.11, http://www.numpy.org/
    * Pandas >= 0.18, http://pandas.pydata.org/
* To build the documentation: Sphinx >= 1.1.3, http://sphinx-doc.org/
* To run the tests: Nose >= 1.3.0, http://nose.readthedocs.org/
* To get code coverage testing: Coverage >= 3.6b3, http://nedbatchelder.com/code/coverage/


Installing
==========

Use ``setup.py``::

   python setup.py install

To install in develop mode::

   python setup.py develop


Reading the docs
================

After installing::

   python setup.py build_sphinx

Then, direct your browser to ``_build/html/index.html``.


Testing
=======

To run the tests, use::

    nosetests


Contact
=======

Please send a mail to elong-wheat@groupes.renater.fr.


Contributing
============

#. Check for open issues or open a fresh issue to start a discussion around a
   feature idea or a bug: XXX
#. If you feel uncomfortable or uncertain about an issue or your changes, feel
   free to email elong-wheat@groupes.renater.fr.
