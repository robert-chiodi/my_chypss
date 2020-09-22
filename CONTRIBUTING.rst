.. _contributing_label:

Contributing to CHyPSS
======================

Overview
---------
This document provides direction on how to contribute to CHyPSS.
It (and the associated work-flow) are in active development. 
Discussion is welcome on it, and it should be considered a living document. 

To file an issue or request a new feature in CHyPSS, please see the `Filing an Issue or Requesting a New Feature`_ Section.

To help in the development of CHyPSS, potentially answering a filed issue or feature request, please see the section on `Submitting New Code`_.
Whenever contributing to CHyPSS, make sure to adhere to the `Style Guidelines`_ discussed at the end of this document.

Filing an Issue or Requesting a New Feature
-------------------------------------------


Submitting New Code
-------------------
Should be corresponding unit tests and potentially regression tests.


Style Guidelines
----------------
CHyPSS is written in C++11 and heavily utilizes the `MFEM <https://mfem.org/>`_ framework.
As such, it borrows some conventions in order to stay consistent with objects and functions used from MFEM.

The major portions of the style guidelines are:

  - `Repository Structure`_
  - `File Layout`_
  - `Variables`_
  - `Class Organization`_
  - `Free-Functions and Class Methods`_
  - `Templates`_
  - `If-Statements`_
  - `Loops`_
  - `Macros`_

Repository Structure
~~~~~~~~~~~~~~~~~~~~
All files should be placed in the ``chypss`` source folder.
No sub-directories are allowed.
Unit tests that test functionalities in CHyPSS are placed in ``tests/unit`` and regression tests should be placed in ``tests/regression``.

Consistent file extensions should be used.
The type of file and extension are listed in the table below. 
The C++ Template File refers to the definition of template functions (and inlined functions), which must be available to be included for instantiation, but should be separated from the public interface of the class.

+-------------------+----------------+
| File Type         | File Extension |
+===================+================+
| C++ Header File   | ``.hpp``       |
+-------------------+----------------+
| C++ Template File | ``.tpp``       |
+-------------------+----------------+
| C++ Source File   | ``.cpp``       |
+-------------------+----------------+
| C   Header File   | ``.h``         |
+-------------------+----------------+
| C   Source File   | ``.c``         |
+-------------------+----------------+


File Layout
~~~~~~~~~~~~~~~
This section discuses the layout of files in CHyPSS.
The individual file components are presented in the order they should be encountered. 
An example file layout is below.::

	================ EXAMPLE FILE LAYOUT ===================
	#ifndef HEADER_GUARD (if Header/Template file)
	#define HEADER_GUARD (if Header/Template file)

	include "chypss/corresponding_header.h" (if source file)

	#include <system_headers> (alphabetical order)

	#include <third_party_headers> (alphabetical order)

	#include "chypss/chypss_headers.h" (alphabetical order)

	namespace chypss{
	
		MAIN CODE BODY HERE (classes, functions, etc.)

	} // namespace chypss


	#endif // HEADER_GUARD (if Header/Template file)
	=======================================================


**License and Copyright**:  A license and copyright notice should be provided at the top of every file associated with CHyPSS. 
The license statement is for the `Mozilla Public License 2.0 (MPL2) <https://www.mozilla.org/en-US/MPL/2.0/>`_. 
The notice should appear as::

	// This file is part of the Coupled Hypersonic Protected System Simulator (CHyPSS).
	//
	// Copyright (C) YEAR AUTHOR NAME <email_address>
	//
	// This Source Code Form is subject to the terms of the Mozilla Public
	// License, v. 2.0. If a copy of the MPL was not distributed with this
	// file, You can obtain one at https://mozilla.org/MPL/2.0/.


**Header Guards**: (For Header/Template Files Only) Header guards should be named according to the path to the file using all capital letters with a underscore appended. 
As an example, for the files ``example.h`` and ``example.tpp`` the header guards would be ``CHYPSS_EXAMPLE_H_`` and ``CHYPSS_EXAMPLE_TPP_``, respectively. 

**Includes**: Include files should be divided into sections, with each section separated by a space and individually sorted in alphabetical order. 
The individual sections are, in order: System Headers, Third-Party Headers, and CHyPSS headers. 
If the file is a source file, the corresponding header file should be given in its own section at the top, before any System Headers.

**namespace**: 
All classes, functions, and variables should be placed in the ``chypss`` namespace. 
The end of the namespace should be marked by the comment ``// namespace chypss``.

**main code body**: The main body of the code in the file, containing any functions, classes, or variables, should appear in the ``chypss`` namespace. 
Style rules governing the text here are discussed in more detail in their own separate sections, such as `Variables`_, `Class Organization`_, `Free-Functions and Class Methods`_, `Templates`_, and others.


Variables
~~~~~~~~~
In general, variable names should be descriptive and of an appropriate length.
Do not fear having ``long_variable_names_if_needed``.
In the case where domain-specific knowledge leads to an obvious variable name choice, it should be used.
As an example, variables ``A``, ``x``, and ``b`` for representing a linear system.

Variable names should be written in snake_case, with underscores ( _ ) used to separate words in the variable name. In addition, member variables should have an appended ``_m`` to them, such as ``member_variable_example_m``, and arguments to a function should have a ``a_`` prepended to them, such as ``a_function_argument_example``. Following this convention makes it readily apparent whether a variable was passed into a function, is a member variable to a class, or was instantiated inside the function/method itself.

Class Organization
~~~~~~~~~~~~~~~~~~
A class should be named following CamelCase, such as ``ClassNameExample``. 
This follows the convention given in MFEM, which is used heavily in CHyPSS. 
Unless two classes are highly related, each class should be placed in its own file of an equivalent snake_case name, such as ``class_name_example.h``.
The sections of a class should be organized as::

	===================== EXAMPLE CLASS LAYOUT =====================

	ClassName {

	 using declarations

	 public:
	  Public Methods
	         - Constructors
		 - Copy Constructor
		 - Move Constructor
		 - Copy Assignment
		 - Move Assignment
		 - Other methods
		 - Destructor

	  Public Members

	 protected:
	  Protected Methods
	         - Constructors
		 - Copy Constructor
		 - Move Constructor
		 - Copy Assignment
		 - Move Assignment
		 - Other methods
		 - Destructor

	  Protected Members


	 private:
	  Private Methods
	         - Constructors
		 - Copy Constructor
		 - Move Constructor
		 - Copy Assignment
		 - Move Assignment
		 - Other methods
		 - Destructor

	  Private Members
	};

	===============================================================

In general, class members should be kept private to ensure invariants are kept. 
For classes that are largely a collection of data expressed through public members, use a ``struct`` instead.

Classes should be limited to two constructors (a default constructor and a non-default one). 
If the non-default constructor consists of a single argument, it should be marked ``explicit`` to prevent unintentional implicit casting. 
The ``explicit`` keyword can be ignored if casting is desired, however the default behavior should be to prevent it.
Other manners of constructing the class type should be through static methods, such as ``static ClassName OtherConstructionApproach(..arguments);``.
The static method should be appropriately named.

Free-Functions and Class Methods
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Templates
~~~~~~~~~


If-Statements
~~~~~~~~~~~~~

Loops
~~~~~


Macros
~~~~~~
