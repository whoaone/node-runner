# Node Runner

*A lightweight, general-purpose pre-processor for Nastran, written in Python.*
---

## 📖 About

Node Runner is a desktop application designed to simplify the creation, editing, and visualization of Finite Element Analysis (FEA) models for Nastran. Built with Python, PySide6, and PyVista, it provides an intuitive graphical interface for common pre-processing tasks.

Whether you're a student learning FEA, an engineer needing a quick tool for simple models, or a hobbyist exploring structural analysis, Node Runner aims to be an accessible and straightforward solution.

## ✨ Key Features

* **Full Lifecycle Editing:** Create, view, edit, and delete core Nastran entities.
* **Entities Supported:** Nodes, Elements (CQUAD4, CTRIA3, CBEAM, CBAR, CROD, RBE2), Properties (PSHELL, PCOMP, PBEAM), and Materials (MAT1, MAT8).
* **Loads & Constraints:** Full support for creating, editing, and deleting Loads (FORCE, MOMENT, PLOAD4) and nodal Constraints (SPC).
* **Advanced Visualization:** An interactive 3D view powered by PyVista with a configurable, relative scaling system for load and constraint vectors.
* **File I/O:** Open and save standard Nastran Bulk Data Files (`.bdf`, `.dat`).
* **Parametric Generator:** Includes a tool to quickly generate a complete fuselage and floor structure model.
* **Customizable UI:** Features both a light and dark theme.

## 💻 Tech Stack

* **Language:** Python 3
* **GUI:** PySide6 (the official Qt for Python library)
* **3D Rendering:** PyVista
* **Nastran Core:** pyNastran

## 🚀 Getting Started

**1. Clone the repository:**
```bash
git clone [https://github.com/your-username/node-runner.git](https://github.com/whoaone/node-runner.git)
cd node-runner
