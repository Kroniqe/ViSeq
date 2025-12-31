Instructions

1. Save this code as aliview_clone.py 

2. Open Command Prompt or Terminal in admin in aliview_clone.py containing folder

3. Copy and paste this exact command into the black terminal window and press Enter. Rebuild your EXE using the special command to include Biopython data/

pyinstaller --noconsole --onefile --name="ViSeq" --collect-all Bio aliview_clone.py

4. Final Setup (Do not skip)
Go to dist: Open the dist folder created inside your project folder.
Move the App: Move ViSeq.exe to your Desktop or a dedicated folder.
Add External Tools: Copy muscle.exe and FastTree.exe (provided) and paste them into the same folder as ViSeq.exe

5. You can now double-click ViSeq.exe to run your fully functional application!

6. This version includes all requested features:
   
UI Overhaul: High-contrast gray list, bold text, ruler scale, and unified scrolling.
Sequence Editing: Drag-and-drop reordering, column/row deletion, and undo/redo (Ctrl+Z).
Primer Design: Strict biological parameters, paired search, visual highlighting, and CSV export.
Phylogeny: FastTree/Neighbor-Joining integration with Tree plotting and Newick export.
Performance: Optimized rendering (no "question marks") and fast search.

📄 License
This project is licensed under the GNU General Public License v3.0 (GPLv3).
This choice is mandated by the project's dependency on PyQt6, which is distributed under GPLv3.

PyQt6: GPLv3
Biopython: Biopython License (Permissive)
NumPy/Matplotlib: BSD/PSF (Permissive)

You are free to use, modify, and distribute this software, provided that any derivative works are also open-source under the same GPLv3 license.
