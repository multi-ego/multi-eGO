# domains.py

`domains.py` is a script used to define intra-region associations for specific energy epsilon parameters in a multi-eGO environment.

## Usage


### Arguments

- **`--intra`**  
  Specifies the path to the reference intra-region matrix file. This file defines the intra-region epsilon parameter relationships.  
  **Example:** `--intra PATH/reference_intramat`

- **`--top`**  
  Specifies the path to the reference topology file. This file contains the necessary topology information for the domain definitions.  
  **Example:** `--top PATH/reference_top`

- **`--dom_res`**  
  Defines the sequence of domain ranges as `start-end` pairs, where each range represents a domain of residues. These ranges must satisfy the following conditions:
  - Be in non-decreasing order.
  - Not overlap with one another.

  **Examples:**  
  Valid domain ranges:  
  - `--dom_res 1-30 45-60`  
  - `--dom_res 4-4 45-60`  

  Invalid domain ranges:  
  - `--dom_res 1-30 25-60` (Overlapping ranges)  
  - `--dom_res 30-20 35-60` (Non-decreasing order)

- **`--out`** *(optional)*  
  Specifies the output directory where the results will be saved. If omitted, the results will be saved in the current working directory.  
  **Example:** `--out output_directory`

- **`--invert`** *(optional)*  
  If this flag is included, the flag assignment is inverted and 0 will be pleaced in the domain ranges instead of 1

---

### Examples of Usage

1. Define domains using two non-overlapping ranges:  
