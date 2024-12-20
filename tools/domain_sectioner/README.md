# domains.py

Use domains.py to define the intra-region to associated a specific energy epsilon parameter in multi-eGO

## Usage:

```
domains.py --intra PATH/reference_intramat --top PATH/reference_top --dom_res r1_start-r1_end r2_start-r2_end (--out ouput_directory)
```

Sequence of domain ranges must be in non-decreasing order and non overlapping
Examples:
V  --dom_res 1-30 45-60 
V  --dom_res 4-4 45-60 
X  --dom_res 1-30 25-60 
X  --dom_res 30-20 35-60 
