```python
def lengthunitconverter(unitfrom:str,unitto:str,value):
    
    ''' conversion ratios of units to m '''
    conratdict={
        'm'         : 1,
        'km'        : 1000,
        'mi'        : 1609.34,
        'Angstrom'  : 1e-10,
        'inch'      : 0.0254,
        'nm'        : 1e-9,
        'micrometer': 1e-6
    }

    return conratdict[unitfrom]/conratdict[unitto] * value
```
