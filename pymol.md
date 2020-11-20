# Code for Pymol images

The following are simply a copy-paste of the commands typed in PyMol and then cleaned in Python with

    import re
    print('\n'.join(re.findall('PyMOL>(.*?)\n', log)))

Followed by `get_view` to get the orientation (`set_view`).
In general they are 100% messy with a lot of indecision in the figure making process.

## common settings
These are present on my machines `.pymolrc` and may or may not be relevant.

    set use_shaders, 0
    set ray_trace_mode, 3
    bg_color white
    alias backless, hide sticks, (name C+N+H+HA+O and not resn pro) or (name C+H+HA+O and resn pro)
    print('cmd "backless" added')
    run ~/.pymol/PyMOLRosettaServer.py
    
## whole
membrane in mustard, protein in teal, ligand in coral


    cmd.alter('chain A', 'resv+=69-1')
    sort
    cmd.alter('chain A and resi 245-9999', 'resv+=283-245')
    sort
    remove resn DUM and byres test-post2 around 1
    remove resn DUM and byres test-post2 around 2
    remove resn DUM and byres test-post2 around 3
    color 0xff7f50, chain A and element C
    color 0x008080, chain A and element C
    show spheres, resn DUM
    color ffdb58, resn DUM
    color 0xffdb58, resn DUM
    color 0xff7f50, chain B and element C
    hide sticks, element H
    set_view (\
         0.785548747,    0.006387103,   -0.618765831,\
        -0.618780136,    0.015995923,   -0.785402596,\
         0.004881511,    0.999851763,    0.016517220,\
         0.000000000,   -0.000000000, -168.780624390,\
        -0.412899017,   -0.417537689,   -2.250850677,\
       125.079566956,  212.481643677,  -20.000000000 )
    set ray_shadows,0
    ray 1280, 807
    
## Poisson-Boltzmann

    set_view (\
        -0.835593224,   -0.010836864,    0.549237430,\
         0.549339592,   -0.020859988,    0.835339010,\
         0.002404293,    0.999723256,    0.023384122,\
         0.000008702,    0.000004761, -227.923904419,\
         1.285539150,   -1.428180814,   -6.046429157,\
       184.222854614,  271.624938965,  -20.000000000 )
       
## Sodium


    color 0x008080, chain A and element C
    color 0xff7f50, chain B and element C
    show sticks, chain B
    hide sticks, element H
    show surface, chain A
    set_view (\
        -0.746522844,    0.664294064,   -0.037595604,\
         0.603245258,    0.651920915,   -0.459442943,\
        -0.280696541,   -0.365665227,   -0.887410283,\
        -0.000034206,   -0.000002660,  -45.226322174,\
        -5.691270828,    0.596085191,   -0.949021637,\
         9.363839149,   81.090744019,  -20.000000000 )
    ray 1280, 914
    
## Variants

    color 0x008080, chain A and element C
    color 0xff7f50, chain B and element C
    hide element H
    hide sticks, element H
    show sticks, chain B
    hide sticks, element H
    show surface, chain A
    hide sticks, element H
    ray 2559, 1615
    ray 2559, 1719
    set h_bond_cutoff_center, 5
    set h_bond_cutoff_edge, 4
    ray 1280, 914
    zoom resi 206
    hide sticks
    show sticks, resi 206
    show sticks, resi 208
    show sticks, resi 296
    375
    show sticks, resi 375
    show sticks, resi 387
    backless
    hide sticks, (name C+N+H+HA+O and not resn pro) or (name C+H+HA+O and resn pro) 
    hide sticks, element H
    select variants, resi 208+296+375+387
    show sticks, variants
    backless
    hide sticks, (name C+N+H+HA+O and not resn pro) or (name C+H+HA+O and resn pro) 
    hide sticks, element H
    color 40e0d0, variants and element C
    color 0x40e0d0, variants and element C
    zoom variants
    ray 2559, 1828
    hide sticks, element H
    ray 2559, 1828
    set_view (\
        -0.558684468,   -0.811633348,   -0.170631483,\
        -0.789526463,    0.583465636,   -0.190271258,\
         0.253988087,    0.028418958,   -0.966787219,\
        -0.000076481,   -0.000045091,  -55.750110626,\
       -12.756890297,   -8.291956902,  -10.460436821,\
        37.776084900,   73.669357300,  -20.000000000 )
    show stick, resn arg and byres resi 391 and resn PO4 around 5
    hide sticks, element H
    color 0x00FFFF, variants and element C
    color 0x7fffd4, variants and element C
    color 0x00FFFF, variants and element C
    color 0x7fffd4, variants and element C
    color 0xf08080, chain A and not variants and element C
    color 0x40e0d0, variants and element C
    color 0x008080, chain A and element C
    color 0x40e0d0, variants and element C
    ray 2559, 1828
    
    
## Redux

    
select targeted, xxx
cmd.alter('chain A and targeted', 'resv+=69-1')
sort
cmd.alter('chain A and resi 245-9999 and targeted', 'resv+=283-245')
sort
color 0xE6E6FA, chain B and element C
show sticks, chain B and not element H
color 0x008080, chain A and element C
remove resn MEM
select variants, resi 208+296+375+387
show sticks, variants
backless
hide sticks, (name C+N+H+HA+O and not resn pro) or (name C+H+HA+O and resn pro) 
hide sticks, element H
show sticks, chain B and not element H
color 0x40e0d0, variants and element C
# turquoise
color 0x00FFFF, variants and element C 
#cyan