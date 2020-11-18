# Code for Pymol images

The following are simply a copy paste of the commands typed in PyMol and then cleaned in Python with

    import re
    print('\n'.join(re.findall('PyMOL>(.*?)\n', log)))

Followed by `get_view` to get the orientation (`set_view`).

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
    hide element H
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